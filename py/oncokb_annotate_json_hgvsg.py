#!/usr/bin/env python

import os,sys
import argparse
import logging
import collections
import json
import requests
import time # -jaq

start_time = time.time() #start timer for elapsed time -jaq

SCRIPT_PATH = os.path.abspath(__file__)
FORMAT = '[%(asctime)s] %(levelname)s %(message)s'
l = logging.getLogger()
lh = logging.StreamHandler()
lh.setFormatter(logging.Formatter(FORMAT))
l.addHandler(lh)
l.setLevel(logging.INFO)
debug = l.debug; info = l.info; warning = l.warning; error = l.error

DESCRIPTION = '''

With a config file and JSON input, iterate over all VARIANTS and add
OncoKB annotation when available. See documentation for details

API calls are to https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange
'''

EPILOG = '''
'''

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG,
  formatter_class=CustomFormatter)

# Config file
parser.add_argument('config')

# Input JSON file
parser.add_argument('json')

# Skip Filtered variants below specified VAF
parser.add_argument('-m', '--min-filtered-vaf', action='store', type=float,
    default=1.0, help='Minimum VAF to annotate Filtered variants')

# Add variant data to oncokb section (for debugging)
parser.add_argument('--include-variant', action='store_true',
    help='Store variant information in oncokb entry')

parser.add_argument('-v', '--verbose', action='store_true',
    help='Set logging level to DEBUG')

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)

# Create an object from Requests response object in the case of 
# an exception-raising or HTTP status not ok API call,
# either containing a string representation of the exception 
# or a status fields from the Requests response object
def get_api_requests(res):
  if 'exception' in res:
    return res
  return { 'status_code': res.status_code, 'reason': res.reason }


# Format GatewaySeq loci and base changes to a comma-separated MAF string 
# compatible with the OncoKB `byGenomicChange` API call.
#
# SNVs, INDELS, and di/tri/quad (maybe more?) nucleotide substitutions
# are handled, but for complex variants with REF,ALT both > 1 and
# len(REF) != len(ALT) this passes through as a complex substitution but
# doesn't appear to work. Looking at 999 variants in our initial set,
# there was one annotation like this picked up by protein change but 
# not genomic. 
#
# TODO for these types of complex variants, querying `byHGVSg` may 
# find more variants
def get_maf_string(variant, columns):
  typ = variant[columns.index('type')]
  chrom = variant[columns.index('chrom')]
  chrom_no_chr = chrom.replace('chr', '')
  pos = int(variant[columns.index('pos')])
  ref = variant[columns.index('ref')]
  alt = variant[columns.index('alt')]
  # SNV
  if (len(ref) == 1) and (len(alt) == 1):
    if typ != 'SNV':
      error('Type mismatch')
      sys.exit(1)
    maf_ref = ref
    maf_alt = alt
    start = pos
    end = pos
  # DELETION
  elif (len(ref) > 1) and (len(alt) == 1):
    if typ != 'INDEL':
      error('Type mismatch')
      sys.exit(1)
    start = pos + 1
    end = pos + len(ref) - 1
    maf_ref = ref[1:]
    maf_alt = '-'
  # INSERTION
  elif (len(ref) == 1) and (len(alt) > 1):
    if typ != 'INDEL':
      error('Type mismatch')
      sys.exit(1)
    start = pos
    end = pos + 1
    maf_ref = '-'
    maf_alt = alt[1:]
  # COMPLEX this might handle DNP/TNP/ONP but probably not 
  # len(ref) != len(alt)
  else: 
    maf_ref = ref
    maf_alt = alt
    start = pos
    end = start + len(alt) - 1

  maf_s = ','.join(
    str(x) for x in (chrom_no_chr, start, end, maf_ref, maf_alt))
  return maf_s
  

def type_mismatch(): # -jaq
  error('Type mismatch')
  sys.exit(1)

  
def get_hgvsg(variant, columns): # -jaq
  variant_type = variant[columns.index('type')]
  chromosome = variant[columns.index('chrom')].replace('chr', '')
  position = int(variant[columns.index('pos')])
  reference = variant[columns.index('ref')]
  alteration = variant[columns.index('alt')]

  hgvsg = f'{chromosome}:g.'

  ### SNV "Substitution"
  if (len(reference) == 1) and (len(alteration) == 1):
    if variant_type != 'SNV':
      type_mismatch()
    hgvsg += f'{position}{reference}>{alteration}'
      
  ### DELETION
  elif (len(reference) > 1) and (len(alteration) == 1):
    if variant_type != 'INDEL':
      type_mismatch()
    position += 1
    if len(reference) == 2: # one nucleotide deletion
      hgvsg += f'{position}del'
    else: # serveral nucleotide deletion
      second_position = position + len(reference) - 2
      hgvsg += f'{position}_{second_position}del'

  ### INSERTION
  elif (len(reference) == 1) and (len(alteration) > 1):
    if variant_type != 'INDEL':
      type_mismatch()
    second_position = position + 1
    hgvsg += f'{position}_{second_position}ins{alteration[1:]}'
    
#  info(f'vartype: {variant_type}, pos: {position}, ref: {reference}, alt: {alteration}, hgvsg: {hgvsg}')
  
  return hgvsg


# returns a tuple of 
#   oncokb_data   OncoKB data or false if request failed
#   res           either str(exception) if exception was raised or 
#                 Requests response object
def get_oncokb(HGVSg, timeout, tumorType):
  url = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg'
  params = {
    'hgvsg': HGVSg,
    'referenceGenome': 'GRCh38'
  }
  if tumorType != None:
    params['tumorType'] = tumorType

  try:
    res = requests.get(url, headers=headers, params=params, timeout=timeout)
#    info(res.url)
  # A timeout or a failed network connection will raise an exception
  # Catch it and return the string version
  except Exception as e:
    return False, { 'exception': str(e) }  

  if not res.ok:
    oncokb_data = False
  else:
    oncokb_data = res.json()
  return oncokb_data, res

# Check for presence of required keys and types in the config file
# Hard exit if anything is missing
def check_gs_config(config, config_p):
  gs_config_constraints = {
    'oncokb_api_key': str,
    'oncokb_api_timeout': int,
    'gs_oncokb_tumor_type_map': dict
  }
  errors = []

  for k in gs_config_constraints.keys():
    if k not in config or \
      type(config[k]) != gs_config_constraints[k]:
      errors.append(f'{config_p}: {k} missing or invalid')

  if errors:
    for e in errors:
      error(e)
    sys.exit(1)


debug('%s begin', SCRIPT_PATH)

info(f'Config file: {args.config}')
with open(args.config) as fh:
  gs_config = json.loads(fh.read())
check_gs_config(gs_config, args.config)

tumor_types = []
for gs_tumor_type_values in gs_config['gs_oncokb_tumor_type_map'].values():
  tumor_types += gs_tumor_type_values
info(f"GS tumor types n: {len(gs_config['gs_oncokb_tumor_type_map'].keys())}")
info(f"OncoKB tumor types n: {len(tumor_types)}")
info('Duplicates are:')
info([i for i, count in collections.Counter(tumor_types).items() if count > 1])
info(f"Unique OncoKB tumor types n: {len(set(tumor_types))}")

tumor_types = list(set(tumor_types))
info(f"Final tumor types n: {len(tumor_types)}")

api_key = gs_config['oncokb_api_key']
headers = {
  'accept': 'application/json',
  'authorization': f'Bearer {api_key}'
}

tiers = ['PASS', 'Filtered']
#tiers = ['PASS']

info(f'tumor_types_n: {len(tumor_types)}')

info(f'Input file: {args.json}')
with open(args.json) as fh:
  gs_data = json.loads(fh.read())

columns = gs_data['VARIANTS']['PASS']['columns']

oncokb_data = {}
annotated_count = {}
skipped_count = {}
total_count = {}
for tier in tiers:
  tier_oncokb_data = []
  annotated_count[tier] = 0
  total_count[tier] = 0
  skipped_count[tier] = 0
  info(f'Fetching OncoKB for {tier}')
  if len(gs_data['VARIANTS'][tier].keys()) == 0:
    continue
  columns = gs_data['VARIANTS'][tier]['columns']
  for variant in gs_data['VARIANTS'][tier]['data']:
    total_count[tier] += 1
    HGVSg = get_hgvsg(variant, columns)
#    info(f'HGVSg String: {HGVSg}')

    variant_oncokb_data = {} 
    if args.include_variant:
      variant_oncokb_data['variant'] = variant

    # Pre-flight check. Ignore Filtered variants below args.min_filtered_vaf
    if tier == 'Filtered':
      vaf = float(variant[columns.index('vaf')].replace('%',''))
      if vaf < args.min_filtered_vaf:
        #info(f'Failed VAF filter: {vaf}')
        variant_oncokb_data['apiStatus'] = 'low_vaf'
        tier_oncokb_data.append(variant_oncokb_data)
        skipped_count[tier] += 1
        continue

    # Pre-flight check. Fetch data without a tumor type.
    preflight_oncokb_data, res = \
      get_oncokb(HGVSg, gs_config['oncokb_api_timeout'], None)
    if not preflight_oncokb_data:
      variant_oncokb_data['apiStatus'] = 'api_failed'
      variant_oncokb_data['apiRequests'] = get_api_requests(res)
      tier_oncokb_data.append(variant_oncokb_data)
      continue
    # If the description entry is empty, was not found
    elif not preflight_oncokb_data['mutationEffect']['description']:
      variant_oncokb_data['apiStatus'] = 'not_found'
      tier_oncokb_data.append(variant_oncokb_data)
      continue
    variant_oncokb_data['apiStatus'] = 'ok'

    # Go ahead and annotate
    annotated_count[tier] += 1

    for tumor_type in tumor_types:
      tumor_type_oncokb_data, res = \
        get_oncokb(HGVSg, gs_config['oncokb_api_timeout'], 
                   tumor_type)
      # If no description, lookup failed, set empty object
      # If the per-tumor lookup returns false, the API call failed
      # add the requests object for debugging
      if not tumor_type_oncokb_data:
        tumor_type_oncokb_data = {}
        tumor_type_oncokb_data['apiStatus'] = 'api_failed'
        tumor_type_oncokb_data['apiRequests'] = get_api_requests(res)
      # If the description entry is empty, was not found
      elif not tumor_type_oncokb_data['mutationEffect']['description']:
        tumor_type_oncokb_data['apiStatus'] = 'not_found'
      else:
        tumor_type_oncokb_data['apiStatus'] = 'ok'
      variant_oncokb_data[tumor_type] = tumor_type_oncokb_data
    tier_oncokb_data.append(variant_oncokb_data)

  if 'REPORTING' not in gs_data:
    gs_data['REPORTING'] = {}
  if 'oncokb' not in gs_data['REPORTING']:
    gs_data['REPORTING']['oncokb'] = {}
  gs_data['REPORTING']['oncokb'][tier] = tier_oncokb_data

json_string = json.dumps(
  gs_data,
  indent=2
)

print(json_string)

total_annotated = 0
for tier in tiers:
  info(f'Total {tier}: {total_count[tier]}')
  info(f'Skipped {tier}: {skipped_count[tier]}')
  info(f'Annotated {tier}: {annotated_count[tier]}')
  total_annotated += annotated_count[tier]
info(f'Total Annotated: {total_annotated}')

### added total elapsed time -jaq
end_time = time.time()
elapsed = end_time - start_time # in seconds

minutes, seconds = divmod(elapsed, 60)
hours, minutes = divmod(minutes, 60)
seconds += elapsed % 1 # adds decimaled seconds

elapsed_str = 'Elapsed:'
if int(hours) > 0:
  elapsed_str += f' {int(hours):01d} hour'
  if int(hours) != 1: # proper english plurals
    elapsed_str += 's'
if int(minutes) > 0:
  elapsed_str += f' {int(minutes):01d} minute'
  if int(minutes) != 1: # proper english plurals
    elapsed_str += 's'
elapsed_str += f' {seconds:.03f} secs'

info(elapsed_str)

debug('%s end', (SCRIPT_PATH))

