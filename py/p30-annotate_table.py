#!/usr/bin/env python

import os,sys
import argparse
import logging
import pandas as pd
import requests
import json

SCRIPT_PATH = os.path.abspath(__file__)
FORMAT = '[%(asctime)s] %(levelname)s %(message)s'
l = logging.getLogger()
lh = logging.StreamHandler()
lh.setFormatter(logging.Formatter(FORMAT))
l.addHandler(lh)
l.setLevel(logging.INFO)
debug = l.debug; info = l.info; warning = l.warning; error = l.error

DESCRIPTION = '''

Given a table of variants, annotate with OncoKB data

'''

EPILOG = '''
'''

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG,
  formatter_class=CustomFormatter)

parser.add_argument('config') # config file with key -jaq
parser.add_argument('variant_table') # inputted table -jaq
parser.add_argument('-g', '--genomic', action='store_true',
    help='Use genomic coordinates/alterations (MAF format)')
parser.add_argument('-v', '--verbose', action='store_true',
    help='Set logging level to DEBUG')

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)
  
# Check for presence of required keys and types in the config file
# Hard exit if anything is missing (ported from oncokb_annotate_json)
# -jaq
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

token = gs_config['oncokb_api_key']

headers = {
  'accept': 'application/json',
  'authorization': f'Bearer {token}'
}

json_out_d = 'json_protein'
if args.genomic:
  json_out_d = 'json_genomic'
os.makedirs(json_out_d)

df = pd.read_csv(args.variant_table, sep='\t')

def add_maf(row):
  typ = row['type']
  pos = row['pos']
  ref = row['ref']
  alt = row['alt']
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
    start = pos + 1
    end = pos + len(ref) - 1
    maf_ref = ref[1:]
    maf_alt = '-'
  # INSERTION
  elif (len(ref) == 1) and (len(alt) > 1):
    start = pos
    end = pos + 1
    maf_ref = '-'
    maf_alt = alt[1:]
  # COMPLEX this might handle DNP/TNP/ONP but what if 
  # len(ref) != len(alt)?
  else: 
    maf_ref = ref
    maf_alt = alt
    start = pos
    end = start + len(alt) - 1

  row['start'] = start
  row['end'] = end
  row['maf_ref'] = maf_ref
  row['maf_alt'] = maf_alt
  return row

def add_oncokb(row):
  chrom = row['chrom']
  chrom_no_chr = chrom.replace('chr', '')
  pos = row['pos']
  gene = row['gene']
  alteration = row['psyntax'].replace('p.', '')
  
  if args.genomic:
    start = row['start']
    end = row['end']
    maf_ref = row['maf_ref']
    maf_alt = row['maf_alt']
    ref = row['ref']
    alt = row['alt']
    genomicLocation = ','.join(str(x) for x in (chrom_no_chr, start, end, maf_ref, maf_alt)) 
  
    url = 'https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange'
    params = {
      'genomicLocation': genomicLocation,
      'referenceGenome': 'GRCh38'
    }
  else:
    url = f'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange'
    params = {
      'hugoSymbol': gene,
      'alteration': alteration,
      'referenceGenome': 'GRCh38'
    }
  r = requests.get(url, headers=headers, params=params)
  if args.genomic:
    info(f'{gene} {alteration} {pos} {ref} {alt} {chrom_no_chr} {start} {end} {maf_ref} {maf_alt} {r.status_code}')
  else:
    info(f'{gene} {alteration} {r.status_code}')
  oncokb_data = r.json()
  json_string = json.dumps(oncokb_data,
    sort_keys=True, indent=2, separators=(',',':'))
  clean_alteration = alteration.replace('*', 'STAR')
  clean_alteration = clean_alteration.replace('>', 'GT')
  clean_alteration = clean_alteration.replace('+', 'PLUS')
  if args.genomic:
    out_p = os.path.join(json_out_d, f'{gene}_{clean_alteration}_{chrom}_{pos}_{ref}_{alt}.json')
  else:
    out_p = os.path.join(json_out_d, f'{gene}_{clean_alteration}.json')
  if args.genomic:
    out_p = os.path.join(json_out_d, f'{gene}_{chrom}_{pos}_{ref}_{alt}.json')

  with open(out_p, 'w') as out_fh:
    out_fh.write(json_string)
  return oncokb_data['mutationEffect']['description']


if args.genomic:
  df = df.apply(add_maf, axis=1)

df['oncokb'] = df.apply(add_oncokb, axis=1)

df.to_csv(sys.stdout, sep='\t', index=None)

debug('%s end', (SCRIPT_PATH))