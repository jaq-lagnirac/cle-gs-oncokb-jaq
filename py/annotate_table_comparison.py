#!/usr/bin/env python

# based off of p30-annotate_table.py

import os,sys
import argparse
import logging
import pandas as pd
import requests
import json

# - jaq
from enum import Enum
import time
import statistics


# global constants
NUM_DECIMALS = 3
VAR_PRECISION = NUM_DECIMALS + 2


# global vars
call_type = 'temp'
elapsed_array = []
successes_x = 0
total_n = 0


# enumerated constants for constants
class Call(Enum): 
  NO_CALL, \
  GENOMIC, \
  PROTEIN, \
  HGVSG \
    = range(4)
  # if adding enum, update range to reflect change

 
call_type = Call.NO_CALL

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

Outputs tsv file to standard output

Outputs statistics to standard error
  -Elapsed time
  -Proportion of successful hits

'''

EPILOG = '''

REQUIRES CONFIG FILE WITH VALID ONCOKB KEY

'''

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG,
  formatter_class=CustomFormatter)

parser.add_argument('config') # config file with key -jaq
parser.add_argument('variant_table') # inputted table -jaq
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
    'oncokb_api_key': str
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

### create JSON directories
json_out_protein = 'json_protein'
os.makedirs(json_out_protein)
json_out_genomic = 'json_genomic'
os.makedirs(json_out_genomic)
json_out_hgvsg = 'json_hgvsg'
os.makedirs(json_out_hgvsg)

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
  ### COMPLEX this might handle DNP/TNP/ONP but what if 
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
  

def add_hgvsg(row): # -jaq
  variant_type = row['type']
  chromosome = row['chrom'].replace('chr', '')
  position = int(row['pos'])
  reference = row['ref']
  alteration = row['alt']
  
  hgvsg = f'{chromosome}:g.'
  
  ### SNV "Substitution"
  if (len(reference) == 1) and (len(alteration) == 1):
    if variant_type != 'SNV':
      error('Type mismatch')
      sys.exit(1)
    hgvsg += f'{position}{reference}>{alteration}'
      
  ### DELETION
  elif (len(reference) > 1) and (len(alteration) == 1):
    if variant_type != 'INDEL':
      error('Type mismatch')
      sys.exit(1)
    position += 1
    if len(reference) == 2: # one nucleotide deletion
      hgvsg += f'{position}del'
    else: # serveral nucleotide deletion
      second_position = position + len(reference) - 2
      hgvsg += f'{position}_{second_position}del'

  ### INSERTION
  elif (len(reference) == 1) and (len(alteration) > 1):
    if variant_type != 'INDEL':
      error('Type mismatch')
      sys.exit(1)
    second_position = position + 1
    hgvsg += f'{position}_{second_position}ins{alteration[1:]}'
    
  ### DELINS (not robust, simply to handle what we're dealing with)
  else: # variant_type == 'SNV': # The HARD one
    if len(reference) == 1:
      hgvsg += f'{position}delins{alteration}'
    else: # len(ref) > 1
      second_position = position + len(reference) - 1
      hgvsg += f'{position}_{second_position}delins{alteration}'
    
  row['hgvsg'] = hgvsg
  return row


def add_oncokb(row):
  ### ties in global vars
  global call_type # not necessary, but clarifies role
  global elapsed_array
  
  ### starts elapsed time
  start_time = time.time()

  ### sets up API vars
  chrom = row['chrom']
  chrom_no_chr = chrom.replace('chr', '')
  pos = row['pos']
  ref = row['ref']
  alt = row['alt']
  gene = row['gene']
  alteration = row['psyntax'].replace('p.', '')
  
  ### generates API calls
  if call_type == Call.GENOMIC:
    start = row['start']
    end = row['end']
    maf_ref = row['maf_ref']
    maf_alt = row['maf_alt']
    genomicLocation = ','.join(str(x) for x in (chrom_no_chr, start, end, maf_ref, maf_alt)) 
  
    url = 'https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange'
    params = {
      'genomicLocation': genomicLocation,
      'referenceGenome': 'GRCh38'
    }
  elif call_type == Call.PROTEIN:
    url = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange'
    params = {
      'hugoSymbol': gene,
      'alteration': alteration,
      'referenceGenome': 'GRCh38'
    }
  elif call_type == Call.HGVSG:
    hgvsg = row['hgvsg']
    
    url = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg'
    params = {
      'hgvsg': hgvsg,
      'referenceGenome': 'GRCh38'
    }
  else: # NO_CALL, safely exit
    info('API mismatch')
    sys.exit(1)
  
  ### requests URL, checks status
  r = requests.get(url, headers=headers, params=params)
  if call_type == Call.GENOMIC:
    info(f'{gene} {alteration} {pos} {ref} {alt} {chrom_no_chr} {start} {end} {maf_ref} {maf_alt} {r.status_code}')
  elif call_type == Call.PROTEIN:
    info(f'{gene} {alteration} {r.status_code}')
  elif call_type == Call.HGVSG:
    info(f'{gene} {alteration} {hgvsg} {r.status_code}')
  
  ### pulls data and creates json filenames
  oncokb_data = r.json()
  json_string = json.dumps(oncokb_data,
    sort_keys=True, indent=2, separators=(',',':'))
  clean_alteration = alteration.replace('*', 'STAR')
  clean_alteration = clean_alteration.replace('>', 'GT')
  clean_alteration = clean_alteration.replace('+', 'PLUS')
  
  ### generates directory paths
  if call_type == Call.GENOMIC:
    out_p = os.path.join(json_out_genomic, f'{gene}_{clean_alteration}_{chrom}_{pos}_{ref}_{alt}.json')
    out_p = os.path.join(json_out_genomic, f'{gene}_{chrom}_{pos}_{ref}_{alt}.json')
  elif call_type == Call.PROTEIN:
    out_p = os.path.join(json_out_protein, f'{gene}_{clean_alteration}.json')
  elif call_type == Call.HGVSG:
    out_p = os.path.join(json_out_hgvsg, f'{gene}_{clean_alteration}_{hgvsg}.json')

  ### writes to file
  with open(out_p, 'w') as out_fh:
    out_fh.write(json_string)
  
  ### ends elapsed time, processes elapsed time
  end_time = time.time()
  elapsed_time = round(end_time - start_time, NUM_DECIMALS)
  elapsed_array.append(elapsed_time)
  
#  return bool(oncokb_data['mutationEffect']['description'])
#  return oncokb_data['mutationEffect']['description']
  try:
    return oncokb_data['mutationEffect']['description']
  except KeyError:
    error(f'Not found in OncoKB: {call_type} {gene} {hgvsg}')
    sys.exit(1)


def check_api_call(df, column):
  global successes_x
  global total_n
  
  new_column = []
  for call in df[column]:
    new_column.append(bool(call))
  type_name = (column.split(' '))[0]
  new_name = f'found {type_name}?'
  df[new_name] = new_column
  
  successes_x = new_column.count(True)
  total_n = len(new_column)


df = df.apply(add_maf, axis=1)
df = df.apply(add_hgvsg, axis=1)

### calls byGenomicChange
info('----------Starts byGenomicChange calls----------')
elapsed_array = [] # resets global var
call_type = Call.GENOMIC
genomic_column = 'genomic oncokb'
df[genomic_column] = df.apply(add_oncokb, axis=1)
check_api_call(df, genomic_column)

# handles elapsed time
genomic_elapsed = 'genomic elapsed (sec)'
df[genomic_elapsed] = elapsed_array
elapsed_avg_genomic = round(statistics.mean(elapsed_array), \
  NUM_DECIMALS)
var_genomic = round(statistics.variance(elapsed_array), \
  VAR_PRECISION)
sd_genomic = round(statistics.stdev(elapsed_array), \
  NUM_DECIMALS)
total_elapsed_genomic = round(sum(elapsed_array), \
  NUM_DECIMALS)

# handles proportions, uses global vars
genomic_successes = successes_x
genomic_total = total_n
genomic_prop = round(genomic_successes  / genomic_total, NUM_DECIMALS)


### calls byProteinChange
info('----------Starts byProteinChange calls----------')
elapsed_array = [] # resets global var
call_type = Call.PROTEIN
protein_column = 'protein oncokb'
df[protein_column] = df.apply(add_oncokb, axis=1)
check_api_call(df, protein_column)

# handles elapsed time
protein_elapsed = 'protein elapsed (sec)'
df[protein_elapsed] = elapsed_array
elapsed_avg_protein = round(statistics.mean(elapsed_array), \
  NUM_DECIMALS)
var_protein = round(statistics.variance(elapsed_array), \
  VAR_PRECISION)
sd_protein = round(statistics.stdev(elapsed_array), \
  NUM_DECIMALS)
total_elapsed_protein = round(sum(elapsed_array), \
  NUM_DECIMALS)
  
# handles proportions, uses global vars
protein_successes = successes_x
protein_total = total_n
protein_prop = round(protein_successes / protein_total, \
  NUM_DECIMALS)


### calls byHGVSg
info('----------Starts byHGVSg calls----------')
elapsed_array = [] # resets global var  
call_type = Call.HGVSG
hgvsg_column = 'hgvsg oncokb'
df[hgvsg_column] = df.apply(add_oncokb, axis=1)
check_api_call(df, hgvsg_column)

# handles elapsed time
hgvsg_elapsed = 'hgvsg elapsed (sec)'
df[hgvsg_elapsed] = elapsed_array
elapsed_avg_hgvsg = round(statistics.mean(elapsed_array), \
  NUM_DECIMALS)
var_hgvsg = round(statistics.variance(elapsed_array), \
  VAR_PRECISION)
sd_hgvsg = round(statistics.stdev(elapsed_array), \
  NUM_DECIMALS)
total_elapsed_hgvsg = round(sum(elapsed_array), \
  NUM_DECIMALS)

# handles proportions, uses global vars
hgvsg_successes = successes_x
hgvsg_total = total_n
hgvsg_prop = round(hgvsg_successes  / hgvsg_total, NUM_DECIMALS)


### prints overview info to .err file
info('----------OVERVIEW OF PERFORMANCE----------')

# prints elapsed time info
info('---ELAPSED TIME---')

info(f'byGenomicChange:')
info(f'\t-Total: {total_elapsed_genomic} secs')
info(f'\t-Mean/Average: {elapsed_avg_genomic} secs')
info(f'\t-Variance: {var_genomic} secs²')
info(f'\t-Standard Deviation: {sd_genomic} secs')

info(f'byProteinChange:')
info(f'\t-Total: {total_elapsed_protein} secs')
info(f'\t-Mean/Average: {elapsed_avg_protein} secs')
info(f'\t-Variance: {var_protein} secs²')
info(f'\t-Standard Deviation: {sd_protein} secs')

info(f'byHGVSg:')
info(f'\t-Total: {total_elapsed_hgvsg} secs')
info(f'\t-Mean/Average: {elapsed_avg_hgvsg} secs')
info(f'\t-Variance: {var_hgvsg} secs²')
info(f'\t-Standard Deviation: {sd_hgvsg} secs')

elapsed_total = round((total_elapsed_genomic \
  + total_elapsed_protein \
  + total_elapsed_hgvsg), \
  NUM_DECIMALS)
seconds = elapsed_total
minutes = 0
# using while loop rather than divmod to prevent rounding discrepancies
while seconds > 60:
  seconds -= 60
  minutes += 1
seconds = round(seconds, NUM_DECIMALS)
info(f'Approx. Total Elapsed Time: {minutes} min, {seconds} secs ({elapsed_total} secs)')

# prints successful hits
info('---SUCCESSFUL HITS---')
info(f'Genomic Proportion of Successful Hits: {genomic_prop} ({genomic_successes}/{genomic_total})')
info(f'Protein Proportion of Successful Hits: {protein_prop} ({protein_successes}/{protein_total})')
info(f'HGVSg Proportion of Successful Hits: {hgvsg_prop} ({hgvsg_successes}/{hgvsg_total})')

df.to_csv(sys.stdout, sep='\t', index=None)

debug('%s end', (SCRIPT_PATH))