#!/usr/bin/env python

# based off of p30-json_to_unique_tier13_variants.py

import os, sys; sys
import argparse
import logging
import pandas as pd
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

Takes input JSON file(s) and converts them into a tsv table to be processed later

Parses both PASS and Filtered variants

'''

EPILOG = '''
'''

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG,
  formatter_class=CustomFormatter)

parser.add_argument('input_files', nargs='+')
parser.add_argument('-v', '--verbose', action='store_true',
    help='Set logging level to DEBUG')

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)

debug('%s begin', SCRIPT_PATH)

sep = '\t'

fields = [
  'type', 'chrom', 'pos', 'ref', 'alt', 'gene', 'transcript', 
  'psyntax', 'csyntax'
]
data = {}
for field in fields:
  data[field] = []
data['filter_type'] = [] # -jaq

filter_types = ['PASS', 'Filtered']

for filter_type in filter_types:
  for in_f in args.input_files:
    with open(in_f) as fh:
      j = json.loads(fh.read())
    if filter_type not in j['VARIANTS']:
      info(f'No {filter_type} in {in_f}')
      continue
    if 'data' not in j['VARIANTS'][filter_type]:
      info(f'No data in {in_f}')
      continue
    tier13_columns = j['VARIANTS'][filter_type]['columns']
    tier13_data = j['VARIANTS'][filter_type]['data']
    for entry in tier13_data:
      for field in fields:
        data[field].append(entry[tier13_columns.index(field)])
      data['filter_type'].append(filter_type) # -jaq

df = pd.DataFrame.from_dict(data)
info(f'Shape before dedup: {df.shape}')
df = df.drop_duplicates()
info(f'Shape after dedup: {df.shape}')
df.to_csv(sys.stdout, sep=sep, index=None)

debug('%s end', (SCRIPT_PATH))


