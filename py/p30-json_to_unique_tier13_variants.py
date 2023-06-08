#!/usr/bin/env python

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
  'psyntax', 'csyntax'#,
  #'start', 'end', 'maf_ref', 'maf_alt' # added this row to test out
]
data = {}
for field in fields:
  data[field] = []

for in_f in args.input_files:
  with open(in_f) as fh:
    j = json.loads(fh.read())
  if 'PASS' not in j['VARIANTS']:
    info(f'No PASS in {in_f}')
    continue
  if 'data' not in j['VARIANTS']['PASS']:
    info(f'No data in {in_f}')
    continue
  tier13_columns = j['VARIANTS']['PASS']['columns']
  tier13_data = j['VARIANTS']['PASS']['data']
  for entry in tier13_data:
    for field in fields:
      data[field].append(entry[tier13_columns.index(field)])

df = pd.DataFrame.from_dict(data)
info(f'Shape before dedup: {df.shape}')
df = df.drop_duplicates()
info(f'Shape after dedup: {df.shape}')
df.to_csv(sys.stdout, sep=sep, index=None)

debug('%s end', (SCRIPT_PATH))


