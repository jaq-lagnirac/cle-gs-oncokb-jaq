#!/usr/bin/env python

import os
import sys
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

Compares gene and transcript ID with a reference table of the genes and transcripts

'''

EPILOG = '''
'''

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
  pass
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG,
  formatter_class=CustomFormatter)

parser.add_argument('generated_table')
parser.add_argument('reference_table')
parser.add_argument('-v', '--verbose', action='store_true',
    help='Set logging level to DEBUG')

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)

debug('%s begin', SCRIPT_PATH)

sep = '\t'

generated = pd.read_csv(args.generated_table, sep=sep)
reference = pd.read_csv(args.reference_table, sep=sep)


def compare_generated(row):
  gene = row['gene']
  id38 = row['transcript']

  data_comparison = reference.loc[(reference['Hugo Symbol'] == gene) \
    & (reference['GRCh38 Isoform'] == id38)] # false if present, true if not
  
  correct_id = not data_comparison.empty # true if present, false if not
  row['Correct ID?'] = correct_id # adds result to row
  
  return row

  
generated.apply(compare_generated, axis=1)
  
gene_test = reference.loc[(reference['Hugo Symbol'] == 'ABL1') & (reference['GRCh38 Isoform'] == 'ENST00000318560')]
print(gene_test)# false if present, true if not

df.to_csv(sys.stdout, sep='\t', index=None)

debug('%s end', (SCRIPT_PATH))