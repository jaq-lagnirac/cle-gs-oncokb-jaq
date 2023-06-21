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

generated = pd.read_csv(args.generated_table, sep='\t')
reference = pd.read_csv(args.reference_table, sep='\t')

def compare_generated(row):
  gene = row['gene']
  id38 = row['transcript']
  
  print(reference.loc[gene])
  
generated.apply(compare_generated, axis=1)