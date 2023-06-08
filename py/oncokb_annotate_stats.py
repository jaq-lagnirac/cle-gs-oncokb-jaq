#!/usr/bin/env python

import os,sys
import re
import argparse
import logging

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

parser.add_argument('err_file')
parser.add_argument('-v', '--verbose', action='store_true',
    help='Set logging level to DEBUG')

args = parser.parse_args()

if args.verbose:
  l.setLevel(logging.DEBUG)

debug('%s begin', SCRIPT_PATH)

all_total_annotated = []
all_elapsed_time = []
with open(args.err_file) as fh:
  for line in fh:
    line = line.strip()
    if 'Total Annotated' in line:
      # INFO Total Annotated: 5
      m = re.search(r'Total Annotated: (\d+)', line)
      if (m):
        total_annotated = m.group(1)
      else:
        error('Could not parse total annotated')
        sys.exit(1)
      all_total_annotated.append(total_annotated)
    elif 'Elapsed' in line:
      reference_point = 'Elapsed: '
      elapsed_index = line.index(reference_point) + len(reference_point)
      m = line[elapsed_index:]
      if (m):
        elapsed_time = m
      else:
        error('Could not parse elapsed time')
        sys.exit(1)
      all_elapsed_time.append(elapsed_time)

for i, total_annotated in enumerate(all_total_annotated):
  print('\t'.join((total_annotated, all_elapsed_time[i])))

debug('%s end', (SCRIPT_PATH))
