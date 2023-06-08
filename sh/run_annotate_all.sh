#!/bin/bash

argn=2

set -o nounset
set -o errexit
set -o pipefail

E_BADARGS=65

ee () { echo -e "$@" 1>&2; }

script=$(basename "$0")

usage=$(cat <<HERE
USAGE: $script [-i] <config> <JSON directory>

  -i  Use --include-variant option


HERE
)

include_flag=''
while getopts ":ih" opt; do
  case $opt in
    i)
      ee 'Adding --include-variant flag'
      include_flag='--include-variant'
      ;;
    h|\?)
      # For unknown options, OPTARG is set to the option
      if [[ -n ${OPTARG=''} ]]; then
        ee "Invalid option: -$OPTARG" >&2
      fi
      ee "$usage"
      exit $E_BADARGS
  esac
done

shift $((OPTIND-1))

if [ $# -ne $argn ]; then
  ee "$usage"
  exit $E_BADARGS
fi

config=$1
json_d=$2

ee $config
ee $json_d

rm -rf annotated
mkdir annotated

for f in $json_d/*.report.json; do 
  echo $f
  /usr/bin/time -v py/oncokb_annotate_json.py $include_flag $config $f > annotated/$(basename $f)
done 2>annotate_all.err
