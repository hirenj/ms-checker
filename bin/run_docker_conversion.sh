#!/bin/bash

wanted_version=$(cat version.json | jq --raw-output '.revision[0]')

input_dir=${1:-/work}
output_parent=${2:-/dist}

output_dir="${output_parent}/msdatas-${wanted_version}"

mkdir -p "$output_dir"

# declare -a msf_folders_to_convert=(
# 	"$HOME/dev/ccg-recipes"
# 	'/Volumes/trs582/From_Sund/msfs'
# #	'/Volumes/trs582/brain_msfs/msfs'
# )

curl --output uniprot.deleted.txt 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/delac_sp.txt'
curl --output uniprot.merged.txt 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt'
curl -O 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt'



for folder in $input_dir/*/; do
	./bin/run_conversion.sh --input "${folder%*/}" --output "$output_dir" --nonetwork
done