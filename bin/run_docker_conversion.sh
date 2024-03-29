#!/bin/bash

#wanted_version=$(cat version.json | jq --raw-output '.revision[0]')

input_dir=${1:-/input}
output_parent=${2:-/dist}

output_dir="${output_parent}"

mkdir -p "$output_dir"

# declare -a msf_folders_to_convert=(
# 	"$HOME/dev/ccg-recipes"
# 	'/Volumes/trs582/From_Sund/msfs'
# #	'/Volumes/trs582/brain_msfs/msfs'
# )


if [[ ! -f "uniprot.deleted.txt" ]]; then
	curl --output uniprot.deleted.txt 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/delac_sp.txt'
fi
if [[ ! -f "uniprot.merged.txt" ]]; then
	curl --output uniprot.merged.txt 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/sec_ac.txt'
fi
if [[ ! -f "cellosaurus.txt" ]]; then
	curl -O 'ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt'
fi


for folder in $input_dir/*/; do
	./bin/run_conversion.sh --input "${folder%*/}" --output "$output_dir" --nonetwork
done