#!/bin/bash

wanted_version=$(cat version.json | jq --raw-output '.revision[0]')

output_dir="/Volumes/trs582/msdatas-out-${wanted_version}"

mkdir -p "$output_dir"

declare -a msf_folders_to_convert=(
	"$HOME/dev/ccg-recipes"
	'/Volumes/trs582/From_Sund/msfs'
#	'/Volumes/trs582/brain_msfs/msfs'
)

for folder in "${msf_folders_to_convert[@]}"; do
	./bin/run_conversion.sh --input "$folder" --output "$output_dir" --nonetwork
done