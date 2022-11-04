#!/bin/bash

(
echo "["
for taxid in 10090 10029 9606 7227 10116 9823; do
	for gene in $(cat glycogenes.txt| tr '[:upper:]' '[:lower:]'); do
		curl --silent "https://mygene.info/v3/query?q=symbol:${gene}&species=${taxid}" | jq -c 'first(.hits[]|select(has("entrezgene")))|{entrez:.entrezgene,name:.symbol,taxid:.taxid}';
		echo ','
	done
done
echo "]"
)
> ref_gene_names.json