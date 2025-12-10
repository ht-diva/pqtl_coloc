#!/bin/bash

source /exchange/healthds/singularity_functions

# parameters from snakemake
OFILE=${snakemake_output[sumstat]}
PREFIX=${snakemake_params[prefix]}
BUFFER=${snakemake_params[tail]}
LOCUS=${snakemake_params[locus]}
LOCUSEQ=${snakemake_wildcards[locuseq]}
QREGION="${LOCUSEQ}_region_pwas.tsv"
QYAML="${LOCUSEQ}_query_pwas.yaml"

echo "Genomic region: $LOCUS"

# extract region bounaries from locus string
chr=$(echo $LOCUS | cut -d'_' -f1)
beg=$(echo $LOCUS | cut -d'_' -f2)
end=$(echo $LOCUS | cut -d'_' -f3)

# extend the boundaries if needed
beg_ext=$((beg - $BUFFER))
end_ext=$((end + $BUFFER))

region=${{chr}}:${{beg_ext}}-${{end_ext}}
echo "Extended region is: $region"

# create region file
printf "%s\t%s\t%s\n" "$chr" "$beg" "$end" > $QREGION

# create YAML file    
# extract "8280.238" then replace dot with dash
seqid=$(echo $LOCUSEQ | awk -F'[._]' '{{print $2"-"$3}}')

# write YAML file
cat <<EOF > $QYAML
project: ${snakemake_params[projA]}
study: ${snakemake_params[studA]}
category: ${snakemake_params[cateA]}

trait:
  - seqid: $seqid

output:
  - build
  - notes.sex
  - notes.source_id
  - trait.desc
  - total.samples
EOF

echo "YAML file created: $QYAML"

# start to query region
gwasstudio export \
  --search-file $QYAML \
  --get-regions-snps $QREGION \
  --output-prefix $PREFIX


# extract directory and the first numeric ID after "seq."
odir=$(dirname "$OFILE")
id=$(basename "$OFILE" | sed -E 's/^seq\.([0-9]+).*/\1/')

# rename output and remove meta, YAML, region, and slurm logs
mv ${snakemake_params[filename]} $OFILE
rm $QYAML $QREGION "$odir/seq.${id}_meta.csv"

echo -e "\nEnd of query."