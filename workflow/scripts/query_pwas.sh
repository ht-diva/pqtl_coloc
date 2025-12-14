#!/bin/bash

set -euo pipefail
# IMPORTANT: explicit redirection
exec >"${snakemake_log[0]}" 2>&1

source /exchange/healthds/singularity_functions

# parameters from snakemake
OFILE=${snakemake_output[sumstat]}
TEMPFILE=${snakemake_params[filename]}
PREFIX=${snakemake_params[prefix]}
BUFFER=${snakemake_params[tail]}
LOCUS=${snakemake_params[locus]}
LOCUSEQ=${snakemake_wildcards[locuseq]}
QREGION="${LOCUSEQ}_region_pwas.tsv"
QYAML="${LOCUSEQ}_query_pwas.yaml"

echo "--- PWAS QUERY START ---"
date

# extract region bounaries from locus string
chr=$(echo $LOCUS | cut -d'_' -f1)
beg=$(echo $LOCUS | cut -d'_' -f2)
end=$(echo $LOCUS | cut -d'_' -f3)

# extend the boundaries if needed
beg_ext=$((beg - BUFFER))
end_ext=$((end + BUFFER))

# make sure beg_ext is not negative
if [ "$beg_ext" -lt 0 ]; then
    beg_ext=0
fi

echo "Original region:  ${chr}:${beg}-${end}"
echo "Buffered region:  ${chr}:${beg_ext}-${end_ext}"
echo "Buffer size:      ${BUFFER} bp"

# create region file
printf "%s\t%s\t%s\n" "$chr" "$beg_ext" "$end_ext" > $QREGION
echo "Region file created: $QREGION"

# for YAML file, extract "8280.238" then replace dot with dash
seqid=$(echo $LOCUSEQ | awk -F'[._]' '{{print $2"-"$3}}')

# write YAML file
cat <<EOF > $QYAML
category: ${snakemake_params[cateA]}
project: ${snakemake_params[projA]}
study: ${snakemake_params[studA]}

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
echo

# start to query region
gwasstudio export \
  --search-file $QYAML \
  --get-regions-snps $QREGION \
  --output-prefix $PREFIX


# extract directory and the first numeric ID after "seq."
odir=$(dirname "$OFILE")
id=$(basename "$OFILE" | sed -E 's/^seq\.([0-9]+).*/\1/')

# atomic move to final output
mv $TEMPFILE $OFILE

# remove metadata, YAML, and region files
rm $QYAML $QREGION "$odir/seq.${id}_meta.csv"

echo "--- PWAS QUERY END ---"
date