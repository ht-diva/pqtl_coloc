#!/bin/bash

set -euo pipefail
# IMPORTANT: explicit redirection
exec >"${snakemake_log[0]}" 2>&1

source /exchange/healthds/singularity_functions

# parameters from snakemake
OFILE=${snakemake_output[sentinel]}
PREFIX=${snakemake_params[prefix]}
BUFFER="${snakemake_params[tail]}"
EXPECTED_N=${snakemake_params[nvarB]}
LOCUS=${snakemake_params[locus]}
LOCUSEQ=${snakemake_wildcards[locuseq]}
QREGION="${LOCUSEQ}_region_gwas.tsv"
QYAML="${LOCUSEQ}_query_gwas.yaml"

echo "--- GWAS QUERY START ---"
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
    beg_ext=1
fi

echo "Original region:  ${chr}:${beg}-${end}"
echo "Buffered region:  ${chr}:${beg_ext}-${end_ext}"
echo "Buffer size:      ${BUFFER} bp"

# create region file
printf "%s\t%s\t%s\n" "$chr" "$beg_ext" "$end_ext" > $QREGION
echo "Region file created: $QREGION"

# write YAML file
cat <<EOF > $QYAML
category: ${snakemake_params[cateB]}
project: ${snakemake_params[projB]}
study: ${snakemake_params[studB]}

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
  --output-prefix $PREFIX \
  || echo "WARNING: gwasstudio exited non-zero (ignored)"

echo
echo "--- VALIDATING OUTPUT FILES ---"

shopt -s nullglob
files=( $PREFIX*.csv.gz )

echo "Found ${#files[@]} CSV.GZ files (expected $EXPECTED_N)"

# HARD FAIL if count is wrong
if [[ ${#files[@]} -lt $EXPECTED_N ]]; then
    echo "ERROR: Expected $EXPECTED_N GWAS outputs, found ${#files[@]}"
    exit 1
fi

valid_files=()

for f in "${files[@]}"; do
    echo
    echo "Checking: $f"

    # must be non-empty
    if [[ ! -s "$f" ]]; then
        echo "  -> empty, skipping"
        continue
    fi

    # must be readable gzip
    if ! gzip -t "$f" 2>/dev/null; then
        echo "  -> gzip corrupted, skipping"
        continue
    fi

    # must have required header columns; SAFE read (no SIGPIPE)
    header=$(gzip -cd "$f" | sed -n '1p')

    if ! echo "$header" | grep -q "SNPID,CHR,POS,BETA,SE"; then
        echo "  -> header invalid"
        echo "     header was: $header"
        continue
    fi

    echo "  -> VALID"
    valid_files+=("$f")
done

echo
echo "Valid files: ${#valid_files[@]}"

if [[ ${#valid_files[@]} -lt $EXPECTED_N ]]; then
    echo "ERROR: Only ${#valid_files[@]} valid GWAS files, expected $EXPECTED_N"
    exit 1
fi


# extract directory and "seq.16300" from "seq.16300.4"
odir=$(dirname "$OFILE")
id=$(basename "$OFILE" | sed -E 's/^seq\.([0-9]+).*/\1/')

# create sentinel file as final output
touch $OFILE

# remove YAML, region, and GWAS metadata files
rm $QYAML $QREGION "$odir/seq.${id}_meta.csv"
#mv slurm*.out $odir

echo "--- GWAS QUERY END ---"
date