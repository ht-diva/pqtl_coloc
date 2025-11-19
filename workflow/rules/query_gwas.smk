
rule gwas_yaml:
    output:
        qregion = temp(ws_path("tmp/gwas/{locuseq}_region.tsv")),
        qyaml = temp(ws_path("tmp/gwas/{locuseq}_query.yaml")),
    params:
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        locuseq = "{locuseq}",
        tail  = config.get("coloc").get("extension"),
        projB = config.get("gwasstudio").get("project_B"),
        studB = config.get("gwasstudio").get("study_B"),
        cateB = config.get("gwasstudio").get("category_B"),
    #conda:
    #    "envs/environment.yml"
    resources:
        runtime=lambda wc, attempt: 5 + attempt * 5,
    shell:
        """
        source /exchange/healthds/singularity_functions

        echo "Genomic region: {params.locus}"
       
        # take region bounaries from locus string
        chr=$(echo {params.locus} | cut -d'_' -f1)
        beg=$(echo {params.locus} | cut -d'_' -f2)
        end=$(echo {params.locus} | cut -d'_' -f3)
        
        # extend the boundaries if needed
        beg_ext=$((beg - {params.tail}))
        end_ext=$((end + {params.tail}))
        
        region=${{chr}}:${{beg_ext}}-${{end_ext}}
        echo "Extended region is: $region"

        # create region file
        echo "$chr\t$beg\t$end" > {output.qregion}

        # Write YAML file
cat <<EOF > {output.qyaml}
project: {params.projB}
study: {params.studB}
category: {params.cateB}

output:
  - build
  - notes.sex
  - notes.source_id
  - trait.desc
  - total.samples
EOF

        echo "YAML file created: {output.qyaml}"
        """


rule query_gwas:
    input:
        qregion = rules.gwas_yaml.output.qregion,
        qyaml = rules.gwas_yaml.output.qyaml,
    output:
        sentinel = temp(ws_path("tmp/gwas/{locuseq}/{locuseq}.sentinel"))
    #conda:
    #    "envs/environment.yml"
    params:
        prefix=lambda wildcards, output: output.sentinel.replace(".sentinel", ""),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        source /exchange/healthds/singularity_functions

        gwasstudio export  --search-file {input.qyaml}  --get-regions-snps {input.qregion}  --output-prefix {params.prefix}

        odir=$(dirname {output.sentinel})
        #mv slurm*.out $odir
        touch {output.sentinel}
        """
