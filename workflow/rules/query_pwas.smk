
rule pwas_yaml:
    output:
        qregion = temp(ws_path("tmp/pwas/{locuseq}_region.tsv")),
        qyaml = temp(ws_path("tmp/pwas/{locuseq}_query.yaml")),
    params:
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        locuseq = "{locuseq}",
        tail = config.get("coloc").get("extension"),
        projA = config.get("gwasstudio").get("project_A"),
        studA = config.get("gwasstudio").get("study_A"),
        cateA = config.get("gwasstudio").get("category_A"),
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

        # create yaml file
        
        # Extract "8280.238" then replace dot with dash
        seqid=$(echo {params.locuseq} | awk -F'[._]' '{{print $2"-"$3}}')
        
        # Write YAML file
cat <<EOF > {output.qyaml}
project: {params.projA}
study: {params.studA}
category: {params.cateA}

trait:
  - seqid: $seqid

output:
  - build
  - notes.sex
  - notes.source_id
  - trait.desc
  - total.samples
EOF

        echo "YAML file created: {output.qyaml}"
        """


rule query_pwas:
    input:
        qregion = rules.pwas_yaml.output.qregion,
        qyaml = rules.pwas_yaml.output.qyaml,
    output:
        sumstat = temp(ws_path("tmp/pwas/{locuseq}_sumstat.csv.gz"))
    #conda:
    #    "envs/environment.yml"
    params:
        filename = lambda wildcards: get_pwasname(wildcards),
        prefix  = lambda wildcards, output: output.sumstat.replace("_sumstat.csv.gz", ""),
        locuseq = "{locuseq}",
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    shell:
        """
        source /exchange/healthds/singularity_functions

        gwasstudio export --search-file {input.qyaml} --get-regions-snps {input.qregion} --output-prefix {params.prefix}
        
        # rename output and remove slurm logs
        mv {params.filename} {output.sumstat}
        #rm slurm*.out
        """
