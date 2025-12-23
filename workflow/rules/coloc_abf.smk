
rule run_coloc:
    input:
        pwas = rules.query_pwas.output.sumstat,
        gwas = rules.query_gwas.output.sentinel
    output:
        result = ws_path("tmp/coloc/{locuseq}_coloc_results.csv")
    conda:
        "../envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 150 + attempt * 60,
    script:
        "../scripts/coloc_abf.R"
