
rule combine_results:
    input:
        expand(rules.run_coloc.output.result, locuseq = data.locuseq)
    output:
        ws_path("combined_coloc_results.tsv")
    conda:
        "../envs/coloc.yml"
    resources:
        runtime=lambda wc, attempt: 30 + attempt * 30,
    script:
        "../scripts/combine_results.R"
