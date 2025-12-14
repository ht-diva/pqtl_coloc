
rule query_gwas:
    output:
        sentinel = ws_path("tmp/gwas/{locuseq}/{locuseq}.sentinel"),
    log:
        ws_path("logs/gwas/{locuseq}.log")
    params:
        prefix=lambda wildcards, output: output.sentinel.replace(".sentinel", ""),
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        tail  = config.get("coloc").get("extension"),
        projB = config.get("gwasstudio").get("project_B"),
        studB = config.get("gwasstudio").get("study_B"),
        cateB = config.get("gwasstudio").get("category_B"),
        nvarB = config.get("gwasstudio").get("ntraits_B"),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    script:
        "../scripts/query_gwas.sh"
