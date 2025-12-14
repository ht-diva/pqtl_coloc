
rule query_pwas:
    output:
        sumstat = ws_path("tmp/pwas/{locuseq}_sumstat.csv.gz"),
    log:
        ws_path("logs/pwas/{locuseq}.log")
    params:
        filename = lambda wildcards: get_pwasname(wildcards),
        prefix = lambda wildcards, output: output.sumstat.replace("_sumstat.csv.gz", ""),
        locus = lambda wildcards: get_locus(wildcards.locuseq),
        tail  = config.get("coloc").get("extension"),
        projA = config.get("gwasstudio").get("project_A"),
        studA = config.get("gwasstudio").get("study_A"),
        cateA = config.get("gwasstudio").get("category_A"),
    resources:
        runtime=lambda wc, attempt: 120 + attempt * 60,
    script:
        "../scripts/query_pwas.sh"
