
from pathlib import Path
import pandas as pd


# Define input for the rules
# read loci list
lb = pd.read_csv(config["path_lb"])

# Create a new column by concatenating 
lb["locus"]  = lb["chr"].astype(str) + "_" + lb["start"].astype(str) + "_" + lb["end"].astype(str)
lb["locuseq"] = lb["seqid"].astype(str) + "_" + lb["locus"].astype(str)

data = (
    pd.DataFrame(lb, columns=["locuseq", "seqid", "chr", "locus", "SNPID"])
    .set_index("locuseq", drop=False)
    .sort_index()
)


def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))

# return locus of locuseq
def get_locus(wildcards):
    return str(data.loc[wildcards, "locus"])

# return GWAS summary results 
def get_gwas(wildcards):
    seqid = data.loc[wildcards, "seqid"]
    file_path = f"{seqid}/{seqid}.gwaslab.tsv.gz"
    return str(Path(config.get("path_gwas"), file_path))

# return filename of tileDB output
def get_pwasname(wildcards):

    odir = config["workspace_path"]
    proj = config["gwasstudio"]["project_A"]
    study = config["gwasstudio"]["study_A"]
    locuseq = wildcards.locuseq     # e.g. seq.16300.4_22_43928847_43998522

    # extract seqid and region parts
    seqid = locuseq.split("_")[0]               # seq.16300.4
    region = "_".join(locuseq.split("_")[1:])   # 22_43928847_43998522

    filename = f"{locuseq}_{proj}_{study}_{seqid}.csv.gz"

    return str(Path(odir, "tmp/pwas", filename))

#results/test/tmp/pwas/seq.16300.4_22_43928847_43998522_hdsc_believe_seq.16300.4.csv.gz
