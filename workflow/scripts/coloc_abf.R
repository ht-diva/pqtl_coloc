
suppressPackageStartupMessages({
  library(R.utils) # to import compressed inputs
  library(glue)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(Rmpfr)
  library(coloc)
  library(purrr)
})


pwas_file <- snakemake@input[["pwas"]]
gwas_file <- snakemake@input[["gwas"]]
resu_file <- snakemake@output[["result"]]

# Handle missings when computing p-value from beta&sd
safe_pnorm <- function(b, se, p=FALSE) {
  
  # Ensure the vectors are of the same length
  if(length(b) != length(se)) {
    stop("Beta and SE must be of the same length")
  }
  
  k  <- length(b)
  b  <- as.numeric(b)
  se <- as.numeric(se)
  
  # Initialize result vector with NA values
  result <- rep(NA, k)
  
  # Identify non-missing and non-zero indices
  i <- which(!is.na(b) & !is.na(se) & se != 0)
  
  # compute z-score and take absolute, raise digits with mpfr, apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * pnorm(z_mpfr)
  mlog10p <- - log10(p_mpfr)
  
  # print p-value in character format and mlog10p in numeric
  if(p==TRUE){
    # reformat to mpfr character, then to numeric (don't set digits for MLOG10P)
    mlog10p_mpfr <- Rmpfr::formatMpfr(p_mpfr, scientific = TRUE, digits = 6)
    result[i] <- mlog10p_mpfr
  } else {
    mlog10p_mpfr <- Rmpfr::formatMpfr(mlog10p, scientific = TRUE)
    result[i] <- as.numeric(mlog10p_mpfr)
  }
  
  return(result)
}


prepare4coloc <- function(data){
  
  temp  <- data |>
    dplyr::rename(position = POS, beta = BETA) |>
    dplyr::distinct(position, .keep_all = TRUE) |> # remove duplicate SNPs
    dplyr::rename_with(~gsub("meta_total_samples", "N", .x)) |>
    dplyr::mutate(
      snp = paste0(CHR, ":", position),
      varbeta = SE^2,
      pvalues = safe_pnorm(beta, SE, p = TRUE),
      MAF = ifelse(EAF < 0.5, EAF, 1- EAF),
      sdY = coloc:::sdY.est(varbeta, MAF, N)
      ) |>
    dplyr::select(position, snp, beta, varbeta, MAF, pvalues, sdY)
  
  temp$type <- "quant"
  odata <- as.list(na.omit(temp))
  odata$type <- unique(odata$type)
  odata$sdY <- unique(odata$sdY)
  
  return(odata)
}

# files with credible sets
sums_lists <- list.files(
  path = dirname(gwas_file),
  pattern = ".csv.gz",
  full.names = TRUE
)


#-------------------------------#
# -----     Run Coloc      -----
#-------------------------------#

run_coloc <- function(pfile, gfile){
  
  pwas <- fread(pfile)
  gwas <- fread(gfile)

  annot_pwas <- prepare4coloc(pwas)
  annot_gwas <- prepare4coloc(gwas)

  # run coloc standard
  res <- coloc::coloc.abf(annot_pwas, annot_gwas)

  res_h4 <- res$summary %>% t() %>% as.data.frame()

  locuseq <- pfile %>% basename() %>% stringr::str_remove("_sumstat.csv.gz")
  seqid   <- unique(pwas$meta_notes_source_id)
  protein <- unique(pwas$meta_trait_desc)
  pheno <- unique(gwas$meta_trait_desc)

  res_final <- data.frame(
    "seqid" = seqid,
    "locus" = locuseq,
    "protein" = protein,
    "phenotype" = pheno
  ) %>%
    cbind(res_h4)

  return(res_final)
}

# all possible protein-phenotype pairs
traits2test <- expand.grid(pwas_file, sums_lists, stringsAsFactors = FALSE)

# run coloc
res <- map2_df(traits2test$Var1, traits2test$Var2, run_coloc)

# save results
write.csv(res, file = resu_file, quote = T, row.names = F)
