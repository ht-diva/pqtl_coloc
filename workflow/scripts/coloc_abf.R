
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

#----------------------------#
# -----      Inputs     -----
#----------------------------#

# Pass parameters defined at Snakemake's rule
pwas_file <- snakemake@input[["pwas"]]
gwas_file <- snakemake@input[["gwas"]]
resu_file <- snakemake@output[["result"]]
bin_pwas  <- snakemake@params[["bin_pwas"]]
bin_gwas  <- snakemake@params[["bin_gwas"]]


#----------------------------#
# -----    Functions    -----
#----------------------------#

# Handle NAs when computing p-value from beta&sd
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
  
  # compute z-score and take absolute, 
  # raise digits with mpfr, 
  # apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * Rmpfr::pnorm(z_mpfr) # dispatch the correct pnorm()
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

#---------------#

# Annotated GWAS sumstat in a format desired for coloc:
#  - for quantitative traits, compute sdY & type = 'quant'
#  - for binary traits, compute case proportion (s) & type = 'cc'
prepare4coloc <- function(data, dichotomous = FALSE){
  
  temp  <- data |>
    dplyr::mutate(
      varbeta = SE^2,
      pvalues = safe_pnorm(BETA, SE, p = TRUE),
      MAF = ifelse(EAF < 0.5, EAF, 1 - EAF)
    ) |>
    dplyr::rename_with(
      ~gsub("meta_total_samples", "N", .x)
      ) |>
    dplyr::rename(
      snp = SNPID, # as coloc requires these names
      beta = BETA
      )
  
  if(dichotomous){
    
    temp <- temp |>
      dplyr::mutate(s = meta_total_cases/N) %>%
      dplyr::select(snp, beta, varbeta, MAF, pvalues, s)
  
  } else {
    
    temp <- temp %>%
      dplyr::mutate(sdY = coloc:::sdY.est(varbeta, MAF, N)) %>%
      dplyr::select(snp, beta, varbeta, MAF, pvalues, sdY)
    }
  
  
  odata <- as.list(na.omit(temp))
  
  if(dichotomous) {
    
    odata$type <- "cc"
    odata$s <- unique(odata$s)
    
    } else {
      
      odata$type <- "quant"
      odata$sdY <- unique(odata$sdY)
  }
  
  return(odata)
}

#---------------#

# Annotate GWAS first, then run coloc ABF
run_coloc <- function(gfile){
  
  gwas <- fread(gfile)
  
  annot_gwas <- prepare4coloc(gwas, dichotomous = bin_gwas)
  
  # run coloc standard
  res <- coloc::coloc.abf(annot_pwas, annot_gwas)
  
  res_h4 <- res$summary %>% t() %>% as.data.frame()
  
  # Define labels helping to combine coloc results later
  locus <- pwas_file %>% basename() %>%
    stringr::str_remove("_sumstat.csv.gz") %>%
    stringr::str_remove("seq.(\\d)+.(\\d)+_")
  
  seqid   <- unique(pwas$meta_notes_source_id)
  protein <- unique(pwas$meta_trait_desc)
  pheno <- unique(gwas$meta_trait_desc)
  
  res_final <- data.frame(
    "seqid" = seqid,
    "locus" = locus,
    "protein" = protein,
    "phenotype" = pheno
  ) %>%
    cbind(res_h4)
  
  return(res_final)
}


#-------------------------------#
# -----     Run Coloc      -----
#-------------------------------#

# List GWASs sumstat extracted by tileDB
sums_lists <- list.files(
  path = dirname(gwas_file),
  pattern = ".csv.gz",
  full.names = TRUE
)

# Read protein GWAS (PWAS)
pwas <- fread(pwas_file)

# Annotate PWAS
annot_pwas <- prepare4coloc(pwas, dichotomous = bin_pwas)

# run coloc
res <- map_df(sums_lists, run_coloc)

# save results
write.csv(res, file = resu_file, quote = T, row.names = F)
