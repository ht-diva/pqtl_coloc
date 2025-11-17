#!/usr/bin/Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

#----------#
# taking variants file as input
sets_path <- snakemake@input
file_path <- snakemake@output

#--------------#
# the path where result files are saved
exmpl_path <- as.character(sets_path[1])
pcmd <- dirname(exmpl_path)
file_path <- as.character(file_path)

#--------------#
# scan the results for all proteins sequence
res_files <- list.files(
  pattern = paste0("seq.(\\d+).(\\d+)_(\\d+)_(\\d+)_(\\d+)_coloc_results.csv"),
  path = pcmd,       # the path where results files live
  recursive = FALSE, # to show the files in subdirectories or subfolders
  full.names = TRUE # to show full path
)

#--------------#
# extract seqid from input sentinel files
input_seqid <- map_dfr(
  sets_path, function(path) {
      base_path = dirname(path)
      file_name = basename(path)
      seqid = gsub("_.*$", "", file_name) # extract the protein sequence id and remove file format
  data.frame(base_path, file_name, seqid)
  }
)

#--------------#
# extract seqid from filenames
seq_list_tbl <- tibble(res_files) %>% mutate(seqid = str_extract(res_files, "seq.\\d+.\\d+"))

# select input seqids from all present output filenames
res_files_input <- res_files[seq_list_tbl$seqid %in% input_seqid$seqid]

#--------------#
# Merge filenames characteristics
combined_results <- data.table::rbindlist(
    fill = TRUE,
    lapply(
      res_files_input, 
      function(x) {
        data.table::fread(x, data.table=F, fill = TRUE) #%>% 
        # mutate(
        #     seqid = stringr::str_split_fixed(basename(x), "_", 2)[,1],
        #     locus = stringr::str_split_fixed(basename(x), "_locus_", 2)[,2],
        #     locus = stringr::str_remove_all(locus, "_conditional_snps.tsv")
        #     )
    }
    )
  )

#combined_results <- combined_results %>% arrange(Chr, bp)

#--------------#
# save the joint results
write.table(combined_results, file = file_path, sep = "\t", quote = T, row.names = F)
