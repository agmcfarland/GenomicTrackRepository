# Downloading NCBI hs1 chromosome lengths 
# https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)


if (!file.exists(file.path(project_paths$data_raw, 'hs1.chrom.sizes.txt'))) {

  system(
    command = paste0(
      'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt', ' ', project_paths$data_raw
    )
  )
  
  }

df_chrom_lengths <- read.csv(file.path(project_paths$data_raw, 'hs1.chrom.sizes.txt'), header = FALSE, sep = '\t')
colnames(df_chrom_lengths) <- c('seqname', 'chr_length')
write.csv(df_chrom_lengths,  file.path(project_paths$data_processed, 'hs1', 'hs1.chrom.sizes.txt'))