# Script to generate cancer gene track
# From COSMIC sanger institute. behind paywall/sign up wall

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

df_cancer <- read.csv(file.path(project_paths$data_external, 'Cosmic_CancerGeneCensus_v101_GRCh38.tsv'), sep = '\t')

df_cancer <- df_cancer %>%
  dplyr::select(GENE_SYMBOL, NAME, CHROMOSOME, GENOME_START, GENOME_STOP) %>%
  dplyr::rename(
    gene_symbol = GENE_SYMBOL,
    description = NAME,
    seqname = CHROMOSOME,
    start = GENOME_START,
    stop = GENOME_STOP,
  ) %>%
  na.omit() %>%
  dplyr::mutate(
    strand = '*',
    seqname = paste0('chr', as.character(seqname))
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(df_cancer, file = file.path(project_paths$data_processed, 'hg38', 'Cosmic_CancerGeneCensus_v101_GRCh38.rds'))

