# Script to generate CpG island track for hg38 genome
# https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=cpgIslandExt

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome_assembly <- 'mm9'

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/cpgIslandExt.txt.gz', ' ', project_paths$data_raw
  )
)

df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, 'cpgIslandExt.txt.gz'), col_names = FALSE, delim = '\t')
colnames(df_raw_table) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")

df_raw_table <- df_raw_table %>%
  dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) %>%
  dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = chromStart,
    end = chromEnd)

grange_raw_table <- df_raw_table %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, genome_assembly, 'cpgIslandExt.rds'))

