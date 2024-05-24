# Script to generate CpG island track for hg38 genome
# https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/

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
    'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes', ' ', project_paths$data_raw
  )
)

df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, 'mm9.chrom.sizes'), col_names = FALSE, delim = '\t')
colnames(df_raw_table) <- c("chromosome", "size")

df_raw_table <- df_raw_table %>%
  dplyr::filter(!stringr::str_detect(chromosome, '_random|_alt|chrUn|_fix')) %>%
  dplyr::mutate(
    strand = '*',
    start = 0,
    end = size)

grange_raw_table <- df_raw_table %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, genome_assembly, 'mm9.chrom.sizes.rds'))

