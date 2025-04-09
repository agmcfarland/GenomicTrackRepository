# Script to generate refseq transcription and transcription start site tracks for hg38 genome

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz', ' ', project_paths$data_raw
  )
)


df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, 'ncbiRefSeqCurated.txt.gz'), col_names = FALSE, delim = '\t')
colnames(df_raw_table) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                            'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')

df_raw_table <- df_raw_table %>%
  dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) %>%
  dplyr::mutate(strand = '*')

# Transcription start and stop
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, name, name2) %>%
  dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = txStart,
    end = txEnd
  ) %>%
  base::unique() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_transcription.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', txStart, txEnd),
    end = start + 1) %>%
  dplyr::select(-c(txStart, txEnd)) %>%
  base::unique() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_transcription_start.rds'))


# Transcription start and stop
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, name, name2) %>%
  dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = txStart,
    end = txEnd
  ) %>%
  base::unique() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  IRanges::reduce()

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_transcription_reduced.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', txStart, txEnd),
    end = start + 1) %>%
  dplyr::select(-c(txStart, txEnd)) %>%
  base::unique() %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_transcription_start_reduced.rds'))



# Make extended list
df_raw_table <- read.csv(file = file.path(project_paths$data_raw, 'ncbiRefSeqCurated.txt.gz'), header = FALSE, sep = '\t')
colnames(df_raw_table) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
                            'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')


df_gene <- df_raw_table %>%
  # as.data.frame() %>%
  # dplyr::rename(
  #   chrom = seqnames,
  #   txStart = start,
  #   txEnd = end
  # ) %>%
  # dplyr::select(-width) %>%
  dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) %>%
  dplyr::mutate(
    cdsStartStat = ifelse(cdsStart == cdsEnd, 'none', 'cmpl'),
    custom_index = seq(1, dplyr::n(), 1),
    name = paste0(name, '_', custom_index)
  ) %>%
  dplyr::mutate(unique_id = paste(chrom, txStart, txEnd, cdsStart, cdsEnd, exonCount, sep = '_')) %>%
  dplyr::group_by(unique_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(cdsStartStat != 'none')

df_expanded_genes <- df_gene %>%
  InvestigateIntegrations::extend_raw_refseq_table()

df_expanded_genes2 <- merge(
  df_expanded_genes,
  df_gene %>% dplyr::select(name, name2) %>%
    base::unique(),
  by = 'name'
)
  
df_expanded_genes2 <- df_expanded_genes2 %>%
  dplyr::mutate(
    feature_new = dplyr::case_when(
      (stringr::str_detect(feature, 'UTR_left') & strand == '+') ~ stringr::str_replace(feature, 'UTR_left', '5UTR'),
      (stringr::str_detect(feature, 'UTR_left') & strand == '-') ~ stringr::str_replace(feature, 'UTR_left', '3UTR'),
      (stringr::str_detect(feature, 'UTR_right') & strand == '+') ~ stringr::str_replace(feature, 'UTR_right', '3UTR'),
      (stringr::str_detect(feature, 'UTR_right') & strand == '-') ~ stringr::str_replace(feature, 'UTR_right', '5UTR'),
      TRUE ~ feature
    ))

df_expanded_genes2 <- df_expanded_genes2 %>%
  dplyr::select(-feature) %>%
  dplyr::rename(feature = feature_new) %>%
  base::unique()

grange_raw_table <- df_expanded_genes2 %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_expanded.rds'))

write.csv(df_expanded_genes2, file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_expanded.csv'), row.names = FALSE)


# df_expanded_genes2 <- read.csv(file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_expanded.csv'))

