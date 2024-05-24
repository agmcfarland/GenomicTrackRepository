# Script to generate refseq transcription and transcription start site tracks for mm9 genome
# https://genome.ucsc.edu/cgi-bin/hgTables?db=mm9&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_doSchema=describe+table+schema

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

genome_assembly <- 'mm9'

project_paths <- build_project_paths(project_dir)

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz', ' ', project_paths$data_raw
  )
)


df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, 'refGene.txt.gz'), col_names = FALSE, delim = '\t')
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
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, genome_assembly, 'ncbiRefSeqCurated_transcription.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', txStart, txEnd),
    end = start + 1) %>%
  dplyr::select(-c(txStart, txEnd)) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, genome_assembly, 'ncbiRefSeqCurated_transcription_start.rds'))


