# Script to generate refseq transcription and transcription start site tracks for hs1 genome
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hub_3671779_hs1&hgta_group=genes&hgta_track=hub_3671779_refSeqComposite&hgta_table=hub_3671779_refSeqComposite&hgta_doSchema=describe+table+schema


rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

conda_path <- '/home/ubuntu/miniconda3/condabin/conda'

genome <- 'hs1'

dir.create(file.path(project_paths$data_raw, genome))

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/gbdb/hs1/ncbiRefSeq/ncbiRefSeq.bb', ' ', file.path(project_paths$data_raw, genome)
  )
)

system(
  paste(conda_path, 'run -n ucsc_track_tools_env bigBedToBed',
        file.path(project_paths$data_raw, genome, 'ncbiRefSeq.bb'), file.path(project_paths$data_raw, genome, 'ncbiRefSeq.bed')
        )
  )


df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, genome, 'ncbiRefSeq.bed'), col_names = FALSE, delim = '\t')
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

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', txStart, txEnd),
    end = start + 1) %>%
  dplyr::select(-c(txStart, txEnd)) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_start.rds'))

# Transcription start and stop
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, name, name2) %>%
  dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = txStart,
    end = txEnd
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  IRanges::reduce()

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_reduced.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(txStart, txEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', txStart, txEnd),
    end = start + 1) %>%
  dplyr::select(-c(txStart, txEnd)) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_start_reduced.rds'))


