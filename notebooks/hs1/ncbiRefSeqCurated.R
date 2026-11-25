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

dir.create(file.path(project_paths$data_raw, genome), showWarnings = F)

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


df_raw_table <- read.csv(file = file.path(project_paths$data_raw, genome, 'ncbiRefSeq.bed'), header = FALSE, sep = '\t') # vroom::vroom(file = file.path(project_paths$data_raw, genome, 'ncbiRefSeq.bed'), col_names = FALSE, delim = '\t')
# colnames(df_raw_table) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount',
#                             'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')
colnames(df_raw_table) <- c("chrom",
                            'chromStart', #"chromStart",
                            'chromEnd', #"chromEnd",
                            "name",
                            "score",
                            "strand",
                            'thickStart', #"thickStart",
                            'thickEnd', #"thickEnd",
                            "reserved",
                            'blockCount', #"blockCount",
                            'blockSizes', #"blockSizes",
                            'chromStarts', #"chromStarts",
                            "name2",
                            "cdsStartStat",
                            "cdsEndStat",
                            "exonFrames",
                            "type",
                            "geneName",
                            "geneName2",
                            "geneType"
                            )
# 
# # View(df_raw_table)
# ## InvestigateIntegrations::extend_raw_refseq_table() work flow input START
# df_raw_table_invint <- df_raw_table %>%
#   dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix'))
# 
# saveRDS(df_raw_table_invint, file.path(project_paths$data_processed, 'hs1', 'filtered_ncbirefseqcurated.rds'))
# 
# 
# # Extend refseq table
# df_processed_table <- InvestigateIntegrations::extend_raw_refseq_table(readRDS(file.path(project_paths$data_processed, 'hs1', 'filtered_ncbirefseqcurated.rds')))
# 
# df_processed_table <- df_processed_table %>%
#   dplyr::left_join( 
#     df_raw_table_invint %>%
#       dplyr::select(chrom, name, name2),
#     by = c('chrom', 'name')
#   )
# 
# # Store processed table results
# saveRDS(df_processed_table, file.path(project_paths$data_processed, 'hs1', 'expanded_ncbirefseqcurated.rds'))
# 
# ## InvestigateIntegrations::extend_raw_refseq_table() work flow input END

# Resume other methods of extracting info

df_raw_table <- df_raw_table %>%
  dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) #%>%
  # dplyr::mutate(strand = '*')

# Transcription start and stop
grange_raw_table <- df_raw_table %>%
  dplyr::select(chromStart, chromEnd, chrom, name, name2, strand) %>%
  # dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = chromStart,
    end = chromEnd
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(chromStart, chromEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', chromStart, chromEnd),
    end = start + 1) %>%
  dplyr::select(-c(chromStart, chromEnd)) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_start.rds'))

# Transcription start and stop
grange_raw_table <- df_raw_table %>%
  dplyr::select(chromStart, chromEnd, chrom, name, name2) %>%
  # dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = chromStart,
    end = chromEnd
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  IRanges::reduce()

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_reduced.rds'))

# Transcription start site
grange_raw_table <- df_raw_table %>%
  dplyr::select(chromStart, chromEnd, chrom, strand, name, name2) %>%
  dplyr::mutate(
    start = ifelse(strand == '+', chromStart, chromEnd),
    end = start + 1) %>%
  dplyr::select(-c(chromStart, chromEnd)) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hs1', 'ncbiRefSeqCurated_transcription_start_reduced.rds'))


