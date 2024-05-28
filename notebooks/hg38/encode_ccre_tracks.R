# Script to generate ENCODE cCRE tracks for hg38 genome
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=encodeCcreCombined&hgta_table=encodeCcreCombined&hgta_doSchema=describe+table+schema
# more info: https://screen.encodeproject.org/

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome_assembly <- 'hg38'

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb', ' ', project_paths$data_raw
  )
)

system(
  command = paste0(
    'gunzip', ' ', project_paths$data_raw, '/mm9.fa.gz'
  )
)

conda_path <- '/home/ubuntu/miniconda3/condabin/conda'

bigbed_filepath <- file.path(project_paths$data_raw, 'encodeCcreCombined.bb')
  
bed_filepath <- file.path(project_paths$data_raw, 'encodeCcreCombined.bed')

system(
  command = paste0(conda_path, ' ', 'run ', '-n ', 'vittoria_env bigBedToBed ', bigbed_filepath, ' ', bed_filepath)
  )

df_raw_table <- vroom::vroom(file.path(bed_filepath), col_names = FALSE, delim = '\t')

colnames(df_raw_table) <-c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "reserved", "ccre", "encodeLabel", "zScore", "ucscLabel", "accessionLabel", "description")

df_raw_table <- df_raw_table %>%
  dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) %>%
  dplyr::mutate(strand = '*') %>%
  dplyr::rename(
    start = chromStart,
    end = chromEnd) %>%
  dplyr::select(-c(thickStart, thickEnd, reserved))

lapply(dplyr::group_split(df_raw_table, ucscLabel), function(df_ucsc_label) {

  # df_ucsc_label <- dplyr::group_split(df_raw_table, ucscLabel)[[1]]
  
  output_filename <- file.path(project_paths$data_processed, genome_assembly, paste0('encode_ccre_', base::unique(df_ucsc_label$ucscLabel), '.rds'))
  
  print(output_filename)
  
  grange_raw_table <- df_ucsc_label %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  saveRDS(grange_raw_table, file = output_filename)
  
})
