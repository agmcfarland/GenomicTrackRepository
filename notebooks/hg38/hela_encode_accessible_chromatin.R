# Script to generate DNA accessible regions from the ENCODE database. Uses filtered peaks.

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome <- 'hg38'


custom_hela_s3_tracks <- data.frame(
  'local_filename' = c(
    'ENCFF024UJN.bed.gz',
    'ENCFF239FBO.bed.gz',
    'ENCFF548ARU.bed.gz'
  ),
  'database_path' = c(
    'https://www.encodeproject.org/files/ENCFF024UJN/@@download/ENCFF024UJN.bed.gz',
    'https://www.encodeproject.org/files/ENCFF239FBO/@@download/ENCFF239FBO.bed.gz',
    'https://www.encodeproject.org/files/ENCFF548ARU/@@download/ENCFF548ARU.bed.gz'
  ),
  'track_name' = c(
    'ENCFF024UJN_encode_filtered_dnaseI.rds',
    'ENCFF239FBO_encode_filtered_RAD21.rds',
    'ENCFF548ARU_encode_filtered_CTCF.rds'
  )
)

overwrite <- TRUE

for (raw_file_ in split(custom_hela_s3_tracks, 1:nrow(custom_hela_s3_tracks))) {
  
  print(raw_file_$local_filename)
  
  if (overwrite & file.exists(file.path(project_paths$data_raw, raw_file_$local_filename))) {
    
    file.remove(file.path(project_paths$data_raw, raw_file_$local_filename))
  }  
  
  system(
    command = paste0(
      'wget ', raw_file_$database_path, ' --directory-prefix ', project_paths$data_raw
    )
  )
  
  df_track <- vroom::vroom(file.path(project_paths$data_raw, raw_file_$local_filename), delim = '\t', col_names = FALSE) %>%
    dplyr::select(X1, X2, X3) %>%
    dplyr::filter(!stringr::str_detect(X1, '_random|_alt|chrUn|_fix')) %>%
    dplyr::rename(
      seqname = X1,
      start = X2, 
      end = X3
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  saveRDS(df_track, file.path(project_paths$data_processed, genome, raw_file_$track_name))
}
