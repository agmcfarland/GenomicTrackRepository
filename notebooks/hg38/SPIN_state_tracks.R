# Script to generate SPIN tracks.
# https://github.com/ma-compbio/SPIN?tab=readme-ov-file

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome <- 'hg38'

url <- "https://drive.google.com/uc?export=download&id=1gdwtrhTctddO9TCBXBaZpZFOAHWCUTli"

download.file(url, file.path(project_paths$data_raw, genome, 'SPIN_tracks.bed'), method = "auto")

df_spin <- read.csv(file.path(project_paths$data_raw, genome, 'SPIN_tracks.bed'), sep ='\t', header = FALSE)
colnames(df_spin) <- c('chromosome', 'start', 'end', 'location')

for (df_spin_state_ in dplyr::group_split(df_spin, location)) {
  
  location_name <- base::unique(df_spin_state_$location)
  
  df_spin_state_ <- df_spin_state_ %>%
    dplyr::filter(!stringr::str_detect(chromosome, '_random|_alt|chrUn|_fix'))
  
  grange_location <- GenomicRanges::makeGRangesFromDataFrame(df_spin_state_, keep.extra.columns = TRUE)
  
  saveRDS(grange_location, file.path(project_paths$data_processed, genome, paste0('SPIN_', location_name, '.rds')))
  
}