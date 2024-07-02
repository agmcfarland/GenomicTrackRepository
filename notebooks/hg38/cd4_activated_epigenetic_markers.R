# Script to generate epigenetic markers for HeLa cells from the ENCODE database. Uses filtered peaks.
# https://www.encodeproject.org/search/?type=Experiment&control_type%21=%2A&status=released&perturbed=false&assay_title=Histone+ChIP-seq&assembly=GRCh38&biosample_ontology.term_name=HeLa-S3&target.investigated_as=histone

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome <- 'hg38'

chromatin_data <- data.frame(
  epigenetic_mark = c("CTCF", "H3K4me3", "H3K27me3", "H3K27ac", "DNase-seq", "H3K4me1", "H3K9me3", "H3K36me3", "H3K36me3", "H3K4me3"),
  database_path = c(
    "https://www.encodeproject.org/files/ENCFF951RPD/@@download/ENCFF951RPD.bed.gz",
    "https://www.encodeproject.org/files/ENCFF116CIR/@@download/ENCFF116CIR.bed.gz",
    "https://www.encodeproject.org/files/ENCFF056VQX/@@download/ENCFF056VQX.bed.gz",
    "https://www.encodeproject.org/files/ENCFF729DVH/@@download/ENCFF729DVH.bed.gz",
    "https://www.encodeproject.org/files/ENCFF849QFO/@@download/ENCFF849QFO.bed.gz",
    "https://www.encodeproject.org/files/ENCFF869FHR/@@download/ENCFF869FHR.bed.gz",
    "https://www.encodeproject.org/files/ENCFF340GQQ/@@download/ENCFF340GQQ.bed.gz",
    "https://www.encodeproject.org/files/ENCFF813XLL/@@download/ENCFF813XLL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF372SOB/@@download/ENCFF372SOB.bed.gz",
    "https://www.encodeproject.org/files/ENCFF283HDN/@@download/ENCFF283HDN.bed.gz"
  )
) %>% 
  dplyr::mutate(
    local_filename = basename(database_path),
    track_name = paste0(stringr::str_replace(local_filename, '.bed.gz', ''), '_act_tcell_encode_filtered_', epigenetic_mark, '.rds'))

overwrite <- TRUE

for (raw_file_ in split(chromatin_data, 1:nrow(chromatin_data))) {
  
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

if (FALSE %in% file.exists(file.path(project_paths$data_processed, genome, chromatin_data$track_name))) {
  print('Not all tracks were created')
}

