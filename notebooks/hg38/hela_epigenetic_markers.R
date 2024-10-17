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
  epigenetic_mark = c("H3K4me3", #"H3K36me3", #"H3K27me3", 
                      "H3K27ac", "H3K9ac", "H3K4me1", "H3K9me3", "H3K4me2", "H3K27me3", "H3K36me3", "H4K20me1", "H2AFZ", "H3K79me2"),
  database_path = c(
    "https://www.encodeproject.org/files/ENCFF903JDG/@@download/ENCFF903JDG.bed.gz", #1
    # "https://www.encodeproject.org/files/ENCFF348FNV/@@download/ENCFF348FNV.bed.gz",#2
    # "https://www.encodeproject.org/files/ENCFF428MSN/@@download/ENCFF428MSN.bed.gz",#3
    "https://www.encodeproject.org/files/ENCFF831XSS/@@download/ENCFF831XSS.bed.gz",#4
    "https://www.encodeproject.org/files/ENCFF021PYM/@@download/ENCFF021PYM.bed.gz",#5
    "https://www.encodeproject.org/files/ENCFF162RSB/@@download/ENCFF162RSB.bed.gz",#6
    "https://www.encodeproject.org/files/ENCFF872YCK/@@download/ENCFF872YCK.bed.gz",#7
    "https://www.encodeproject.org/files/ENCFF429OQI/@@download/ENCFF429OQI.bed.gz",
    "https://www.encodeproject.org/files/ENCFF584RYA/@@download/ENCFF584RYA.bed.gz",
    "https://www.encodeproject.org/files/ENCFF582MED/@@download/ENCFF582MED.bed.gz",
    "https://www.encodeproject.org/files/ENCFF406DAM/@@download/ENCFF406DAM.bed.gz",
    "https://www.encodeproject.org/files/ENCFF094MFL/@@download/ENCFF094MFL.bed.gz",
    "https://www.encodeproject.org/files/ENCFF238XWI/@@download/ENCFF238XWI.bed.gz"
  )
) %>% 
  dplyr::mutate(
    local_filename = basename(database_path),
    track_name = paste0(stringr::str_replace(local_filename, '.bed.gz', ''), '_encode_filtered_', epigenetic_mark, '.rds'))

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

