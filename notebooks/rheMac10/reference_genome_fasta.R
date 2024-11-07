# Downloading NCBI hg38 assembly from hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

dir.create(file.path(project_paths$data_processed, 'rheMac10'))


if (!file.exists(file.path(project_paths$data_raw, 'rheMac10.fa.gz'))) {

  system(
    command = paste0(
      'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz', ' ', project_paths$data_raw
    )
  )
  
  }

genomic_sequence <- Biostrings::readDNAStringSet(file.path(project_paths$data_raw, 'rheMac10.fa.gz'))

chromosomes_to_keep <- names(genomic_sequence)[!stringr::str_detect(names(genomic_sequence), '_random|_alt|chrUn|_fix')]
  
genomic_sequence <- genomic_sequence[names(genomic_sequence) %in% chromosomes_to_keep]
  
Biostrings::writeXStringSet(x = genomic_sequence, filepath = file.path(project_paths$data_processed, 'rheMac10', 'rheMac10.fa.gz'), compress = TRUE)


if (!file.exists(file.path(project_paths$data_raw, 'rheMac10.2bit'))) {
  system(
    command = paste0(
      'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.2bit', ' ', file.path(project_paths$data_processed, 'rheMac10')
    )
  )
}


