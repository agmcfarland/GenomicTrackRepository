# Downloading NCBI hs1 assembly from hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)


if (!file.exists(file.path(project_paths$data_raw, 'hs1.fa.gz'))) {

  system(
    command = paste0(
      'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz', ' ', project_paths$data_raw
    )
  )
  
  }

genomic_sequence <- Biostrings::readDNAStringSet(file.path(project_paths$data_raw, 'hs1.fa.gz'))

chromosomes_to_keep <- names(genomic_sequence)[!stringr::str_detect(names(genomic_sequence), '_random|_alt|chrUn|_fix')]
  
genomic_sequence <- genomic_sequence[names(genomic_sequence) %in% chromosomes_to_keep]
  
Biostrings::writeXStringSet(x = genomic_sequence, filepath = file.path(project_paths$data_processed, 'hs1', 'hs1.fasta.gz'), compress = TRUE)
