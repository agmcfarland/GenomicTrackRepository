

if (!file.exists(file.path(working_dir, 'hg38.fa.gz'))) {

  system(
    command = paste0(
      'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz', ' ', working_dir
    )
  )
  
  }


genomic_sequence <- Biostrings::readDNAStringSet(file.path(working_dir, 'hg38.fa.gz'))

chromosomes_to_keep <- names(genomic_sequence)[!stringr::str_detect(names(genomic_sequence), '_random|_alt|chrUn|_fix')]
  
genomic_sequence <- genomic_sequence[names(genomic_sequence) %in% chromosomes_to_keep]
  