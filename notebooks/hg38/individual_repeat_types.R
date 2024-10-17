# Script to extract hg38 integration sites from published gene therapy samples stored in github database
# https://github.com/helixscript/published_gt_sample_data/blob/main/published_gt_samples.tsv
# https://raw.githubusercontent.com/helixscript/published_gt_sample_data/refs/heads/main/published_gt_samples.tsv
# "all_intSites.tsv": {
#   // "type": "file",
#   // "description": "All integration sites in microb120. Taken by finding catting all intSites.",
#   // "source": ["John Everett", "VMHSC"],
#   // "link": "",
#   // "purpose": "Comparison of relative integration site frequency",
#   // "relative_location": "data/external"
# },


rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
genome <- 'hg38'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

project_paths$data_interim_analysis <- file.path(project_paths$data_interim, genome)
project_paths$data_processed_analysis <- file.path(project_paths$data_processed, genome)

dir.create(project_paths$data_interim_analysis, showWarnings = FALSE)
dir.create(project_paths$data_processed_analysis, showWarnings = FALSE)

relabel_repeat_categories <- function(df) {
  df <- df %>%
    dplyr::mutate(
      repeat_simplified = repeat_class,
      repeat_simplified = ifelse(is.na(repeat_simplified), 'not_a_repeat', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'LINE'), 'LINE', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'SINE'), 'SINE', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'LTR'), 'LTR', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'atellite'), 'satellite', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'etrop'), 'retrotransposon', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'DNA'), 'DNA', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'known'), 'unknown', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'rRNA|snRNA|scRNA|tRNA|srpRNA'), 'RNA', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'Low_complexity'), 'low_complexity', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'RC/Helitron'), 'other', repeat_simplified),
      repeat_simplified = ifelse(stringr::str_detect(repeat_simplified, 'Simple_repeat'), 'simple_repeat', repeat_simplified)
      ) 
  return(df)
}

overwrite <- TRUE

path_to_aavenger_repeat_table <- '/data/AAVengeR/data/genomeAnnotations/hg38.repeatTable.gz'

df_repeat_table <- vroom::vroom(path_to_aavenger_repeat_table)

df_repeat_table <- df_repeat_table %>%
  relabel_repeat_categories() %>%
  dplyr::select(-repeat_class) %>%
  dplyr::rename(repeat_class = repeat_simplified)


repeat_classes_to_process <- c(base::unique(df_repeat_table$repeat_class), 'all')

if (overwrite) {
    
  for (repeat_class_ in repeat_classes_to_process) {
    
    print(repeat_class_)
      
    if (repeat_class_ != 'all') {
      df_track <- df_repeat_table %>%
        dplyr::filter(repeat_class == repeat_class_)
    } else {
      df_track <- df_repeat_table
    }
    
    df_track <- df_track %>%
      dplyr::select(query_seq, query_start, query_end, repeat_class) %>%
      dplyr::filter(!stringr::str_detect(query_seq, '_random|_alt|chrUn|_fix')) %>%
      dplyr::rename(
        seqname = query_seq,
        start = query_start, 
        end = query_end
      ) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
    filename_ <- paste0('repeatclass_', repeat_class_, '.rds')
    
    saveRDS(df_track, file.path(project_paths$data_processed, genome, filename_))
  }
}
