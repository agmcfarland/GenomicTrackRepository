# Script to add finishing touches to the Vincent Wu HIV compartmentalization dataset.

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))


project_paths <- build_project_paths(project_dir)

df_all_sites <- read.csv(file.path(project_paths$data_raw, 'reference_integration_sites/wu_hiv/wu_hiv_hs1.csv'))

df_all_sites <- df_all_sites %>%
  relabel_repeat_categories()

sequencing_run_ids <- base::unique(df_all_sites$run_ID)
sequencing_run_ids <- paste0(file.path(project_paths$data_processed, 'reference_integration_sites', 'wu_hiv/'), sequencing_run_ids)

skip_mhc_calculation_ <- FALSE

if (!skip_mhc_calculation_) {
  
  df_repeat_table <- vroom::vroom(file.path('/data/AAVengeR/data/genomeAnnotations/hs1.repeatTable.gz'))
  
  df_top_mhc_repeats <- do.call(rbind, lapply(sequencing_run_ids, function(run_path) {
    
    # run_path <- sequencing_run_ids[[1]]
    
    df_mhc <- readRDS(file.path(run_path, 'core', 'multiHitClusters.rds'))
    
    overlaps <- find_overlaps_in_repeats(aavenger_mhc = df_mhc, repeat_table = df_repeat_table)
    
    df_mhc_repeats <- map_repeat_overlaps_to_mhc(aavenger_mhc = df_mhc, repeat_table = df_repeat_table, overlaps = overlaps)
    
    df_mhc_repeat_counts <- count_repeat_class_per_cluster(aavenger_mhc_repeats = df_mhc_repeats)
    
    df_mhc_repeat_counts_top <- extract_most_abundant_repeat_class_per_cluster(mhc_repeats_counts_per_cluster = df_mhc_repeat_counts)
    
    return(df_mhc_repeat_counts_top)
  }))
  
  
  df_top_mhc_repeats_with_metadata <- merge(
    df_top_mhc_repeats,
    df_metadata,
    by.x = 'sample',
    by.y = 'SpecimenAccNum',
    all.x = TRUE
  )
  
  saveRDS(df_top_mhc_repeats_with_metadata, file.path(project_paths$data_interim_analysis, 'majority_repeat_per_cluster.rds'))
  
}
# 
# df_raw_table <- vroom::vroom(file = file.path(project_paths$data_raw, 'cpgIslandExt.txt.gz'), col_names = FALSE, delim = '\t')
# colnames(df_raw_table) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
# 
# df_raw_table <- df_raw_table %>%
#   dplyr::filter(!stringr::str_detect(chrom, '_random|_alt|chrUn|_fix')) %>%
#   dplyr::mutate(strand = '*') %>%
#   dplyr::rename(
#     start = chromStart,
#     end = chromEnd)
# 
# grange_raw_table <- df_raw_table %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
# 
# saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, 'hg38', 'cpgIslandExt.rds'))
# 
# 
