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
notebook_name <- 'reference_integration_sites'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

project_paths$data_interim_analysis <- file.path(project_paths$data_interim, notebook_name)
project_paths$data_processed_analysis <- file.path(project_paths$data_processed, notebook_name)

dir.create(project_paths$data_interim_analysis, showWarnings = FALSE)
dir.create(project_paths$data_processed_analysis, showWarnings = FALSE)

overwrite <- TRUE

system(
  command = paste0(
    'wget ', 'https://raw.githubusercontent.com/helixscript/published_gt_sample_data/refs/heads/main/published_gt_samples.tsv', ' --directory-prefix ', project_paths$data_interim_analysis
  )
)

df_samples <- vroom::vroom(file.path(project_paths$data_interim_analysis, 'published_gt_samples.tsv'), delim = '\t')

df_intsites <- vroom::vroom(file.path(project_paths$data_external, 'all_intSites.tsv.gz')) %>%
  dplyr::select(-c(nearestFeature, inFeature, nearestFeatureStrand, inFeatureExon, nearestFeatureDist, nearestOncoFeature, nearestOncoFeatureDist))

ref_genome_ <- 'hg38'

for (vector_type_ in c('gamma', 'lenti')) {

  print(vector_type_)
  
  df_temp_samples <- df_samples %>%
    dplyr::mutate(filter_vector = tolower(Vector)) %>%
    dplyr::filter(stringr::str_detect(filter_vector, vector_type_))
  
  unique_gtsps <- base::unique(df_temp_samples$GTSP)
  
  df_temp_intsites <- df_intsites %>%
    dplyr::filter(refGenome == ref_genome_) %>%
    dplyr::filter(internalSampleID %in% unique_gtsps)
  
  df_temp_intsites <- df_temp_intsites %>%
    dplyr::mutate(timepoint_in_days = convert_to_days(timePoint))
  
  write.csv(df_temp_intsites, file.path(project_paths$data_processed_analysis, paste0(vector_type_, '_hg38_published_gt_set.csv')), row.names = FALSE)
  
}
