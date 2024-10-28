library(dplyr)


#' List Non-Global Variables
#'
#' This function takes a list of all variables and a list of global variables to keep. It returns a list of non-global variables by finding the set difference between all variables and the global variables to keep.
#'
#' @param all_variables A character vector containing all variables.
#'
#' @param global_variables_to_keep A character vector containing global variables that should be retained.
#'
#' @return A character vector of non-global variables that are not in the list of global variables to keep.
#'
#' @examples
#' # List of all variables
#' all_vars <- c("var1", "var2", "var3", "var4")
#'
#' # List of global variables to keep
#' global_vars <- c("var1", "var3")
#'
#' # List non-global variables
#' non_global_vars <- list_non_global_variables(all_vars, global_vars)
#'
#' @export
list_non_global_variables <- function(all_variables, global_variables_to_keep) {
  return(setdiff(all_variables, global_variables_to_keep))
}

#' Check if a file or directory exists
#'
#' This function checks whether a file or directory exists at the specified path.
#'
#' @param x A character string specifying the path to the file or directory to be checked.
#' @param type A character string specifying the type of path to check, either 'file' (default) or 'directory'.
#'
#' @return NULL
#'
#' @examples
#' # Check if a file exists
#' check_exists('path/to/myfile.txt')
#'
#' # Check if a directory exists
#' check_exists('path/to/mydirectory', type = 'directory')
#'
#' @seealso \code{\link{file.exists}}, \code{\link{dir.exists}}
#'
#' @export
check_exists <- function(x, type = 'file'){
  if ( type == 'file'){
    if ( file.exists(x) == FALSE ){
      stop(paste0(x, ' does not exist'))
    }
  } 
  if ( type == 'directory'){
    if ( dir.exists(x) == FALSE ){
      stop(paste0(x, ' does not exist'))
    }
  } 
}

#' Build a set of project paths based on a root directory
#'
#' This function constructs a data frame containing various project paths by concatenating the specified 'root_directory' with subdirectory names. The resulting paths can be used to organize and access different project-related directories.
#'
#' @param root_directory A character string specifying the root directory where the project-related subdirectories will be created.
#'
#' @return A data frame containing project paths as character strings.
#'
#' @examples
#' # Define project root directory
#' root_dir <- 'path/to/project'
#'
#' # Build project paths
#' project_paths <- build_project_paths(root_dir)
#'
#' # Access specific paths
#' data_raw_path <- project_paths$data_raw
#'
#' @seealso \code{\link{check_exists}} - This function depends on 'check_exists' for verifying the existence of the 'root_directory'.
#'
#' @export
build_project_paths <- function(root_directory){
  
  check_exists(root_directory, type = 'directory')
  
  project_paths <-  data.frame(
    'data' = file.path(root_directory, 'data'),
    'data_processed' = file.path(root_directory, 'data', 'processed'),
    'data_raw' = file.path(root_directory, 'data', 'raw'),
    'data_interim' = file.path(root_directory, 'data', 'interim'),
    'data_external' = file.path(root_directory, 'data', 'external'),
    'references' = file.path(root_directory, 'references'),
    'src' = file.path(root_directory, 'src'),
    'reports' = file.path(root_directory, 'reports'),
    'integration_site_full' = file.path(root_directory, 'data', 'interim', 'aavenger', 'combined_runs')
  )
  return(project_paths)
}

#' List Installed R Packages
#'
#' This function retrieves a list of installed R packages and their versions.
#'
#' @param list_namespace_only Logical. If TRUE, the function lists only the packages that are loaded in the current namespace.
#'
#' @return A tibble with two columns, "Package" and "Version," listing the names and versions of installed R packages.
#'
#' @examples
#' list_installed_R_packages()
# list_installed_R_packages(list_namespace_only = TRUE)
#'
#' @import tibble
#' @import dplyr
#' @import utils
#'
#' @export
list_installed_R_packages <- function(list_namespace_only = FALSE) {
  df_installed_packages <- tibble::tibble(
    Package = names(installed.packages()[,3]),
    Version = unname(installed.packages()[,3])
  )
  
  if (list_namespace_only) {
    df_installed_packages <- df_installed_packages %>%
      dplyr::filter(Package %in% loadedNamespaces())
  }
  
  return(df_installed_packages)
}


convert_to_days <- function(time_vector) {
  # Apply the conversion to each element in the vector
  sapply(time_vector, function(time_string) {
    # Extract the time unit (last character) and numeric value
    unit <- substring(time_string, 1, 1)
    value <- as.numeric(substring(time_string, 2))
    
    # Convert based on the unit
    if (unit == "M") {
      return(value * 30)        # Convert months to days
    } else if (unit == "Y") {
      return(value * 365)       # Convert years to days
    } else if (unit == "D") {
      return(value)             # Days remain as is
    } else if (unit == "W") {
      return(value * 7)         # Convert weeks to days
    } else {
      return(NA)                # Return NA if unrecognized format
    }
  })
}


### REPEAT ANALYSIS FUNCTIONS TEMP START ###

#' Format AAVengeR sites
#' 
#' Applies formatting to standard AAVengeR sites table
#' 
format_aavenger_sites <- function(df) {
  df <- df %>%
    split_posid_into_chromosome_position() %>%
    ordered_aavenger_posid()
  
  return(df)
}

split_posid_into_chromosome_position <- function(df) {
  return(df %>%
           tidyr::separate(col = posid, into = c('chromosome', 'position', 'extra'), remove = FALSE) %>%
           dplyr::mutate(
             position = as.numeric(position),
             strand = gsub("[^+-]", "", posid)))
}
  

ordered_aavenger_posid <- function(df, remove_chromosome_number = TRUE) {
  df <- df %>%
    dplyr::mutate(
      chromosome_number = stringr::str_replace(chromosome, 'chr', ''),
      chromosome_number = stringr::str_replace(chromosome_number, 'X', '23'),
      chromosome_number = stringr::str_replace(chromosome_number, 'Y', '24'),
      chromosome_number = as.numeric(chromosome_number)
    ) %>%
    dplyr::arrange(chromosome_number, as.numeric(position))
  
  df$posid_factor <- factor(df$posid, levels = base::unique(df$posid))
  
  if (remove_chromosome_number) {
    df <- dplyr::select(df, -chromosome_number)
  }
  
  return(df)
}


relabel_repeat_categories <- function(df) {
  df <- df %>%
    dplyr::mutate(
      repeat_simplified = ifelse(is.na(repeat_class), 'not_a_repeat', repeat_class),
      repeat_simplified = ifelse(repeat_class == '', 'not_a_repeat', repeat_class),
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

unnest_mhc_table <- function(aavenger_mhc) {
  aavenger_mhc <- aavenger_mhc %>%
    tidyr::unnest(posids) %>%
    dplyr::rename(posid = posids)
  return(aavenger_mhc)
}

split_posid_into_separate_values <- function(df) {
  df <- df %>%
    tidyr::separate(posid, into = c('chromosome', 'start', 'extra'), remove = FALSE)
  return(df)
}

relabel_NA_repeat_as_not_a_repeat <- function(df) {
  df <- df %>%
    dplyr::mutate(repeat_class = ifelse(is.na(repeat_class), 'not_a_repeat', repeat_class)) %>%
  return(df)
}

find_overlaps_in_repeats <- function(aavenger_mhc, repeat_table) {

  grange_repeat_table <- GenomicRanges::makeGRangesFromDataFrame(repeat_table %>%
                                                                   dplyr::select(query_seq, query_start, query_end, repeat_name, repeat_class),
                           keep.extra.columns = TRUE,
                           ignore.strand = TRUE,
                           seqinfo = NULL,
                           seqnames.field = 'query_seq',
                           start.field = 'query_start',
                           end.field = 'query_end',
                           strand.field = 'strand',
                           starts.in.df.are.0based = FALSE)
  
  grange_mhc <-  GenomicRanges::makeGRangesFromDataFrame(aavenger_mhc %>%
                                                           unnest_mhc_table() %>%
                                                           split_posid_into_separate_values() %>%
                                                           dplyr::mutate(end = start),
                           keep.extra.columns = TRUE,
                           ignore.strand = TRUE,
                           seqinfo = NULL,
                           seqnames.field = 'chromosome',
                           start.field = 'start',
                           end.field = 'end',
                           strand.field = 'strand',
                           starts.in.df.are.0based = FALSE)
  
  
  overlaps <- GenomicRanges::findOverlaps(grange_mhc, grange_repeat_table) %>%
    as.data.frame()
    
  return(overlaps)                        
}

extract_overlap_repeat_by_indices <- function(overlaps, repeat_table) {
  
  extracted_repeat_indices <- lapply(unlist(overlaps$subjectHits), function(repeat_index){
    return(repeat_table[repeat_index, ])
    })
  
  extracted_repeat_indices <- do.call(rbind, extracted_repeat_indices)
  
  overlaps <- cbind(overlaps, extracted_repeat_indices)
  
  return(overlaps)
}

map_repeat_overlaps_to_mhc <- function(aavenger_mhc, repeat_table, overlaps) {

  overlaps <- extract_overlap_repeat_by_indices(repeat_table = repeat_table, overlaps = overlaps)
  
  aavenger_mhc <- aavenger_mhc %>%
    unnest_mhc_table() %>%
    split_posid_into_separate_values() %>%
    dplyr::mutate(row_index = 1:nrow(.))
  
  aavenger_mhc_repeats <- merge(
    aavenger_mhc,
    overlaps,
    by.x = 'row_index',
    by.y = 'queryHits',
    all.x = TRUE
  )
  
  return(aavenger_mhc_repeats)
}

count_repeat_class_per_cluster <- function(aavenger_mhc_repeats) {
  
  df_mhc_repeats_counts <- aavenger_mhc_repeats %>%
    # dplyr::mutate(repeat_class = ifelse(is.na(repeat_class), 'not_a_repeat', repeat_class)) %>%
    # relabel_repeat_categories() %>%
    relabel_repeat_categories() %>%
    dplyr::select(-repeat_class) %>%
    dplyr::rename(repeat_class = repeat_simplified) %>%
    dplyr::group_by(trial, subject, sample, clusterID, repeat_class) %>%
    dplyr::mutate(repeat_class_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(trial, subject, sample, clusterID, repeat_class, repeat_class_count) %>%
    base::unique() %>%
    dplyr::group_by(trial, subject, sample, clusterID) %>%
    dplyr::mutate(
      total_repeat_detections = sum(repeat_class_count),
      percentage_of_repeat_detections = 100 * (repeat_class_count/total_repeat_detections)
      )
  
  return(df_mhc_repeats_counts)
  
}

extract_most_abundant_repeat_class_per_cluster <- function(mhc_repeats_counts_per_cluster) {
  mhc_repeats_counts_per_cluster <- mhc_repeats_counts_per_cluster %>%
    dplyr::group_by(trial, subject, sample, clusterID) %>%
    dplyr::arrange(dplyr::desc(percentage_of_repeat_detections)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  return(mhc_repeats_counts_per_cluster)
}

### REPEAT ANALYSIS TEMP END ###