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


