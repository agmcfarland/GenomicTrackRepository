# Script to generate GC content track for hg38 genome
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# more info: https://www.biostars.org/p/121751/

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

genome_assembly <- 'mm9'

system(
  command = paste0(
    'rsync', ' -avzP', ' rsync://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz', ' ', project_paths$data_raw
  )
)

system(
  command = paste0(
    'gunzip', ' ', project_paths$data_raw, '/mm9.fa.gz'
  )
)


python_script_path <- file.path(project_dir, 'src', 'Utils', 'calculate_gc_content_by_window.py')

genome_fasta_path <- file.path(project_paths$data_raw, 'mm9.fa')

window_size <- rev(c(100, 1000, 10000, 100000, 1000000))

overwrite = FALSE
for (window_ in window_size) {
  output_file_name <- file.path(project_paths$data_raw, paste0('gc_content_', window_, '.csv'))
  if (!file.exists(output_file_name)|overwrite){
    print(window_)
    system(paste0('/home/ubuntu/miniconda3/condabin/conda', ' run', ' -n', ' repeat_project', ' python ', python_script_path, ' ', genome_fasta_path, ' ', output_file_name, ' ', window_))
    }
  }


gc_content_threshold <- 0 #45

for (window_ in window_size) {
  print(window_)
  grange_raw_table <- vroom::vroom(file.path(project_paths$data_raw, paste0('gc_content_', window_, '.csv')), delim = ',') %>%
    dplyr::filter(gc_percentage >= gc_content_threshold) %>%
    dplyr::mutate(strand = '*') %>%
    dplyr::rename(
      start = start_window,
      end = end_window
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  saveRDS(grange_raw_table, file = file.path(project_paths$data_processed, genome_assembly, paste0('gc_content_', window_, '.rds')))
  
  }
