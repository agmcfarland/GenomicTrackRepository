# Script to convert hg38 genomic tracks to HS1 (T2T) genomic coordinates


rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

conda_path <- '/home/ubuntu/miniconda3/condabin/conda'

source_genome <- 'hg38'
target_genome <- 'hs1'

encode_tracks <- c(
  "ENCFF024UJN_encode_filtered_dnaseI.rds",
  "ENCFF239FBO_encode_filtered_RAD21.rds",
  "ENCFF548ARU_encode_filtered_CTCF.rds",
  "ENCFF903JDG_encode_filtered_H3K4me3.rds",
  "ENCFF831XSS_encode_filtered_H3K27ac.rds",
  "ENCFF021PYM_encode_filtered_H3K9ac.rds",
  "ENCFF162RSB_encode_filtered_H3K4me1.rds",
  "ENCFF872YCK_encode_filtered_H3K9me3.rds",
  "ENCFF429OQI_encode_filtered_H3K4me2.rds",
  "ENCFF584RYA_encode_filtered_H3K27me3.rds",
  "ENCFF582MED_encode_filtered_H3K36me3.rds",
  "ENCFF406DAM_encode_filtered_H4K20me1.rds",
  "ENCFF094MFL_encode_filtered_H2AFZ.rds",
  "ENCFF238XWI_encode_filtered_H3K79me2.rds",
  "ENCFF951RPD_act_tcell_encode_filtered_CTCF.rds",
  "ENCFF116CIR_act_tcell_encode_filtered_H3K4me3.rds",
  "ENCFF056VQX_act_tcell_encode_filtered_H3K27me3.rds",
  "ENCFF729DVH_act_tcell_encode_filtered_H3K27ac.rds",
  "ENCFF849QFO_act_tcell_encode_filtered_DNase-seq.rds",
  "ENCFF869FHR_act_tcell_encode_filtered_H3K4me1.rds",
  "ENCFF340GQQ_act_tcell_encode_filtered_H3K9me3.rds",
  "ENCFF813XLL_act_tcell_encode_filtered_H3K36me3.rds",
  "ENCFF372SOB_act_tcell_encode_filtered_H3K36me3.rds",
  "ENCFF283HDN_act_tcell_encode_filtered_H3K4me3.rds"
  )


chain_file_path <- file.path(project_paths$data_processed, 'liftover_chains', 'hg38ToHs1.over.chain.gz')

overwrite <- FALSE

for (track_ in encode_tracks) {

output_file_path <- file.path(project_paths$data_processed, target_genome, track_)

  if (!file.exists(output_file_path) | overwrite) {
    
    print(paste0('processing: ', output_file_path))
  
    bedfile_name <- stringr::str_replace(track_, '.rds', '')
    
    all_bedfile_paths <- data.frame(
      'raw' = file.path(project_paths$data_interim, target_genome, paste0(bedfile_name, '.bed')),
      'lifted'= file.path(project_paths$data_interim, target_genome, paste0(bedfile_name, '_lifted.bed')),
      'unlifted'= file.path(project_paths$data_interim, target_genome, paste0(bedfile_name, '_unlifted.bed'))
      )
  
    df_track <- readRDS(file.path(project_paths$data_processed, source_genome, track_))
    
    df_bed_track <- data.frame(seqnames=seqnames(df_track), starts=format(start(df_track)-1, scientific=F), ends=format(end(df_track), scientific=F))
    
    write.table(df_bed_track, all_bedfile_paths$raw, quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    system(
      paste(conda_path, 'run -n ucsc_track_tools_env liftOver',
            all_bedfile_paths$raw, chain_file_path, all_bedfile_paths$lifted,  all_bedfile_paths$unlifted)
      )
    
    df_bed_track_lifted <- read.csv(all_bedfile_paths$lifted, sep = '\t', col.names = c('seqnames', 'start', 'stop'))
    
    df_bed_track_lifted_granges <- df_bed_track_lifted %>% 
      GenomicRanges::makeGRangesFromDataFrame()
    
    saveRDS(df_bed_track_lifted_granges, output_file_path)
    
  }
}
  
  
