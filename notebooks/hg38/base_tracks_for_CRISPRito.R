
## Based on an already developed track at ncbiRefSeqCurated_track.R
# These are formatted for pyranges used by CRISPRito. Will not work with InvestigateIntegrations heatmap functions

rm(list = ls())
project_dir <- '/data/GenomicTrackRepository'
# Libraries
libraries <- list('dplyr')
invisible(suppressMessages(do.call(library, libraries)))
invisible(suppressMessages(sapply(list.files(file.path(project_dir, 'src', 'Rlib'), pattern = "utils.R", full.names = TRUE), source)))

project_paths <- build_project_paths(project_dir)

df_refseq <- readRDS(file.path(project_paths$data_processed, 'hg38', 'ncbiRefSeqCurated_expanded.rds')) %>%
  as.data.frame() %>%
  dplyr::rename(
    chromosome = seqnames
    )

# Presence absence features
for (feature_ in c('exon', 'intron', '3UTR', '5UTR')) {
  df_feature_ <- df_refseq %>%
    dplyr::filter(feature == feature_) %>%
    dplyr::select(chromosome, start, end, feature)
  
  write.csv(
    df_feature_,
    file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', feature_, '.csv')),
    row.names = F
    )
}

# Named features
df_feature_ <- df_refseq %>%
  dplyr::select(chromosome, start, end, name2) %>%
  dplyr::rename(annotation = name2)

write.csv(
  df_feature_, 
  file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', 'gene', '.csv')),
  row.names = F
  )

df_feature_ <- readRDS(file.path(project_paths$data_processed, 'hg38', 'Cosmic_CancerGeneCensus_v101_GRCh38.rds')) %>%
  as.data.frame() %>%
  dplyr::rename(
    chromosome = seqnames,
    annotation = gene_symbol
  ) %>%
  dplyr::select(chromosome, start, end, annotation)

write.csv(
  df_feature_, 
  file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', 'oncogene', '.csv')),
  row.names = F
)
