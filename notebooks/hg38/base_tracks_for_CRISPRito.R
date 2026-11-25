
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

# Gene component annotations
for (feature_ in c('exon', 'intron', '3UTR', '5UTR')) {
  df_feature_ <- df_refseq %>%
    dplyr::filter(feature == feature_) %>%
    dplyr::rename(annotation = feature) %>% 
    dplyr::select(chromosome, start, end, annotation) %>% 
    base::unique()
  
  write.csv(
    df_feature_,
    file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', feature_, '.csv')),
    row.names = F
    )
}

# Named features
df_feature_ <- df_refseq %>%
  dplyr::select(chromosome, start, end, name2) %>%
  dplyr::rename(annotation = name2) %>% 
  base::unique()

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
  dplyr::select(chromosome, start, end, annotation) %>% 
  base::unique()

write.csv(
  df_feature_, 
  file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', 'oncogene', '.csv')),
  row.names = F
)

# cCRE
df_ccre <- readRDS(file.path(project_paths$data_processed, 'hg38', 'encode_ccre_all.rds')) %>%
  as.data.frame() %>%
  dplyr::rename(
    chromosome = seqnames,
    annotation = ucscLabel
  ) %>%
  dplyr::select(chromosome, start, end, annotation) %>% 
  base::unique()

write.csv(
  df_ccre, 
  file.path(project_paths$data_processed, 'hg38', paste0('crisprito_', 'encode_ccre_all', '.csv')),
  row.names = F
)
