#!/usr/bin/env Rscript
#
# annotate_integration_sites.R
#
# Annotate the integration sites in a SINGLE AAVengeR sites.rds table.
#
# Standalone script: it does not source any other project's utils.R or use
# a project_dir/build_project_paths() layout. The only paths it needs are
# passed on the command line or set as constants near the top of main()
# (REFERENCE DATA PATHS section) — edit those in place if they move.
#
#   1. Loads the sites table (checking for `posid` and `sample` columns).
#   2. Runs each annotation module (GC content, gene location, orientation,
#      expression, nuclear localization, repeats, ENCODE/genomic features,
#      Reactome pathways, oncogene lists, ZNF/KRAB-ZNF genes, COSMIC).
#   3. Combines everything into one wide table.
#   4. Writes all output under <output-dir>/<run-name>/, with per-feature
#      intermediate files in an intermediate/ subfolder that can optionally
#      be deleted once the final table has been written (--clean).
#
# USAGE
#   Rscript annotate_integration_sites.R --sites-file SITES_RDS --output-dir OUTPUT_DIR [options]
#
# --sites-file and --output-dir are required; the script errors out if
# either is missing.
#
# EXAMPLES
#   # Simplest case: run_name is taken from the file name
#   Rscript annotate_integration_sites.R \
#       --sites-file /data/aavenger_runs/GTSP1234/sites.rds \
#       --output-dir /data/annotated_runs
#
#   # Explicit run name, delete intermediates when done
#   Rscript annotate_integration_sites.R \
#       --sites-file /data/aavenger_runs/GTSP1234/sites.rds \
#       --run-name GTSP1234 \
#       --output-dir /data/annotated_runs \
#       --clean
#
# NOTE: sample-level metadata (emi_experimental_mouse_sample_metadata.csv)
# is intentionally no longer required. Any of the original metadata-derived
# columns (Emi_StudyID, CD8_present, etc.) are gone; everything downstream
# now runs off of whatever columns AAVengeR itself puts in the sites table.

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Install with install.packages('optparse').")
  }
  library(optparse)
  library(dplyr)
  library(magrittr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

build_option_parser <- function() {
  option_list <- list(
    make_option(c("-s", "--sites-file"), type = "character", default = NULL,
                help = "REQUIRED. Path to the AAVengeR sites.rds table."),
    make_option(c("-o", "--output-dir"), type = "character", default = NULL,
                help = "REQUIRED. Base directory under which a <run-name>/ folder is created."),
    make_option(c("-r", "--run-name"), type = "character", default = NULL,
                help = "Identifier for this run; used as the output folder name. [default: basename of --sites-file, minus extension]"),
    make_option(c("--reference-cache-dir"), type = "character", default = NULL,
                help = paste0("Shared, cross-run cache for expensive reference downloads/conversions ",
                               "(refseq table, Reactome table). Not touched by --clean. ",
                               "[default: <output-dir>/_reference_cache]")),
    make_option(c("--clean"), action = "store_true", default = FALSE,
                help = "Delete this run's intermediate/ folder once the final tables have been written."),
    make_option(c("--write-csv"), action = "store_true", default = FALSE,
                help = "Also write CSV copies of the final sites/annotated-features tables."),
    make_option(c("--keep-going"), action = "store_true", default = TRUE,
                help = "Continue running remaining annotation modules if one module fails. [default: TRUE]")
  )
  OptionParser(
    usage = "%prog --sites-file SITES_RDS --output-dir OUTPUT_DIR [options]",
    option_list = option_list,
    description = "Annotate a single AAVengeR sites.rds table."
  )
}

# Run one annotation module, logging and (optionally) tolerating failure so a
# single missing reference file doesn't take down the whole pipeline.
run_step <- function(step_name, fn, keep_going = TRUE) {
  message(sprintf("[%s] %s ...", format(Sys.time(), "%H:%M:%S"), step_name))
  if (keep_going) {
    result <- tryCatch(
      fn(),
      error = function(e) {
        warning(sprintf("Step '%s' failed: %s -- skipping.", step_name, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
  } else {
    result <- fn()
  }
  invisible(result)
}

# ---------------------------------------------------------------------------
# Loading a single AAVengeR sites table (replaces pull_aavenger_run_inputs /
# pull_aavenger_runs / make_full_sites_table, which were batch + metadata
# driven).
# ---------------------------------------------------------------------------

load_sites_table <- function(sites_file) {
  df_sites <- readRDS(sites_file)

  if (!"posid" %in% colnames(df_sites)) {
    if (!all(c("chromosome", "strand", "position") %in% colnames(df_sites))) {
      stop("sites table has no 'posid' column and is missing chromosome/strand/position to build one")
    }
    df_sites$posid <- paste0(df_sites$chromosome, df_sites$strand, df_sites$position)
  }

  if (!"sample" %in% colnames(df_sites)) {
    stop("sites table has no 'sample' column")
  }

  df_sites %>% dplyr::mutate(unique_integration_site_id = seq(1, dplyr::n()))
}

# ---------------------------------------------------------------------------
# Annotation modules
#
# Each module takes what it needs, saves its own completed_*_feature.rds
# into `intermediate_dir`, and returns the path it wrote (invisibly).
# ---------------------------------------------------------------------------

annotate_gc_content <- function(df_sites, intermediate_dir, genomic_tracks_path) {
  genomic_sequence <- Biostrings::readDNAStringSet(file.path(genomic_tracks_path, "hg38.fasta.gz"))

  chromosome_lengths <- read.csv(file.path(genomic_tracks_path, "hg38.chrom.sizes.txt"), sep = ",")
  chromosome_lengths <- chromosome_lengths[colnames(chromosome_lengths) %in% c("seqname", "chr_length")]
  chromosome_lengths <- chromosome_lengths %>%
    dplyr::filter(!stringr::str_detect(seqname, "_"))

  testthat::expect_equal(nrow(chromosome_lengths), 25)

  result <- InvestigateIntegrations::GC_percentage_per_insertion(
    df_sites %>%
      dplyr::mutate(
        strand = "*",
        start = position,
        end = position + 1,
        mid = position,
        seqname = chromosome
      ) %>%
      dplyr::select(unique_integration_site_id, seqname, start, mid, end, strand),
    genomic_sequence = genomic_sequence,
    window_sizes_to_test = c(50, 100, 1000),
    chromosome_lengths = chromosome_lengths,
    fasta_name = "genome"
  )

  result <- result %>%
    dplyr::select(unique_integration_site_id, genome_gc.fasta.50, genome_gc.fasta.100, genome_gc.fasta.1000)
  colnames(result) <- stringr::str_replace_all(colnames(result), "genome_gc.fasta.", "gc_bp_")

  out <- file.path(intermediate_dir, "completed_gc_content_per_int_site_feature.rds")
  saveRDS(result, out)
  invisible(out)
}

# Downloads/builds (and caches, once, in reference_cache_dir) the expanded
# NCBI refseq table, then annotates intron/exon/UTR occurrence for this run.
annotate_gene_location <- function(df_sites, df_hits_grange, intermediate_dir, reference_cache_dir) {
  refseq_cache <- file.path(reference_cache_dir, "expanded_ncbirefseqcurated.rds")

  if (!file.exists(refseq_cache)) {
    dir.create(reference_cache_dir, showWarnings = FALSE, recursive = TRUE)
    raw_gz <- file.path(reference_cache_dir, "ncbiRefSeqCurated.txt.gz")
    system(paste0(
      "rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz ",
      reference_cache_dir
    ))

    df_raw_table <- read.csv(file = raw_gz, header = FALSE, sep = "\t")
    colnames(df_raw_table) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                                 "exonCount", "exonStarts", "exonEnds", "score", "name2",
                                 "cdsStartStat", "cdsEndStat", "exonFrames")
    df_raw_table <- df_raw_table %>%
      dplyr::filter(!stringr::str_detect(chrom, "_random|_alt|chrUn|_fix"))

    df_processed_table <- InvestigateIntegrations::extend_raw_refseq_table(df_raw_table)
    saveRDS(df_processed_table, refseq_cache)
  }

  df_processed_table <- readRDS(refseq_cache) %>%
    dplyr::group_by(name, feature) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(order = ifelse(strand == "+", seq(1, dplyr::n()), rev(seq(1, dplyr::n())))) %>%
    dplyr::ungroup()

  subject_grange <- df_processed_table
  rownames(subject_grange) <- seq_along(1:nrow(subject_grange))
  subject_grange <- subject_grange %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  overlap_results <- InvestigateIntegrations::extract_overlaps(query = df_hits_grange, subject = subject_grange)

  overlap_results <- overlap_results %>%
    dplyr::mutate(
      is_rRNA = ifelse(stringr::str_detect(name_subject, "NR"), TRUE, FALSE),
      is_mRNA = ifelse(stringr::str_detect(name_subject, "NM"), TRUE, FALSE)
    ) %>%
    dplyr::mutate(
      feature_subject = dplyr::case_when(
        (feature_subject == "UTR_left"  & strand_subject == "+") ~ "UTR5",
        (feature_subject == "UTR_right" & strand_subject == "+") ~ "UTR3",
        (feature_subject == "UTR_left"  & strand_subject == "-") ~ "UTR3",
        (feature_subject == "UTR_right" & strand_subject == "-") ~ "UTR5",
        TRUE ~ feature_subject
      )
    ) %>%
    dplyr::mutate(
      in_UTR5_mRNA   = ifelse(feature_subject == "UTR5"   & is_mRNA, TRUE, FALSE),
      in_UTR3_mRNA   = ifelse(feature_subject == "UTR3"   & is_mRNA, TRUE, FALSE),
      in_exon_mRNA   = ifelse(feature_subject == "exon"   & is_mRNA, TRUE, FALSE),
      in_intron_mRNA = ifelse(feature_subject == "intron" & is_mRNA, TRUE, FALSE),
      in_UTR5_rRNA   = ifelse(feature_subject == "UTR5"   & is_rRNA, TRUE, FALSE),
      in_UTR3_rRNA   = ifelse(feature_subject == "UTR3"   & is_rRNA, TRUE, FALSE),
      in_exon_rRNA   = ifelse(feature_subject == "exon"   & is_rRNA, TRUE, FALSE),
      in_intron_rRNA = ifelse(feature_subject == "intron" & is_rRNA, TRUE, FALSE)
    ) %>%
    dplyr::group_by(unique_integration_site_id_query) %>%
    dplyr::mutate(
      first_intron_rRNA = ifelse(order_subject == 1 & feature_subject == "intron" & is_rRNA, TRUE, FALSE),
      first_intron_mRNA = ifelse(order_subject == 1 & feature_subject == "intron" & is_mRNA, TRUE, FALSE),
      rRNA_coding_isoforms = sum(is_rRNA),
      mRNA_coding_isoforms = sum(is_mRNA)
    ) %>%
    dplyr::ungroup()

  overlap_results_summarized <- overlap_results %>%
    dplyr::group_by(unique_integration_site_id_query) %>%
    dplyr::mutate(
      rRNA_coding = max(is_rRNA),
      mRNA_coding = max(is_mRNA),
      in_first_intron_rRNA = max(first_intron_rRNA),
      in_UTR5_rRNA = max(in_UTR5_rRNA),
      in_UTR3_rRNA = max(in_UTR3_rRNA),
      in_exon_rRNA = max(in_exon_rRNA),
      in_intron_rRNA = max(in_intron_rRNA),
      in_first_intron_mRNA = max(first_intron_mRNA),
      in_UTR5_mRNA = max(in_UTR5_mRNA),
      in_UTR3_mRNA = max(in_UTR3_mRNA),
      in_exon_mRNA = max(in_exon_mRNA),
      in_intron_mRNA = max(in_intron_mRNA)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_integration_site_id_query,
                  in_UTR5_rRNA, in_UTR3_rRNA, in_exon_rRNA, in_intron_rRNA, in_first_intron_rRNA,
                  in_UTR5_mRNA, in_UTR3_mRNA, in_exon_mRNA, in_intron_mRNA, in_first_intron_mRNA,
                  rRNA_coding_isoforms, mRNA_coding_isoforms) %>%
    base::unique() %>%
    dplyr::mutate(is_intergenic = FALSE) %>%
    dplyr::rename(unique_integration_site_id = unique_integration_site_id_query)

  testthat::expect_equal(length(base::unique(overlap_results_summarized$unique_integration_site_id)),
                          length(overlap_results_summarized$unique_integration_site_id))

  df_intergenic <- df_sites %>%
    dplyr::filter(!unique_integration_site_id %in% overlap_results_summarized$unique_integration_site_id) %>%
    dplyr::select(unique_integration_site_id) %>%
    dplyr::mutate(is_intergenic = TRUE)

  overlap_results_summarized <- dplyr::bind_rows(overlap_results_summarized, df_intergenic) %>%
    base::replace(is.na(.), FALSE)

  testthat::expect_equal(nrow(overlap_results_summarized), nrow(df_sites))

  out <- file.path(intermediate_dir, "completed_genic_integration_feature.rds")
  saveRDS(overlap_results_summarized, out)
  invisible(out)
}

annotate_orientation <- function(df_sites, intermediate_dir) {
  df_sites_orientation <- df_sites %>%
    dplyr::select(unique_integration_site_id, strand, nearestGeneStrand) %>%
    dplyr::mutate(integration_in_same_orientation_as_nearestGene = ifelse(strand == nearestGeneStrand, TRUE, FALSE)) %>%
    dplyr::select(unique_integration_site_id, integration_in_same_orientation_as_nearestGene)

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_orientation))

  out <- file.path(intermediate_dir, "completed_integration_orientation_feature.rds")
  saveRDS(df_sites_orientation, out)
  invisible(out)
}

annotate_expression <- function(df_sites, intermediate_dir, rna_seq_abundance_path) {
  if (!file.exists(rna_seq_abundance_path)) {
    stop("RNA-seq abundance file not found: ", rna_seq_abundance_path)
  }

  df_rna_seq <- read.csv(rna_seq_abundance_path)

  df_expressed <- df_sites %>%
    dplyr::select(unique_integration_site_id, nearestGene) %>%
    tidyr::separate_rows(nearestGene, sep = ",", convert = TRUE) %>%
    dplyr::mutate(nearestGene = stringr::str_replace_all(nearestGene, " ", ""))

  df_expressed <- merge(
    df_expressed,
    df_rna_seq %>% dplyr::select(-ENSEMBL_ID),
    by.x = "nearestGene", by.y = "gene_symbol",
    all.x = TRUE
  )

  df_expressed_summarized <- df_expressed %>%
    tidyr::pivot_longer(cols = -c(nearestGene, unique_integration_site_id)) %>%
    dplyr::mutate(
      name = dplyr::case_when(
        name == "X19_S19_R1_001" ~ "CD4_uninfected_expression",
        name == "X20_S20_R1_001" ~ "CD4_HIV_infected_expression",
        TRUE ~ "lookatme"
      )
    ) %>%
    dplyr::group_by(unique_integration_site_id, name) %>%
    dplyr::mutate(
      max_expression_value = max(value),
      min_expression_value = min(value),
      max_expression_value = ifelse(is.na(max_expression_value), 0, max_expression_value),
      min_expression_value = ifelse(is.na(min_expression_value), 0, min_expression_value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-value, -nearestGene) %>%
    base::unique() %>%
    tidyr::pivot_wider(names_from = "name", values_from = c("max_expression_value", "min_expression_value"), values_fill = 0)

  testthat::expect_equal(nrow(df_sites), nrow(df_expressed_summarized))

  out <- file.path(intermediate_dir, "completed_expression_profiles_feature.rds")
  saveRDS(df_expressed_summarized, out)
  invisible(out)
}

annotate_nuclear_localization <- function(df_sites, df_hits_grange, intermediate_dir, spin_state_bed_path) {
  if (!file.exists(spin_state_bed_path)) {
    stop("Nuclear localization (SPIN state) file not found: ", spin_state_bed_path)
  }

  df_nuclear_localization <- read.table(spin_state_bed_path)
  colnames(df_nuclear_localization) <- c("seqname", "start", "end", "description")
  df_nuclear_localization <- df_nuclear_localization %>%
    GenomicRanges::makeGRangesFromDataFrame(ignore.strand = TRUE, keep.extra.columns = TRUE)

  overlap_results <- InvestigateIntegrations::extract_overlaps(query = df_hits_grange, subject = df_nuclear_localization)

  overlap_results <- overlap_results %>%
    dplyr::select(unique_integration_site_id_query, description_subject) %>%
    dplyr::rename(unique_integration_site_id = unique_integration_site_id_query, description = description_subject)

  missing_localization <- df_sites %>%
    dplyr::filter(!unique_integration_site_id %in% overlap_results$unique_integration_site_id) %>%
    dplyr::mutate(description = "missing_nuclear_localization") %>%
    dplyr::select(unique_integration_site_id, description)

  overlap_results_summarized <- rbind(overlap_results, missing_localization) %>%
    dplyr::mutate(
      spin_state = description,
      in_interior_active_spin_state = ifelse(spin_state %in% c("Speckle", "Interior_Act1", "Interior_Act2", "Interior_Act3"), TRUE, FALSE)
    ) %>%
    dplyr::mutate(count = 1) %>%
    tidyr::pivot_wider(names_from = "spin_state", values_from = "count") %>%
    dplyr::select(-c(description)) %>%
    base::replace(is.na(.), FALSE)

  testthat::expect_equal(nrow(df_sites), nrow(overlap_results_summarized))

  out <- file.path(intermediate_dir, "completed_nuclear_localization_feature.rds")
  saveRDS(overlap_results_summarized, out)
  invisible(out)
}

annotate_repeats <- function(df_sites, intermediate_dir) {
  # NOTE: this used to call a project-specific relabel_repeat_categories()
  # helper before checking for NA. That helper isn't available in this
  # standalone script, so this treats any blank/NA repeat_class as "not in
  # a repeat". If richer repeat-class relabeling is needed, add it here.
  df_sites_repeats <- df_sites %>%
    dplyr::select(unique_integration_site_id, repeat_class) %>%
    dplyr::mutate(in_repeat = ifelse(is.na(repeat_class) | repeat_class == "", FALSE, TRUE)) %>%
    dplyr::select(-repeat_class)

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_repeats))

  out <- file.path(intermediate_dir, "completed_in_repeat_class_feature.rds")
  saveRDS(df_sites_repeats, out)
  invisible(out)
}

annotate_genomic_features <- function(df_sites, df_hits_grange, intermediate_dir, genomic_tracks_path) {
  feature_files <- c(
    "cpgIslandExt.rds",
    "encode_ccre_CTCF.rds",
    "encode_ccre_enhD.rds",
    "encode_ccre_K4m3.rds",
    "encode_ccre_enhP.rds",
    "encode_ccre_prom.rds",
    "ncbiRefSeqCurated_transcription_start_reduced.rds",
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

  df_hits_grange_larger_window <- df_hits_grange %>%
    as.data.frame() %>%
    dplyr::select(-position, -strand, -width) %>%
    dplyr::mutate(start = start - 50, end = end + 50) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

  df_all_features <- data.frame()
  for (f in feature_files) {
    feature_path <- file.path(genomic_tracks_path, f)
    if (!file.exists(feature_path)) {
      warning("Genomic feature track not found, skipping: ", feature_path, call. = FALSE)
      next
    }

    feature_ranges <- readRDS(feature_path)

    overlap_results <- InvestigateIntegrations::extract_overlaps(query = df_hits_grange_larger_window, subject = feature_ranges) %>%
      dplyr::rename(unique_integration_site_id = unique_integration_site_id_query) %>%
      dplyr::select(unique_integration_site_id) %>%
      dplyr::mutate(feature_name = f, in_feature = TRUE) %>%
      base::unique()

    df_not_in_overlap <- df_sites %>%
      dplyr::select(unique_integration_site_id) %>%
      dplyr::filter(!unique_integration_site_id %in% overlap_results$unique_integration_site_id) %>%
      dplyr::mutate(feature_name = f, in_feature = FALSE)

    overlap_results <- rbind(overlap_results, df_not_in_overlap)

    testthat::expect_equal(nrow(df_sites), nrow(overlap_results))

    df_all_features <- rbind(df_all_features, overlap_results)
  }

  df_all_features_summarized <- df_all_features %>%
    dplyr::mutate(feature_name = stringr::str_replace(feature_name, ".rds", "")) %>%
    tidyr::pivot_wider(names_from = "feature_name", values_from = "in_feature") %>%
    base::replace(is.na(.), FALSE)

  testthat::expect_equal(nrow(df_sites), nrow(df_all_features_summarized))

  out <- file.path(intermediate_dir, "completed_various_genomic_features_feature.rds")
  saveRDS(df_all_features_summarized, out)
  invisible(out)
}

annotate_pathways <- function(df_sites, intermediate_dir, reference_cache_dir) {
  reactome_cache <- file.path(reference_cache_dir, "modified_NCBI2Reactome_PE_All_Levels.rds")

  if (!file.exists(reactome_cache)) {
    dir.create(reference_cache_dir, showWarnings = FALSE, recursive = TRUE)
    raw_txt <- file.path(reference_cache_dir, "NCBI2Reactome_PE_All_Levels.txt")
    download.file("https://reactome.org/download/current/NCBI2Reactome_PE_All_Levels.txt", raw_txt)

    df_raw_table <- read.csv(file = raw_txt, header = FALSE, sep = "\t")
    colnames(df_raw_table) <- c("V1_remove", "pathway", "gene_symbol_and_localization", "connecting_pathway",
                                 "internet_link", "annotation", "unknown", "species")

    df_raw_table <- df_raw_table %>%
      dplyr::select(-c(V1_remove, pathway, connecting_pathway, internet_link, unknown)) %>%
      dplyr::filter(species == "Homo sapiens", !stringr::str_detect(gene_symbol_and_localization, "\\(")) %>%
      dplyr::mutate(
        localization = sapply(gene_symbol_and_localization, function(x) stringr::str_split(x, "\\[")[[1]][2]),
        symbol = sapply(gene_symbol_and_localization, function(x) stringr::str_split(x, "\\[")[[1]][1]),
        localization = stringr::str_replace(localization, "\\]", "")
      )

    df_raw_table_annotation <- df_raw_table %>%
      dplyr::select(-c(gene_symbol_and_localization, species, localization)) %>%
      base::unique() %>%
      dplyr::mutate(count = 1, annotation = paste0("uniprot_pathway_", annotation)) %>%
      tidyr::pivot_wider(names_from = "annotation", values_from = "count")

    df_raw_table_localization <- df_raw_table %>%
      dplyr::select(-c(gene_symbol_and_localization, species, annotation)) %>%
      base::unique() %>%
      dplyr::mutate(count = 1, localization = paste0("uniprot_localization_", localization)) %>%
      tidyr::pivot_wider(names_from = "localization", values_from = "count")

    df_processed_table <- merge(df_raw_table_annotation, df_raw_table_localization, by = "symbol")

    df_processed_table <- df_processed_table %>%
      dplyr::filter(
        !stringr::str_detect(symbol, "\\*"),
        !stringr::str_detect(symbol, "\\;"),
        !stringr::str_detect(symbol, "\\'"),
        !stringr::str_detect(symbol, "\\_")
      ) %>%
      dplyr::mutate(symbol = stringr::str_replace_all(symbol, " ", ""))

    colnames(df_processed_table) <- stringr::str_replace_all(colnames(df_processed_table), " ", "_")

    saveRDS(df_processed_table, reactome_cache)
  }

  df_sites_annotated <- merge(
    df_sites %>% dplyr::select(unique_integration_site_id, nearestGene),
    readRDS(reactome_cache),
    by.x = "nearestGene", by.y = "symbol", all.x = TRUE
  ) %>%
    dplyr::select(-nearestGene) %>%
    base::replace(is.na(.), FALSE)

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_annotated))

  out <- file.path(intermediate_dir, "completed_reactome_pathway_and_localization_feature.rds")
  saveRDS(df_sites_annotated, out)
  invisible(out)
}

annotate_onco_gene_list <- function(df_sites, intermediate_dir, onco_gene_list_path) {
  if (!file.exists(onco_gene_list_path)) {
    stop("Oncogene list not found: ", onco_gene_list_path)
  }
  onco_gene_list <- read.csv(onco_gene_list_path, sep = "\t")

  df_sites_onco <- df_sites %>%
    dplyr::select(unique_integration_site_id, nearestGene) %>%
    tidyr::separate_rows(nearestGene, sep = ",", convert = TRUE) %>%
    dplyr::mutate(
      nearestGene = stringr::str_replace_all(nearestGene, " ", ""),
      in_bushmanOncoGeneList = ifelse(nearestGene %in% onco_gene_list, TRUE, FALSE)
    ) %>%
    dplyr::group_by(unique_integration_site_id) %>%
    dplyr::mutate(in_bushmanOncoGeneList = max(in_bushmanOncoGeneList)) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_integration_site_id, in_bushmanOncoGeneList) %>%
    base::unique()

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_onco))

  out <- file.path(intermediate_dir, "completed_in_bushman_onco_gene_list_feature.rds")
  saveRDS(df_sites_onco, out)
  invisible(out)
}

annotate_znf_genes <- function(df_sites, intermediate_dir) {
  df_sites_znf <- df_sites %>%
    dplyr::select(unique_integration_site_id, nearestGene) %>%
    tidyr::separate_rows(nearestGene, sep = ",", convert = TRUE) %>%
    dplyr::mutate(
      nearestGene = stringr::str_replace_all(nearestGene, " ", ""),
      in_znf_gene = stringr::str_starts(nearestGene, "ZNF")
    ) %>%
    dplyr::group_by(unique_integration_site_id) %>%
    dplyr::mutate(in_znf_gene = max(in_znf_gene)) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_integration_site_id, in_znf_gene) %>%
    base::unique()

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_znf))

  out <- file.path(intermediate_dir, "completed_in_znf_genes_feature.rds")
  saveRDS(df_sites_znf, out)
  invisible(out)
}

annotate_krab_znf_genes <- function(df_sites, intermediate_dir, krab_znf_csv_path) {
  if (!file.exists(krab_znf_csv_path)) {
    stop("KRAB-ZNF gene list not found: ", krab_znf_csv_path)
  }
  df_krab_znf <- read.csv(krab_znf_csv_path)

  df_sites_znf <- df_sites %>%
    dplyr::select(unique_integration_site_id, nearestGene) %>%
    tidyr::separate_rows(nearestGene, sep = ",", convert = TRUE) %>%
    dplyr::mutate(in_krab_znf_gene = sapply(nearestGene, function(x) x %in% df_krab_znf$gene)) %>%
    dplyr::group_by(unique_integration_site_id) %>%
    dplyr::mutate(in_krab_znf_gene = max(in_krab_znf_gene)) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_integration_site_id, in_krab_znf_gene) %>%
    base::unique()

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_znf))

  out <- file.path(intermediate_dir, "completed_in_krab_znf_genes_feature.rds")
  saveRDS(df_sites_znf, out)
  invisible(out)
}

annotate_cosmic <- function(df_sites, intermediate_dir, cosmic_tsv_gz_path) {
  if (!file.exists(cosmic_tsv_gz_path)) {
    stop("COSMIC cancer gene census file not found: ", cosmic_tsv_gz_path)
  }
  df_cosmic <- vroom::vroom(cosmic_tsv_gz_path)

  df_sites_onco <- df_sites %>%
    dplyr::select(unique_integration_site_id, nearestGene) %>%
    tidyr::separate_rows(nearestGene, sep = ",", convert = TRUE) %>%
    dplyr::mutate(
      nearestGene = stringr::str_replace_all(nearestGene, " ", ""),
      in_cosmic_cancergenelist = ifelse(nearestGene %in% df_cosmic$GENE_SYMBOL, TRUE, FALSE)
    ) %>%
    dplyr::group_by(unique_integration_site_id) %>%
    dplyr::mutate(in_cosmic_cancergenelist = max(in_cosmic_cancergenelist)) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_integration_site_id, in_cosmic_cancergenelist) %>%
    base::unique()

  testthat::expect_equal(nrow(df_sites), nrow(df_sites_onco))

  out <- file.path(intermediate_dir, "completed_in_cosmic_onco_gene_list_feature.rds")
  saveRDS(df_sites_onco, out)
  invisible(out)
}

# Merge every completed_*_feature.rds in intermediate_dir onto the core
# AAVengeR columns of df_sites, and sanitize column names for downstream use.
combine_annotations <- function(df_sites, intermediate_dir) {
  list_of_feature_files <- Sys.glob(file.path(intermediate_dir, "completed*feature.rds"))

  df_all <- df_sites %>% dplyr::select(unique_integration_site_id)

  for (f in list_of_feature_files) {
    df_all <- merge(df_all, readRDS(f), by = "unique_integration_site_id", all = TRUE)
  }

  core_cols <- intersect(
    c("unique_integration_site_id", "trial", "subject", "sample", "posid", "position", "chromosome",
      "strand", "inGene", "inExon", "repeat_name", "repeat_class", "nearestGeneStrand",
      "sonicLengths", "reads", "nRepsObs", "nearestGene", "nearestGeneDist"),
    colnames(df_sites)
  )

  df_all <- merge(df_sites %>% dplyr::select(dplyr::all_of(core_cols)), df_all,
                   by = "unique_integration_site_id", all = TRUE)

  message(sprintf("  combined table: %d rows x %d columns", nrow(df_all), ncol(df_all)))

  # sanitize column names
  bad_chars <- c(",", ";", "\\+", "\\-", "\\.", "\\/", "\\:", "\\(", "\\)", "'", '"', " ")
  for (ch in bad_chars) {
    colnames(df_all) <- stringr::str_replace_all(colnames(df_all), ch, "_")
  }
  colnames(df_all) <- stringr::str_replace_all(colnames(df_all), "__", "_")

  testthat::expect_equal(nrow(df_sites), nrow(df_all))

  df_all
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# REFERENCE DATA PATHS
#
# Standalone, hardcoded absolute paths — edit these in place if the data
# ever moves. Nothing here depends on a project_dir/build_project_paths()
# layout.
# ---------------------------------------------------------------------------

GENOMIC_TRACKS_PATH   <- "/data/GenomicTrackRepository/data/processed/hg38"
RNA_SEQ_ABUNDANCE_PATH <- "/data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/interim/rna_seq/salmon_to_deseq_tpm_output/abundance.csv"
SPIN_STATE_BED_PATH   <- "/data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/external/K562_SPIN_25kb_hg38.bed"
KRAB_ZNF_CSV_PATH     <- "/data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/external/krab_znf_S_HUNTLEY_SupplS6.csv"
COSMIC_TSV_GZ_PATH    <- "/data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/external/Cosmic_CancerGeneCensus_v100_GRCh38.tsv.gz"
ONCO_GENE_LIST_PATH   <- "/data/AAVengeR/data/geneLists/COSMIC_v10_genes"

main <- function() {
  parser <- build_option_parser()
  opts <- parse_args(parser)

  if (is.null(opts$`sites-file`)) {
    stop("--sites-file is required.", call. = FALSE)
  }
  sites_file <- normalizePath(opts$`sites-file`, mustWork = TRUE)

  if (is.null(opts$`output-dir`)) {
    stop("--output-dir is required.", call. = FALSE)
  }
  output_base <- opts$`output-dir`

  run_name <- opts$`run-name` %||% tools::file_path_sans_ext(basename(sites_file))
  reference_cache_dir <- opts$`reference-cache-dir` %||% file.path(output_base, "_reference_cache")

  run_output_dir <- file.path(output_base, run_name)
  intermediate_dir <- file.path(run_output_dir, "intermediate")
  dir.create(intermediate_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(reference_cache_dir, showWarnings = FALSE, recursive = TRUE)

  message("Run name:             ", run_name)
  message("Sites file:           ", sites_file)
  message("Output directory:     ", run_output_dir)
  message("Reference cache:      ", reference_cache_dir)

  # --- load this run's sites table --------------------------------------
  df_sites <- load_sites_table(sites_file)

  df_hits_grange <- df_sites %>%
    dplyr::select(unique_integration_site_id, chromosome, position, strand) %>%
    dplyr::mutate(start = position, end = position) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  # --- reference data paths (see constants above) -----------------------
  genomic_tracks_path <- GENOMIC_TRACKS_PATH
  rna_seq_abundance_path <- RNA_SEQ_ABUNDANCE_PATH
  spin_state_bed_path <- SPIN_STATE_BED_PATH
  krab_znf_csv_path <- KRAB_ZNF_CSV_PATH
  cosmic_tsv_gz_path <- COSMIC_TSV_GZ_PATH
  onco_gene_list_path <- ONCO_GENE_LIST_PATH

  keep_going <- isTRUE(opts$`keep-going`)

  # --- run every annotation module --------------------------------------
  run_step("GC percentage per insertion",
           function() annotate_gc_content(df_sites, intermediate_dir, genomic_tracks_path), keep_going)

  run_step("Gene location (intron/exon/UTR)",
           function() annotate_gene_location(df_sites, df_hits_grange, intermediate_dir, reference_cache_dir), keep_going)

  run_step("Orientation vs. nearest gene",
           function() annotate_orientation(df_sites, intermediate_dir), keep_going)

  run_step("Inferred expression levels",
           function() annotate_expression(df_sites, intermediate_dir, rna_seq_abundance_path), keep_going)

  run_step("Nuclear localization (SPIN state)",
           function() annotate_nuclear_localization(df_sites, df_hits_grange, intermediate_dir, spin_state_bed_path), keep_going)

  run_step("In repeats",
           function() annotate_repeats(df_sites, intermediate_dir), keep_going)

  run_step("Genomic features (ENCODE/CpG/etc.)",
           function() annotate_genomic_features(df_sites, df_hits_grange, intermediate_dir, genomic_tracks_path), keep_going)

  run_step("Reactome pathway/localization",
           function() annotate_pathways(df_sites, intermediate_dir, reference_cache_dir), keep_going)

  run_step("Bushman oncogene list",
           function() annotate_onco_gene_list(df_sites, intermediate_dir, onco_gene_list_path), keep_going)

  run_step("ZNF genes",
           function() annotate_znf_genes(df_sites, intermediate_dir), keep_going)

  run_step("KRAB-ZNF genes",
           function() annotate_krab_znf_genes(df_sites, intermediate_dir, krab_znf_csv_path), keep_going)

  run_step("COSMIC cancer gene census",
           function() annotate_cosmic(df_sites, intermediate_dir, cosmic_tsv_gz_path), keep_going)

  # --- combine everything -------------------------------------------------
  message("Combining annotations ...")
  df_all <- combine_annotations(df_sites, intermediate_dir)

  saveRDS(df_sites, file.path(run_output_dir, "sites.rds"))
  saveRDS(df_all, file.path(run_output_dir, "annotated_features.rds"))

  if (isTRUE(opts$`write-csv`)) {
    write.csv(df_sites, file.path(run_output_dir, "sites.csv"), row.names = FALSE)
    write.csv(df_all, file.path(run_output_dir, "annotated_features.csv"), row.names = FALSE)
  }

  if (isTRUE(opts$clean)) {
    message("Removing intermediate directory: ", intermediate_dir)
    unlink(intermediate_dir, recursive = TRUE)
  }

  message("Done. Final tables written to: ", run_output_dir)
}

if (sys.nframe() == 0) {
  main()
}