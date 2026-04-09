feature_theme_map <- tibble::tribble(
  ~feature_description, ~theme, ~evidence

  "CTCF", "Chromatin structure",
  "CTCF_activated_t_cell", "Chromatin structure",
  "RAD21", "Chromatin structure",

  "DNase-seq_activated_t_cell", "Chromatin accessibility",
  "dnaseI", "Chromatin accessibility",

  "H2AFZ", "Histone variant",

  "H3K27ac", "Active regulatory",
  "H3K27ac_activated_t_cell", "Active regulatory",

  "H3K27me3", "Repressive",
  "H3K27me3_activated_t_cell", "Repressive",

  "H3K36me3", "Transcription",
  "H3K36me3_activated_t_cell", "Transcription",

  "H3K4me1", "Enhancer",
  "H3K4me1_activated_t_cell", "Enhancer",

  "H3K4me2", "Enhancer",

  "H3K4me3", "Promoter",
  "H3K4me3_activated_t_cell", "Promoter",

  "H3K79me2", "Transcription",

  "H3K9ac", "Active regulatory",

  "H3K9me3", "Repressive",
  "H3K9me3_activated_t_cell", "Repressive",

  "H4K20me1", "Histone variant",

  "DNA", "Sequence/Repeats",
  "satellite", "Sequence/Repeats",
  "gc", "Sequence/Repeats",

  "transcription", "Genome annotation",
  "transcription-start", "Genome annotation"
)