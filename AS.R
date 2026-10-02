#AS 
# ============================================================
# STEP 6 — Download ALL OV exon-skipping events
# No patient filtering yet
# ============================================================

library(httr)

download_url <- "https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload"

res_test <- POST(
  download_url,
  body = list(
    tissue     = "OV",
    genes      = "",
    samples    = "",
    splicetype = "ES",
    pctwithval = 0,
    avgexppct  = 0,
    psirange   = 0,
    psistd     = 0
  ),
  encode = "form",
  timeout(300)
)

cat("HTTP:", status_code(res_test), "\n")
cat("Type:", headers(res_test)[["content-type"]], "\n")
cat(
  "Size:",
  length(content(res_test, as = "raw")),
  "bytes\n"
)

# Save ZIP
zip_test <- tempfile(fileext = ".zip")

writeBin(
  content(res_test, as = "raw"),
  zip_test
)

# What's inside?
print(
  unzip(zip_test, list = TRUE)
)
# ============================================================
# STEP 7 — Extract and inspect OV ES PSI data
# ============================================================

# Extract ZIP
es_dir <- tempfile()
dir.create(es_dir)

unzip(
  zip_test,
  exdir = es_dir
)

es_file <- list.files(
  es_dir,
  pattern = "PSI_download_OV",
  full.names = TRUE
)

cat("File:", es_file, "\n")
cat(
  "Size:",
  round(file.info(es_file)$size / 1024^2, 2),
  "MB\n"
)


# ------------------------------------------------------------
# Read header + first few rows only
# ------------------------------------------------------------

library(data.table)

es_preview <- fread(
  es_file,
  nrows = 5
)

cat("\nDimensions of preview:\n")
print(dim(es_preview))

cat("\nColumn names:\n")
print(colnames(es_preview))

cat("\nFirst rows:\n")
print(es_preview)
# ============================================================
# STEP 8 — Download OV ES + clinical/sample information
# ============================================================

library(httr)

res_es_clinical <- POST(
  download_url,
  
  body = list(
    tissue       = "OV",
    genes        = "",
    samples      = "",
    splicetype   = "ES",
    pctwithval   = 0,
    avgexppct    = 0,
    psirange     = 0,
    psistd       = 0,
    annotation   = "true",
    clinicalData = "true"
  ),
  
  encode = "form",
  timeout(300)
)

cat("HTTP:", status_code(res_es_clinical), "\n")
cat(
  "Size:",
  round(length(content(res_es_clinical, as = "raw")) / 1024^2, 2),
  "MB\n"
)

# Save
zip_es_clinical <- tempfile(fileext = ".zip")

writeBin(
  content(res_es_clinical, as = "raw"),
  zip_es_clinical
)

# Show ZIP contents
zip_contents <- unzip(
  zip_es_clinical,
  list = TRUE
)

print(zip_contents)
# ============================================================
# STEP 9 — Inspect raw beginning of PSI + clinical file
# ============================================================

inspect_dir <- tempfile()
dir.create(inspect_dir)

unzip(
  zip_es_clinical,
  exdir = inspect_dir
)

psi_clinical_file <- list.files(
  inspect_dir,
  pattern = "PSI_download_OV",
  full.names = TRUE
)

cat("File:", psi_clinical_file, "\n")
cat(
  "Size:",
  round(file.info(psi_clinical_file)$size / 1024^2, 2),
  "MB\n\n"
)

# IMPORTANT: read raw lines, not fread()
raw_lines <- readLines(
  psi_clinical_file,
  n = 15,
  warn = FALSE
)

for (i in seq_along(raw_lines)) {
  cat("\n========== LINE", i, "==========\n")
  cat(substr(raw_lines[i], 1, 3000))
  cat("\n")
}
# ============================================================
# STEP 10 — Match SpliceSeq patients with our 277-state cohort
# ============================================================

library(data.table)
library(dplyr)

# ------------------------------------------------------------
# 1. Read ONLY header
# ------------------------------------------------------------

header <- fread(
  psi_clinical_file,
  nrows = 0
)

all_cols <- colnames(header)

cat("Total columns:", length(all_cols), "\n")


# ------------------------------------------------------------
# 2. Identify TCGA patient columns
# ------------------------------------------------------------

splice_samples <- all_cols[
  grepl("^TCGA_", all_cols)
]

cat("SpliceSeq OV patients:", length(splice_samples), "\n")


# ------------------------------------------------------------
# 3. Convert SpliceSeq IDs
#
# TCGA_13_0913
# ->
# TCGA-13-0913
# ------------------------------------------------------------

splice_patient_ids <- gsub(
  "_",
  "-",
  splice_samples
)


# ------------------------------------------------------------
# 4. Our 277 molecular-state patients
# ------------------------------------------------------------

our_ids <- robust_state_df$submitter_id

cat("Our molecular-state patients:", length(our_ids), "\n")


# ------------------------------------------------------------
# 5. Overlap
# ------------------------------------------------------------

overlap_ids <- intersect(
  our_ids,
  splice_patient_ids
)

cat("\n====================================\n")
cat("OVERLAP:", length(overlap_ids), "/", length(our_ids), "\n")
cat("====================================\n")


# ------------------------------------------------------------
# 6. State distribution in overlapping patients
# ------------------------------------------------------------

overlap_states <- robust_state_df %>%
  filter(submitter_id %in% overlap_ids)

print(
  table(overlap_states$Robust_State)
)


# ------------------------------------------------------------
# 7. Missing patients
# ------------------------------------------------------------

missing_ids <- setdiff(
  our_ids,
  splice_patient_ids
)

cat("\nPatients without SpliceSeq data:",
    length(missing_ids), "\n")

if(length(missing_ids) > 0) {
  print(missing_ids)
}
# ============================================================
# STEP 10 FIX — Read sample IDs directly from raw header
# ============================================================

# Read ONLY first raw line
header_line <- readLines(
  psi_clinical_file,
  n = 1,
  warn = FALSE
)

# Split by TAB
header_fields <- strsplit(
  header_line,
  "\t",
  fixed = TRUE
)[[1]]

cat("Total header fields:", length(header_fields), "\n")


# ------------------------------------------------------------
# Identify TCGA columns directly
# ------------------------------------------------------------

splice_samples_raw <- header_fields[
  grepl("^TCGA[_-]", header_fields)
]

cat("SpliceSeq TCGA samples:", length(splice_samples_raw), "\n")

cat("\nFirst 10 raw IDs:\n")
print(head(splice_samples_raw, 10))


# ------------------------------------------------------------
# Standardize IDs
# TCGA_04_1331 -> TCGA-04-1331
# ------------------------------------------------------------

splice_patient_ids <- gsub(
  "_",
  "-",
  splice_samples_raw,
  fixed = TRUE
)

cat("\nFirst standardized IDs:\n")
print(head(splice_patient_ids, 10))


# ------------------------------------------------------------
# Standardize our IDs too, just in case
# ------------------------------------------------------------

our_ids <- toupper(
  trimws(
    as.character(robust_state_df$submitter_id)
  )
)

splice_patient_ids <- toupper(
  trimws(splice_patient_ids)
)


# ------------------------------------------------------------
# Calculate overlap
# ------------------------------------------------------------

overlap_ids <- intersect(
  our_ids,
  splice_patient_ids
)

cat("\n====================================\n")
cat("OUR PATIENTS:", length(unique(our_ids)), "\n")
cat("SPLICESEQ PATIENTS:", length(unique(splice_patient_ids)), "\n")
cat("OVERLAP:", length(overlap_ids), "\n")
cat("====================================\n\n")

cat("First overlapping patients:\n")
print(head(overlap_ids, 20))


# ------------------------------------------------------------
# State distribution
# ------------------------------------------------------------

overlap_states <- robust_state_df[
  robust_state_df$submitter_id %in% overlap_ids,
]

cat("\nState distribution:\n")
print(
  table(overlap_states$Robust_State)
)


# ------------------------------------------------------------
# Missing from SpliceSeq
# ------------------------------------------------------------

missing_ids <- setdiff(
  our_ids,
  splice_patient_ids
)

cat(
  "\nMissing from SpliceSeq:",
  length(missing_ids),
  "\n"
)
# ============================================================
# STEP 11 — Differential Exon Skipping
# State1 vs State2
# ============================================================

library(data.table)
library(dplyr)

# ------------------------------------------------------------
# 1. Read PSI table
# ------------------------------------------------------------

es <- fread(
  psi_clinical_file,
  skip = 14,
  header = FALSE,
  na.strings = c("null", "NA", "")
)

cat("Raw ES dimensions:", dim(es), "\n")


# ------------------------------------------------------------
# 2. Assign correct column names from original header
# ------------------------------------------------------------

colnames(es) <- header_fields

cat("ES dimensions after naming:", dim(es), "\n")


# ------------------------------------------------------------
# 3. Keep metadata + our overlapping patients
# ------------------------------------------------------------

meta_cols <- c(
  "symbol",
  "as_id",
  "splice_type",
  "exons",
  "from_exon",
  "to_exon",
  "novel_splice",
  "pct_with_values",
  "psi_range",
  "std_psi"
)

# Convert TCGA IDs back to SpliceSeq format
overlap_splice_names <- gsub(
  "-",
  "_",
  overlap_ids,
  fixed = TRUE
)

es_sub <- es[
  ,
  c(meta_cols, overlap_splice_names),
  with = FALSE
]

cat("Patients retained:", length(overlap_splice_names), "\n")
cat("ES events:", nrow(es_sub), "\n")


# ------------------------------------------------------------
# 4. State membership
# ------------------------------------------------------------

state_map <- robust_state_df %>%
  filter(submitter_id %in% overlap_ids) %>%
  select(submitter_id, Robust_State)

state1_ids <- state_map$submitter_id[
  state_map$Robust_State == 1
]

state2_ids <- state_map$submitter_id[
  state_map$Robust_State == 2
]

state1_cols <- gsub("-", "_", state1_ids, fixed = TRUE)
state2_cols <- gsub("-", "_", state2_ids, fixed = TRUE)

cat("\nState1:", length(state1_cols), "\n")
cat("State2:", length(state2_cols), "\n")


# ------------------------------------------------------------
# 5. Convert PSI columns to numeric
# ------------------------------------------------------------

patient_cols <- c(state1_cols, state2_cols)

es_sub[
  ,
  (patient_cols) := lapply(.SD, as.numeric),
  .SDcols = patient_cols
]


# ------------------------------------------------------------
# 6. Differential PSI
# ------------------------------------------------------------

results <- vector(
  "list",
  nrow(es_sub)
)

for (i in seq_len(nrow(es_sub))) {
  
  x1 <- as.numeric(
    unlist(es_sub[i, ..state1_cols])
  )
  
  x2 <- as.numeric(
    unlist(es_sub[i, ..state2_cols])
  )
  
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  med1 <- if(n1 > 0) median(x1) else NA_real_
  med2 <- if(n2 > 0) median(x2) else NA_real_
  
  delta <- med2 - med1
  
  # Require enough observations in both states
  if(n1 >= 10 && n2 >= 10) {
    
    p <- tryCatch(
      wilcox.test(
        x2,
        x1,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    
  } else {
    
    p <- NA_real_
    
  }
  
  results[[i]] <- data.frame(
    N_State1   = n1,
    N_State2   = n2,
    Median_State1 = med1,
    Median_State2 = med2,
    Delta_PSI  = delta,
    P_value    = p
  )
}


# ------------------------------------------------------------
# 7. Combine results
# ------------------------------------------------------------

stats <- bind_rows(results)

es_results <- bind_cols(
  as.data.frame(es_sub[, ..meta_cols]),
  stats
)

es_results$FDR <- p.adjust(
  es_results$P_value,
  method = "BH"
)


# ------------------------------------------------------------
# 8. Direction
# ------------------------------------------------------------

es_results <- es_results %>%
  mutate(
    Direction = case_when(
      Delta_PSI > 0 ~ "State2_higher_PSI",
      Delta_PSI < 0 ~ "State1_higher_PSI",
      TRUE ~ "No_difference"
    )
  ) %>%
  arrange(FDR, desc(abs(Delta_PSI)))


# ------------------------------------------------------------
# 9. Significant differential AS events
# ------------------------------------------------------------

es_significant <- es_results %>%
  filter(
    FDR < 0.05,
    abs(Delta_PSI) >= 0.10
  )


# ------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DIFFERENTIAL EXON SKIPPING\n")
cat("========================================\n")

cat("Total ES events:", nrow(es_results), "\n")

cat(
  "Tested events:",
  sum(!is.na(es_results$P_value)),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(es_results$FDR < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "FDR < 0.05 & |Delta PSI| >= 0.10:",
  nrow(es_significant),
  "\n"
)

cat("\nDirection:\n")

print(
  table(es_significant$Direction)
)


# ------------------------------------------------------------
# 11. Top events
# ------------------------------------------------------------

print(
  es_significant %>%
    select(
      symbol,
      as_id,
      exons,
      from_exon,
      to_exon,
      N_State1,
      N_State2,
      Median_State1,
      Median_State2,
      Delta_PSI,
      FDR,
      Direction
    ) %>%
    head(30)
)
# ------------------------------------------------------------
# 4. State membership
# ------------------------------------------------------------

state_map <- robust_state_df %>%
  dplyr::filter(submitter_id %in% overlap_ids) %>%
  dplyr::select(submitter_id, Robust_State)

state1_ids <- state_map$submitter_id[
  state_map$Robust_State == 1
]

state2_ids <- state_map$submitter_id[
  state_map$Robust_State == 2
]

state1_cols <- gsub("-", "_", state1_ids, fixed = TRUE)
state2_cols <- gsub("-", "_", state2_ids, fixed = TRUE)

cat("\nState1:", length(state1_cols), "\n")
cat("State2:", length(state2_cols), "\n")
# ============================================================
# STEP 11B — Differential Exon Skipping
# ============================================================

library(data.table)
library(dplyr)

patient_cols <- c(state1_cols, state2_cols)

# PSI columns numeric
es_sub[
  ,
  (patient_cols) := lapply(.SD, as.numeric),
  .SDcols = patient_cols
]

# ------------------------------------------------------------
# Differential PSI
# ------------------------------------------------------------

results <- vector("list", nrow(es_sub))

for (i in seq_len(nrow(es_sub))) {
  
  x1 <- as.numeric(
    unlist(es_sub[i, ..state1_cols])
  )
  
  x2 <- as.numeric(
    unlist(es_sub[i, ..state2_cols])
  )
  
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  med1 <- if (n1 > 0) median(x1) else NA_real_
  med2 <- if (n2 > 0) median(x2) else NA_real_
  
  delta <- med2 - med1
  
  if (n1 >= 10 && n2 >= 10) {
    
    p <- tryCatch(
      wilcox.test(
        x2,
        x1,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    
  } else {
    
    p <- NA_real_
  }
  
  results[[i]] <- data.frame(
    N_State1 = n1,
    N_State2 = n2,
    Median_State1 = med1,
    Median_State2 = med2,
    Delta_PSI = delta,
    P_value = p
  )
}

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

stats <- dplyr::bind_rows(results)

es_results <- dplyr::bind_cols(
  as.data.frame(es_sub[, ..meta_cols]),
  stats
)

es_results$FDR <- p.adjust(
  es_results$P_value,
  method = "BH"
)

es_results <- es_results %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      Delta_PSI > 0 ~ "State2_higher_PSI",
      Delta_PSI < 0 ~ "State1_higher_PSI",
      TRUE ~ "No_difference"
    )
  ) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(Delta_PSI)))

# ------------------------------------------------------------
# Significant events
# ------------------------------------------------------------

es_significant <- es_results %>%
  dplyr::filter(
    FDR < 0.05,
    abs(Delta_PSI) >= 0.10
  )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("STATE1 vs STATE2 — EXON SKIPPING\n")
cat("========================================\n")

cat("State1 patients:", length(state1_cols), "\n")
cat("State2 patients:", length(state2_cols), "\n")

cat("Total ES events:", nrow(es_results), "\n")

cat(
  "Tested events:",
  sum(!is.na(es_results$P_value)),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(es_results$FDR < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "FDR < 0.05 & |Delta PSI| >= 0.10:",
  nrow(es_significant),
  "\n\n"
)

print(table(es_significant$Direction))

# ------------------------------------------------------------
# Top 30
# ------------------------------------------------------------

top30_ES <- es_significant %>%
  dplyr::select(
    symbol,
    as_id,
    exons,
    from_exon,
    to_exon,
    N_State1,
    N_State2,
    Median_State1,
    Median_State2,
    Delta_PSI,
    P_value,
    FDR,
    Direction
  ) %>%
  head(30)

print(top30_ES)
cat("Total ES events:", nrow(es_results), "\n")

cat(
  "FDR < 0.05:",
  sum(es_results$FDR < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "FDR < 0.05 & |Delta PSI| >= 0.10:",
  nrow(es_significant),
  "\n"
)
# ============================================================
# STEP 12 — Download remaining OV alternative-splicing types
# AA, AD, RI, AP, AT, ME
# ============================================================

library(httr)

download_url <- "https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload"

as_types <- c("AA", "AD", "RI", "AP", "AT", "ME")

output_dir <- "~/Downloads/TCGA_OV_SpliceSeq"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

for (type in as_types) {
  
  cat("\n====================================\n")
  cat("Downloading:", type, "\n")
  cat("====================================\n")
  
  res <- POST(
    download_url,
    
    body = list(
      tissue       = "OV",
      genes        = "",
      samples      = "",
      splicetype   = type,
      pctwithval   = 0,
      avgexppct    = 0,
      psirange     = 0,
      psistd       = 0,
      clinicalData = "true"
    ),
    
    encode = "form",
    timeout(600)
  )
  
  cat("HTTP:", status_code(res), "\n")
  
  zip_file <- file.path(
    output_dir,
    paste0("OV_", type, ".zip")
  )
  
  writeBin(
    content(res, as = "raw"),
    zip_file
  )
  
  extract_dir <- file.path(
    output_dir,
    type
  )
  
  dir.create(
    extract_dir,
    showWarnings = FALSE
  )
  
  unzip(
    zip_file,
    exdir = extract_dir
  )
  
  extracted <- list.files(
    extract_dir,
    full.names = TRUE
  )
  
  cat("Extracted file:\n")
  print(extracted)
  
  if(length(extracted) > 0) {
    
    cat(
      "Size:",
      round(file.info(extracted[1])$size / 1024^2, 2),
      "MB\n"
    )
  }
}

cat("\n====================================\n")
cat("ALL DOWNLOADS COMPLETE\n")
cat("Location:", output_dir, "\n")
cat("====================================\n")
library(ggplot2)
library(dplyr)

plot_df <- es_results %>%
  mutate(
    Significant = case_when(
      FDR < 0.05 & Delta_PSI >= 0.10  ~ "State2 higher",
      FDR < 0.05 & Delta_PSI <= -0.10 ~ "State1 higher",
      TRUE ~ "Not significant"
    )
  )

ggplot(
  plot_df,
  aes(
    x = Delta_PSI,
    y = -log10(FDR),
    color = Significant
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.4
  ) +
  
  geom_vline(
    xintercept = c(-0.10, 0.10),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "State1 higher" = "#2166AC",
      "Not significant" = "grey75",
      "State2 higher" = "#B2182B"
    )
  ) +
  
  labs(
    x = expression(Delta*"PSI (State2 - State1)"),
    y = expression(-log[10]*"(FDR)"),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  )
library(ggplot2)
library(dplyr)

plot_df <- es_results %>%
  mutate(
    Category = case_when(
      
      # Statistically significant + large effect
      FDR < 0.05 & Delta_PSI <= -0.10 ~
        "State1 higher (|ΔPSI| ≥ 0.10)",
      
      FDR < 0.05 & Delta_PSI >= 0.10 ~
        "State2 higher (|ΔPSI| ≥ 0.10)",
      
      # Significant but smaller effect
      FDR < 0.05 & abs(Delta_PSI) < 0.10 ~
        "FDR < 0.05, |ΔPSI| < 0.10",
      
      # Everything else
      TRUE ~ "Not significant"
    )
  )


ggplot(
  plot_df,
  aes(
    x = Delta_PSI,
    y = -log10(FDR),
    color = Category
  )
) +
  
  geom_point(
    alpha = 0.60,
    size = 1.5
  ) +
  
  # Effect-size thresholds
  geom_vline(
    xintercept = c(-0.10, 0.10),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Statistical significance threshold
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey75",
      "FDR < 0.05, |ΔPSI| < 0.10" = "darkgrey",
      "State1 higher (|ΔPSI| ≥ 0.10)" = "#2166AC",
      "State2 higher (|ΔPSI| ≥ 0.10)" = "#B2182B"
    )
  ) +
  
  labs(
    x = expression(Delta*"PSI (State2 - State1)"),
    y = expression(-log[10]*"(FDR)"),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10)
  )
library(ggplot2)
library(dplyr)

plot_df <- es_results %>%
  mutate(
    Category = case_when(
      
      FDR < 0.05 & Delta_PSI <= -0.10 ~
        "State1 higher |ΔPSI| ≥ 0.10",
      
      FDR < 0.05 & Delta_PSI >= 0.10 ~
        "State2 higher |ΔPSI| ≥ 0.10",
      
      FDR < 0.05 ~
        "FDR < 0.05",
      
      TRUE ~
        "Not significant"
    )
  )

ggplot(
  plot_df,
  aes(
    x = Delta_PSI,
    y = -log10(FDR),
    color = Category
  )
) +
  
  geom_point(
    alpha = 0.60,
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = c(-0.10, 0.10),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "Not significant" = "grey80",
      "FDR < 0.05" = "#E69F00",
      "State1 higher |ΔPSI| ≥ 0.10" = "#2166AC",
      "State2 higher |ΔPSI| ≥ 0.10" = "#B2182B"
    )
  ) +
  
  labs(
    x = expression(Delta*"PSI (State2 - State1)"),
    y = expression(-log[10]*"(FDR)"),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10)
  )
# ============================================================
# STEP 13 — GLOBAL ES PROFILE
# Unsupervised PCA of exon-skipping PSI
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. PSI matrix
# rows    = ES events
# columns = 273 patients
# ------------------------------------------------------------

psi_matrix <- as.matrix(
  es_sub[, ..patient_cols]
)

storage.mode(psi_matrix) <- "numeric"

rownames(psi_matrix) <- paste0(
  es_sub$symbol,
  "_ES_",
  es_sub$as_id
)

cat("Initial ES matrix:", dim(psi_matrix), "\n")


# ------------------------------------------------------------
# 2. Coverage filter
# Require PSI in >=75% of patients
# ------------------------------------------------------------

coverage <- rowMeans(!is.na(psi_matrix))

psi_filtered <- psi_matrix[
  coverage >= 0.75,
  ,
  drop = FALSE
]

cat(
  "After >=75% coverage filter:",
  dim(psi_filtered),
  "\n"
)


# ------------------------------------------------------------
# 3. Variance of each ES event
# State information is NOT used here
# ------------------------------------------------------------

event_variance <- apply(
  psi_filtered,
  1,
  var,
  na.rm = TRUE
)

summary(event_variance)


# ------------------------------------------------------------
# 4. Keep top 10% most variable ES events
# UNSUPERVISED feature selection
# ------------------------------------------------------------

variance_cutoff <- quantile(
  event_variance,
  0.90,
  na.rm = TRUE
)

psi_variable <- psi_filtered[
  event_variance >= variance_cutoff &
    event_variance > 0,
  ,
  drop = FALSE
]

cat(
  "Top 10% variable ES events:",
  nrow(psi_variable),
  "\n"
)


# ------------------------------------------------------------
# 5. Median imputation
# Event-wise
# ------------------------------------------------------------

psi_imputed <- psi_variable

for(i in seq_len(nrow(psi_imputed))) {
  
  med <- median(
    psi_imputed[i, ],
    na.rm = TRUE
  )
  
  psi_imputed[
    i,
    is.na(psi_imputed[i, ])
  ] <- med
}

cat(
  "Remaining NA:",
  sum(is.na(psi_imputed)),
  "\n"
)


# ------------------------------------------------------------
# 6. PCA
# Patients must be rows
# ------------------------------------------------------------

pca_es <- prcomp(
  t(psi_imputed),
  center = TRUE,
  scale. = TRUE
)

variance_explained <- (
  pca_es$sdev^2 /
    sum(pca_es$sdev^2)
) * 100

cat(
  "PC1 variance:",
  round(variance_explained[1], 2),
  "%\n"
)

cat(
  "PC2 variance:",
  round(variance_explained[2], 2),
  "%\n"
)


# ------------------------------------------------------------
# 7. PCA dataframe
# ------------------------------------------------------------

pca_df <- data.frame(
  submitter_id = gsub(
    "_",
    "-",
    rownames(pca_es$x),
    fixed = TRUE
  ),
  PC1 = pca_es$x[, 1],
  PC2 = pca_es$x[, 2]
)


# ------------------------------------------------------------
# 8. Add molecular state AFTER PCA
# ------------------------------------------------------------

pca_df <- pca_df %>%
  left_join(
    robust_state_df %>%
      dplyr::select(
        submitter_id,
        Robust_State
      ),
    by = "submitter_id"
  )

pca_df$Robust_State <- factor(
  pca_df$Robust_State,
  levels = c(1, 2),
  labels = c("State1", "State2")
)

print(table(pca_df$Robust_State))


# ------------------------------------------------------------
# 9. PCA plot
# ------------------------------------------------------------

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Robust_State
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.70
  ) +
  
  stat_ellipse(
    level = 0.95,
    linewidth = 0.8
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top"
  )
# ============================================================
# STEP 14 — Global ES difference between State1 and State2
# PERMANOVA
# ============================================================

# Install only if necessary
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}

library(vegan)

# ------------------------------------------------------------
# 1. Patients x ES events
# ------------------------------------------------------------

X_es <- t(psi_imputed)

cat("Matrix:", dim(X_es), "\n")


# ------------------------------------------------------------
# 2. State vector in exactly the same patient order
# ------------------------------------------------------------

patient_ids_es <- gsub(
  "_",
  "-",
  rownames(X_es),
  fixed = TRUE
)

state_test_df <- data.frame(
  submitter_id = patient_ids_es
)

state_test_df$State <- robust_state_df$Robust_State[
  match(
    state_test_df$submitter_id,
    robust_state_df$submitter_id
  )
]

state_test_df$State <- factor(
  state_test_df$State,
  levels = c(1, 2),
  labels = c("State1", "State2")
)

print(table(state_test_df$State))

cat(
  "Missing state assignments:",
  sum(is.na(state_test_df$State)),
  "\n"
)


# ------------------------------------------------------------
# 3. Euclidean distance on standardized PSI
#
# Same conceptual scaling used for PCA:
# each ES event contributes on comparable scale
# ------------------------------------------------------------

X_es_scaled <- scale(X_es)

es_distance <- dist(
  X_es_scaled,
  method = "euclidean"
)


# ------------------------------------------------------------
# 4. PERMANOVA
# ------------------------------------------------------------

set.seed(123)

permanova_es <- vegan::adonis2(
  es_distance ~ State,
  data = state_test_df,
  permutations = 9999
)

print(permanova_es)


# ------------------------------------------------------------
# 5. Also test dispersion
# Important for interpreting PERMANOVA
# ------------------------------------------------------------

dispersion_es <- vegan::betadisper(
  es_distance,
  state_test_df$State
)

dispersion_test <- permutest(
  dispersion_es,
  permutations = 9999
)

cat("\n====================================\n")
cat("DISPERSION TEST\n")
cat("====================================\n")

print(dispersion_test)
print(permanova_es)
# ============================================================
# STEP 13 — Check downloaded AS files
# ============================================================

library(data.table)

as_types <- c("AA", "AD", "RI", "AP", "AT", "ME")

as_files <- list()

for (type in as_types) {
  
  type_dir <- file.path(
    output_dir,
    type
  )
  
  files <- list.files(
    type_dir,
    full.names = TRUE
  )
  
  if(length(files) == 0) {
    
    cat("\n", type, ": NO FILE FOUND\n")
    
  } else {
    
    as_files[[type]] <- files[1]
    
    cat("\n====================================\n")
    cat("TYPE:", type, "\n")
    cat("FILE:", basename(files[1]), "\n")
    
    cat(
      "SIZE:",
      round(file.info(files[1])$size / 1024^2, 2),
      "MB\n"
    )
    
    # Read raw header
    h <- readLines(
      files[1],
      n = 1,
      warn = FALSE
    )
    
    fields <- strsplit(
      h,
      "\t",
      fixed = TRUE
    )[[1]]
    
    tcga_cols <- fields[
      grepl("^TCGA[_-]", fields)
    ]
    
    cat(
      "TOTAL COLUMNS:",
      length(fields),
      "\n"
    )
    
    cat(
      "TCGA PATIENTS:",
      length(tcga_cols),
      "\n"
    )
    
    # Standardize IDs
    ids <- gsub(
      "_",
      "-",
      tcga_cols,
      fixed = TRUE
    )
    
    ov <- intersect(
      robust_state_df$submitter_id,
      ids
    )
    
    cat(
      "OVERLAP WITH STATES:",
      length(ov),
      "/ 277\n"
    )
  }
}

cat("\n====================================\n")
cat("FILE CHECK COMPLETE\n")
cat("====================================\n")

# ============================================================
# STEP 14 — Differential AS analysis
# AA, AD, RI, AP, AT, ME
# ============================================================

library(data.table)
library(dplyr)

as_types <- c("AA", "AD", "RI", "AP", "AT", "ME")

all_as_results <- list()
all_as_significant <- list()
as_summary <- list()


# ============================================================
# FUNCTION
# ============================================================

analyse_AS_type <- function(type, file) {
  
  cat("\n\n")
  cat("====================================================\n")
  cat("ANALYSING:", type, "\n")
  cat("====================================================\n")
  
  # --------------------------------------------------------
  # 1. Read header
  # --------------------------------------------------------
  
  header_line <- readLines(
    file,
    n = 1,
    warn = FALSE
  )
  
  header_fields <- strsplit(
    header_line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  # --------------------------------------------------------
  # 2. Read actual AS events
  # first 14 lines = header + clinical information
  # --------------------------------------------------------
  
  dat <- fread(
    file,
    skip = 14,
    header = FALSE,
    na.strings = c("null", "NA", "")
  )
  
  colnames(dat) <- header_fields
  
  cat("Total events:", nrow(dat), "\n")
  
  # --------------------------------------------------------
  # 3. Metadata
  # --------------------------------------------------------
  
  meta_cols <- c(
    "symbol",
    "as_id",
    "splice_type",
    "exons",
    "from_exon",
    "to_exon",
    "novel_splice",
    "pct_with_values",
    "psi_range",
    "std_psi"
  )
  
  # --------------------------------------------------------
  # 4. Patients overlapping our states
  # --------------------------------------------------------
  
  splice_cols <- header_fields[
    grepl("^TCGA[_-]", header_fields)
  ]
  
  splice_ids <- gsub(
    "_",
    "-",
    splice_cols,
    fixed = TRUE
  )
  
  overlap <- intersect(
    robust_state_df$submitter_id,
    splice_ids
  )
  
  state_map <- robust_state_df[
    robust_state_df$submitter_id %in% overlap,
    c("submitter_id", "Robust_State")
  ]
  
  state1_ids <- state_map$submitter_id[
    state_map$Robust_State == 1
  ]
  
  state2_ids <- state_map$submitter_id[
    state_map$Robust_State == 2
  ]
  
  state1_cols <- gsub(
    "-",
    "_",
    state1_ids,
    fixed = TRUE
  )
  
  state2_cols <- gsub(
    "-",
    "_",
    state2_ids,
    fixed = TRUE
  )
  
  patient_cols <- c(
    state1_cols,
    state2_cols
  )
  
  cat(
    "State1:", length(state1_cols),
    "| State2:", length(state2_cols),
    "\n"
  )
  
  # --------------------------------------------------------
  # 5. Subset
  # --------------------------------------------------------
  
  dat_sub <- dat[
    ,
    c(meta_cols, patient_cols),
    with = FALSE
  ]
  
  dat_sub[
    ,
    (patient_cols) := lapply(.SD, as.numeric),
    .SDcols = patient_cols
  ]
  
  # --------------------------------------------------------
  # 6. Coverage
  # Require >=75% observations WITHIN EACH STATE
  # --------------------------------------------------------
  
  psi1 <- as.matrix(
    dat_sub[, ..state1_cols]
  )
  
  psi2 <- as.matrix(
    dat_sub[, ..state2_cols]
  )
  
  storage.mode(psi1) <- "numeric"
  storage.mode(psi2) <- "numeric"
  
  coverage1 <- rowMeans(!is.na(psi1))
  coverage2 <- rowMeans(!is.na(psi2))
  
  keep <- (
    coverage1 >= 0.75 &
      coverage2 >= 0.75
  )
  
  dat_test <- dat_sub[keep]
  
  psi1 <- psi1[keep, , drop = FALSE]
  psi2 <- psi2[keep, , drop = FALSE]
  
  cat(
    "Events after >=75% coverage:",
    nrow(dat_test),
    "\n"
  )
  
  # --------------------------------------------------------
  # 7. Event-level statistics
  # --------------------------------------------------------
  
  N <- nrow(dat_test)
  
  n1 <- rowSums(!is.na(psi1))
  n2 <- rowSums(!is.na(psi2))
  
  median1 <- apply(
    psi1,
    1,
    median,
    na.rm = TRUE
  )
  
  median2 <- apply(
    psi2,
    1,
    median,
    na.rm = TRUE
  )
  
  delta <- median2 - median1
  
  p_values <- numeric(N)
  
  cat("Running Wilcoxon tests...\n")
  
  for(i in seq_len(N)) {
    
    x1 <- psi1[i, ]
    x2 <- psi2[i, ]
    
    x1 <- x1[is.finite(x1)]
    x2 <- x2[is.finite(x2)]
    
    p_values[i] <- tryCatch(
      
      wilcox.test(
        x2,
        x1,
        exact = FALSE
      )$p.value,
      
      error = function(e) NA_real_
    )
  }
  
  # --------------------------------------------------------
  # 8. Results
  # --------------------------------------------------------
  
  result <- data.frame(
    AS_Type = type,
    
    symbol = dat_test$symbol,
    as_id = dat_test$as_id,
    exons = dat_test$exons,
    from_exon = dat_test$from_exon,
    to_exon = dat_test$to_exon,
    
    N_State1 = n1,
    N_State2 = n2,
    
    Median_State1 = median1,
    Median_State2 = median2,
    
    Delta_PSI = delta,
    P_value = p_values
  )
  
  result$FDR <- p.adjust(
    result$P_value,
    method = "BH"
  )
  
  result$Direction <- dplyr::case_when(
    
    result$Delta_PSI > 0 ~
      "State2_higher_PSI",
    
    result$Delta_PSI < 0 ~
      "State1_higher_PSI",
    
    TRUE ~
      "No_difference"
  )
  
  result <- result[
    order(result$FDR),
  ]
  
  # --------------------------------------------------------
  # 9. Significant
  # --------------------------------------------------------
  
  sig_FDR <- result[
    !is.na(result$FDR) &
      result$FDR < 0.05,
  ]
  
  sig_strong <- result[
    !is.na(result$FDR) &
      result$FDR < 0.05 &
      abs(result$Delta_PSI) >= 0.10,
  ]
  
  cat("\nRESULTS\n")
  
  cat(
    "FDR < 0.05:",
    nrow(sig_FDR),
    "\n"
  )
  
  cat(
    "FDR < 0.05 & |Delta PSI| >= 0.10:",
    nrow(sig_strong),
    "\n"
  )
  
  if(nrow(sig_strong) > 0) {
    
    cat("\nStrong events:\n")
    
    print(
      sig_strong[
        ,
        c(
          "symbol",
          "as_id",
          "Median_State1",
          "Median_State2",
          "Delta_PSI",
          "FDR",
          "Direction"
        )
      ]
    )
  }
  
  # --------------------------------------------------------
  # 10. Summary
  # --------------------------------------------------------
  
  summary_row <- data.frame(
    
    AS_Type = type,
    
    Total_events = nrow(dat),
    
    Coverage_filtered_events =
      nrow(dat_test),
    
    FDR_005 =
      nrow(sig_FDR),
    
    Strong_DPSI_010 =
      nrow(sig_strong),
    
    State1_strong =
      sum(
        sig_strong$Delta_PSI <= -0.10,
        na.rm = TRUE
      ),
    
    State2_strong =
      sum(
        sig_strong$Delta_PSI >= 0.10,
        na.rm = TRUE
      )
  )
  
  return(
    list(
      results = result,
      significant = sig_strong,
      summary = summary_row
    )
  )
}


# ============================================================
# RUN ALL SIX AS TYPES
# ============================================================

for(type in as_types) {
  
  tmp <- analyse_AS_type(
    type = type,
    file = as_files[[type]]
  )
  
  all_as_results[[type]] <-
    tmp$results
  
  all_as_significant[[type]] <-
    tmp$significant
  
  as_summary[[type]] <-
    tmp$summary
}


# ============================================================
# COMBINE SUMMARY
# ============================================================

as_summary_df <- dplyr::bind_rows(
  as_summary
)

cat("\n\n")
cat("====================================================\n")
cat("FINAL AS SUMMARY\n")
cat("====================================================\n")

print(as_summary_df)
# ============================================================
# STEP 15 — Reanalyse ES using identical pipeline
# ============================================================

es_same_pipeline <- analyse_AS_type(
  type = "ES",
  file = psi_clinical_file
)

all_as_results[["ES"]] <- es_same_pipeline$results
all_as_significant[["ES"]] <- es_same_pipeline$significant

# Add ES to summary
as_summary_complete <- dplyr::bind_rows(
  es_same_pipeline$summary,
  as_summary_df
)

# Order AS types
as_summary_complete$AS_Type <- factor(
  as_summary_complete$AS_Type,
  levels = c("ES", "AA", "AD", "RI", "AP", "AT", "ME")
)

as_summary_complete <- as_summary_complete %>%
  dplyr::arrange(AS_Type)

print(as_summary_complete)
# ============================================================
# STEP 16 — AS event-level comparison plot
# ============================================================

library(dplyr)
library(ggplot2)

as_plot <- as_summary_complete %>%
  mutate(
    FDR_percent =
      100 * FDR_005 / Coverage_filtered_events,
    
    Strong_percent =
      100 * Strong_DPSI_010 / Coverage_filtered_events
  )

print(
  as_plot %>%
    dplyr::select(
      AS_Type,
      Coverage_filtered_events,
      FDR_005,
      FDR_percent,
      Strong_DPSI_010,
      Strong_percent
    )
)


ggplot(
  as_plot,
  aes(
    x = AS_Type,
    y = Strong_percent
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_text(
    aes(
      label = Strong_DPSI_010
    ),
    vjust = -0.4,
    size = 4
  ) +
  
  labs(
    x = "Alternative splicing type",
    y = "Strong differential events (%)"
  ) +
  
  theme_classic(base_size = 14)
# ============================================================
# STEP 17 — GLOBAL AS PROFILE ANALYSIS
# PCA-related variance + PERMANOVA + dispersion
# ES, AA, AD, RI, AP, AT, ME
# ============================================================

library(data.table)
library(dplyr)
library(vegan)

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

all_as_files <- as_files
all_as_files[["ES"]] <- psi_clinical_file

types_global <- c(
  "ES", "AA", "AD", "RI", "AP", "AT", "ME"
)

global_results <- list()


# ============================================================
# FUNCTION
# ============================================================

global_AS_test <- function(type, file) {
  
  cat("\n============================================\n")
  cat("GLOBAL ANALYSIS:", type, "\n")
  cat("============================================\n")
  
  # --------------------------------------------------------
  # 1. Header
  # --------------------------------------------------------
  
  header_line <- readLines(
    file,
    n = 1,
    warn = FALSE
  )
  
  header_fields <- strsplit(
    header_line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  # --------------------------------------------------------
  # 2. Read data
  # --------------------------------------------------------
  
  dat <- fread(
    file,
    skip = 14,
    header = FALSE,
    na.strings = c("null", "NA", "")
  )
  
  colnames(dat) <- header_fields
  
  # --------------------------------------------------------
  # 3. Match our 273 patients
  # --------------------------------------------------------
  
  splice_cols <- header_fields[
    grepl("^TCGA[_-]", header_fields)
  ]
  
  splice_ids <- gsub(
    "_",
    "-",
    splice_cols,
    fixed = TRUE
  )
  
  overlap <- intersect(
    robust_state_df$submitter_id,
    splice_ids
  )
  
  patient_cols <- gsub(
    "-",
    "_",
    overlap,
    fixed = TRUE
  )
  
  # --------------------------------------------------------
  # 4. PSI matrix
  # --------------------------------------------------------
  
  psi <- as.matrix(
    dat[, ..patient_cols]
  )
  
  storage.mode(psi) <- "numeric"
  
  cat(
    "Initial events:",
    nrow(psi),
    "\n"
  )
  
  # --------------------------------------------------------
  # 5. >=75% overall coverage
  # --------------------------------------------------------
  
  coverage <- rowMeans(
    !is.na(psi)
  )
  
  psi <- psi[
    coverage >= 0.75,
    ,
    drop = FALSE
  ]
  
  cat(
    "Coverage filtered:",
    nrow(psi),
    "\n"
  )
  
  # --------------------------------------------------------
  # 6. Variance
  # IMPORTANT:
  # state labels NOT used
  # --------------------------------------------------------
  
  event_var <- apply(
    psi,
    1,
    var,
    na.rm = TRUE
  )
  
  keep_var <- (
    is.finite(event_var) &
      event_var > 0
  )
  
  psi <- psi[
    keep_var,
    ,
    drop = FALSE
  ]
  
  event_var <- event_var[
    keep_var
  ]
  
  # --------------------------------------------------------
  # 7. Top 10% variable events
  # --------------------------------------------------------
  
  cutoff <- quantile(
    event_var,
    0.90,
    na.rm = TRUE
  )
  
  psi_var <- psi[
    event_var >= cutoff,
    ,
    drop = FALSE
  ]
  
  cat(
    "Top 10% variable:",
    nrow(psi_var),
    "\n"
  )
  
  # --------------------------------------------------------
  # 8. Median imputation
  # --------------------------------------------------------
  
  for(i in seq_len(nrow(psi_var))) {
    
    med <- median(
      psi_var[i, ],
      na.rm = TRUE
    )
    
    psi_var[
      i,
      is.na(psi_var[i, ])
    ] <- med
  }
  
  # --------------------------------------------------------
  # 9. Patients x events
  # --------------------------------------------------------
  
  X <- t(psi_var)
  
  # Standardize events
  X_scaled <- scale(X)
  
  # Remove any problematic columns
  good_cols <- apply(
    X_scaled,
    2,
    function(x) all(is.finite(x))
  )
  
  X_scaled <- X_scaled[
    ,
    good_cols,
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # 10. State labels
  # --------------------------------------------------------
  
  patient_ids <- gsub(
    "_",
    "-",
    rownames(X_scaled),
    fixed = TRUE
  )
  
  state <- robust_state_df$Robust_State[
    match(
      patient_ids,
      robust_state_df$submitter_id
    )
  ]
  
  state <- factor(
    state,
    levels = c(1, 2),
    labels = c("State1", "State2")
  )
  
  cat("\nStates:\n")
  print(table(state))
  
  # --------------------------------------------------------
  # 11. PCA
  # --------------------------------------------------------
  
  pca <- prcomp(
    X_scaled,
    center = FALSE,
    scale. = FALSE
  )
  
  var_exp <- (
    pca$sdev^2 /
      sum(pca$sdev^2)
  ) * 100
  
  # --------------------------------------------------------
  # 12. Distance
  # --------------------------------------------------------
  
  d <- dist(
    X_scaled,
    method = "euclidean"
  )
  
  state_df <- data.frame(
    State = state
  )
  
  # --------------------------------------------------------
  # 13. PERMANOVA
  # --------------------------------------------------------
  
  set.seed(123)
  
  perm <- vegan::adonis2(
    d ~ State,
    data = state_df,
    permutations = 9999
  )
  
  # --------------------------------------------------------
  # 14. Dispersion
  # --------------------------------------------------------
  
  disp <- vegan::betadisper(
    d,
    state
  )
  
  set.seed(123)
  
  disp_test <- permutest(
    disp,
    permutations = 9999
  )
  
  # --------------------------------------------------------
  # 15. Extract results
  # --------------------------------------------------------
  
  result <- data.frame(
    
    AS_Type = type,
    
    Events_Coverage =
      nrow(psi),
    
    Events_Top10 =
      ncol(X_scaled),
    
    PC1_percent =
      var_exp[1],
    
    PC2_percent =
      var_exp[2],
    
    PERMANOVA_R2 =
      perm$R2[1],
    
    PERMANOVA_F =
      perm$F[1],
    
    PERMANOVA_P =
      perm$`Pr(>F)`[1],
    
    Dispersion_P =
      disp_test$tab$`Pr(>F)`[1]
  )
  
  print(result)
  
  return(result)
}


# ============================================================
# RUN ALL 7 TYPES
# ============================================================

for(type in types_global) {
  
  global_results[[type]] <- global_AS_test(
    type,
    all_as_files[[type]]
  )
}


# ============================================================
# FINAL TABLE
# ============================================================

global_AS_summary <- dplyr::bind_rows(
  global_results
)

global_AS_summary <- global_AS_summary %>%
  mutate(
    PERMANOVA_R2_percent =
      PERMANOVA_R2 * 100
  )

cat("\n\n============================================\n")
cat("GLOBAL AS RESULTS\n")
cat("============================================\n")

print(global_AS_summary)
# ============================================================
# STEP 18 — Multiple-testing correction across AS types
# ============================================================

global_AS_summary$PERMANOVA_FDR <- p.adjust(
  global_AS_summary$PERMANOVA_P,
  method = "BH"
)

global_AS_summary <- global_AS_summary %>%
  dplyr::arrange(PERMANOVA_FDR)

print(
  global_AS_summary %>%
    dplyr::select(
      AS_Type,
      PERMANOVA_R2_percent,
      PERMANOVA_P,
      PERMANOVA_FDR,
      Dispersion_P
    )
)# ============================================================
# STEP 19 — Global AS association with molecular states
# PERMANOVA R2 + FDR
# ============================================================

library(dplyr)
library(ggplot2)

plot_global <- global_AS_summary %>%
    mutate(
        AS_Type = factor(
            AS_Type,
            levels = c("ES", "AA", "AD", "RI", "AP", "AT", "ME")
        ),
        FDR_label = ifelse(
            PERMANOVA_FDR < 0.001,
            "FDR < 0.001",
            paste0("FDR = ", sprintf("%.3f", PERMANOVA_FDR))
        )
    )

ggplot(
    plot_global,
    aes(
        x = AS_Type,
        y = PERMANOVA_R2_percent
    )
) +
    geom_col(
        width = 0.70
    ) +
    geom_text(
        aes(label = FDR_label),
        vjust = -0.5,
        size = 3.7
    ) +
    labs(
        x = "Alternative splicing type",
        y = expression("Variance explained by molecular state ("*R^2*" %)"),
        title = "Global association of alternative splicing profiles with molecular states"
    ) +
    coord_cartesian(
        ylim = c(0, max(plot_global$PERMANOVA_R2_percent) * 1.20)
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(
            hjust = 0.5,
            face = "bold"
        )
    )
# ============================================================
# STEP 19 — Global AS association with molecular states
# PERMANOVA R2 + FDR
# ============================================================

library(dplyr)
library(ggplot2)

plot_global <- global_AS_summary %>%
  mutate(
    AS_Type = factor(
      AS_Type,
      levels = c("ES", "AA", "AD", "RI", "AP", "AT", "ME")
    ),
    FDR_label = ifelse(
      PERMANOVA_FDR < 0.001,
      "FDR < 0.001",
      paste0("FDR = ", sprintf("%.3f", PERMANOVA_FDR))
    )
  )

ggplot(
  plot_global,
  aes(
    x = AS_Type,
    y = PERMANOVA_R2_percent
  )
) +
  geom_col(
    width = 0.70
  ) +
  geom_text(
    aes(label = FDR_label),
    vjust = -0.5,
    size = 3.7
  ) +
  labs(
    x = "Alternative splicing type",
    y = expression("Variance explained by molecular state ("*R^2*" %)"),
    title = "Global association of alternative splicing profiles with molecular states"
  ) +
  coord_cartesian(
    ylim = c(0, max(plot_global$PERMANOVA_R2_percent) * 1.20)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
# ============================================================
# STEP 20 — EVENT-WEIGHTED OVERALL ALTERNATIVE SPLICING
# ES + AA + AD + RI + AP + AT + ME
# ============================================================

library(data.table)
library(dplyr)
library(vegan)
library(ggplot2)

types_all <- c("ES", "AA", "AD", "RI", "AP", "AT", "ME")

all_as_files <- as_files
all_as_files[["ES"]] <- psi_clinical_file


# ============================================================
# 1. Function to extract coverage-filtered PSI matrix
# ============================================================

extract_AS_matrix <- function(type, file) {
  
  cat("\nReading:", type, "\n")
  
  # Header
  header_line <- readLines(
    file,
    n = 1,
    warn = FALSE
  )
  
  header_fields <- strsplit(
    header_line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  # Full data
  dat <- fread(
    file,
    skip = 14,
    header = FALSE,
    na.strings = c("null", "NA", "")
  )
  
  colnames(dat) <- header_fields
  
  # Available TCGA patients
  splice_cols <- header_fields[
    grepl("^TCGA[_-]", header_fields)
  ]
  
  splice_ids <- gsub(
    "_",
    "-",
    splice_cols,
    fixed = TRUE
  )
  
  # Same molecular-state patients
  overlap <- intersect(
    robust_state_df$submitter_id,
    splice_ids
  )
  
  patient_cols <- gsub(
    "-",
    "_",
    overlap,
    fixed = TRUE
  )
  
  # PSI matrix
  psi <- as.matrix(
    dat[, ..patient_cols]
  )
  
  storage.mode(psi) <- "numeric"
  
  # --------------------------------------------------------
  # >=75% overall coverage
  # --------------------------------------------------------
  
  coverage <- rowMeans(
    !is.na(psi)
  )
  
  psi <- psi[
    coverage >= 0.75,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # Remove zero/undefined variance
  # --------------------------------------------------------
  
  event_var <- apply(
    psi,
    1,
    var,
    na.rm = TRUE
  )
  
  keep <- is.finite(event_var) &
    event_var > 0
  
  psi <- psi[
    keep,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # Unique event names
  # --------------------------------------------------------
  
  rownames(psi) <- paste0(
    type,
    "_",
    dat$symbol[coverage >= 0.75][keep],
    "_",
    dat$as_id[coverage >= 0.75][keep]
  )
  
  cat(
    type,
    "events retained:",
    nrow(psi),
    "\n"
  )
  
  return(psi)
}


# ============================================================
# 2. Extract all 7 AS matrices
# ============================================================

as_matrices <- list()

for(type in types_all) {
  
  as_matrices[[type]] <- extract_AS_matrix(
    type,
    all_as_files[[type]]
  )
}


# ============================================================
# 3. Check identical patient order
# ============================================================

reference_patients <- colnames(
  as_matrices[[1]]
)

for(type in types_all) {
  
  as_matrices[[type]] <-
    as_matrices[[type]][
      ,
      reference_patients,
      drop = FALSE
    ]
}

cat(
  "\nPatients:",
  length(reference_patients),
  "\n"
)


# ============================================================
# 4. Combine ALL AS events
# ============================================================

overall_AS <- do.call(
  rbind,
  as_matrices
)

cat("\n====================================\n")
cat("COMBINED AS MATRIX\n")
cat("====================================\n")

cat(
  "Total events:",
  nrow(overall_AS),
  "\n"
)

cat(
  "Patients:",
  ncol(overall_AS),
  "\n"
)


# ============================================================
# 5. Variance across patients
# State information is NOT used
# ============================================================

overall_variance <- apply(
  overall_AS,
  1,
  var,
  na.rm = TRUE
)

summary(overall_variance)


# ============================================================
# 6. Top 10% most variable events
# ============================================================

variance_cutoff <- quantile(
  overall_variance,
  0.90,
  na.rm = TRUE
)

overall_AS_variable <- overall_AS[
  overall_variance >= variance_cutoff,
  ,
  drop = FALSE
]

cat(
  "\nTop 10% variable AS events:",
  nrow(overall_AS_variable),
  "\n"
)


# ============================================================
# 7. Median imputation
# ============================================================

overall_AS_imputed <- overall_AS_variable

for(i in seq_len(nrow(overall_AS_imputed))) {
  
  med <- median(
    overall_AS_imputed[i, ],
    na.rm = TRUE
  )
  
  overall_AS_imputed[
    i,
    is.na(overall_AS_imputed[i, ])
  ] <- med
}

cat(
  "Remaining NA:",
  sum(is.na(overall_AS_imputed)),
  "\n"
)


# ============================================================
# 8. Patients x AS events
# ============================================================

X_overall <- t(
  overall_AS_imputed
)

# Standardize every AS event
X_overall_scaled <- scale(
  X_overall
)

good_cols <- apply(
  X_overall_scaled,
  2,
  function(x) all(is.finite(x))
)

X_overall_scaled <- X_overall_scaled[
  ,
  good_cols,
  drop = FALSE
]

cat(
  "Final features:",
  ncol(X_overall_scaled),
  "\n"
)


# ============================================================
# 9. Molecular state
# ============================================================

patient_ids <- gsub(
  "_",
  "-",
  rownames(X_overall_scaled),
  fixed = TRUE
)

overall_state <- robust_state_df$Robust_State[
  match(
    patient_ids,
    robust_state_df$submitter_id
  )
]

overall_state <- factor(
  overall_state,
  levels = c(1, 2),
  labels = c("State1", "State2")
)

print(
  table(overall_state)
)


# ============================================================
# 10. PCA
# ============================================================

pca_overall <- prcomp(
  X_overall_scaled,
  center = FALSE,
  scale. = FALSE
)

variance_explained <- (
  pca_overall$sdev^2 /
    sum(pca_overall$sdev^2)
) * 100

cat(
  "\nPC1:",
  round(variance_explained[1], 3),
  "%\n"
)

cat(
  "PC2:",
  round(variance_explained[2], 3),
  "%\n"
)


# ============================================================
# 11. PCA dataframe
# ============================================================

pca_overall_df <- data.frame(
  
  Patient = patient_ids,
  
  PC1 = pca_overall$x[, 1],
  
  PC2 = pca_overall$x[, 2],
  
  State = overall_state
)


# ============================================================
# 12. PCA plot
# ============================================================

ggplot(
  pca_overall_df,
  aes(
    x = PC1,
    y = PC2,
    color = State
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.70
  ) +
  
  stat_ellipse(
    level = 0.95,
    linewidth = 0.8
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    
    color = NULL,
    
    title =
      "Overall alternative-splicing landscape"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


# ============================================================
# 13. Global distance
# ============================================================

overall_distance <- dist(
  X_overall_scaled,
  method = "euclidean"
)

overall_state_df <- data.frame(
  State = overall_state
)


# ============================================================
# 14. PERMANOVA
# ============================================================

set.seed(123)

overall_permanova <- vegan::adonis2(
  
  overall_distance ~ State,
  
  data = overall_state_df,
  
  permutations = 9999
)

cat("\n====================================\n")
cat("OVERALL AS — PERMANOVA\n")
cat("====================================\n")

print(overall_permanova)


# ============================================================
# 15. Dispersion test
# ============================================================

overall_dispersion <- vegan::betadisper(
  overall_distance,
  overall_state
)

set.seed(123)

overall_dispersion_test <- permutest(
  overall_dispersion,
  permutations = 9999
)

cat("\n====================================\n")
cat("OVERALL AS — DISPERSION\n")
cat("====================================\n")

print(overall_dispersion_test)


# ============================================================
# 16. Extract key statistics
# ============================================================

overall_R2 <- overall_permanova$R2[1]

overall_P <- overall_permanova$`Pr(>F)`[1]

overall_disp_P <-
  overall_dispersion_test$tab$`Pr(>F)`[1]


cat("\n====================================\n")
cat("FINAL OVERALL AS RESULT\n")
cat("====================================\n")

cat(
  "PERMANOVA R2:",
  overall_R2,
  "\n"
)

cat(
  "Variance explained:",
  round(overall_R2 * 100, 3),
  "%\n"
)

cat(
  "PERMANOVA P:",
  overall_P,
  "\n"
)

cat(
  "Dispersion P:",
  overall_disp_P,
  "\n"
)
# ============================================================
# GLOBAL AS — PCA
# ============================================================

library(ggplot2)

ggplot(
  pca_overall_df,
  aes(
    x = PC1,
    y = PC2,
    color = State
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.70
  ) +
  stat_ellipse(
    level = 0.95,
    linewidth = 0.9
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    color = NULL,
    title = "Overall alternative-splicing landscape"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
# ============================================================
# GLOBAL AS — variance explained plot
# ============================================================

library(dplyr)
library(ggplot2)

overall_row <- data.frame(
  AS_Type = "Overall AS",
  PERMANOVA_R2_percent = overall_R2 * 100,
  PERMANOVA_P = overall_P,
  Dispersion_P = overall_disp_P
)

plot_r2 <- global_AS_summary %>%
  dplyr::select(
    AS_Type,
    PERMANOVA_R2_percent,
    PERMANOVA_P,
    Dispersion_P
  ) %>%
  dplyr::bind_rows(overall_row)

plot_r2$AS_Type <- factor(
  plot_r2$AS_Type,
  levels = c(
    "ES", "AA", "AD", "RI",
    "AP", "AT", "ME",
    "Overall AS"
  )
)

plot_r2 <- plot_r2 %>%
  mutate(
    P_label = case_when(
      PERMANOVA_P < 0.001 ~
        paste0(
          "P = ",
          format(
            PERMANOVA_P,
            scientific = TRUE,
            digits = 1
          )
        ),
      
      TRUE ~
        paste0(
          "P = ",
          sprintf("%.3f", PERMANOVA_P)
        )
    )
  )

ggplot(
  plot_r2,
  aes(
    x = AS_Type,
    y = PERMANOVA_R2_percent
  )
) +
  geom_col(
    width = 0.70
  ) +
  geom_text(
    aes(label = P_label),
    vjust = -0.45,
    size = 3.4
  ) +
  labs(
    x = "Alternative-splicing profile",
    y = expression(
      "Variance explained by molecular state (" *
        R^2 * " %)"
    ),
    title = "Association between alternative splicing and molecular state"
  ) +
  coord_cartesian(
    ylim = c(
      0,
      max(plot_r2$PERMANOVA_R2_percent) * 1.20
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
# ============================================================
# OVERALL AS — PCoA visualization of PERMANOVA distance
# ============================================================

library(ggplot2)
library(dplyr)

# PCoA from the SAME distance matrix used for PERMANOVA
pcoa_as <- cmdscale(
  overall_distance,
  k = 2,
  eig = TRUE
)

# Variance represented by PCoA axes
positive_eig <- pcoa_as$eig[pcoa_as$eig > 0]

pcoa_var1 <- 100 * pcoa_as$eig[1] / sum(positive_eig)
pcoa_var2 <- 100 * pcoa_as$eig[2] / sum(positive_eig)

# Data frame
pcoa_df <- data.frame(
  Patient = patient_ids,
  PCoA1 = pcoa_as$points[, 1],
  PCoA2 = pcoa_as$points[, 2],
  State = overall_state
)

# State centroids
centroids <- pcoa_df %>%
  group_by(State) %>%
  summarise(
    PCoA1_centroid = mean(PCoA1),
    PCoA2_centroid = mean(PCoA2),
    .groups = "drop"
  )

# Add centroid coordinates to every patient
pcoa_df <- pcoa_df %>%
  left_join(
    centroids,
    by = "State"
  )

# Plot
ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = State
  )
) +
  
  # connect each patient to its state centroid
  geom_segment(
    aes(
      xend = PCoA1_centroid,
      yend = PCoA2_centroid
    ),
    alpha = 0.12,
    linewidth = 0.35
  ) +
  
  # patients
  geom_point(
    size = 2.4,
    alpha = 0.65
  ) +
  
  # 95% ellipse
  stat_ellipse(
    level = 0.95,
    linewidth = 1
  ) +
  
  # centroids
  geom_point(
    data = centroids,
    aes(
      x = PCoA1_centroid,
      y = PCoA2_centroid,
      color = State
    ),
    shape = 4,
    size = 6,
    stroke = 2,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = paste0(
      "PCoA1 (",
      round(pcoa_var1, 1),
      "%)"
    ),
    
    y = paste0(
      "PCoA2 (",
      round(pcoa_var2, 1),
      "%)"
    ),
    
    color = NULL,
    
    title =
      "Overall alternative-splicing landscape",
    
    subtitle =
      "PERMANOVA: R² = 0.0066, P = 3 × 10⁻⁴; dispersion P = 0.443"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "top",
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )
# ============================================================
# STEP 21 — Patient-level global AS distance contrast
# ============================================================

library(dplyr)
library(ggplot2)

D <- as.matrix(overall_distance)

state_vec <- overall_state
names(state_vec) <- rownames(D)

contrast_df <- data.frame(
  Patient = rownames(D),
  State = state_vec,
  AS_distance_contrast = NA_real_
)

for(i in seq_len(nrow(D))) {
  
  patient <- rownames(D)[i]
  
  # Exclude self-distance
  other_patients <- setdiff(
    rownames(D),
    patient
  )
  
  state1_patients <- other_patients[
    state_vec[other_patients] == "State1"
  ]
  
  state2_patients <- other_patients[
    state_vec[other_patients] == "State2"
  ]
  
  mean_to_state1 <- mean(
    D[patient, state1_patients]
  )
  
  mean_to_state2 <- mean(
    D[patient, state2_patients]
  )
  
  # Positive = closer to State1
  # Negative = closer to State2
  contrast_df$AS_distance_contrast[i] <-
    mean_to_state2 - mean_to_state1
}


# ============================================================
# Plot
# ============================================================

ggplot(
  contrast_df,
  aes(
    x = State,
    y = AS_distance_contrast,
    fill = State
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.30,
    width = 0.8
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.10,
    size = 1.4,
    alpha = 0.45
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = NULL,
    y = "Relative AS distance\n(State2 distance - State1 distance)",
    title = "Patient-level alternative-splicing state similarity",
    subtitle = "Positive values indicate greater similarity to State1; negative values to State2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )
# ============================================================
# STEP 22 — Strong AS events
# ============================================================

library(dplyr)

strong_list <- lapply(
  all_as_results,
  function(x) {
    
    x <- as.data.frame(x)
    
    x <- dplyr::select(
      x,
      AS_Type,
      symbol,
      as_id,
      Median_State1,
      Median_State2,
      Delta_PSI,
      P_value,
      FDR
    )
    
    # Make identifiers consistent
    x$AS_Type <- as.character(x$AS_Type)
    x$symbol  <- as.character(x$symbol)
    x$as_id   <- as.character(x$as_id)
    
    return(x)
  }
)

strong_AS_events <- dplyr::bind_rows(strong_list) %>%
  dplyr::filter(
    FDR < 0.05,
    abs(Delta_PSI) >= 0.10
  ) %>%
  dplyr::mutate(
    Direction = ifelse(
      Delta_PSI > 0,
      "State2_high",
      "State1_high"
    )
  ) %>%
  dplyr::arrange(AS_Type, FDR)


# ============================================================
# RESULTS
# ============================================================

cat(
  "Total strong AS events:",
  nrow(strong_AS_events),
  "\n\n"
)

print(
  strong_AS_events %>%
    dplyr::count(AS_Type, Direction)
)

print(
  strong_AS_events,
  n = Inf
)
# ============================================================
# Print full 34-event table safely
# ============================================================

strong_AS_print <- as.data.frame(strong_AS_events)

strong_AS_print <- strong_AS_print[
  order(
    strong_AS_print$AS_Type,
    strong_AS_print$FDR
  ),
]

row.names(strong_AS_print) <- NULL

write.table(
  strong_AS_print,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
# ============================================================
# STEP 23 — Visualize 34 strong differential AS events
# ============================================================

library(dplyr)
library(ggplot2)

as_plot <- strong_AS_events %>%
  mutate(
    Event = paste0(symbol, " (", AS_Type, ":", as_id, ")"),
    Event = reorder(Event, Delta_PSI)
  )

ggplot(
  as_plot,
  aes(
    x = Delta_PSI,
    y = Event,
    color = AS_Type
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Delta_PSI,
      y = Event,
      yend = Event
    ),
    linewidth = 0.6,
    alpha = 0.6
  ) +
  geom_point(
    size = 3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-0.10, 0.10),
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    x = expression(Delta*"PSI (State2 - State1)"),
    y = NULL,
    color = "AS type",
    title = "Strong differential alternative-splicing events",
    subtitle = "FDR < 0.05 and |ΔPSI| ≥ 0.10"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold"
    )
  )
# ============================================================
# STEP 24 — GO enrichment of strong AS genes
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# Strong AS genes
strong_AS_genes <- unique(strong_AS_events$symbol)

cat("Unique strong AS genes:", length(strong_AS_genes), "\n")
print(strong_AS_genes)


# ============================================================
# Background = genes represented in all tested AS events
# ============================================================

background_AS_genes <- unique(
  unlist(
    lapply(
      all_as_results,
      function(x) as.character(x$symbol)
    )
  )
)

background_AS_genes <- background_AS_genes[
  !is.na(background_AS_genes) &
    background_AS_genes != ""
]

cat(
  "AS background genes:",
  length(background_AS_genes),
  "\n"
)


# ============================================================
# Convert SYMBOL -> ENTREZ
# ============================================================

strong_map <- bitr(
  strong_AS_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_map <- bitr(
  background_AS_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# ============================================================
# GO Biological Process
# ============================================================

go_AS <- enrichGO(
  gene = unique(strong_map$ENTREZID),
  universe = unique(background_map$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

go_AS_df <- as.data.frame(go_AS)

go_AS_df <- go_AS_df %>%
  arrange(p.adjust)

print(
  go_AS_df[
    ,
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "geneID"
    )
  ],
  row.names = FALSE
)
# ============================================================
# Plot significant GO terms
# ============================================================

go_sig <- go_AS_df %>%
  filter(p.adjust < 0.05) %>%
  slice_head(n = 15)

if(nrow(go_sig) > 0) {
  
  go_sig$Description <- factor(
    go_sig$Description,
    levels = rev(go_sig$Description)
  )
  
  ggplot(
    go_sig,
    aes(
      x = -log10(p.adjust),
      y = Description
    )
  ) +
    geom_col() +
    labs(
      x = expression(-log[10]("FDR")),
      y = NULL,
      title = "Biological processes associated with differential AS genes"
    ) +
    theme_classic(base_size = 13)
  
} else {
  
  cat(
    "\nNo GO Biological Process terms reached FDR < 0.05.\n"
  )
}
# ============================================================
# STEP 25 — PSI matrix for the 34 strong AS events
# ============================================================

library(dplyr)

strong_event_names <- paste0(
  strong_AS_events$AS_Type,
  "_",
  strong_AS_events$symbol,
  "_",
  strong_AS_events$as_id
)

strong_PSI_list <- list()

for(type in unique(strong_AS_events$AS_Type)) {
  
  wanted <- strong_event_names[
    strong_AS_events$AS_Type == type
  ]
  
  available <- intersect(
    wanted,
    rownames(as_matrices[[type]])
  )
  
  strong_PSI_list[[type]] <-
    as_matrices[[type]][
      available,
      ,
      drop = FALSE
    ]
}

strong_PSI <- do.call(
  rbind,
  strong_PSI_list
)

cat("Strong AS events recovered:", nrow(strong_PSI), "\n")
cat("Patients:", ncol(strong_PSI), "\n")

print(rownames(strong_PSI))
# ============================================================
# STEP 26 — 34-event AS signature
# Heatmap + PCA + PERMANOVA
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggplot2)

# ============================================================
# 1. Median imputation
# ============================================================

strong_PSI_imp <- strong_PSI

for(i in seq_len(nrow(strong_PSI_imp))) {
  
  med <- median(
    strong_PSI_imp[i, ],
    na.rm = TRUE
  )
  
  strong_PSI_imp[
    i,
    is.na(strong_PSI_imp[i, ])
  ] <- med
}


# ============================================================
# 2. State annotation
# ============================================================

strong_patient_ids <- gsub(
  "_",
  "-",
  colnames(strong_PSI_imp),
  fixed = TRUE
)

strong_state <- robust_state_df$Robust_State[
  match(
    strong_patient_ids,
    robust_state_df$submitter_id
  )
]

strong_state <- factor(
  strong_state,
  levels = c(1, 2),
  labels = c("State1", "State2")
)

print(table(strong_state))


# ============================================================
# 3. Z-score PSI for HEATMAP only
# ============================================================

strong_PSI_z <- t(
  scale(
    t(strong_PSI_imp)
  )
)

# Order patients by State
patient_order_heat <- order(strong_state)

strong_PSI_z_heat <- strong_PSI_z[
  ,
  patient_order_heat,
  drop = FALSE
]

state_heat <- strong_state[
  patient_order_heat
]


# ============================================================
# 4. Heatmap annotation
# ============================================================

ha <- HeatmapAnnotation(
  State = state_heat
)

Heatmap(
  strong_PSI_z_heat,
  
  name = "PSI\nZ-score",
  
  top_annotation = ha,
  
  show_column_names = FALSE,
  
  show_row_names = TRUE,
  
  row_names_gp = grid::gpar(
    fontsize = 8
  ),
  
  cluster_rows = TRUE,
  
  cluster_columns = FALSE,
  
  column_split = state_heat,
  
  column_title = NULL,
  
  row_title =
    "Strong differential AS events",
  
  border = FALSE
)


# ============================================================
# 5. PCA
# ============================================================

X_strong <- t(strong_PSI_imp)

# Standardize AS events
X_strong_scaled <- scale(X_strong)

good_cols <- apply(
  X_strong_scaled,
  2,
  function(x) all(is.finite(x))
)

X_strong_scaled <- X_strong_scaled[
  ,
  good_cols,
  drop = FALSE
]

pca_strong <- prcomp(
  X_strong_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_var <- (
  pca_strong$sdev^2 /
    sum(pca_strong$sdev^2)
) * 100

pca_strong_df <- data.frame(
  Patient = strong_patient_ids,
  PC1 = pca_strong$x[,1],
  PC2 = pca_strong$x[,2],
  State = strong_state
)


# ============================================================
# 6. PCA plot
# ============================================================

ggplot(
  pca_strong_df,
  aes(
    x = PC1,
    y = PC2,
    color = State
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  stat_ellipse(
    level = 0.95,
    linewidth = 0.9
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    ),
    
    color = NULL,
    
    title =
      "Differential AS signature across molecular states"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )


# ============================================================
# 7. PERMANOVA
# ============================================================

strong_distance <- dist(
  X_strong_scaled,
  method = "euclidean"
)

strong_state_df <- data.frame(
  State = strong_state
)

set.seed(123)

strong_permanova <- adonis2(
  strong_distance ~ State,
  data = strong_state_df,
  permutations = 9999
)

cat("\n====================================\n")
cat("34-EVENT AS SIGNATURE — PERMANOVA\n")
cat("====================================\n")

print(strong_permanova)


# ============================================================
# 8. Dispersion
# ============================================================

strong_disp <- betadisper(
  strong_distance,
  strong_state
)

set.seed(123)

strong_disp_test <- permutest(
  strong_disp,
  permutations = 9999
)

cat("\n====================================\n")
cat("34-EVENT AS SIGNATURE — DISPERSION\n")
cat("====================================\n")

print(strong_disp_test)


# ============================================================
# 9. Final numbers
# ============================================================

cat(
  "\nVariance explained by State:",
  round(
    strong_permanova$R2[1] * 100,
    2
  ),
  "%\n"
)

cat(
  "PERMANOVA P:",
  strong_permanova$`Pr(>F)`[1],
  "\n"
)

cat(
  "Dispersion P:",
  strong_disp_test$tab$`Pr(>F)`[1],
  "\n"
)
# ============================================================
# STEP 27 — Align age with AS patients
# ============================================================

# AS patient IDs
as_patient_ids <- gsub(
  "_",
  "-",
  colnames(strong_PSI),
  fixed = TRUE
)

# Match clinical data
age_AS <- clinical_final$age_years[
  match(
    as_patient_ids,
    clinical_final$submitter_id
  )
]

state_AS <- robust_state_df$Robust_State[
  match(
    as_patient_ids,
    robust_state_df$submitter_id
  )
]

# Combined metadata
AS_metadata <- data.frame(
  Patient = as_patient_ids,
  State = factor(
    state_AS,
    levels = c(1, 2),
    labels = c("State1", "State2")
  ),
  Age = age_AS
)

rownames(AS_metadata) <- AS_metadata$Patient


# ============================================================
# Check
# ============================================================

cat("Total AS patients:", nrow(AS_metadata), "\n")
cat("Missing age:", sum(is.na(AS_metadata$Age)), "\n\n")

print(table(AS_metadata$State, useNA = "ifany"))

cat("\nAge by state:\n")

print(
  aggregate(
    Age ~ State,
    data = AS_metadata,
    FUN = function(x)
      c(
        N = length(x),
        Median = median(x, na.rm = TRUE),
        Mean = mean(x, na.rm = TRUE)
      )
  )
)
# ============================================================
# STEP 28 — Age-adjusted differential AS analysis
# PSI ~ State + Age
# ============================================================

library(limma)
library(dplyr)

# Complete clinical cases
complete_patients <- AS_metadata %>%
  filter(
    !is.na(Age),
    !is.na(State)
  )

# Design matrix
design_AS <- model.matrix(
  ~ State + Age,
  data = complete_patients
)

colnames(design_AS)

# Expected:
# "(Intercept)" "StateState2" "Age"


# ============================================================
# Function for each AS type
# ============================================================

age_adjust_AS <- function(type) {
  
  mat <- as_matrices[[type]]
  
  # Convert column IDs
  mat_ids <- gsub(
    "_",
    "-",
    colnames(mat),
    fixed = TRUE
  )
  
  # Patients available in both matrix + metadata
  keep_patients <- complete_patients$Patient[
    complete_patients$Patient %in% mat_ids
  ]
  
  # Match matrix columns
  idx <- match(
    keep_patients,
    mat_ids
  )
  
  mat_sub <- mat[, idx, drop = FALSE]
  
  # Force same order as metadata
  meta_sub <- complete_patients[
    match(
      keep_patients,
      complete_patients$Patient
    ),
  ]
  
  # --------------------------------------------------------
  # Coverage >=75% within BOTH states
  # --------------------------------------------------------
  
  state1_idx <- meta_sub$State == "State1"
  state2_idx <- meta_sub$State == "State2"
  
  cov1 <- rowMeans(
    !is.na(mat_sub[, state1_idx, drop = FALSE])
  )
  
  cov2 <- rowMeans(
    !is.na(mat_sub[, state2_idx, drop = FALSE])
  )
  
  keep_events <- cov1 >= 0.75 & cov2 >= 0.75
  
  mat_sub <- mat_sub[
    keep_events,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------
  # Median imputation per event
  # --------------------------------------------------------
  
  for(i in seq_len(nrow(mat_sub))) {
    
    med <- median(
      mat_sub[i, ],
      na.rm = TRUE
    )
    
    mat_sub[
      i,
      is.na(mat_sub[i, ])
    ] <- med
  }
  
  # --------------------------------------------------------
  # Design
  # --------------------------------------------------------
  
  design <- model.matrix(
    ~ State + Age,
    data = meta_sub
  )
  
  # --------------------------------------------------------
  # limma
  # --------------------------------------------------------
  
  fit <- lmFit(
    mat_sub,
    design
  )
  
  fit <- eBayes(fit)
  
  res <- topTable(
    fit,
    coef = "StateState2",
    number = Inf,
    adjust.method = "none",
    sort.by = "none"
  )
  
  # --------------------------------------------------------
  # Event information
  # --------------------------------------------------------
  
  event_names <- rownames(mat_sub)
  
  res$Event <- event_names
  res$AS_Type <- type
  
  res$Delta_PSI_adjusted <- res$logFC
  
  return(res)
}


# ============================================================
# Run all 7 AS types
# ============================================================

AS_types <- c(
  "ES", "AA", "AD",
  "RI", "AP", "AT", "ME"
)

age_adjusted_results <- lapply(
  AS_types,
  age_adjust_AS
)

names(age_adjusted_results) <- AS_types


# ============================================================
# Combine ALL events
# ============================================================

age_adjusted_all <- bind_rows(
  age_adjusted_results
)


# ============================================================
# GLOBAL BH across all AS events
# ============================================================

age_adjusted_all$Global_FDR <- p.adjust(
  age_adjusted_all$P.Value,
  method = "BH"
)


# ============================================================
# Summary
# ============================================================

cat(
  "Total age-adjusted AS events tested:",
  nrow(age_adjusted_all),
  "\n"
)

cat(
  "Global FDR < 0.05:",
  sum(
    age_adjusted_all$Global_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Global FDR < 0.05 & |adjusted ΔPSI| >= 0.10:",
  sum(
    age_adjusted_all$Global_FDR < 0.05 &
      abs(age_adjusted_all$Delta_PSI_adjusted) >= 0.10,
    na.rm = TRUE
  ),
  "\n"
)
# ============================================================
# STEP 29 — Final age-adjusted strong AS events
# ============================================================

final_AS18 <- age_adjusted_all %>%
  dplyr::filter(
    Global_FDR < 0.05,
    abs(Delta_PSI_adjusted) >= 0.10
  ) %>%
  dplyr::mutate(
    Direction = ifelse(
      Delta_PSI_adjusted > 0,
      "State2_high",
      "State1_high"
    )
  ) %>%
  dplyr::arrange(Global_FDR)


cat(
  "Final strong age-adjusted AS events:",
  nrow(final_AS18),
  "\n\n"
)

# Number by AS type
print(
  final_AS18 %>%
    dplyr::count(AS_Type, Direction)
)


# Full table
write.table(
  as.data.frame(
    final_AS18[
      ,
      c(
        "AS_Type",
        "Event",
        "Delta_PSI_adjusted",
        "P.Value",
        "Global_FDR",
        "Direction"
      )
    ]
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
# Extract gene name from Event:
# e.g. AP_ADPRHL1_26374 -> ADPRHL1

final_AS18$Gene <- sub(
  "^[^_]+_([^_]+)_.*$",
  "\\1",
  final_AS18$Event
)

cat(
  "\nUnique genes:",
  length(unique(final_AS18$Gene)),
  "\n\n"
)

cat(
  paste(
    unique(final_AS18$Gene),
    collapse = "\n"
  )
)
# ============================================================
# STEP 29 — Final age-adjusted strong AS events
# Global FDR < 0.05 & |adjusted Delta PSI| >= 0.10
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Select final strong events
# ------------------------------------------------------------

final_AS18 <- age_adjusted_all %>%
  dplyr::filter(
    Global_FDR < 0.05,
    abs(Delta_PSI_adjusted) >= 0.10
  ) %>%
  dplyr::mutate(
    Direction = ifelse(
      Delta_PSI_adjusted > 0,
      "State2_high",
      "State1_high"
    )
  ) %>%
  dplyr::arrange(Global_FDR)


# ------------------------------------------------------------
# 2. Extract gene name from Event
# Example: AP_ADPRHL1_26374 -> ADPRHL1
# ------------------------------------------------------------

final_AS18$Gene <- sub(
  "^[^_]+_([^_]+)_.*$",
  "\\1",
  final_AS18$Event
)


# ------------------------------------------------------------
# 3. Basic summary
# ------------------------------------------------------------

cat("\n====================================\n")
cat("FINAL AGE-ADJUSTED AS SIGNATURE\n")
cat("====================================\n\n")

cat(
  "Strong AS events:",
  nrow(final_AS18),
  "\n"
)

cat(
  "Unique genes:",
  length(unique(final_AS18$Gene)),
  "\n\n"
)


# ------------------------------------------------------------
# 4. AS type + direction
# ------------------------------------------------------------

cat("Events by AS type and direction:\n\n")

print(
  final_AS18 %>%
    dplyr::count(
      AS_Type,
      Direction
    )
)


# ------------------------------------------------------------
# 5. Unique genes
# ------------------------------------------------------------

cat("\n====================================\n")
cat("UNIQUE GENES\n")
cat("====================================\n\n")

cat(
  paste(
    unique(final_AS18$Gene),
    collapse = "\n"
  )
)


# ------------------------------------------------------------
# 6. Full final table
# ------------------------------------------------------------

cat("\n\n====================================\n")
cat("FINAL EVENT TABLE\n")
cat("====================================\n\n")

final_table <- final_AS18 %>%
  dplyr::select(
    AS_Type,
    Gene,
    Event,
    Delta_PSI_adjusted,
    P.Value,
    Global_FDR,
    Direction
  )

write.table(
  as.data.frame(final_table),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)



# ============================================================
# STEP 30 — Age-adjusted AS signature plot
# ============================================================

library(dplyr)
library(ggplot2)

plot_AS18 <- final_AS18 %>%
  mutate(
    Label = paste0(
      Gene,
      " (", AS_Type, ":", sub(".*_", "", Event), ")"
    ),
    Label = reorder(
      Label,
      Delta_PSI_adjusted
    )
  )

ggplot(
  plot_AS18,
  aes(
    x = Delta_PSI_adjusted,
    y = Label,
    color = AS_Type
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Delta_PSI_adjusted,
      y = Label,
      yend = Label
    ),
    linewidth = 0.7,
    alpha = 0.65
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-0.10, 0.10),
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    x = expression(
      "Age-adjusted " * Delta * "PSI (State2 - State1)"
    ),
    y = NULL,
    color = "AS type",
    title = "Age-adjusted state-associated alternative splicing",
    subtitle = "Global FDR < 0.05 and |adjusted ΔPSI| ≥ 0.10"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )
# ============================================================
# STEP 31 — Draw heatmap without rasterization
# ============================================================

Heatmap(
  final18_z,
  
  name = "PSI\nZ-score",
  
  top_annotation = ha18,
  
  column_split = state18_ordered,
  
  cluster_columns = FALSE,
  
  cluster_rows = TRUE,
  
  show_column_names = FALSE,
  show_row_names = TRUE,
  
  row_names_gp = grid::gpar(
    fontsize = 9
  ),
  
  row_title = "Age-adjusted AS events",
  
  column_title = NULL,
  
  border = FALSE,
  
  use_raster = FALSE
)
# ============================================================
# STEP 32 — PCA + PERMANOVA
# Final 18 age-adjusted AS events
# ============================================================

library(vegan)
library(ggplot2)

# ------------------------------------------------------------
# 1. Patient x AS-event matrix
# ------------------------------------------------------------

X18 <- t(final18_imp)

# Standardize events
X18_scaled <- scale(X18)

# Remove problematic columns if any
good_cols <- apply(
  X18_scaled,
  2,
  function(x) all(is.finite(x))
)

X18_scaled <- X18_scaled[
  ,
  good_cols,
  drop = FALSE
]

cat(
  "Events used:",
  ncol(X18_scaled),
  "\n"
)


# ------------------------------------------------------------
# 2. PCA
# ------------------------------------------------------------

pca18 <- prcomp(
  X18_scaled,
  center = FALSE,
  scale. = FALSE
)

pca18_var <- (
  pca18$sdev^2 /
    sum(pca18$sdev^2)
) * 100


pca18_df <- data.frame(
  Patient = matrix_patient_ids,
  PC1 = pca18$x[,1],
  PC2 = pca18$x[,2],
  State = state18
)


# ------------------------------------------------------------
# 3. PCA plot
# ------------------------------------------------------------

ggplot(
  pca18_df,
  aes(
    x = PC1,
    y = PC2,
    color = State
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  stat_ellipse(
    level = 0.95,
    linewidth = 0.9
  ) +
  labs(
    x = paste0(
      "PC1 (",
      round(pca18_var[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca18_var[2], 1),
      "%)"
    ),
    color = NULL,
    title =
      "Age-adjusted AS signature across molecular states"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )


# ------------------------------------------------------------
# 4. Distance matrix
# ------------------------------------------------------------

dist18 <- dist(
  X18_scaled,
  method = "euclidean"
)


# ------------------------------------------------------------
# 5. PERMANOVA
# ------------------------------------------------------------

meta18 <- data.frame(
  State = state18
)

set.seed(123)

permanova18 <- adonis2(
  dist18 ~ State,
  data = meta18,
  permutations = 9999
)

cat("\n====================================\n")
cat("FINAL 18 AS EVENTS — PERMANOVA\n")
cat("====================================\n")

print(permanova18)


# ------------------------------------------------------------
# 6. Dispersion test
# ------------------------------------------------------------

disp18 <- betadisper(
  dist18,
  state18
)

set.seed(123)

disp18_test <- permutest(
  disp18,
  permutations = 9999
)

cat("\n====================================\n")
cat("FINAL 18 AS EVENTS — DISPERSION\n")
cat("====================================\n")

print(disp18_test)


# ------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------

cat(
  "\nVariance explained by State:",
  round(
    permanova18$R2[1] * 100,
    2
  ),
  "%\n"
)

cat(
  "PERMANOVA P:",
  permanova18$`Pr(>F)`[1],
  "\n"
)

cat(
  "Dispersion P:",
  disp18_test$tab$`Pr(>F)`[1],
  "\n"
)

cat(
  "PC1 variance:",
  round(pca18_var[1], 2),
  "%\n"
)

cat(
  "PC2 variance:",
  round(pca18_var[2], 2),
  "%\n"
)
# ============================================================
# MOLECULAR STATE → RNA + AS WORKFLOW
# ============================================================

library(DiagrammeR)

grViz("
digraph molecular_state_workflow {

graph [
  layout = dot,
  rankdir = TB,
  nodesep = 0.45,
  ranksep = 0.55
]

node [
  shape = box,
  style = 'rounded,filled',
  fontname = Helvetica,
  fontsize = 10,
  color = '#455A64',
  penwidth = 1.1,
  margin = 0.12
]

edge [
  color = '#607D8B',
  penwidth = 1.2,
  arrowsize = 0.7
]


# ============================================================
# DISCOVERY
# ============================================================

TCGA [
  label = 'TCGA-OV PRIMARY TUMORS\\n277 patients\\nMatched RNA + CNV + mutation data',
  fillcolor = '#DCEAF5'
]

MOFA [
  label = 'MOFA2 MULTI-OMICS INTEGRATION\\n1,722 RNA + 1,902 CNV + 107 mutation features\\n20 latent factors → Factors 1–10 retained',
  fillcolor = '#DDEFD8'
]

EXPLORE [
  label = 'EXPLORATORY 5-STATE SOLUTION\\nState 1–5\\nApparently distinct molecular substructure',
  fillcolor = '#F8DDDD'
]

ROBUSTNESS [
  label = 'ROBUSTNESS TESTING\\nFactor-number sensitivity + bootstrap validation\\n5-state solution unstable\\nMedian bootstrap ARI ≈ 0.35',
  fillcolor = '#FFF0CF'
]

STATE2 [
  label = 'ROBUST 2-STATE SOLUTION\\nState 1: n = 226   |   State 2: n = 51\\nK=2 bootstrap median ARI ≈ 0.80',
  fillcolor = '#D6EEDB'
]


TCGA -> MOFA
MOFA -> EXPLORE
EXPLORE -> ROBUSTNESS
ROBUSTNESS -> STATE2


# ============================================================
# RNA BRANCH
# ============================================================

RNA [
  label = 'STATE-ASSOCIATED RNA SIGNAL\\n110 genes at FDR < 0.05\\n34 State 1-high   |   76 State 2-high',
  fillcolor = '#E7E0F3'
]

STRING [
  label = 'STRING NETWORK FILTERING\\n37 state-associated genes connected\\nwithin protein-interaction networks',
  fillcolor = '#E7E0F3'
]

LASSO [
  label = 'NETWORK-GUIDED NESTED CV LASSO\\nOuter-fold DE selection + STRING filtering + inner LASSO\\nAUC = 0.763 (95% CI: 0.688–0.838)\\nBalanced accuracy = 0.729',
  fillcolor = '#DCE7F4'
]

CORE9 [
  label = 'STABLE 9-GENE RNA CORE\\nSelected in ≥80% of outer folds\\nADGRL3 · PIWIL1 · SLC15A1 · ZIC5 · GATA2\\nFOXG1 · MAGEA11 · PALM3 · SLCO1A2',
  fillcolor = '#D4ECDF'
]

CORE4 [
  label = 'HIGHEST-STABILITY RNA CORE\\nSelected in 10/10 outer folds\\nADGRL3 · PIWIL1 · SLC15A1 · ZIC5',
  fillcolor = '#BFE2CF'
]


# ============================================================
# ALTERNATIVE SPLICING BRANCH
# ============================================================

SPLICE [
  label = 'TCGA SPLICESEQ INTEGRATION\\n273/277 patients matched\\nES · AA · AD · RI · AP · AT · ME',
  fillcolor = '#E5EAF4'
]

GLOBALAS [
  label = 'GLOBAL AS LANDSCAPE\\nModest but significant state association\\nPERMANOVA R² = 0.0066, P = 3 × 10⁻⁴\\nDispersion P = 0.443',
  fillcolor = '#E8E4F3'
]

ASTYPES [
  label = 'AS-TYPE-SPECIFIC ANALYSIS\\nState-associated shifts across multiple AS classes\\nAP and AT contained the largest numbers\\nof strong differential events',
  fillcolor = '#E8E4F3'
]

AS34 [
  label = 'INITIAL STRONG AS SIGNAL\\n34 events · 21 unique genes\\nFDR < 0.05 and |ΔPSI| ≥ 0.10',
  fillcolor = '#F6E5CF'
]

AGE [
  label = 'AGE-ADJUSTED PAN-AS ANALYSIS\\nPSI ~ Molecular State + Age\\nState 1 median age = 57.8 y\\nState 2 median age = 66.1 y\\nGlobal BH correction across AS events',
  fillcolor = '#FFE4B8'
]

AS18 [
  label = 'FINAL AGE-ADJUSTED AS SIGNAL\\n126 events at global FDR < 0.05\\n18 strong events · 13 unique genes\\nGlobal FDR < 0.05 and |adjusted ΔPSI| ≥ 0.10',
  fillcolor = '#D5EDDF'
]

ASGENES [
  label = '13-GENE AS SET\\nCASK · TMEM8B · FLOT2 · ZNF576 · TCF4\\nSLMAP · SIMC1 · ADPRHL1 · SMARCD3\\nFAM86B1 · ZNF559 · DMKN · SPEG',
  fillcolor = '#C5E5D4'
]

ASCHAR [
  label = 'AS SIGNATURE CHARACTERIZATION\\nPC1 = 23.4%   |   PC2 = 10.6%\\nPERMANOVA R² = 0.072, P = 1 × 10⁻⁴\\nDispersion P = 9 × 10⁻⁴\\nState-associated pattern with within-state heterogeneity',
  fillcolor = '#DCE7F4'
]


# ============================================================
# BRANCHES
# ============================================================

STATE2 -> RNA
RNA -> STRING
STRING -> LASSO
LASSO -> CORE9
CORE9 -> CORE4

STATE2 -> SPLICE
SPLICE -> GLOBALAS
GLOBALAS -> ASTYPES
ASTYPES -> AS34
AS34 -> AGE
AGE -> AS18
AS18 -> ASGENES
ASGENES -> ASCHAR


# Keep parallel branches aligned
{rank=same; RNA; SPLICE}
{rank=same; STRING; GLOBALAS}
{rank=same; LASSO; ASTYPES}
{rank=same; CORE9; AS34}
{rank=same; CORE4; AGE}

}
")
# ============================================================
# COMPLETE MOLECULAR STATE WORKFLOW
# MULTI-OMICS + RNA + ALTERNATIVE SPLICING + DRUG RESPONSE
# ============================================================

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <- "~/Desktop"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# CREATE COMPLETE WORKFLOW
# ============================================================

workflow_complete <- grViz("
digraph complete_workflow {

graph [
  layout = dot,
  rankdir = TB,
  nodesep = 0.45,
  ranksep = 0.55,
  bgcolor = white
]

node [
  shape = box,
  style = 'rounded,filled',
  fontname = Helvetica,
  fontsize = 10,
  color = '#455A64',
  fontcolor = '#263238',
  penwidth = 1.1,
  margin = 0.12
]

edge [
  color = '#607D8B',
  penwidth = 1.2,
  arrowsize = 0.7
]


# ============================================================
# MULTI-OMICS STATE DISCOVERY
# ============================================================

TCGA [
  label = 'TCGA-OV PRIMARY TUMORS\\n277 patients\\nMatched RNA + CNV + mutation data',
  fillcolor = '#DCEAF5'
]

MOFA [
  label = 'MOFA2 MULTI-OMICS INTEGRATION\\n1,722 RNA + 1,902 CNV + 107 mutation features\\n20 latent factors → Factors 1–10 retained',
  fillcolor = '#DDEFD8'
]

EXPLORE [
  label = 'EXPLORATORY 5-STATE SOLUTION\\nState 1–5\\nApparently distinct molecular substructure',
  fillcolor = '#F8DDDD'
]

ROBUSTNESS [
  label = 'ROBUSTNESS TESTING\\nFactor-number sensitivity + bootstrap validation\\n5-state solution unstable\\nMedian bootstrap ARI ≈ 0.35',
  fillcolor = '#FFF0CF'
]

STATE2 [
  label = 'ROBUST 2-STATE SOLUTION\\nState 1: n = 226   |   State 2: n = 51\\nK=2 bootstrap median ARI ≈ 0.80',
  fillcolor = '#D6EEDB'
]


TCGA -> MOFA
MOFA -> EXPLORE
EXPLORE -> ROBUSTNESS
ROBUSTNESS -> STATE2


# ============================================================
# RNA BRANCH
# ============================================================

RNA [
  label = 'STATE-ASSOCIATED RNA SIGNAL\\n110 genes at FDR < 0.05\\n34 State 1-high   |   76 State 2-high',
  fillcolor = '#E7E0F3'
]

STRING [
  label = 'STRING NETWORK FILTERING\\n37 state-associated genes connected\\nwithin protein-interaction networks',
  fillcolor = '#E7E0F3'
]

LASSO [
  label = 'NETWORK-GUIDED NESTED CV LASSO\\nOuter-fold DE selection + STRING filtering + inner LASSO\\nAUC = 0.763 (95% CI: 0.688–0.838)\\nBalanced accuracy = 0.729',
  fillcolor = '#DCE7F4'
]

CORE9 [
  label = 'STABLE 9-GENE RNA CORE\\nSelected in ≥80% of outer folds\\nADGRL3 · PIWIL1 · SLC15A1 · ZIC5 · GATA2\\nFOXG1 · MAGEA11 · PALM3 · SLCO1A2',
  fillcolor = '#D4ECDF'
]

CORE4 [
  label = 'HIGHEST-STABILITY RNA CORE\\nSelected in 10/10 outer folds\\nADGRL3 · PIWIL1 · SLC15A1 · ZIC5',
  fillcolor = '#BFE2CF'
]


STATE2 -> RNA
RNA -> STRING
STRING -> LASSO
LASSO -> CORE9
CORE9 -> CORE4


# ============================================================
# ALTERNATIVE SPLICING BRANCH
# ============================================================

SPLICE [
  label = 'TCGA SPLICESEQ PSI INTEGRATION\\n273/277 patients matched\\nES · AA · AD · RI · AP · AT · ME',
  fillcolor = '#E5EAF4'
]

GLOBALAS [
  label = 'GLOBAL AS LANDSCAPE\\nModest but significant state association\\nPERMANOVA R² = 0.0066, P = 3 × 10⁻⁴\\nDispersion P = 0.443',
  fillcolor = '#E8E4F3'
]

ASTYPES [
  label = 'AS-TYPE-SPECIFIC ANALYSIS\\nSignificant global associations for\\nES · AA · AD · RI · AP · AT\\nME not significant after BH correction',
  fillcolor = '#E8E4F3'
]

AS34 [
  label = 'INITIAL STRONG AS SIGNAL\\n34 events · 21 unique genes\\nFDR < 0.05 and |ΔPSI| ≥ 0.10',
  fillcolor = '#F6E5CF'
]

AGE [
  label = 'AGE-ADJUSTED PAN-AS ANALYSIS\\nPSI ~ Molecular State + Age\\nState 1 median age = 57.8 y\\nState 2 median age = 66.1 y\\nGlobal BH correction across all AS events',
  fillcolor = '#FFE4B8'
]

AS18 [
  label = 'FINAL AGE-ADJUSTED AS SIGNAL\\n126 events at global FDR < 0.05\\n18 strong events · 13 unique genes\\nGlobal FDR < 0.05 and |adjusted ΔPSI| ≥ 0.10',
  fillcolor = '#D5EDDF'
]

ASGENES [
  label = '13-GENE AS SET\\nCASK · TMEM8B · FLOT2 · ZNF576 · TCF4\\nSLMAP · SIMC1 · ADPRHL1 · SMARCD3\\nFAM86B1 · ZNF559 · DMKN · SPEG',
  fillcolor = '#C5E5D4'
]

ASCHAR [
  label = 'AS SIGNATURE CHARACTERIZATION\\nPC1 = 23.4%   |   PC2 = 10.6%\\nPERMANOVA R² = 0.072, P = 1 × 10⁻⁴\\nDispersion P = 9 × 10⁻⁴\\nState-associated pattern with within-state heterogeneity',
  fillcolor = '#DCE7F4'
]


STATE2 -> SPLICE
SPLICE -> GLOBALAS
GLOBALAS -> ASTYPES
ASTYPES -> AS34
AS34 -> AGE
AGE -> AS18
AS18 -> ASGENES
ASGENES -> ASCHAR


# ============================================================
# DRUG RESPONSE BRANCH
# ============================================================

DEPMAP [
  label = 'DEPMAP HGSOC PROJECTION\\n21 HGSOC cell lines with expression data\\n110 TCGA state-associated genes matched',
  fillcolor = '#F3E5D8'
]

STATESCORE [
  label = 'MOLECULAR STATE PROJECTION\\nContinuous State2-likeness score\\nTCGA RNA state signature projected\\nto HGSOC cell lines',
  fillcolor = '#F3E5D8'
]

DRUGSCREEN [
  label = 'PHARMACOLOGICAL SCREENING\\nGDSC2 + PRISM\\nDrug sensitivity correlated with\\ncontinuous State2-likeness',
  fillcolor = '#F8E7CF'
]

CONCORDANCE [
  label = 'CROSS-PLATFORM CONCORDANCE\\n18 drugs represented in both platforms\\nNo drug significant after multiple-testing correction\\nDirectional concordance used for prioritization',
  fillcolor = '#FFF0CF'
]

GSK [
  label = 'GSK-1904529A\\nCross-platform replicated directional candidate\\nState2-associated sensitivity direction\\nin both GDSC2 and PRISM\\nNot independently validated',
  fillcolor = '#D6EEDB'
]

NONOVERLAP [
  label = 'NON-OVERLAPPING CELL-LINE CHECK\\nGDSC-only: n = 5, rho = -0.70\\nPRISM-only: n = 6, rho ≈ -0.66\\nDirection retained after removing shared models',
  fillcolor = '#DCE7F4'
]

TARGET [
  label = 'TARGET-LEVEL FOLLOW-UP\\nGSK-1904529A → IGF1R / INSR axis\\nNo significant state difference in\\nIGF1R or INSR RNA or CNV',
  fillcolor = '#E7E0F3'
]

GSKGENES [
  label = 'STATE–DRUG RESPONSE FOLLOW-UP\\nExploratory genes with concordant directional\\nrelationships with GSK response\\nGATA2 · ZIC5 · ZIC2 · POU6F2\\nSAMD11 · LIN28B\\nSmall cell-line sample sizes',
  fillcolor = '#E7E0F3'
]


# Drug projection derives from the RNA state signal

RNA -> DEPMAP
DEPMAP -> STATESCORE
STATESCORE -> DRUGSCREEN
DRUGSCREEN -> CONCORDANCE
CONCORDANCE -> GSK
GSK -> NONOVERLAP
NONOVERLAP -> TARGET
TARGET -> GSKGENES


# ============================================================
# ALIGN MAJOR BRANCHES
# ============================================================

{rank=same; RNA; SPLICE}
{rank=same; STRING; GLOBALAS; DEPMAP}
{rank=same; LASSO; ASTYPES; STATESCORE}
{rank=same; CORE9; AS34; DRUGSCREEN}
{rank=same; CORE4; AGE; CONCORDANCE}

}
")


# ============================================================
# DISPLAY IN RSTUDIO
# ============================================================

workflow_complete


# ============================================================
# EXPORT SVG
# ============================================================

svg_code <- export_svg(workflow_complete)

svg_file <- file.path(
  output_dir,
  'Molecular_State_RNA_AS_Drug_Workflow.svg'
)

writeLines(
  svg_code,
  svg_file
)


# ============================================================
# EXPORT HIGH-RESOLUTION PNG
# ============================================================

png_file <- file.path(
  output_dir,
  'Molecular_State_RNA_AS_Drug_Workflow.png'
)

rsvg_png(
  charToRaw(svg_code),
  file = png_file,
  width = 5000,
  height = 5200
)


# ============================================================
# EXPORT PDF
# ============================================================

pdf_file <- file.path(
  output_dir,
  'Molecular_State_RNA_AS_Drug_Workflow.pdf'
)

rsvg_pdf(
  charToRaw(svg_code),
  file = pdf_file
)


# ============================================================
# CONFIRM FILES
# ============================================================

cat(
  '\n============================================\n',
  'COMPLETE WORKFLOW SAVED\n',
  '============================================\n\n',
  'SVG:\n', svg_file, '\n\n',
  'PNG:\n', png_file, '\n\n',
  'PDF:\n', pdf_file, '\n\n',
  '============================================\n'
)
# ============================================================
# CHECK AVAILABLE OBJECTS FOR AS / DRUG ANALYSIS
# ============================================================

# Objects currently in R environment
objs <- ls()

# Search for objects potentially related to:
# splicing / PSI / DepMap / GDSC / PRISM / GSK

grep(
  "splice|psi|depmap|gdsc|prism|gsk|drug|sens|auc|ic50",
  objs,
  ignore.case = TRUE,
  value = TRUE
)
# ============================================================
# FINAL AGE-ADJUSTED AS GENES
# ============================================================

AS_genes_13 <- c(
  "CASK",
  "TMEM8B",
  "FLOT2",
  "ZNF576",
  "TCF4",
  "SLMAP",
  "SIMC1",
  "ADPRHL1",
  "SMARCD3",
  "FAM86B1",
  "ZNF559",
  "DMKN",
  "SPEG"
)

AS_genes_13
# Find data frames containing AS event-related columns

candidate_objects <- ls()

for (x in candidate_objects) {
  
  obj <- try(get(x), silent = TRUE)
  
  if (is.data.frame(obj) || is.matrix(obj)) {
    
    cn <- colnames(obj)
    
    if (!is.null(cn) &&
        any(grepl(
          "PSI|splice|as_id|adjusted|Global_FDR|delta",
          cn,
          ignore.case = TRUE
        ))) {
      
      cat("\n============================\n")
      cat("OBJECT:", x, "\n")
      cat("DIM:", paste(dim(obj), collapse = " x "), "\n")
      cat("COLUMNS:\n")
      print(cn)
    }
  }
}
# ============================================================
# CHECK GSK DATA STRUCTURE
# ============================================================

cat("\n===== GDSC GSK =====\n")
print(dim(gsk_gdsc_unique))
print(colnames(gsk_gdsc_unique))
print(head(gsk_gdsc_unique))

cat("\n===== PRISM GSK =====\n")
print(dim(gsk_prism_unique))
print(colnames(gsk_prism_unique))
print(head(gsk_prism_unique))

cat("\n===== DEPMAP HGSOC =====\n")
print(dim(depmap_hgsoc))
print(colnames(depmap_hgsoc))

cat("\n===== FINAL 18 AS EVENTS =====\n")
print(dim(event_info))
print(
  event_info[, c(
    "Event",
    "AS_Type",
    "Gene",
    "Delta_PSI_adjusted",
    "Global_FDR",
    "Direction"
  )]
)

cat("\n===== FINAL 13 AS GENES =====\n")

AS_genes_13 <- unique(event_info$Gene)

print(AS_genes_13)
cat("\nN genes:", length(AS_genes_13), "\n")
cat("\n===== GDSC GSK =====\n")
print(dim(gsk_gdsc_unique))
print(colnames(gsk_gdsc_unique))
print(gsk_gdsc_unique)

cat("\n===== PRISM GSK =====\n")
print(dim(gsk_prism_unique))
print(colnames(gsk_prism_unique))
print(gsk_prism_unique)
AS_genes_13 <- c(
  "CASK",
  "FLOT2",
  "SLMAP",
  "FAM86B1",
  "DMKN",
  "TMEM8B",
  "SIMC1",
  "ZNF576",
  "TCF4",
  "ADPRHL1",
  "SMARCD3",
  "ZNF559",
  "SPEG"
)
cat("\n===== GDSC GSK =====\n")
print(dim(gsk_gdsc_unique))
print(colnames(gsk_gdsc_unique))
print(gsk_gdsc_unique)

cat("\n===== PRISM GSK =====\n")
print(dim(gsk_prism_unique))
print(colnames(gsk_prism_unique))
print(gsk_prism_unique)
# ============================================================
# 13 AS-ASSOCIATED GENES vs GSK-1904529A RESPONSE
# GDSC + PRISM
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Final 13 AS-associated genes
# ------------------------------------------------------------

AS_genes_13 <- c(
  "CASK",
  "FLOT2",
  "SLMAP",
  "FAM86B1",
  "DMKN",
  "TMEM8B",
  "SIMC1",
  "ZNF576",
  "TCF4",
  "ADPRHL1",
  "SMARCD3",
  "ZNF559",
  "SPEG"
)


# ------------------------------------------------------------
# 2. Find corresponding DepMap expression columns
# Example column format: CASK (8573)
# ------------------------------------------------------------

gene_col_map <- data.frame(
  Gene = AS_genes_13,
  Expression_column = sapply(
    AS_genes_13,
    function(g) {
      
      hit <- grep(
        paste0("^", g, " \\("),
        colnames(expr),
        value = TRUE
      )
      
      if (length(hit) == 0) {
        return(NA_character_)
      }
      
      hit[1]
    }
  ),
  stringsAsFactors = FALSE
)

cat("\n===== GENE MATCHING =====\n")
print(gene_col_map)

cat(
  "\nMatched:",
  sum(!is.na(gene_col_map$Expression_column)),
  "/",
  length(AS_genes_13),
  "\n"
)


# ------------------------------------------------------------
# 3. Keep genes available in DepMap
# ------------------------------------------------------------

gene_col_map <- gene_col_map[
  !is.na(gene_col_map$Expression_column),
]

genes_available <- gene_col_map$Gene


# ------------------------------------------------------------
# 4. Create clean expression table
# ------------------------------------------------------------

expr_AS <- expr[
  ,
  c(
    "ModelID",
    gene_col_map$Expression_column
  ),
  drop = FALSE
]

# Rename expression columns to gene symbols

colnames(expr_AS) <- c(
  "ModelID",
  gene_col_map$Gene
)


# ------------------------------------------------------------
# 5. GDSC: merge expression with GSK sensitivity
# ------------------------------------------------------------

gdsc_AS <- merge(
  gsk_gdsc_unique,
  expr_AS,
  by = "ModelID"
)

cat("\n===== GDSC MATCHED CELL LINES =====\n")
print(gdsc_AS$ModelID)

cat("\nN =", nrow(gdsc_AS), "\n")


# ------------------------------------------------------------
# 6. PRISM: merge expression with GSK sensitivity
# ------------------------------------------------------------

prism_AS <- merge(
  gsk_prism_unique,
  expr_AS,
  by.x = "depmap_id",
  by.y = "ModelID"
)

cat("\n===== PRISM MATCHED CELL LINES =====\n")
print(prism_AS$depmap_id)

cat("\nN =", nrow(prism_AS), "\n")


# ------------------------------------------------------------
# 7. GDSC correlations
#
# Lower IC50 = greater sensitivity
# Therefore negative rho:
# higher gene expression -> greater sensitivity
# ------------------------------------------------------------

gdsc_gene_results <- lapply(
  genes_available,
  function(g) {
    
    x <- gdsc_AS[[g]]
    y <- gdsc_AS$IC50
    
    ok <- complete.cases(x, y)
    
    if (sum(ok) < 4) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        x[ok],
        y[ok],
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
      Gene = g,
      GDSC_N = sum(ok),
      GDSC_rho = unname(test$estimate),
      GDSC_P = test$p.value
    )
  }
)

gdsc_gene_results <- bind_rows(gdsc_gene_results)

gdsc_gene_results$GDSC_FDR <- p.adjust(
  gdsc_gene_results$GDSC_P,
  method = "BH"
)


# ------------------------------------------------------------
# 8. PRISM correlations
#
# Lower AUC = greater sensitivity
# Therefore negative rho:
# higher expression -> greater sensitivity
# ------------------------------------------------------------

prism_gene_results <- lapply(
  genes_available,
  function(g) {
    
    x <- prism_AS[[g]]
    y <- prism_AS$auc
    
    ok <- complete.cases(x, y)
    
    if (sum(ok) < 4) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        x[ok],
        y[ok],
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
      Gene = g,
      PRISM_N = sum(ok),
      PRISM_rho = unname(test$estimate),
      PRISM_P = test$p.value
    )
  }
)

prism_gene_results <- bind_rows(prism_gene_results)

prism_gene_results$PRISM_FDR <- p.adjust(
  prism_gene_results$PRISM_P,
  method = "BH"
)


# ------------------------------------------------------------
# 9. Combine GDSC + PRISM
# ------------------------------------------------------------

AS_GSK_results <- merge(
  gdsc_gene_results,
  prism_gene_results,
  by = "Gene",
  all = TRUE
)


# ------------------------------------------------------------
# 10. Determine cross-platform direction
# ------------------------------------------------------------

AS_GSK_results$Same_direction <- with(
  AS_GSK_results,
  sign(GDSC_rho) == sign(PRISM_rho)
)


# Mean absolute correlation strength

AS_GSK_results$Mean_abs_rho <- rowMeans(
  abs(
    AS_GSK_results[
      ,
      c("GDSC_rho", "PRISM_rho")
    ]
  ),
  na.rm = TRUE
)


# Biological interpretation
# Since lower IC50/AUC = more sensitivity:

AS_GSK_results$Association <- ifelse(
  AS_GSK_results$Same_direction &
    AS_GSK_results$GDSC_rho < 0 &
    AS_GSK_results$PRISM_rho < 0,
  
  "Higher expression -> GSK sensitivity",
  
  ifelse(
    AS_GSK_results$Same_direction &
      AS_GSK_results$GDSC_rho > 0 &
      AS_GSK_results$PRISM_rho > 0,
    
    "Higher expression -> GSK resistance",
    
    "Discordant"
  )
)


# ------------------------------------------------------------
# 11. Sort strongest cross-platform associations first
# ------------------------------------------------------------

AS_GSK_results <- AS_GSK_results[
  order(
    !AS_GSK_results$Same_direction,
    -AS_GSK_results$Mean_abs_rho
  ),
]


# ------------------------------------------------------------
# 12. Print final results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL AS-GENE / GSK RESULTS\n")
cat("============================================\n")

print(
  AS_GSK_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Concordant genes only
# ------------------------------------------------------------

AS_GSK_concordant <- AS_GSK_results[
  AS_GSK_results$Same_direction == TRUE,
]

cat("\n============================================\n")
cat("CROSS-PLATFORM CONCORDANT GENES\n")
cat("============================================\n")

print(
  AS_GSK_concordant,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Specifically sensitivity-associated direction
# ------------------------------------------------------------

AS_GSK_sensitivity <- AS_GSK_results[
  AS_GSK_results$Association ==
    "Higher expression -> GSK sensitivity",
]

cat("\n============================================\n")
cat("CONCORDANT GSK-SENSITIVITY GENES\n")
cat("============================================\n")

print(
  AS_GSK_sensitivity,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Save results
# ------------------------------------------------------------

write.csv(
  AS_GSK_results,
  "~/Desktop/AS13_GSK1904529A_crossplatform_results.csv",
  row.names = FALSE
)

write.csv(
  AS_GSK_concordant,
  "~/Desktop/AS13_GSK1904529A_concordant_genes.csv",
  row.names = FALSE
)

cat(
  "\nSaved to Desktop:\n",
  "AS13_GSK1904529A_crossplatform_results.csv\n",
  "AS13_GSK1904529A_concordant_genes.csv\n"
)

# ============================================================
# 13 AS-ASSOCIATED GENE EXPRESSION PROGRAM
# vs GSK-1904529A RESPONSE
# ============================================================

library(dplyr)
library(ggplot2)

AS_genes_13 <- c(
  "CASK", "FLOT2", "SLMAP", "FAM86B1", "DMKN",
  "TMEM8B", "SIMC1", "ZNF576", "TCF4", "ADPRHL1",
  "SMARCD3", "ZNF559", "SPEG"
)

# ------------------------------------------------------------
# 1. Find DepMap expression columns
# ------------------------------------------------------------

gene_cols <- sapply(
  AS_genes_13,
  function(g) {
    hit <- grep(
      paste0("^", g, " \\("),
      colnames(expr),
      value = TRUE
    )
    
    if (length(hit) == 0) NA_character_ else hit[1]
  }
)

gene_cols <- gene_cols[!is.na(gene_cols)]

cat("Matched genes:", length(gene_cols), "/", length(AS_genes_13), "\n")
print(names(gene_cols))


# ------------------------------------------------------------
# 2. Expression matrix
# ------------------------------------------------------------

expr13 <- expr[
  ,
  c("ModelID", unname(gene_cols)),
  drop = FALSE
]

colnames(expr13) <- c(
  "ModelID",
  names(gene_cols)
)


# ------------------------------------------------------------
# 3. Restrict PCA definition to HGSOC cell lines
#
# Use all available HGSOC models, NOT only the drug-tested lines.
# ------------------------------------------------------------

hgsoc_ids <- unique(depmap_hgsoc$ModelID)

expr13_hgsoc <- expr13[
  expr13$ModelID %in% hgsoc_ids,
]

cat("\nHGSOC models used for program definition:",
    nrow(expr13_hgsoc), "\n")


# ------------------------------------------------------------
# 4. Prepare expression matrix
# ------------------------------------------------------------

X <- as.matrix(
  expr13_hgsoc[, names(gene_cols), drop = FALSE]
)

rownames(X) <- expr13_hgsoc$ModelID

# Remove genes with zero variance
keep_var <- apply(X, 2, sd, na.rm = TRUE) > 0
X <- X[, keep_var, drop = FALSE]

# Median imputation if necessary
for (j in seq_len(ncol(X))) {
  
  X[is.na(X[, j]), j] <- median(
    X[, j],
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# 5. PCA-based AS-gene expression program
# ------------------------------------------------------------

pca_AS13 <- prcomp(
  X,
  center = TRUE,
  scale. = TRUE
)

AS13_score <- pca_AS13$x[, 1]

variance_PC1 <- (
  pca_AS13$sdev[1]^2 /
    sum(pca_AS13$sdev^2)
) * 100

cat("\nPC1 variance explained:",
    round(variance_PC1, 2), "%\n")


# ------------------------------------------------------------
# 6. Inspect PC1 loadings
# ------------------------------------------------------------

AS13_loadings <- data.frame(
  Gene = rownames(pca_AS13$rotation),
  Loading = pca_AS13$rotation[, 1]
)

AS13_loadings <- AS13_loadings[
  order(abs(AS13_loadings$Loading),
        decreasing = TRUE),
]

cat("\n===== PC1 LOADINGS =====\n")
print(AS13_loadings)


# ------------------------------------------------------------
# 7. Score table
# ------------------------------------------------------------

AS13_score_df <- data.frame(
  ModelID = names(AS13_score),
  AS13_PC1 = as.numeric(AS13_score)
)


# ------------------------------------------------------------
# 8. Merge with GDSC
# ------------------------------------------------------------

gdsc_score <- merge(
  gsk_gdsc_unique,
  AS13_score_df,
  by = "ModelID"
)

cat("\n===== GDSC =====\n")
print(gdsc_score)


# ------------------------------------------------------------
# 9. Merge with PRISM
# ------------------------------------------------------------

prism_score <- merge(
  gsk_prism_unique,
  AS13_score_df,
  by.x = "depmap_id",
  by.y = "ModelID"
)

cat("\n===== PRISM =====\n")
print(prism_score)


# ------------------------------------------------------------
# 10. Spearman: AS13 program vs GSK response
# ------------------------------------------------------------

gdsc_test <- cor.test(
  gdsc_score$AS13_PC1,
  gdsc_score$IC50,
  method = "spearman",
  exact = FALSE
)

prism_test <- cor.test(
  prism_score$AS13_PC1,
  prism_score$auc,
  method = "spearman",
  exact = FALSE
)


cat("\n====================================\n")
cat("AS13 PROGRAM vs GSK-1904529A\n")
cat("====================================\n")

cat(
  "\nGDSC\n",
  "N =", nrow(gdsc_score),
  "\nrho =", round(unname(gdsc_test$estimate), 3),
  "\nP =", signif(gdsc_test$p.value, 4),
  "\n"
)

cat(
  "\nPRISM\n",
  "N =", nrow(prism_score),
  "\nrho =", round(unname(prism_test$estimate), 3),
  "\nP =", signif(prism_test$p.value, 4),
  "\n"
)


# ------------------------------------------------------------
# 11. Cross-platform direction
# ------------------------------------------------------------

same_direction <-
  sign(unname(gdsc_test$estimate)) ==
  sign(unname(prism_test$estimate))

cat(
  "\nSame direction across platforms:",
  same_direction,
  "\n"
)


# ------------------------------------------------------------
# 12. Plots
# ------------------------------------------------------------

p_gdsc <- ggplot(
  gdsc_score,
  aes(x = AS13_PC1, y = IC50)
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "GDSC: AS-associated expression program",
    subtitle = paste0(
      "Spearman rho = ",
      round(unname(gdsc_test$estimate), 2),
      ", P = ",
      signif(gdsc_test$p.value, 3)
    ),
    x = "13-gene AS-associated expression program (PC1)",
    y = "GSK-1904529A IC50"
  )

print(p_gdsc)


p_prism <- ggplot(
  prism_score,
  aes(x = AS13_PC1, y = auc)
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "PRISM: AS-associated expression program",
    subtitle = paste0(
      "Spearman rho = ",
      round(unname(prism_test$estimate), 2),
      ", P = ",
      signif(prism_test$p.value, 3)
    ),
    x = "13-gene AS-associated expression program (PC1)",
    y = "GSK-1904529A AUC"
  )

print(p_prism)


# ------------------------------------------------------------
# 13. Save
# ------------------------------------------------------------

write.csv(
  AS13_loadings,
  "~/Desktop/AS13_expression_program_loadings.csv",
  row.names = FALSE
)

write.csv(
  gdsc_score,
  "~/Desktop/AS13_GSK_GDSC_score.csv",
  row.names = FALSE
)

write.csv(
  prism_score,
  "~/Desktop/AS13_GSK_PRISM_score.csv",
  row.names = FALSE
)
cat(
  "GDSC rho =", unname(gdsc_test$estimate),
  "\nGDSC P =", gdsc_test$p.value,
  "\nPRISM rho =", unname(prism_test$estimate),
  "\nPRISM P =", prism_test$p.value,
  "\n"
)


