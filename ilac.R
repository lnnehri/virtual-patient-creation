
> setwd("/Users/lemannur/Downloads")
> model <- read.csv(
  +     "Model.csv",
  +     check.names = FALSE
  + )
> 
  > dim(model)
[1] 2154   49
> expr <- read.csv(
  +     "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  +     check.names = FALSE
  + )
> 
  > dim(expr)
[1]  1775 19221
> library(dplyr)
> 
  > # ============================================================
> # 1. OVARIAN / FALLOPIAN MODELLERI SEC
  > # ============================================================
> 
  > ov_model <- model %>%
    +     dplyr::filter(
      +         grepl(
        +             "ovar|fallopian",
        +             paste(
          +                 OncotreeLineage,
          +                 OncotreePrimaryDisease,
          +                 OncotreeSubtype
          +             ),
        +             ignore.case = TRUE
        +         )
      +     )
  > 
    > cat("Ovarian/Fallopian models:", nrow(ov_model), "\n\n")
  Ovarian/Fallopian models: 76 
  
  > 
    > 
    > # ============================================================
  > # 2. HANGI SUBTYPE'LAR VAR?
    > # ============================================================
  > 
    > print(
      +     ov_model %>%
        +         dplyr::count(
          +             OncotreePrimaryDisease,
          +             OncotreeSubtype,
          +             sort = TRUE
          +         ),
      +     n = 100
      + )
  Error in print.default(m, ..., quote = quote, right = right, max = max) : 
    invalid 'na.print' specification
  
  > subtype_table <- ov_model %>%
    +     dplyr::count(
      +         OncotreePrimaryDisease,
      +         OncotreeSubtype,
      +         sort = TRUE
      +     )
  > 
    > subtype_table
  OncotreePrimaryDisease                   OncotreeSubtype  n
  1  Ovarian Epithelial Tumor  High-Grade Serous Ovarian Cancer 23
  2  Ovarian Epithelial Tumor             Serous Ovarian Cancer 23
  3  Ovarian Epithelial Tumor         Clear Cell Ovarian Cancer  9
  4  Ovarian Epithelial Tumor       Endometrioid Ovarian Cancer  6
  5  Ovarian Epithelial Tumor           Mucinous Ovarian Cancer  6
  6  Ovarian Epithelial Tumor Small Cell Carcinoma of the Ovary  3
  7             Non-Cancerous        Immortalized Ovarian Cells  1
  8     Ovarian Cancer, Other             Ovarian Cancer, Other  1
  9  Ovarian Epithelial Tumor          Brenner Tumor, Malignant  1
  10 Ovarian Epithelial Tumor           Mixed Ovarian Carcinoma  1
  11  Ovarian Germ Cell Tumor             Mixed Germ Cell Tumor  1
  12   Sex Cord Stromal Tumor              Granulosa Cell Tumor  1
  > ov_model_small <- ov_model %>%
    +     dplyr::select(
      +         ModelID,
      +         CellLineName,
      +         OncotreeLineage,
      +         OncotreePrimaryDisease,
      +         OncotreeSubtype,
      +         OncotreeCode
      +     )
  > 
    > expr_ov <- expr %>%
    +     dplyr::filter(ModelID %in% ov_model$ModelID)
  Error in `dplyr::filter()`:
    ! Can't transform a data frame with `NA` or `""` names.
Run `rlang::last_trace()` to see where the error occurred.

> # Boş veya NA sütun isimlerini kontrol et
> which(is.na(colnames(expr)) | colnames(expr) == "")
[1] 1
> # İsimsiz sütunları çıkar
> expr <- expr[
+     ,
+     !is.na(colnames(expr)) & colnames(expr) != "",
+     drop = FALSE
+ ]
> 
> # Kontrol
> dim(expr)
[1]  1775 19220
> any(is.na(colnames(expr)) | colnames(expr) == "")
[1] FALSE
> expr_ov <- expr[
+     expr$ModelID %in% ov_model$ModelID,
+     ,
+     drop = FALSE
+ ]
> 
> cat("Ovarian/Fallopian models:", nrow(ov_model), "\n")
Ovarian/Fallopian models: 76 
> cat("Ovarian expression rows:", nrow(expr_ov), "\n")
Ovarian expression rows: 67 
> cat("Unique ovarian ModelIDs:", length(unique(expr_ov$ModelID)), "\n")
Unique ovarian ModelIDs: 67 
> hgsoc_ids <- ov_model$ModelID[
+     ov_model$OncotreeSubtype == "High-Grade Serous Ovarian Cancer"
+ ]
> 
> cat("HGSOC metadata models:", length(hgsoc_ids), "\n")
HGSOC metadata models: 23 
> cat(
+     "HGSOC with expression:",
+     length(unique(expr_ov$ModelID[expr_ov$ModelID %in% hgsoc_ids])),
+     "\n"
+ )
HGSOC with expression: 21 

###sonrasi
# ============================================================
# 1. HGSOC EXPRESSION DATA
# ============================================================

expr_hgsoc <- expr_ov[
  expr_ov$ModelID %in% hgsoc_ids,
  ,
  drop = FALSE
]

cat("HGSOC expression profiles:", nrow(expr_hgsoc), "\n")


# ============================================================
# 2. EXPRESSION GENE COLUMNS
# ============================================================

metadata_cols <- intersect(
  c("ProfileID", "is_default_entry", "ModelID"),
  colnames(expr_hgsoc)
)

gene_cols <- setdiff(
  colnames(expr_hgsoc),
  metadata_cols
)

# DepMap gene names may be:
# TP53 (7157)
# Remove Entrez ID in parentheses
depmap_gene_symbols <- sub(
  " \\([0-9]+\\)$",
  "",
  gene_cols
)

cat("DepMap protein-coding genes:", length(depmap_gene_symbols), "\n")


# ============================================================
# 3. TCGA 110-GENE SIGNATURE
# ============================================================

# gsea_stats was generated previously from our TCGA state analysis

tcga_signature <- gsea_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(
    Gene,
    Difference,
    Z,
    FDR
  ) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

cat("TCGA significant genes:", nrow(tcga_signature), "\n")

cat(
  "State1-high:",
  sum(tcga_signature$Difference < 0),
  "\n"
)

cat(
  "State2-high:",
  sum(tcga_signature$Difference > 0),
  "\n"
)


# ============================================================
# 4. MATCH WITH DEPMAP
# ============================================================

common_genes <- intersect(
  tcga_signature$Gene,
  depmap_gene_symbols
)

cat(
  "\nSignature genes found in DepMap:",
  length(common_genes),
  "/",
  nrow(tcga_signature),
  "\n"
)


# Missing genes
missing_genes <- setdiff(
  tcga_signature$Gene,
  depmap_gene_symbols
)

cat("\nMissing genes:\n")
print(missing_genes)


# ============================================================
# 5. SHOW MATCHED GENES
# ============================================================

matched_signature <- tcga_signature %>%
  dplyr::filter(Gene %in% common_genes) %>%
  dplyr::arrange(Difference)

print(matched_signature)
####
# ============================================================
# 1. DEPMAP HGSOC EXPRESSION MATRIX
# ============================================================

# Gene columns
metadata_cols <- intersect(
  c("ProfileID", "is_default_entry", "ModelID"),
  colnames(expr_hgsoc)
)

gene_cols <- setdiff(
  colnames(expr_hgsoc),
  metadata_cols
)

# Clean gene names
depmap_gene_symbols <- sub(
  " \\([0-9]+\\)$",
  "",
  gene_cols
)

# expression matrix:
# rows = genes
# columns = HGSOC cell lines

hgsoc_expr <- t(
  as.matrix(
    expr_hgsoc[, gene_cols, drop = FALSE]
  )
)

storage.mode(hgsoc_expr) <- "numeric"

rownames(hgsoc_expr) <- depmap_gene_symbols
colnames(hgsoc_expr) <- expr_hgsoc$ModelID


# ============================================================
# 2. KEEP THE 110 TCGA STATE GENES
# ============================================================

common_genes <- intersect(
  tcga_signature$Gene,
  rownames(hgsoc_expr)
)

hgsoc_sig <- hgsoc_expr[
  common_genes,
  ,
  drop = FALSE
]

cat("Signature genes used:", nrow(hgsoc_sig), "\n")
cat("HGSOC cell lines:", ncol(hgsoc_sig), "\n")


# ============================================================
# 3. Z-SCORE EACH GENE ACROSS HGSOC CELL LINES
# ============================================================

hgsoc_sig_z <- t(
  scale(
    t(hgsoc_sig)
  )
)

# Remove genes with zero variance
keep <- apply(
  hgsoc_sig_z,
  1,
  function(x) all(is.finite(x))
)

hgsoc_sig_z <- hgsoc_sig_z[
  keep,
  ,
  drop = FALSE
]

cat(
  "Genes after variance filtering:",
  nrow(hgsoc_sig_z),
  "\n"
)


# ============================================================
# 4. TCGA WEIGHTS
#
# Positive Z = State2
# Negative Z = State1
# ============================================================

weights <- tcga_signature$Z[
  match(
    rownames(hgsoc_sig_z),
    tcga_signature$Gene
  )
]

names(weights) <- rownames(hgsoc_sig_z)

# normalize weights
weights <- weights / sum(abs(weights))


# ============================================================
# 5. CALCULATE STATE2-LIKENESS SCORE
# ============================================================

state2_score <- as.numeric(
  crossprod(
    weights,
    hgsoc_sig_z
  )
)

names(state2_score) <- colnames(hgsoc_sig_z)


# ============================================================
# 6. CREATE RESULT TABLE
# ============================================================

hgsoc_state <- data.frame(
  ModelID = names(state2_score),
  State2_Likeness = state2_score,
  stringsAsFactors = FALSE
)

hgsoc_state <- hgsoc_state %>%
  dplyr::left_join(
    ov_model %>%
      dplyr::select(
        ModelID,
        CellLineName,
        OncotreeSubtype
      ),
    by = "ModelID"
  ) %>%
  dplyr::arrange(State2_Likeness)


# ============================================================
# 7. DEFINE EXTREME GROUPS
#
# Bottom 30% = State1-like
# Top 30%    = State2-like
# Middle     = Intermediate
# ============================================================

low_cut <- quantile(
  hgsoc_state$State2_Likeness,
  0.30,
  na.rm = TRUE
)

high_cut <- quantile(
  hgsoc_state$State2_Likeness,
  0.70,
  na.rm = TRUE
)

hgsoc_state <- hgsoc_state %>%
  dplyr::mutate(
    State_Like = dplyr::case_when(
      State2_Likeness <= low_cut ~ "State1-like",
      State2_Likeness >= high_cut ~ "State2-like",
      TRUE ~ "Intermediate"
    )
  )


# ============================================================
# 8. RESULTS
# ============================================================

cat("\n===== GROUP SIZES =====\n")
print(table(hgsoc_state$State_Like))

cat("\n===== CELL LINES =====\n")
print(
  hgsoc_state %>%
    dplyr::select(
      CellLineName,
      ModelID,
      State2_Likeness,
      State_Like
    )
)


# ============================================================
# 9. VISUALIZATION
# ============================================================

library(ggplot2)

plot_df <- hgsoc_state %>%
  dplyr::arrange(State2_Likeness) %>%
  dplyr::mutate(
    CellLineName = factor(
      CellLineName,
      levels = CellLineName
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = CellLineName,
    y = State2_Likeness,
    fill = State_Like
  )
) +
  geom_col(width = 0.75) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +

  coord_flip() +

  labs(
    title = "Projection of TCGA-OV Molecular States into HGSOC Cell Lines",
    subtitle = "TCGA-derived 110-gene transcriptional signature",
    x = NULL,
    y = "State2-likeness score",
    fill = "Projected state"
  ) +

  theme_classic(base_size = 12)

print(p)


# ============================================================
# 10. SAVE
# ============================================================

write.csv(
  hgsoc_state,
  "HGSOC_TCGA_State_Projection.csv",
  row.names = FALSE
)

ggsave(
  "HGSOC_TCGA_State_Projection.pdf",
  p,
  width = 8,
  height = 7
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PharmacoGx")

library(PharmacoGx)

availablePSets()
gdsc2 <- PharmacoGx::downloadPSet(
  name = "GDSC_2020(v2-8.2)"
)
gdsc2
class(gdsc2)
PharmacoGx::cellNames(gdsc2)[1:20]

length(PharmacoGx::cellNames(gdsc2))

PharmacoGx::drugNames(gdsc2)[1:20]

length(PharmacoGx::drugNames(gdsc2))
hgsoc_names <- hgsoc_state$CellLineName

hgsoc_names
# ============================================================
# GDSC2 SAMPLE METADATA
# ============================================================

gdsc_samples <- PharmacoGx::sampleInfo(gdsc2)

dim(gdsc_samples)
colnames(gdsc_samples)
head(gdsc_samples)
# ============================================================
# GDSC2 DRUG METADATA
# ============================================================

gdsc_drugs <- PharmacoGx::treatmentInfo(gdsc2)

dim(gdsc_drugs)
colnames(gdsc_drugs)
head(gdsc_drugs)
#matching
library(dplyr)
library(stringr)

# ============================================================
# 1. NAME NORMALIZATION FUNCTION
# ============================================================

clean_cell_name <- function(x) {
  x <- toupper(x)
  x <- gsub("NIH:", "", x)
  x <- gsub("NCI-", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}


# ============================================================
# 2. DEPMAP HGSOC NAMES
# ============================================================

depmap_hgsoc <- hgsoc_state %>%
  dplyr::mutate(
    CleanName = clean_cell_name(CellLineName)
  )

# ============================================================
# 3. GDSC NAMES
# ============================================================

gdsc_match <- gdsc_samples %>%
  dplyr::mutate(
    GDSC_Name = Sample.Name,
    CleanName = clean_cell_name(Sample.Name)
  )


# ============================================================
# 4. MATCH
# ============================================================

cell_match <- depmap_hgsoc %>%
  dplyr::left_join(
    gdsc_match %>%
      dplyr::select(
        GDSC_Name,
        sampleid,
        CleanName,
        COSMIC.identifier,
        GDSC..Tissue.descriptor.1,
        GDSC..Tissue..descriptor.2,
        Cancer.Type...matching.TCGA.label.
      ),
    by = "CleanName"
  )


# ============================================================
# 5. RESULTS
# ============================================================

cat(
  "Total HGSOC lines:",
  nrow(depmap_hgsoc),
  "\n"
)

cat(
  "Matched to GDSC:",
  sum(!is.na(cell_match$GDSC_Name)),
  "\n\n"
)

cat("===== MATCHED =====\n")

print(
  cell_match %>%
    dplyr::filter(!is.na(GDSC_Name)) %>%
    dplyr::select(
      CellLineName,
      GDSC_Name,
      State2_Likeness,
      State_Like,
      COSMIC.identifier
    )
)


cat("\n===== NOT MATCHED =====\n")

print(
  cell_match %>%
    dplyr::filter(is.na(GDSC_Name)) %>%
    dplyr::select(
      CellLineName,
      ModelID,
      State2_Likeness,
      State_Like
    )
)
# ============================================================
# CHECK AVAILABLE GDSC2 SENSITIVITY MEASURES
# ============================================================

sens_measures <- PharmacoGx::sensitivityMeasures(gdsc2)

sens_measures
# Sensitivity experiment structure
sn <- PharmacoGx::sensNumber(gdsc2)

dim(sn)
sn[1:5, 1:5]
grep(
  "sens",
  getNamespaceExports("PharmacoGx"),
  value = TRUE,
  ignore.case = TRUE
)
sens <- PharmacoGx::sensitivityProfiles(gdsc2)

class(sens)
dim(sens)
dimnames(sens)
# GDSC sensitivity experiments -> cell line x drug summary

gdsc2_summary <- PharmacoGx::summarizeSensitivityProfiles(
    gdsc2,
    sensitivity.measure = "ic50_recomputed"
)

class(gdsc2_summary)
dim(gdsc2_summary)
library(dplyr)

# ============================================================
# 1. CHECK CELL LINE MATCH
# ============================================================

common_lines <- intersect(
    cell_match$GDSC_Name[!is.na(cell_match$GDSC_Name)],
    colnames(gdsc2_summary)
)

cat("HGSOC lines found in IC50 matrix:", length(common_lines), "\n")
print(common_lines)


# ============================================================
# 2. STATE2 LIKENESS VECTOR
# ============================================================

matched_hgsoc <- cell_match %>%
    dplyr::filter(!is.na(GDSC_Name)) %>%
    dplyr::distinct(GDSC_Name, .keep_all = TRUE)

state_vector <- matched_hgsoc$State2_Likeness[
    match(common_lines, matched_hgsoc$GDSC_Name)
]

names(state_vector) <- common_lines

print(state_vector)


# ============================================================
# 3. SPEARMAN: STATE2 LIKENESS vs IC50
# ============================================================

drug_results <- lapply(
    rownames(gdsc2_summary),
    function(drug) {

        ic50_values <- gdsc2_summary[
            drug,
            common_lines
        ]

        ok <- is.finite(ic50_values) &
              is.finite(state_vector)

        n <- sum(ok)

        # Require >= 6 cell lines
        if (n < 6)
            return(NULL)

        test <- suppressWarnings(
            cor.test(
                state_vector[ok],
                ic50_values[ok],
                method = "spearman",
                exact = FALSE
            )
        )

        data.frame(
            Drug = drug,
            N = n,
            Rho = unname(test$estimate),
            P_value = test$p.value
        )
    }
)

drug_results <- dplyr::bind_rows(drug_results)


# ============================================================
# 4. FDR
# ============================================================

drug_results$FDR <- p.adjust(
    drug_results$P_value,
    method = "BH"
)


# ============================================================
# 5. INTERPRET DIRECTION
#
# lower IC50 = greater sensitivity
#
# rho < 0:
# State2_Likeness increases + IC50 decreases
# -> State2-sensitive
#
# rho > 0:
# State2_Likeness increases + IC50 increases
# -> State1-sensitive
# ============================================================

drug_results <- drug_results %>%
    dplyr::mutate(
        Direction = dplyr::case_when(
            Rho < 0 ~ "State2-sensitive",
            Rho > 0 ~ "State1-sensitive",
            TRUE ~ "No direction"
        )
    ) %>%
    dplyr::arrange(P_value)


# ============================================================
# 6. SUMMARY
# ============================================================

cat("\nDrugs tested:", nrow(drug_results), "\n")
cat("FDR < 0.05:", sum(drug_results$FDR < 0.05), "\n")
cat("Nominal P < 0.05:", sum(drug_results$P_value < 0.05), "\n")


cat("\n===== TOP STATE2-SENSITIVE =====\n")

print(
    drug_results %>%
        dplyr::filter(Rho < 0) %>%
        dplyr::arrange(P_value) %>%
        head(15)
)


cat("\n===== TOP STATE1-SENSITIVE =====\n")

print(
    drug_results %>%
        dplyr::filter(Rho > 0) %>%
        dplyr::arrange(P_value) %>%
        head(15)
)


# ============================================================
# 7. SAVE
# ============================================================

write.csv(
    drug_results,
    "HGSOC_State2Likeness_GDSC2_IC50.csv",
    row.names = FALSE
)
library(dplyr)
library(ggplot2)
library(ggrepel)

plot_drugs <- drug_results %>%
    mutate(
        neglog10P = -log10(P_value),

        Highlight = case_when(
            P_value < 0.05 & Rho < 0 ~ "State2-sensitive",
            P_value < 0.05 & Rho > 0 ~ "State1-sensitive",
            TRUE ~ "NS"
        ),

        Label = ifelse(
            P_value < 0.05 | abs(Rho) >= 0.60,
            Drug,
            NA
        )
    )

p_drug <- ggplot(
    plot_drugs,
    aes(
        x = Rho,
        y = neglog10P,
        color = Highlight
    )
) +
    geom_point(
        aes(size = N),
        alpha = 0.75
    ) +

    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +

    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dotted"
    ) +

    ggrepel::geom_text_repel(
        aes(label = Label),
        size = 3.5,
        max.overlaps = 20
    ) +

    scale_color_manual(
        values = c(
            "State1-sensitive" = "#377EB8",
            "State2-sensitive" = "#E41A1C",
            "NS" = "grey70"
        )
    ) +

    labs(
        title = "Association of TCGA-derived Molecular State with Drug Sensitivity",
        subtitle = "GDSC2 HGSOC cell lines",
        x = "Spearman rho: State2-likeness vs IC50",
        y = "-log10(P-value)",
        color = NULL,
        size = "Cell lines"
    ) +

    theme_classic(base_size = 12)

print(p_drug)

ggsave(
    "HGSOC_State_DrugSensitivity_GDSC2.pdf",
    p_drug,
    width = 8,
    height = 6
)
##bundan bi sey cikmadi prism e bakiyoruz
# Eğer depmap kurulu değilse:
if (!requireNamespace("depmap", quietly = TRUE)) {
    BiocManager::install("depmap")
}

library(depmap)

# Paket sürümü
packageVersion("depmap")

# depmap içindeki fonksiyonlar
grep(
    "file|download|get",
    getNamespaceExports("depmap"),
    value = TRUE,
    ignore.case = TRUE
)
library(depmap)

files <- dmfiles()

class(files)
dim(files)
colnames(files)

# PRISM geçen kayıtları bul
prism_files <- files[
    apply(files, 1, function(x)
        any(grepl("PRISM|prism", x, ignore.case = TRUE))
    ),
]

prism_files
library(depmap)

# Aradığımız dosyanın kaydını seç
prism_row <- prism_files[
    prism_files$name ==
    "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
]

prism_row

# Gerçek download URL
prism_url <- prism_row$download_url

# Kaydedilecek yer
prism_file <- "/Users/lemannur/Downloads/PRISM_secondary_20Q2.csv"

# İndir
download.file(
    url = prism_url,
    destfile = prism_file,
    mode = "wb"
)

# Dosya gerçekten yaklaşık 290 MB mı kontrol et
file.info(prism_file)$size / 1024^2
library(data.table)

prism <- fread(
    "/Users/lemannur/Downloads/PRISM_secondary_20Q2.csv"
)

dim(prism)
colnames(prism)

# İlk birkaç satır
head(prism)
library(dplyr)

prism_match <- state_projection %>%
    dplyr::filter(ModelID %in% prism$depmap_id) %>%
    dplyr::select(
        ModelID,
        CellLineName,
        State2_Likeness,
        State_Like
    )

cat("HGSOC models:", nrow(state_projection), "\n")
cat("Matched to PRISM:", nrow(prism_match), "\n\n")

prism_match %>%
    dplyr::arrange(State2_Likeness)
    
library(dplyr)

# ============================================================
# 1. PRISM'i 13 HGSOC modeline indir
# ============================================================

prism_hgsoc <- prism %>%
  dplyr::filter(depmap_id %in% prism_match$ModelID) %>%
  dplyr::left_join(
    prism_match %>%
      dplyr::select(
        ModelID,
        State2_Likeness,
        State_Like
      ),
    by = c("depmap_id" = "ModelID")
  )

cat("PRISM rows:", nrow(prism_hgsoc), "\n")
cat("Cell lines:", dplyr::n_distinct(prism_hgsoc$depmap_id), "\n")
cat("Drugs:", dplyr::n_distinct(prism_hgsoc$name), "\n")


# ============================================================
# 2. Her drug için State2_Likeness vs AUC
#
# düşük AUC = daha fazla sensitivity
#
# rho < 0 = State2-sensitive
# rho > 0 = State1-sensitive
# ============================================================

prism_results <- prism_hgsoc %>%
  dplyr::filter(
    is.finite(auc),
    is.finite(State2_Likeness)
  ) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(dplyr::n_distinct(depmap_id) >= 6) %>%
  dplyr::summarise(
    N = dplyr::n_distinct(depmap_id),
    
    Rho = suppressWarnings(
      cor(
        State2_Likeness,
        auc,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    
    P_value = suppressWarnings(
      cor.test(
        State2_Likeness,
        auc,
        method = "spearman",
        exact = FALSE
      )$p.value
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(Rho), is.finite(P_value)) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH"),
    
    Direction = dplyr::case_when(
      Rho < 0 ~ "State2-sensitive",
      Rho > 0 ~ "State1-sensitive",
      TRUE ~ "No direction"
    )
  ) %>%
  dplyr::arrange(P_value)


# ============================================================
# 3. Özet
# ============================================================

cat("\nDrugs tested:", nrow(prism_results), "\n")
cat("P < 0.05:", sum(prism_results$P_value < 0.05), "\n")
cat("FDR < 0.05:", sum(prism_results$FDR < 0.05), "\n\n")

head(prism_results, 30)


# ============================================================
# 4. State2 yönündeki en güçlü ilaçlar
# ============================================================

prism_results %>%
  dplyr::filter(Rho < 0) %>%
  dplyr::arrange(P_value) %>%
  head(20)


# ============================================================
# 5. State1 yönündeki en güçlü ilaçlar
# ============================================================

prism_results %>%
  dplyr::filter(Rho > 0) %>%
  dplyr::arrange(P_value) %>%
  head(20)


# ============================================================
# 6. Kaydet
# ============================================================

write.csv(
  prism_results,
  "/Users/lemannur/Downloads/HGSOC_State2Likeness_PRISM_AUC.csv",
  row.names = FALSE
)
sig_prism <- prism_results %>%
  dplyr::filter(FDR < 0.05)

sig_prism

prism_hgsoc %>%
  dplyr::filter(name %in% sig_prism$name) %>%
  dplyr::select(
    name,
    depmap_id,
    State2_Likeness,
    auc
  ) %>%
  dplyr::arrange(name, State2_Likeness)
##bastan
library(dplyr)

# ============================================================
# 1. Drug × cell line duplicate ölçümleri collapse et
# ============================================================

prism_collapsed <- prism_hgsoc %>%
  dplyr::filter(
    is.finite(auc),
    is.finite(State2_Likeness)
  ) %>%
  dplyr::group_by(
    name,
    depmap_id,
    State2_Likeness
  ) %>%
  dplyr::summarise(
    auc = median(auc, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 2. Spearman'ı yeniden hesapla
# ============================================================

prism_results_v2 <- prism_collapsed %>%
  dplyr::group_by(name) %>%
  dplyr::filter(dplyr::n() >= 6) %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    Rho = cor(
      State2_Likeness,
      auc,
      method = "spearman"
    ),
    
    P_value = cor.test(
      State2_Likeness,
      auc,
      method = "spearman",
      exact = TRUE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH"),
    
    Direction = dplyr::case_when(
      Rho < 0 ~ "State2-sensitive",
      Rho > 0 ~ "State1-sensitive",
      TRUE ~ "No direction"
    )
  ) %>%
  dplyr::arrange(P_value)


# ============================================================
# 3. Sonuç
# ============================================================

cat("Drugs tested:", nrow(prism_results_v2), "\n")
cat("P < 0.05:",
    sum(prism_results_v2$P_value < 0.05), "\n")
cat("FDR < 0.05:",
    sum(prism_results_v2$FDR < 0.05), "\n\n")

print(
  prism_results_v2 %>%
    dplyr::filter(FDR < 0.05),
  n = Inf
)

# İlk 20
print(
  prism_results_v2 %>%
    dplyr::arrange(P_value) %>%
    head(20),
  n = 20
)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

# ============================================================
# FIGURE A — Global PRISM drug sensitivity landscape
# ============================================================

label_drugs <- c(
  "NSC-697923",
  "cholecalciferol",
  "pilaralisib",
  "dolutegravir",
  "satraplatin",
  "tepoxalin",
  "niraparib",
  "linsitinib",
  "resatorvid"
)

plot_df <- prism_results_v2 %>%
  dplyr::mutate(
    neglog10P = -log10(P_value),
    
    Significance = dplyr::case_when(
      FDR < 0.05 & Rho < 0 ~ "FDR significant: State2",
      FDR < 0.05 & Rho > 0 ~ "FDR significant: State1",
      P_value < 0.05 & Rho < 0 ~ "Nominal: State2",
      P_value < 0.05 & Rho > 0 ~ "Nominal: State1",
      TRUE ~ "Not significant"
    ),
    
    Label = ifelse(
      name %in% label_drugs,
      name,
      NA
    )
  )

p1 <- ggplot(
  plot_df,
  aes(
    x = Rho,
    y = neglog10P
  )
) +
  
  geom_point(
    aes(
      color = Significance,
      size = N
    ),
    alpha = 0.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted"
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3.4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  
  scale_color_manual(
    values = c(
      "FDR significant: State2" = "#B2182B",
      "FDR significant: State1" = "#2166AC",
      "Nominal: State2" = "#EF8A62",
      "Nominal: State1" = "#67A9CF",
      "Not significant" = "grey75"
    )
  ) +
  
  scale_size_continuous(
    range = c(2, 5)
  ) +
  
  labs(
    title = "PRISM Drug Sensitivity across the TCGA-derived HGSOC State Axis",
    subtitle = "Spearman association between State2-likeness and drug-response AUC",
    x = "Spearman rho",
    y = expression(-log[10](P)),
    color = NULL,
    size = "Cell lines"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p1)

ggsave(
  "/Users/lemannur/Downloads/PRISM_State_Drug_Landscape.pdf",
  p1,
  width = 9,
  height = 6.5
)


# ============================================================
# FIGURE B — Individual cell-line responses
# ============================================================

selected_drugs <- c(
  "NSC-697923",
  "cholecalciferol",
  "pilaralisib",
  "tepoxalin",
  "niraparib",
  "linsitinib"
)

selected_df <- prism_collapsed %>%
  dplyr::filter(name %in% selected_drugs) %>%
  dplyr::left_join(
    prism_match %>%
      dplyr::select(
        ModelID,
        CellLineName,
        State_Like
      ),
    by = c("depmap_id" = "ModelID")
  )

p2 <- ggplot(
  selected_df,
  aes(
    x = State2_Likeness,
    y = auc
  )
) +
  
  geom_point(
    aes(shape = State_Like),
    size = 3
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7,
    linetype = "dashed"
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = CellLineName),
    size = 2.5,
    max.overlaps = 30
  ) +
  
  facet_wrap(
    ~ name,
    scales = "free_y",
    ncol = 3
  ) +
  
  labs(
    title = "Selected PRISM Drug Responses along the HGSOC Molecular State Axis",
    subtitle = "Lower AUC indicates greater drug sensitivity",
    x = "TCGA-derived State2-likeness score",
    y = "PRISM AUC",
    shape = "Molecular state"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p2)

ggsave(
  "/Users/lemannur/Downloads/PRISM_Selected_Drug_Response.pdf",
  p2,
  width = 10,
  height = 7
)
library(dplyr)
library(ggplot2)
library(ggrepel)

selected_drugs <- c(
  "NSC-697923",
  "cholecalciferol",
  "pilaralisib",
  "tepoxalin",
  "niraparib",
  "linsitinib"
)

# ============================================================
# Panel statistics
# ============================================================

drug_stats <- prism_results_v2 %>%
  dplyr::filter(name %in% selected_drugs) %>%
  dplyr::mutate(
    stat_label = paste0(
      "rho = ", sprintf("%.2f", Rho),
      "\np = ", format.pval(P_value, digits = 2),
      "\nFDR = ", sprintf("%.3f", FDR),
      "\nN = ", N
    )
  ) %>%
  dplyr::select(name, stat_label)

# ============================================================
# Plot data
# ============================================================

selected_df <- prism_collapsed %>%
  dplyr::filter(name %in% selected_drugs) %>%
  dplyr::left_join(
    prism_match %>%
      dplyr::select(
        ModelID,
        CellLineName,
        State_Like
      ),
    by = c("depmap_id" = "ModelID")
  )

# ============================================================
# Figure
# ============================================================

p_prism <- ggplot(
  selected_df,
  aes(
    x = State2_Likeness,
    y = auc
  )
) +
  
  geom_point(
    aes(shape = State_Like),
    size = 3.2
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.65,
    linetype = "dashed"
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = CellLineName),
    size = 2.5,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15
  ) +
  
  geom_text(
    data = drug_stats,
    aes(
      x = -Inf,
      y = Inf,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.1,
    size = 3.1
  ) +
  
  facet_wrap(
    ~ name,
    scales = "free_y",
    ncol = 3
  ) +
  
  labs(
    title = "PRISM Drug Responses along the HGSOC Molecular State Axis",
    subtitle = "Lower AUC indicates greater drug sensitivity",
    x = "TCGA-derived State2-likeness score",
    y = "PRISM AUC",
    shape = "Molecular state"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

print(p_prism)

ggsave(
  "/Users/lemannur/Downloads/PRISM_Selected_Drugs_Final.pdf",
  p_prism,
  width = 11,
  height = 7.5
)

ggsave(
  "/Users/lemannur/Downloads/PRISM_Selected_Drugs_Final.png",
  p_prism,
  width = 11,
  height = 7.5,
  dpi = 400
)
library(dplyr)

# ============================================================
# GDSC2 × PRISM CONCORDANCE
# ============================================================

clean_drug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]", "", x)
  return(x)
}

# -----------------------------
# GDSC
# -----------------------------

gdsc_compare <- drug_results %>%
  dplyr::mutate(
    Drug_clean = clean_drug(Drug)
  ) %>%
  dplyr::select(
    Drug_GDSC = Drug,
    Drug_clean,
    N_GDSC = N,
    Rho_GDSC = Rho,
    P_GDSC = P_value,
    FDR_GDSC = FDR
  )

# -----------------------------
# PRISM
# -----------------------------

prism_compare <- prism_results_v2 %>%
  dplyr::mutate(
    Drug_clean = clean_drug(name)
  ) %>%
  dplyr::select(
    Drug_PRISM = name,
    Drug_clean,
    N_PRISM = N,
    Rho_PRISM = Rho,
    P_PRISM = P_value,
    FDR_PRISM = FDR
  )

# -----------------------------
# Match common drugs
# -----------------------------

drug_concordance <- dplyr::inner_join(
  gdsc_compare,
  prism_compare,
  by = "Drug_clean"
) %>%
  dplyr::mutate(
    
    Concordance = dplyr::case_when(
      
      Rho_GDSC < 0 & Rho_PRISM < 0 ~
        "Concordant State2",
      
      Rho_GDSC > 0 & Rho_PRISM > 0 ~
        "Concordant State1",
      
      TRUE ~
        "Discordant"
    ),
    
    Mean_abs_Rho =
      (abs(Rho_GDSC) + abs(Rho_PRISM)) / 2
  )

# -----------------------------
# State2 concordant
# -----------------------------

state2_concordant <- drug_concordance %>%
  dplyr::filter(
    Concordance == "Concordant State2"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_abs_Rho)
  )

# -----------------------------
# State1 concordant
# -----------------------------

state1_concordant <- drug_concordance %>%
  dplyr::filter(
    Concordance == "Concordant State1"
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_abs_Rho)
  )

# ============================================================
# RESULTS
# ============================================================

cat(
  "Common drugs:",
  nrow(drug_concordance),
  "\n\n"
)

print(
  table(drug_concordance$Concordance)
)

cat("\n--- STATE2 CONCORDANT ---\n")

print(
  state2_concordant,
  n = Inf
)

cat("\n--- STATE1 CONCORDANT ---\n")

print(
  state1_concordant,
  n = Inf
)
# State2 concordant ilaçlar
as.data.frame(state2_concordant)

# State1 concordant ilaçlar
as.data.frame(state1_concordant)
library(ggplot2)
library(ggrepel)
library(dplyr)

plot_concordance <- drug_concordance %>%
  dplyr::mutate(
    Label = ifelse(
      Mean_abs_Rho >= 0.30,
      Drug_PRISM,
      NA
    )
  )

p_concordance <- ggplot(
  plot_concordance,
  aes(
    x = Rho_GDSC,
    y = Rho_PRISM
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      color = Concordance,
      size = Mean_abs_Rho
    ),
    alpha = 0.85
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3.5,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Concordant State2" = "#B2182B",
      "Concordant State1" = "#2166AC",
      "Discordant" = "grey65"
    )
  ) +
  
  coord_cartesian(
    xlim = c(-1, 1),
    ylim = c(-1, 1)
  ) +
  
  labs(
    title = "Cross-platform Concordance of Drug Sensitivity",
    subtitle = "GDSC2 and PRISM associations with TCGA-derived HGSOC State2-likeness",
    x = "GDSC2 Spearman rho",
    y = "PRISM Spearman rho",
    color = NULL,
    size = "Mean |rho|"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_concordance)

ggsave(
  "/Users/lemannur/Downloads/GDSC_PRISM_Concordance.pdf",
  p_concordance,
  width = 8,
  height = 7
)

ggsave(
  "/Users/lemannur/Downloads/GDSC_PRISM_Concordance.png",
  p_concordance,
  width = 8,
  height = 7,
  dpi = 400
)# ============================================================
# GDSC × PRISM CELL-LINE OVERLAP
# ============================================================

library(dplyr)

# GDSC'deki 13 eşleşmiş HGSOC line
gdsc_hgsoc_names <- gdsc_match$Sample.Name

# Önce objenin mevcut olduğunu kontrol
cat("GDSC matched lines:", length(gdsc_hgsoc_names), "\n")
print(gdsc_hgsoc_names)

# Concordant State2 drugs
candidate_drugs <- state2_concordant$Drug_GDSC

overlap_results <- lapply(candidate_drugs, function(drug) {
  
  # -------------------------
  # GDSC available lines
  # -------------------------
  
  gdsc_values <- gdsc2_summary[
    drug,
    intersect(
      colnames(gdsc2_summary),
      gdsc_hgsoc_names
    )
  ]
  
  gdsc_lines <- names(gdsc_values)[
    is.finite(gdsc_values)
  ]
  
  
  # -------------------------
  # PRISM available lines
  # -------------------------
  
  prism_drug <- state2_concordant %>%
    dplyr::filter(Drug_GDSC == drug) %>%
    dplyr::pull(Drug_PRISM)
  
  prism_lines <- prism_collapsed %>%
    dplyr::filter(name == prism_drug) %>%
    dplyr::pull(depmap_id)
  
  
  # Convert GDSC names to DepMap IDs
  gdsc_depmap <- gdsc_match %>%
    dplyr::filter(Sample.Name %in% gdsc_lines) %>%
    dplyr::pull(ModelID)
  
  
  common <- intersect(
    gdsc_depmap,
    prism_lines
  )
  
  
  data.frame(
    Drug = drug,
    
    GDSC_N = length(gdsc_depmap),
    
    PRISM_N = length(unique(prism_lines)),
    
    Shared_N = length(common),
    
    Shared_fraction =
      length(common) /
      min(
        length(gdsc_depmap),
        length(unique(prism_lines))
      ),
    
    Shared_lines =
      paste(common, collapse = ", ")
  )
})

overlap_results <- bind_rows(overlap_results)

overlap_results %>%
  arrange(desc(Shared_fraction))
# ============================================================
# REBUILD CORRECT GDSC HGSOC MATCH TABLE
# ============================================================

library(dplyr)

clean_cell_name <- function(x) {
  x <- toupper(x)
  x <- gsub("NIH:", "", x)
  x <- gsub("NCI-", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

# DepMap HGSOC projection
depmap_hgsoc <- state_projection %>%
  dplyr::mutate(
    clean_name = clean_cell_name(CellLineName)
  )

# GDSC metadata
gdsc_meta <- PharmacoGx::sampleInfo(gdsc2) %>%
  as.data.frame() %>%
  dplyr::mutate(
    clean_name = clean_cell_name(Sample.Name)
  )

# Correct HGSOC ↔ GDSC mapping
gdsc_hgsoc_match <- depmap_hgsoc %>%
  dplyr::inner_join(
    gdsc_meta %>%
      dplyr::select(
        Sample.Name,
        COSMIC.identifier,
        clean_name
      ),
    by = "clean_name"
  ) %>%
  dplyr::select(
    ModelID,
    CellLineName,
    State2_Likeness,
    State_Like,
    Sample.Name,
    COSMIC.identifier
  ) %>%
  dplyr::distinct()

cat("Matched HGSOC models:", nrow(gdsc_hgsoc_match), "\n")

gdsc_hgsoc_match %>%
  dplyr::arrange(State2_Likeness)
library(dplyr)

# ============================================================
# NORMALIZE GDSC SUMMARY COLUMN NAMES
# ============================================================

clean_cell_name <- function(x) {
  x <- toupper(x)
  x <- gsub("NIH:", "", x)
  x <- gsub("NCI-", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

# GDSC summary column mapping
gdsc_col_map <- data.frame(
  GDSC_col = colnames(gdsc2_summary),
  clean_name = clean_cell_name(colnames(gdsc2_summary)),
  stringsAsFactors = FALSE
)

# Add actual gdsc2_summary column name to our 13 HGSOC mapping
gdsc_hgsoc_map <- gdsc_hgsoc_match %>%
  dplyr::mutate(
    clean_name = clean_cell_name(Sample.Name)
  ) %>%
  dplyr::left_join(
    gdsc_col_map,
    by = "clean_name"
  )

cat(
  "Mapped to GDSC sensitivity matrix:",
  sum(!is.na(gdsc_hgsoc_map$GDSC_col)),
  "/",
  nrow(gdsc_hgsoc_map),
  "\n"
)

gdsc_hgsoc_map %>%
  dplyr::select(
    ModelID,
    CellLineName,
    Sample.Name,
    GDSC_col,
    State2_Likeness
  )
# ============================================================
# CORRECT GDSC × PRISM OVERLAP
# ============================================================

candidate_drugs <- state2_concordant$Drug_GDSC

overlap_results <- lapply(candidate_drugs, function(drug) {
  
  # ---------------------------
  # GDSC
  # ---------------------------
  
  valid_map <- gdsc_hgsoc_map %>%
    dplyr::filter(!is.na(GDSC_col))
  
  gdsc_values <- as.numeric(
    gdsc2_summary[
      drug,
      valid_map$GDSC_col
    ]
  )
  
  names(gdsc_values) <- valid_map$ModelID
  
  gdsc_depmap <- names(gdsc_values)[
    is.finite(gdsc_values)
  ]
  
  
  # ---------------------------
  # PRISM
  # ---------------------------
  
  prism_drug <- state2_concordant %>%
    dplyr::filter(Drug_GDSC == drug) %>%
    dplyr::pull(Drug_PRISM)
  
  prism_depmap <- prism_collapsed %>%
    dplyr::filter(
      name == prism_drug,
      is.finite(auc)
    ) %>%
    dplyr::pull(depmap_id) %>%
    unique()
  
  
  # ---------------------------
  # OVERLAP
  # ---------------------------
  
  shared <- intersect(
    gdsc_depmap,
    prism_depmap
  )
  
  data.frame(
    Drug = drug,
    
    GDSC_N = length(gdsc_depmap),
    PRISM_N = length(prism_depmap),
    
    Shared_N = length(shared),
    
    GDSC_only_N =
      length(setdiff(gdsc_depmap, prism_depmap)),
    
    PRISM_only_N =
      length(setdiff(prism_depmap, gdsc_depmap)),
    
    Shared_fraction =
      ifelse(
        min(length(gdsc_depmap),
            length(prism_depmap)) > 0,
        
        length(shared) /
          min(
            length(gdsc_depmap),
            length(prism_depmap)
          ),
        
        NA
      ),
    
    Shared_lines =
      paste(shared, collapse = ", ")
  )
})

overlap_results <- dplyr::bind_rows(
  overlap_results
)

print(
  as.data.frame(overlap_results),
  row.names = FALSE
)
library(dplyr)

# ============================================================
# NON-OVERLAPPING CELL-LINE ANALYSIS
# ============================================================

nonoverlap_results <- lapply(
  state2_concordant$Drug_GDSC,
  function(drug) {
    
    # --------------------------------
    # Drug names
    # --------------------------------
    
    prism_drug <- state2_concordant %>%
      dplyr::filter(Drug_GDSC == drug) %>%
      dplyr::pull(Drug_PRISM)
    
    
    # --------------------------------
    # GDSC data
    # --------------------------------
    
    valid_map <- gdsc_hgsoc_map %>%
      dplyr::filter(!is.na(GDSC_col))
    
    gdsc_values <- as.numeric(
      gdsc2_summary[
        drug,
        valid_map$GDSC_col
      ]
    )
    
    gdsc_df <- data.frame(
      ModelID = valid_map$ModelID,
      State2_Likeness = valid_map$State2_Likeness,
      Response = gdsc_values
    ) %>%
      dplyr::filter(is.finite(Response))
    
    
    # --------------------------------
    # PRISM data
    # --------------------------------
    
    prism_df <- prism_collapsed %>%
      dplyr::filter(
        name == prism_drug,
        is.finite(auc)
      ) %>%
      dplyr::transmute(
        ModelID = depmap_id,
        State2_Likeness,
        Response = auc
      )
    
    
    # --------------------------------
    # Remove shared cell lines
    # --------------------------------
    
    shared <- intersect(
      gdsc_df$ModelID,
      prism_df$ModelID
    )
    
    gdsc_unique <- gdsc_df %>%
      dplyr::filter(
        !ModelID %in% shared
      )
    
    prism_unique <- prism_df %>%
      dplyr::filter(
        !ModelID %in% shared
      )
    
    
    # --------------------------------
    # Spearman correlations
    # --------------------------------
    
    rho_gdsc <- if(nrow(gdsc_unique) >= 4) {
      cor(
        gdsc_unique$State2_Likeness,
        gdsc_unique$Response,
        method = "spearman"
      )
    } else NA
    
    rho_prism <- if(nrow(prism_unique) >= 4) {
      cor(
        prism_unique$State2_Likeness,
        prism_unique$Response,
        method = "spearman"
      )
    } else NA
    
    
    data.frame(
      Drug = drug,
      
      Shared_removed = length(shared),
      
      GDSC_unique_N = nrow(gdsc_unique),
      Rho_GDSC_unique = rho_gdsc,
      
      PRISM_unique_N = nrow(prism_unique),
      Rho_PRISM_unique = rho_prism,
      
      Direction_preserved =
        !is.na(rho_gdsc) &
        !is.na(rho_prism) &
        rho_gdsc < 0 &
        rho_prism < 0
    )
  }
)

nonoverlap_results <- dplyr::bind_rows(
  nonoverlap_results
)

nonoverlap_results <- nonoverlap_results %>%
  dplyr::arrange(
    dplyr::desc(Direction_preserved),
    Rho_GDSC_unique,
    Rho_PRISM_unique
  )

print(
  as.data.frame(nonoverlap_results),
  row.names = FALSE
)
# PRISM annotation
prism %>%
  dplyr::filter(
    grepl("GSK.?1904529A", name, ignore.case = TRUE)
  ) %>%
  dplyr::select(
    name,
    moa,
    target,
    indication,
    phase
  ) %>%
  dplyr::distinct()

# GDSC annotation
gdsc_drugs %>%
  dplyr::filter(
    grepl("GSK.?1904529A", DRUG_NAME, ignore.case = TRUE)
  ) %>%
  dplyr::select(
    DRUG_NAME,
    TARGET,
    TARGET_PATHWAY,
    FDA
  )
# ============================================================
# IGF1R / INSR EXPRESSION IN TCGA ROBUST STATES
# ============================================================

targets <- c("IGF1R", "INSR")

target_expr <- rna_vst_pc[
  intersect(targets, rownames(rna_vst_pc)),
  rownames(factor_scores),
  drop = FALSE
]

target_expr_df <- data.frame(
  Patient = rep(colnames(target_expr), each = nrow(target_expr)),
  Gene = rep(rownames(target_expr), times = ncol(target_expr)),
  Expression = as.vector(target_expr)
)

target_expr_df <- target_expr_df %>%
  dplyr::left_join(
    robust_state_df %>%
      dplyr::select(
        submitter_id,
        Robust_State
      ),
    by = c("Patient" = "submitter_id")
  )

# State comparison
target_expr_results <- target_expr_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Median_State1 = median(
      Expression[Robust_State == 1],
      na.rm = TRUE
    ),
    Median_State2 = median(
      Expression[Robust_State == 2],
      na.rm = TRUE
    ),
    Difference = Median_State2 - Median_State1,
    
    P_value = wilcox.test(
      Expression[Robust_State == 2],
      Expression[Robust_State == 1]
    )$p.value,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  )

target_expr_results
# ============================================================
# IGF1R / INSR CNV IN STATE1 vs STATE2
# ============================================================

targets <- c("IGF1R", "INSR")

target_cnv <- cnv_centered[
  intersect(targets, rownames(cnv_centered)),
  rownames(factor_scores),
  drop = FALSE
]

target_cnv_df <- data.frame(
  Patient = rep(colnames(target_cnv), each = nrow(target_cnv)),
  Gene = rep(rownames(target_cnv), times = ncol(target_cnv)),
  CNV = as.vector(target_cnv)
) %>%
  dplyr::left_join(
    robust_state_df %>%
      dplyr::select(submitter_id, Robust_State),
    by = c("Patient" = "submitter_id")
  )

target_cnv_results <- target_cnv_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Median_State1 = median(
      CNV[Robust_State == 1],
      na.rm = TRUE
    ),
    Median_State2 = median(
      CNV[Robust_State == 2],
      na.rm = TRUE
    ),
    Difference = Median_State2 - Median_State1,
    P_value = wilcox.test(
      CNV[Robust_State == 2],
      CNV[Robust_State == 1]
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  )

target_cnv_results
# ============================================================
# ARE SPLICING PROGRAMS ENRICHED BETWEEN STATE1 AND STATE2?
# ============================================================

library(dplyr)

splicing_gsea <- as.data.frame(gsea_go) %>%
  dplyr::filter(
    grepl(
      "splic|spliceosome|mRNA processing|RNA processing",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    ID,
    Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    core_enrichment
  ) %>%
  dplyr::arrange(p.adjust)

print(
  splicing_gsea,
  row.names = FALSE
)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# ============================================================
# 1. STATE1 ↔ STATE2 RANKED GENE LIST
# ============================================================

rank_df <- gsea_stats %>%
  dplyr::filter(
    !is.na(Z),
    !is.na(Gene)
  )

# SYMBOL -> ENTREZ
gene_map <- bitr(
  rank_df$Gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

rank_df2 <- rank_df %>%
  inner_join(
    gene_map,
    by = c("Gene" = "SYMBOL")
  ) %>%
  arrange(desc(abs(Z))) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- rank_df2$Z
names(gene_list) <- rank_df2$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

# ============================================================
# 2. GSEA WITHOUT SIGNIFICANCE FILTER
# ============================================================

gsea_splicing_test <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

all_go <- as.data.frame(gsea_splicing_test)

# ============================================================
# 3. EXTRACT SPLICING / RNA PROCESSING TERMS
# ============================================================

splicing_results <- all_go %>%
  filter(
    grepl(
      "splic|spliceosom|RNA processing|RNA splicing|mRNA processing",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    ID,
    Description,
    setSize,
    NES,
    pvalue,
    p.adjust,
    core_enrichment
  ) %>%
  arrange(p.adjust)

print(
  splicing_results,
  row.names = FALSE
)
splicing_results <- all_go %>%
  dplyr::filter(
    grepl(
      "splic|spliceosom|RNA processing|RNA splicing|mRNA processing",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    ID,
    Description,
    setSize,
    NES,
    pvalue,
    p.adjust,
    core_enrichment
  ) %>%
  dplyr::arrange(p.adjust)

print(
  as.data.frame(splicing_results),
  row.names = FALSE
)
# ============================================================
# GSK-1904529A:
# TCGA STATE GENES vs PRISM DRUG RESPONSE
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. GSK PRISM response
# ------------------------------------------------------------

gsk_prism <- prism_collapsed %>%
  dplyr::filter(
    name == "GSK1904529A",
    is.finite(auc)
  ) %>%
  dplyr::select(
    depmap_id,
    State2_Likeness,
    auc
  ) %>%
  dplyr::distinct()

cat("GSK PRISM HGSOC lines:", nrow(gsk_prism), "\n")

print(gsk_prism)


# ------------------------------------------------------------
# 2. DepMap expression column names
# ------------------------------------------------------------

# Gene names in DepMap expression:
depmap_gene_names <- gsub(
  " \\([0-9]+\\)$",
  "",
  colnames(expr)
)

# TCGA 110 differential genes
state_genes <- gsea_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(
    Gene,
    Difference,
    Z,
    FDR
  ) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

cat("TCGA state genes:", nrow(state_genes), "\n")


# ------------------------------------------------------------
# 3. Match TCGA genes to DepMap expression
# ------------------------------------------------------------

gene_idx <- match(
  state_genes$Gene,
  depmap_gene_names
)

state_genes$DepMap_col <- gene_idx

state_genes_matched <- state_genes %>%
  dplyr::filter(!is.na(DepMap_col))

cat(
  "Genes matched to DepMap:",
  nrow(state_genes_matched),
  "/",
  nrow(state_genes),
  "\n"
)


# ------------------------------------------------------------
# 4. Match PRISM cell lines to DepMap expression rows
# ------------------------------------------------------------

gsk_expr_idx <- match(
  gsk_prism$depmap_id,
  expr$ModelID
)

keep_lines <- !is.na(gsk_expr_idx)

gsk_prism2 <- gsk_prism[keep_lines, ]

gsk_expr_idx <- gsk_expr_idx[keep_lines]

cat(
  "GSK lines with expression:",
  length(gsk_expr_idx),
  "\n"
)


# ------------------------------------------------------------
# 5. Gene-by-gene correlation with GSK AUC
# ------------------------------------------------------------

gsk_gene_results <- lapply(
  seq_len(nrow(state_genes_matched)),
  function(i) {
    
    gene <- state_genes_matched$Gene[i]
    
    col_idx <- state_genes_matched$DepMap_col[i]
    
    x <- as.numeric(
      expr[gsk_expr_idx, col_idx]
    )
    
    y <- gsk_prism2$auc
    
    keep <- is.finite(x) & is.finite(y)
    
    if(sum(keep) < 5) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        x[keep],
        y[keep],
        method = "spearman",
        exact = TRUE
      )
    )
    
    data.frame(
      Gene = gene,
      
      TCGA_Difference =
        state_genes_matched$Difference[i],
      
      TCGA_Z =
        state_genes_matched$Z[i],
      
      TCGA_FDR =
        state_genes_matched$FDR[i],
      
      N = sum(keep),
      
      Drug_Rho =
        unname(test$estimate),
      
      Drug_P =
        test$p.value
    )
  }
)

gsk_gene_results <- dplyr::bind_rows(
  gsk_gene_results
)


# ------------------------------------------------------------
# 6. Multiple-testing correction
# ------------------------------------------------------------

gsk_gene_results <- gsk_gene_results %>%
  dplyr::mutate(
    
    Drug_FDR =
      p.adjust(
        Drug_P,
        method = "BH"
      ),
    
    TCGA_direction =
      ifelse(
        TCGA_Difference > 0,
        "State2-high",
        "State1-high"
      ),
    
    # AUC lower = more sensitive
    #
    # State2-high gene + rho < 0:
    # higher expression -> lower AUC
    # = compatible with State2 sensitivity
    #
    # State1-high gene + rho > 0:
    # higher expression -> higher AUC
    # = also compatible with State2 sensitivity
    
    Supports_State2_Sensitivity =
      (
        TCGA_Difference > 0 &
          Drug_Rho < 0
      ) |
      (
        TCGA_Difference < 0 &
          Drug_Rho > 0
      )
  )


# ------------------------------------------------------------
# 7. Show strongest candidates
# ------------------------------------------------------------

gsk_gene_results %>%
  dplyr::arrange(
    dplyr::desc(abs(Drug_Rho))
  ) %>%
  dplyr::select(
    Gene,
    TCGA_direction,
    TCGA_Difference,
    Drug_Rho,
    Drug_P,
    Drug_FDR,
    Supports_State2_Sensitivity
  ) %>%
  head(30) %>%
  print(row.names = FALSE)
# ============================================================
# GSK-1904529A:
# TCGA STATE GENES vs GDSC DRUG RESPONSE
# ============================================================

library(dplyr)

drug <- "GSK-1904529A"

# ------------------------------------------------------------
# 1. GDSC GSK response
# ------------------------------------------------------------

valid_map <- gdsc_hgsoc_map %>%
  dplyr::filter(!is.na(GDSC_col))

gsk_values <- as.numeric(
  gdsc2_summary[
    drug,
    valid_map$GDSC_col
  ]
)

gsk_gdsc <- data.frame(
  ModelID = valid_map$ModelID,
  State2_Likeness = valid_map$State2_Likeness,
  IC50 = gsk_values
) %>%
  dplyr::filter(is.finite(IC50))

cat("GSK GDSC lines:", nrow(gsk_gdsc), "\n")


# ------------------------------------------------------------
# 2. Match GDSC lines to DepMap expression
# ------------------------------------------------------------

gdsc_expr_idx <- match(
  gsk_gdsc$ModelID,
  expr$ModelID
)

keep_lines <- !is.na(gdsc_expr_idx)

gsk_gdsc2 <- gsk_gdsc[keep_lines, ]
gdsc_expr_idx <- gdsc_expr_idx[keep_lines]

cat(
  "GSK GDSC lines with expression:",
  length(gdsc_expr_idx),
  "\n"
)


# ------------------------------------------------------------
# 3. Correlate 110 state genes with GSK IC50
# ------------------------------------------------------------

gsk_gdsc_gene_results <- lapply(
  seq_len(nrow(state_genes_matched)),
  function(i) {
    
    gene <- state_genes_matched$Gene[i]
    col_idx <- state_genes_matched$DepMap_col[i]
    
    x <- as.numeric(
      expr[gdsc_expr_idx, col_idx]
    )
    
    y <- gsk_gdsc2$IC50
    
    keep <- is.finite(x) & is.finite(y)
    
    if(sum(keep) < 5)
      return(NULL)
    
    test <- suppressWarnings(
      cor.test(
        x[keep],
        y[keep],
        method = "spearman",
        exact = TRUE
      )
    )
    
    data.frame(
      Gene = gene,
      N_GDSC = sum(keep),
      GDSC_Rho = unname(test$estimate),
      GDSC_P = test$p.value
    )
  }
)

gsk_gdsc_gene_results <- bind_rows(
  gsk_gdsc_gene_results
)


# ------------------------------------------------------------
# 4. Join PRISM + GDSC results
# ------------------------------------------------------------

gsk_crossplatform_genes <- gsk_gene_results %>%
  dplyr::select(
    Gene,
    TCGA_Difference,
    TCGA_Z,
    TCGA_FDR,
    Drug_Rho,
    Drug_P
  ) %>%
  dplyr::rename(
    PRISM_Rho = Drug_Rho,
    PRISM_P = Drug_P
  ) %>%
  dplyr::inner_join(
    gsk_gdsc_gene_results,
    by = "Gene"
  ) %>%
  dplyr::mutate(
    
    Same_Direction =
      sign(PRISM_Rho) == sign(GDSC_Rho),
    
    Supports_State2_Both =
      Same_Direction &
      (
        (TCGA_Difference > 0 &
           PRISM_Rho < 0 &
           GDSC_Rho < 0) |
          
          (TCGA_Difference < 0 &
             PRISM_Rho > 0 &
             GDSC_Rho > 0)
      ),
    
    Mean_abs_rho =
      (
        abs(PRISM_Rho) +
          abs(GDSC_Rho)
      ) / 2
  ) %>%
  dplyr::arrange(
    dplyr::desc(Supports_State2_Both),
    dplyr::desc(Mean_abs_rho)
  )


# ------------------------------------------------------------
# 5. Results
# ------------------------------------------------------------

gsk_crossplatform_genes %>%
  dplyr::select(
    Gene,
    TCGA_Difference,
    PRISM_Rho,
    GDSC_Rho,
    Same_Direction,
    Supports_State2_Both,
    Mean_abs_rho
  ) %>%
  head(30) %>%
  print(row.names = FALSE)


gsk_state2_program %>%
  dplyr::select(
    Gene,
    TCGA_Difference,
    PRISM_Rho,
    GDSC_Rho,
    Mean_abs_rho
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
# ============================================================
# ENRICHMENT OF 41 CROSS-PLATFORM GSK/STATE GENES
# Background = 110 TCGA State1/State2 genes
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# 41 cross-platform genes
genes_41 <- unique(gsk_state2_program$Gene)

# Background = 110 TCGA state-differential genes
background_110 <- unique(state_genes$Gene)

# ------------------------------------------------------------
# SYMBOL -> ENTREZ
# ------------------------------------------------------------

map_41 <- clusterProfiler::bitr(
  genes_41,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

map_bg <- clusterProfiler::bitr(
  background_110,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("41 genes mapped:", nrow(map_41), "\n")
cat("Background genes mapped:", nrow(map_bg), "\n")

# ------------------------------------------------------------
# GO Biological Process enrichment
# ------------------------------------------------------------

gsk_go <- clusterProfiler::enrichGO(
  gene          = unique(map_41$ENTREZID),
  universe      = unique(map_bg$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

gsk_go_df <- as.data.frame(gsk_go)

# ------------------------------------------------------------
# Show results
# ------------------------------------------------------------

gsk_go_df %>%
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    geneID,
    Count
  ) %>%
  dplyr::arrange(p.adjust) %>%
  head(30) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
# ============================================================
# CONTINUOUS CROSS-PLATFORM GSK / STATE SCORE
# ============================================================

gsk_rank <- gsk_crossplatform_genes %>%
  dplyr::mutate(
    
    # TCGA direction:
    # +1 = State2-high
    # -1 = State1-high
    TCGA_sign = sign(TCGA_Difference),
    
    # Drug response direction
    # PRISM AUC and GDSC IC50:
    # negative rho = greater sensitivity with higher expression
    Mean_drug_rho =
      (PRISM_Rho + GDSC_Rho) / 2,
    
    # Positive score =
    # association compatible with State2 GSK sensitivity
    GSK_State2_Score =
      -TCGA_sign * Mean_drug_rho
    
  ) %>%
  dplyr::arrange(
    dplyr::desc(GSK_State2_Score)
  )

gsk_rank %>%
  dplyr::select(
    Gene,
    TCGA_Difference,
    PRISM_Rho,
    GDSC_Rho,
    GSK_State2_Score
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
##
# ============================================================
# GSK-1904529A
# INDEPENDENT CELL-LINE ANALYSIS
# Remove all shared GDSC–PRISM cell lines
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Identify shared cell lines
# ------------------------------------------------------------

shared_ids <- intersect(
  gsk_gdsc2$ModelID,
  gsk_prism2$depmap_id
)

cat("Shared cell lines:\n")
print(shared_ids)

# ------------------------------------------------------------
# 2. Remove shared lines
# ------------------------------------------------------------

gsk_gdsc_unique <- gsk_gdsc2 %>%
  dplyr::filter(!ModelID %in% shared_ids)

gsk_prism_unique <- gsk_prism2 %>%
  dplyr::filter(!depmap_id %in% shared_ids)

cat(
  "\nGDSC-only N:",
  nrow(gsk_gdsc_unique),
  "\n"
)

cat(
  "PRISM-only N:",
  nrow(gsk_prism_unique),
  "\n"
)

# Verify zero overlap
cat(
  "Remaining overlap:",
  length(intersect(
    gsk_gdsc_unique$ModelID,
    gsk_prism_unique$depmap_id
  )),
  "\n"
)


# ============================================================
# 3. Match unique lines to DepMap expression
# ============================================================

gdsc_unique_idx <- match(
  gsk_gdsc_unique$ModelID,
  expr$ModelID
)

prism_unique_idx <- match(
  gsk_prism_unique$depmap_id,
  expr$ModelID
)

# Safety checks
stopifnot(!anyNA(gdsc_unique_idx))
stopifnot(!anyNA(prism_unique_idx))


# ============================================================
# 4. Recalculate correlations for all 110 TCGA state genes
# ============================================================

independent_gene_results <- lapply(
  seq_len(nrow(state_genes_matched)),
  function(i) {
    
    gene <- state_genes_matched$Gene[i]
    col_idx <- state_genes_matched$DepMap_col[i]
    
    # ------------------------
    # GDSC-only
    # ------------------------
    
    x_gdsc <- as.numeric(
      expr[gdsc_unique_idx, col_idx]
    )
    
    y_gdsc <- gsk_gdsc_unique$IC50
    
    keep_gdsc <-
      is.finite(x_gdsc) &
      is.finite(y_gdsc)
    
    if(sum(keep_gdsc) >= 4) {
      
      gdsc_test <- suppressWarnings(
        cor.test(
          x_gdsc[keep_gdsc],
          y_gdsc[keep_gdsc],
          method = "spearman",
          exact = TRUE
        )
      )
      
      gdsc_rho <- unname(gdsc_test$estimate)
      gdsc_p   <- gdsc_test$p.value
      
    } else {
      
      gdsc_rho <- NA
      gdsc_p   <- NA
    }
    
    
    # ------------------------
    # PRISM-only
    # ------------------------
    
    x_prism <- as.numeric(
      expr[prism_unique_idx, col_idx]
    )
    
    y_prism <- gsk_prism_unique$auc
    
    keep_prism <-
      is.finite(x_prism) &
      is.finite(y_prism)
    
    if(sum(keep_prism) >= 4) {
      
      prism_test <- suppressWarnings(
        cor.test(
          x_prism[keep_prism],
          y_prism[keep_prism],
          method = "spearman",
          exact = TRUE
        )
      )
      
      prism_rho <- unname(prism_test$estimate)
      prism_p   <- prism_test$p.value
      
    } else {
      
      prism_rho <- NA
      prism_p   <- NA
    }
    
    
    # ------------------------
    # Output
    # ------------------------
    
    data.frame(
      
      Gene = gene,
      
      TCGA_Difference =
        state_genes_matched$Difference[i],
      
      TCGA_Z =
        state_genes_matched$Z[i],
      
      TCGA_FDR =
        state_genes_matched$FDR[i],
      
      GDSC_N =
        sum(keep_gdsc),
      
      GDSC_Rho =
        gdsc_rho,
      
      GDSC_P =
        gdsc_p,
      
      PRISM_N =
        sum(keep_prism),
      
      PRISM_Rho =
        prism_rho,
      
      PRISM_P =
        prism_p
    )
  }
)

independent_gene_results <- dplyr::bind_rows(
  independent_gene_results
)


# ============================================================
# 5. Determine directional concordance
# ============================================================

independent_gene_results <- independent_gene_results %>%
  dplyr::mutate(
    
    Same_Direction =
      !is.na(GDSC_Rho) &
      !is.na(PRISM_Rho) &
      sign(GDSC_Rho) ==
      sign(PRISM_Rho),
    
    Supports_State2_Both =
      Same_Direction &
      (
        (
          TCGA_Difference > 0 &
            GDSC_Rho < 0 &
            PRISM_Rho < 0
        ) |
          (
            TCGA_Difference < 0 &
              GDSC_Rho > 0 &
              PRISM_Rho > 0
          )
      ),
    
    Mean_abs_rho =
      (
        abs(GDSC_Rho) +
          abs(PRISM_Rho)
      ) / 2
  )


# ============================================================
# 6. Multiple-testing correction separately
# ============================================================

independent_gene_results <- independent_gene_results %>%
  dplyr::mutate(
    
    GDSC_FDR =
      p.adjust(
        GDSC_P,
        method = "BH"
      ),
    
    PRISM_FDR =
      p.adjust(
        PRISM_P,
        method = "BH"
      )
  )


# ============================================================
# 7. Show State2-compatible genes
# ============================================================

independent_candidates <- independent_gene_results %>%
  dplyr::filter(
    Supports_State2_Both == TRUE
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_abs_rho)
  )

cat(
  "\nIndependent cross-platform candidates:",
  nrow(independent_candidates),
  "\n\n"
)

independent_candidates %>%
  dplyr::select(
    Gene,
    TCGA_Difference,
    GDSC_N,
    GDSC_Rho,
    GDSC_P,
    PRISM_N,
    PRISM_Rho,
    PRISM_P,
    Mean_abs_rho
  ) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# ============================================================
# CROSS-PLATFORM GSK CONCORDANCE PLOT
# Non-overlapping GDSC vs PRISM cell lines
# ============================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

plot_df <- independent_gene_results %>%
  dplyr::filter(
    !is.na(GDSC_Rho),
    !is.na(PRISM_Rho)
  ) %>%
  dplyr::mutate(
    TCGA_State = ifelse(
      TCGA_Difference > 0,
      "State2-high",
      "State1-high"
    ),
    
    Concordance = dplyr::case_when(
      Supports_State2_Both ~
        "Supports State2 sensitivity",
      
      Same_Direction ~
        "Same drug-response direction",
      
      TRUE ~
        "Discordant"
    )
  )

# Genes to label
label_genes <- c(
  "IGFN1",
  "GATA2",
  "ZIC5",
  "ZIC2",
  "POU6F2",
  "SAMD11",
  "LIN28B"
)

p <- ggplot(
  plot_df,
  aes(
    x = GDSC_Rho,
    y = PRISM_Rho
  )
) +
  
  # zero lines
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  # all genes
  geom_point(
    aes(
      color = Concordance,
      shape = TCGA_State
    ),
    size = 3,
    alpha = 0.8
  ) +
  
  # labels for strongest candidates
  ggrepel::geom_text_repel(
    data = plot_df %>%
      dplyr::filter(Gene %in% label_genes),
    
    aes(label = Gene),
    
    size = 4,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  
  coord_fixed() +
  
  labs(
    title =
      "Cross-platform molecular correlates of GSK-1904529A sensitivity",
    
    subtitle =
      "Independent HGSOC cell-line sets: GDSC (n=5) and PRISM (n=6)",
    
    x =
      "GDSC: Spearman correlation with GSK IC50",
    
    y =
      "PRISM: Spearman correlation with GSK AUC",
    
    color =
      "Cross-platform relationship",
    
    shape =
      "TCGA molecular state"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title =
      element_text(face = "bold"),
    
    legend.position =
      "right"
  )

p
# ============================================================
# TOP INDEPENDENT CROSS-PLATFORM CANDIDATES
# ============================================================

top_candidates <- independent_candidates %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(
    Gene,
    GDSC_Rho,
    PRISM_Rho
  ) %>%
  tidyr::pivot_longer(
    cols = c(GDSC_Rho, PRISM_Rho),
    names_to = "Dataset",
    values_to = "Rho"
  ) %>%
  dplyr::mutate(
    Dataset = dplyr::recode(
      Dataset,
      GDSC_Rho = "GDSC",
      PRISM_Rho = "PRISM"
    ),
    Gene = factor(
      Gene,
      levels = rev(
        unique(independent_candidates$Gene[1:15])
      )
    )
  )

ggplot(
  top_candidates,
  aes(
    x = Rho,
    y = Gene,
    shape = Dataset
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3.5,
    position = position_dodge(width = 0.5)
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  labs(
    title =
      "State-associated genes reproducibly linked to GSK sensitivity",
    
    subtitle =
      "Non-overlapping HGSOC cell-line sets",
    
    x =
      "Spearman correlation with drug response",
    
    y = NULL,
    
    shape =
      "Drug-response dataset"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title =
      element_text(face = "bold")
  )
# ============================================================
# TOP INDEPENDENT CROSS-PLATFORM CANDIDATES
# ============================================================

top_candidates <- independent_candidates %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::select(
    Gene,
    GDSC_Rho,
    PRISM_Rho
  ) %>%
  tidyr::pivot_longer(
    cols = c(GDSC_Rho, PRISM_Rho),
    names_to = "Dataset",
    values_to = "Rho"
  ) %>%
  dplyr::mutate(
    Dataset = dplyr::recode(
      Dataset,
      GDSC_Rho = "GDSC",
      PRISM_Rho = "PRISM"
    ),
    Gene = factor(
      Gene,
      levels = rev(
        unique(independent_candidates$Gene[1:15])
      )
    )
  )

ggplot(
  top_candidates,
  aes(
    x = Rho,
    y = Gene,
    shape = Dataset
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3.5,
    position = position_dodge(width = 0.5)
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  labs(
    title =
      "State-associated genes reproducibly linked to GSK sensitivity",
    
    subtitle =
      "Non-overlapping HGSOC cell-line sets",
    
    x =
      "Spearman correlation with drug response",
    
    y = NULL,
    
    shape =
      "Drug-response dataset"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title =
      element_text(face = "bold")
  )
# ============================================================
# TCGA-OV ALTERNATIVE SPLICING DATA
# Step 1 — Download PSI data
# ============================================================

library(httr)
library(readr)
library(dplyr)

# TCGA SpliceSeq download page
base_url <- "https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIdownload.jsp"

# First confirm that the server is reachable
res <- GET(base_url)

cat("HTTP status:", status_code(res), "\n")
cat("Content type:", headers(res)[["content-type"]], "\n")

if (status_code(res) == 200) {
  cat("TCGA SpliceSeq server reachable.\n")
} else {
  cat("Server returned:", status_code(res), "\n")
}
robust_state_df
# ============================================================
# STEP 2 — Inspect TCGA SpliceSeq download form
# ============================================================

library(rvest)
library(httr)

page <- read_html(base_url)

forms <- html_elements(page, "form")

cat("Number of forms:", length(forms), "\n\n")

for (i in seq_along(forms)) {
  
  cat("========== FORM", i, "==========\n")
  
  cat(
    "Action:",
    html_attr(forms[[i]], "action"),
    "\n"
  )
  
  cat(
    "Method:",
    html_attr(forms[[i]], "method"),
    "\n"
  )
  
  inputs <- html_elements(
    forms[[i]],
    "input, select, textarea"
  )
  
  info <- data.frame(
    tag   = html_name(inputs),
    name  = html_attr(inputs, "name"),
    type  = html_attr(inputs, "type"),
    value = html_attr(inputs, "value")
  )
  
  print(info)
  cat("\n")
}
# ============================================================
# STEP 3 — Download TCGA-OV ES PSI for our 277 patients
# ============================================================

library(httr)
library(dplyr)

# ------------------------------------------------------------
# 1. Our patient IDs
# ------------------------------------------------------------

patient_ids <- unique(robust_state_df$submitter_id)

cat("Our patients:", length(patient_ids), "\n")

# SpliceSeq textarea: one sample per line
sample_string <- paste(patient_ids, collapse = "\n")


# ------------------------------------------------------------
# 2. POST request
# ------------------------------------------------------------

download_url <- "https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload"

res_es <- POST(
  download_url,
  
  body = list(
    tissue       = "OV",
    genes        = "",
    samples      = sample_string,
    splicetype   = "ES",
    
    # Keep events with PSI available in >=75% samples
    pctwithval   = 75,
    
    # Default expression filter
    avgexppct    = 50,
    
    # Don't impose PSI range/SD filtering yet
    psirange     = 0,
    psistd       = 0,
    
    annotation   = "true"
  ),
  
  encode = "form",
  timeout(300)
)


# ------------------------------------------------------------
# 3. Check response
# ------------------------------------------------------------

cat("HTTP status:", status_code(res_es), "\n")
cat(
  "Content type:",
  headers(res_es)[["content-type"]],
  "\n"
)

cat(
  "Content length:",
  length(content(res_es, as = "raw")),
  "bytes\n"
)
# ============================================================
# STEP 4 — Save and inspect returned ZIP
# ============================================================

zip_file <- tempfile(fileext = ".zip")

writeBin(
  content(res_es, as = "raw"),
  zip_file
)

cat("ZIP saved:", zip_file, "\n")
cat("ZIP size:", file.info(zip_file)$size, "bytes\n\n")

# List files inside ZIP
zip_info <- unzip(zip_file, list = TRUE)

print(zip_info)
# ============================================================
# STEP 5 — Extract ZIP
# ============================================================

extract_dir <- tempfile()
dir.create(extract_dir)

unzip(
  zip_file,
  exdir = extract_dir
)

files_inside <- list.files(
  extract_dir,
  recursive = TRUE,
  full.names = TRUE
)

print(files_inside)

# Print contents of returned files
for (f in files_inside) {
  
  cat("\n==============================\n")
  cat("FILE:", basename(f), "\n")
  cat("==============================\n")
  
  x <- readLines(
    f,
    warn = FALSE
  )
  
  cat(
    paste(head(x, 30), collapse = "\n"),
    "\n"
  )
}
\
