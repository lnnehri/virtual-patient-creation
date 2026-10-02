#23 eylul - virtual patient creation
#first get TCGA data
# Sadece ilk kez
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

clinical_tcga <- GDCquery_clinic(
  project = "TCGA-OV",
  type = "clinical"
)
dim(clinical_tcga)

colnames(clinical_tcga)

View(clinical_tcga)
library(TCGAbiolinks)

query_mut <- GDCquery(
  project = "TCGA-OV",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query_mut)

mutation_tcga <- GDCprepare(query_mut)
dim(mutation_tcga)

colnames(mutation_tcga)

head(mutation_tcga[, c(
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Variant_Classification",
  "Tumor_Sample_Barcode"
)])
mutation_tcga$patient_id <- substr(
  mutation_tcga$Tumor_Sample_Barcode,
  1, 12
)

length(unique(mutation_tcga$patient_id))
length(unique(clinical_tcga$submitter_id))

length(intersect(
  unique(clinical_tcga$submitter_id),
  unique(mutation_tcga$patient_id)
))
patient_cohort <- clinical_tcga[
  clinical_tcga$submitter_id %in% unique(mutation_tcga$patient_id),
]

dim(patient_cohort)
#kacinda RNA-seq var?
query_rna <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

rna_cases <- getResults(query_rna)

dim(rna_cases)
colnames(rna_cases)

length(unique(substr(rna_cases$cases, 1, 12)))

#somatitc + rna-seq ve clinical datasi olanlar
# RNA hasta ID'leri
rna_ids <- unique(substr(rna_cases$cases, 1, 12))

# Ortak hastalar
common_ids <- Reduce(intersect, list(
  unique(clinical_tcga$submitter_id),
  unique(mutation_tcga$patient_id),
  rna_ids
))

length(common_ids)
master_cohort <- clinical_tcga[
  clinical_tcga$submitter_id %in% common_ids,
]

dim(master_cohort)
#CNV verisi
query_cnv <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  access = "open"
)

cnv_cases <- getResults(query_cnv)

dim(cnv_cases)
colnames(cnv_cases)

# Kaç farklı hasta?
cnv_ids <- unique(substr(cnv_cases$cases, 1, 12))
length(cnv_ids)
#kac hastada CNV bilgisi var?
common_4layer <- intersect(common_ids, cnv_ids)

length(common_4layer)
final_ids <- common_4layer

final_cohort <- clinical_tcga[
  clinical_tcga$submitter_id %in% final_ids,
]

dim(final_cohort)
# Sadece final cohort mutasyonları
mut_final <- mutation_tcga[
  mutation_tcga$patient_id %in% final_ids,
]

# BRCA1 / BRCA2 mutation status
brca_status <- data.frame(
  submitter_id = final_ids,
  BRCA1_mut = final_ids %in% unique(
    mut_final$patient_id[mut_final$Hugo_Symbol == "BRCA1"]
  ),
  BRCA2_mut = final_ids %in% unique(
    mut_final$patient_id[mut_final$Hugo_Symbol == "BRCA2"]
  )
)

brca_status$BRCA_any_mut <- 
  brca_status$BRCA1_mut | brca_status$BRCA2_mut

table(brca_status$BRCA1_mut)
table(brca_status$BRCA2_mut)
table(brca_status$BRCA_any_mut)
final_cohort <- merge(
  final_cohort,
  brca_status,
  by = "submitter_id",
  all.x = TRUE
)

dim(final_cohort)
final_cohort <- merge(
  final_cohort,
  brca_status,
  by = "submitter_id",
  all.x = TRUE
)

dim(final_cohort)
#RNA-sew indiriyoruz
query_rna_final <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = final_ids
)

getResults(query_rna_final) |> dim()
GDCdownload(
  query_rna_final,
  method = "api",
  files.per.chunk = 20
)
rna_tcga <- GDCprepare(query_rna_final)

rna_counts <- assay(rna_tcga, "unstranded")

dim(rna_counts)

#error from library
library(SummarizedExperiment)

rna_counts <- assay(rna_tcga, "unstranded")

dim(rna_counts)

#
rna_patient_ids <- substr(colnames(rna_counts), 1, 12)

length(unique(rna_patient_ids))
table(duplicated(rna_patient_ids))

table(colData(rna_tcga)$sample_type)

rna_patient_ids[duplicated(rna_patient_ids)]

#nukseden tumoru olan hastalar
sample_info$patient_id[
  
  sample_info$sample_type == "Recurrent Tumor"
  
]

unique(sample_info$patient_id[
  
  sample_info$sample_type == "Recurrent Tumor"
  
])

recurrent_idx <- which(
  colData(rna_tcga)$sample_type == "Recurrent Tumor"
)

colnames(rna_counts)[recurrent_idx]

substr(colnames(rna_counts)[recurrent_idx], 1, 12)

#6 hastada tekrar eden tumor var. bunlari sonraki asama icin sakliyorum metastatik
recurrent_counts <- rna_counts[, recurrent_idx]

recurrent_patient_ids <- substr(
  colnames(recurrent_counts), 1, 12
)
#ana analiz icin primary tumoru olan hastalari tuutorum, recurrient olanlarin da primary bilgisini tutuyorum
primary_idx <- which(
  colData(rna_tcga)$sample_type == "Primary Tumor"
)

primary_counts <- rna_counts[, primary_idx]

primary_patient_ids <- substr(
  colnames(primary_counts), 1, 12
)

dim(primary_counts)
length(unique(primary_patient_ids))
#
colnames(primary_counts) <- primary_patient_ids
final_cohort_primary <- final_cohort[
  final_cohort$submitter_id %in% primary_patient_ids,
]

dim(final_cohort_primary)
query_cnv_final <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  access = "open",
  barcode = primary_patient_ids
)

dim(getResults(query_cnv_final))
#tekrarli CNV sonucu var
cnv_res <- getResults(query_cnv_final)

table(cnv_res$sample_type)
table(cnv_res$analysis_workflow_type)
cnv_res$patient_id <- substr(cnv_res$cases, 1, 12)

aggregate(
  patient_id ~ analysis_workflow_type,
  data = cnv_res,
  FUN = function(x) length(unique(x))
)
#Ascat2 ile en fazla CNV analiz edilmis onlari aliyorum ortak olsun diye
query_cnv_ascat2 <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  workflow.type = "ASCAT2",
  access = "open",
  barcode = primary_patient_ids
)

cnv_ascat2_res <- getResults(query_cnv_ascat2)

length(unique(substr(cnv_ascat2_res$cases, 1, 12)))
GDCdownload(
  query_cnv_ascat2,
  method = "api",
  files.per.chunk = 20
)
cnv_tcga <- GDCprepare(query_cnv_ascat2)

class(cnv_tcga)

dim(cnv_tcga)
#bunun icin de dublicate veriyor
cnv_ascat2_res <- getResults(query_cnv_ascat2)

table(cnv_ascat2_res$sample_type)


cnv_primary_res <- cnv_ascat2_res[
  grepl("Primary Tumor", cnv_ascat2_res$sample_type) &
    !grepl("Recurrent Tumor", cnv_ascat2_res$sample_type),
]

length(unique(substr(cnv_primary_res$cases, 1, 12)))
dim(cnv_primary_res)

cnv_primary_ids <- unique(
  substr(cnv_primary_res$cases, 1, 12)
)

query_cnv_primary <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  workflow.type = "ASCAT2",
  access = "open",
  barcode = cnv_primary_ids,
  sample.type = "Primary Tumor"
)

dim(getResults(query_cnv_primary))
#error

table(cnv_primary_res$sample_type)

query_cnv_primary <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  workflow.type = "ASCAT2",
  access = "open",
  barcode = cnv_primary_ids,
  sample.type = c(
    "Blood Derived Normal;Primary Tumor",
    "Primary Tumor;Blood Derived Normal",
    "Primary Tumor;Solid Tissue Normal",
    "Solid Tissue Normal;Primary Tumor"
  )
)

dim(getResults(query_cnv_primary))
head(cnv_primary_res[, c("id", "file_name", "cases", "sample_type")])
list.files(
  "GDCdata",
  pattern = "gene_level_copy_number",
  recursive = TRUE,
  full.names = TRUE
)[1:10]
all_cnv_files <- list.files(
  "GDCdata",
  pattern = "gene_level_copy_number",
  recursive = TRUE,
  full.names = TRUE
)

primary_cnv_files <- all_cnv_files[
  sapply(cnv_primary_res$id, function(x)
    any(grepl(x, all_cnv_files))
  )
]

length(primary_cnv_files)
sum(cnv_primary_res$id %in% basename(dirname(all_cnv_files)))
file_ids <- basename(dirname(all_cnv_files))

primary_cnv_files <- all_cnv_files[
  file_ids %in% cnv_primary_res$id
]

length(primary_cnv_files)
cnv_test <- read.delim(primary_cnv_files[1])

dim(cnv_test)
colnames(cnv_test)
head(cnv_test)
sum(!is.na(cnv_test$copy_number))
summary(cnv_test$copy_number)

cnv_file_ids <- basename(dirname(primary_cnv_files))

cnv_patient_ids <- cnv_primary_res$patient_id[
  match(cnv_file_ids, cnv_primary_res$id)
]

length(unique(cnv_patient_ids))
#ustteki eror verdi
cnv_patient_ids <- substr(
  cnv_primary_res$cases[
    match(cnv_file_ids, cnv_primary_res$id)
  ],
  1, 12
)

length(unique(cnv_patient_ids))
sum(is.na(cnv_patient_ids))

cnv_matrix <- sapply(
  primary_cnv_files,
  function(f) read.delim(f)$copy_number
)

rownames(cnv_matrix) <- cnv_test$gene_name
colnames(cnv_matrix) <- cnv_patient_ids

dim(cnv_matrix)
final_multiomics_ids <- intersect(
  primary_patient_ids,
  colnames(cnv_matrix)
)

length(final_multiomics_ids)
final_ids <- final_multiomics_ids

clinical_final <- final_cohort_primary[
  match(final_ids, final_cohort_primary$submitter_id),
]

rna_final <- primary_counts[, final_ids]

cnv_final <- cnv_matrix[, final_ids]

mutation_final <- mutation_tcga[
  mutation_tcga$patient_id %in% final_ids,
]
nrow(clinical_final)
ncol(rna_final)
ncol(cnv_final)
length(unique(mutation_final$patient_id))
clinical_missing <- data.frame(
  variable = colnames(clinical_final),
  missing_n = sapply(clinical_final, function(x)
    sum(is.na(x) | x == "")
  )
)

clinical_missing$missing_percent <-
  round(100 * clinical_missing$missing_n / nrow(clinical_final), 1)

clinical_missing <- clinical_missing[
  order(clinical_missing$missing_percent),
]

View(clinical_missing)
library(openxlsx)

write.xlsx(
  clinical_missing,
  file = "TCGA_OV_277_clinical_missingness.xlsx",
  rowNames = FALSE
)
clinical_missing[
  clinical_missing$missing_percent <= 30,
]
clinical_final$BRCA1_mut <- clinical_final$BRCA1_mut.x
clinical_final$BRCA2_mut <- clinical_final$BRCA2_mut.x
clinical_final$BRCA_any_mut <- clinical_final$BRCA_any_mut.x

clinical_final <- clinical_final[
  , !grepl("BRCA1_mut\\.|BRCA2_mut\\.|BRCA_any_mut\\.",
           colnames(clinical_final))
]

dim(clinical_final)
table(clinical_final$figo_stage, useNA = "ifany")
table(clinical_final$tumor_grade, useNA = "ifany")
table(clinical_final$prior_treatment, useNA = "ifany")
table(clinical_final$vital_status, useNA = "ifany")

summary(clinical_final$age_at_diagnosis)
clinical_final$age_years <- 
  clinical_final$age_at_diagnosis / 365.25

summary(clinical_final$age_years)
mutation_freq <- sort(
  table(mutation_final$Hugo_Symbol),
  decreasing = TRUE
)

head(mutation_freq, 30)
mutation_patient_freq <- aggregate(
  patient_id ~ Hugo_Symbol,
  data = mutation_final,
  FUN = function(x) length(unique(x))
)

mutation_patient_freq <- mutation_patient_freq[
  order(-mutation_patient_freq$patient_id),
]

head(mutation_patient_freq, 30)
library(ComplexHeatmap)

top_genes <- mutation_patient_freq$Hugo_Symbol[1:20]

mut_mat <- matrix(
  "",
  nrow = length(top_genes),
  ncol = length(final_ids),
  dimnames = list(top_genes, final_ids)
)

for (g in top_genes) {
  ids <- unique(
    mutation_final$patient_id[
      mutation_final$Hugo_Symbol == g
    ]
  )
  
  mut_mat[g, colnames(mut_mat) %in% ids] <- "Mutation"
}

alter_fun <- list(
  background = function(x, y, w, h)
    grid::grid.rect(x, y, w, h,
                    gp = grid::gpar(fill = "#F0F0F0", col = NA)),
  
  Mutation = function(x, y, w, h)
    grid::grid.rect(x, y, w * 0.9, h * 0.9,
                    gp = grid::gpar(fill = "#D73027", col = NA))
)

oncoPrint(
  mut_mat,
  alter_fun = alter_fun,
  col = c(Mutation = "#D73027"),
  remove_empty_columns = FALSE,
  column_title = "TCGA-OV – 277 Primary Tumors",
  row_title = "Somatic mutations",
  show_column_names = FALSE
)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


dim(cnv_final)

summary(as.vector(cnv_final))

head(rownames(cnv_final))

identical(
  colnames(cnv_final),
  clinical_final$submitter_id
)
#CNV genlerini secmek icin varyans hazirliyoruz
cnv_variance <- apply(cnv_final, 1, var, na.rm = TRUE)

cnv_variance <- sort(cnv_variance, decreasing = TRUE)

head(cnv_variance, 30)
#diploid 2 den farkli olana bakiyoruz
cnv_changed_n <- apply(
  cnv_final,
  1,
  function(x) sum(!is.na(x) & x != 2)
)

cnv_changed_n <- sort(cnv_changed_n, decreasing = TRUE)

head(cnv_changed_n, 30)
colnames(rowData(rna_tcga))
protein_genes <- rowData(rna_tcga)$gene_name[
  rowData(rna_tcga)$gene_type == "protein_coding"
]

cnv_pc <- cnv_final[
  rownames(cnv_final) %in% protein_genes,
]

dim(cnv_pc)
gene_chr <- cnv_test[, c("gene_name", "chromosome")]

autosomal_genes <- gene_chr$gene_name[
  gene_chr$chromosome %in% as.character(1:22)
]

cnv_pc_auto <- cnv_pc[
  rownames(cnv_pc) %in% autosomal_genes,
]

dim(cnv_pc_auto)
unique(cnv_test$chromosome)
autosomal_genes <- gene_chr$gene_name[
  gene_chr$chromosome %in% paste0("chr", 1:22)
]

cnv_pc_auto <- cnv_pc[
  rownames(cnv_pc) %in% autosomal_genes,
]

dim(cnv_pc_auto)
cnv_variance_auto <- apply(
  cnv_pc_auto,
  1,
  var,
  na.rm = TRUE
)

cnv_variance_auto <- sort(
  cnv_variance_auto,
  decreasing = TRUE
)

head(cnv_variance_auto, 30)
top_cnv_genes <- names(cnv_variance_auto)[1:30]

cnv_top <- cnv_pc_auto[top_cnv_genes, ]

Heatmap(
  cnv_top,
  name = "Copy number",
  show_column_names = FALSE,
  column_title = "TCGA-OV – 277 Primary Tumors",
  row_title = "Top variable CNV genes",
  cluster_columns = TRUE,
  cluster_rows = TRUE
)
#gorsel kotu oldu
library(circlize)

col_fun <- colorRamp2(
  c(0, 2, 4, 6, 10),
  c("navy", "blue", "white", "orange", "red")
)

Heatmap(
  cnv_top,
  name = "Copy number",
  col = col_fun,
  show_column_names = FALSE,
  column_title = "TCGA-OV – 277 Primary Tumors",
  row_title = "Top variable CNV genes",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  use_raster = TRUE
)
library(circlize)

col_fun <- colorRamp2(
  c(0, 2, 4, 6, 10),
  c("navy", "blue", "white", "orange", "red")
)

Heatmap(
  cnv_top,
  name = "Copy number",
  col = col_fun,
  show_column_names = FALSE,
  column_title = "TCGA-OV – 277 Primary Tumors",
  row_title = "Top variable CNV genes",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  use_raster = FALSE
)

#RNA seq katmani
rna_log <- log2(rna_final + 1)

rna_gene_names <- rowData(rna_tcga)$gene_name
rna_gene_types <- rowData(rna_tcga)$gene_type

rna_pc <- rna_log[
  rna_gene_types == "protein_coding",
]

rownames(rna_pc) <- rna_gene_names[
  rna_gene_types == "protein_coding"
]

rna_variance <- apply(rna_pc, 1, var, na.rm = TRUE)
rna_variance <- sort(rna_variance, decreasing = TRUE)

head(rna_variance, 30)
#RNA-seq expression heterogeneity
top_rna_genes <- names(rna_variance)[1:30]

rna_top <- rna_pc[top_rna_genes, ]

# Her gen için hastalar arasında Z-score
rna_top_z <- t(scale(t(rna_top)))

Heatmap(
  rna_top_z,
  name = "Expression\nZ-score",
  show_column_names = FALSE,
  column_title = "TCGA-OV – 277 Primary Tumors",
  row_title = "Top variable RNA genes",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  use_raster = FALSE
)
#combined
# RNA clustering ile hasta sırasını belirle
rna_dist <- dist(t(rna_top_z))
rna_hc <- hclust(rna_dist)

patient_order <- colnames(rna_top_z)[rna_hc$order]

length(patient_order)
head(patient_order)

# Tüm katmanların aynı hastaları içerdiğini doğrula
all(patient_order %in% colnames(rna_final))
all(patient_order %in% colnames(cnv_final))
all(patient_order %in% clinical_final$submitter_id)
all(patient_order %in% mutation_final$patient_id)
clinical_plot <- clinical_final[
  match(patient_order, clinical_final$submitter_id),
]

clinical_anno <- HeatmapAnnotation(
  FIGO = clinical_plot$figo_stage,
  Grade = clinical_plot$tumor_grade,
  BRCA = ifelse(clinical_plot$BRCA_any_mut, "Mut", "WT"),
  Age = clinical_plot$age_years,
  annotation_name_side = "left"
)

clinical_anno
# Mutation
mut_integrated <- mut_mat[, patient_order]

# CNV
cnv_integrated <- cnv_top[, patient_order]

# RNA
rna_integrated <- rna_top_z[, patient_order]

# CNV renk skalası
cnv_col <- circlize::colorRamp2(
  c(0, 2, 4, 6, 10),
  c("navy", "blue", "white", "orange", "red")
)

# RNA renk skalası
rna_col <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

ht_mut <- oncoPrint(
  mut_integrated,
  alter_fun = alter_fun,
  col = c(Mutation = "#D73027"),
  show_column_names = FALSE,
  show_pct = FALSE,
  top_annotation = clinical_anno,
  column_order = patient_order,
  cluster_columns = FALSE,
  row_title = "Mutation"
)

ht_cnv <- Heatmap(
  cnv_integrated,
  name = "CNV",
  col = cnv_col,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_order = patient_order,
  row_title = "CNV"
)

ht_rna <- Heatmap(
  rna_integrated,
  name = "RNA Z-score",
  col = rna_col,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_order = patient_order,
  row_title = "RNA expression"
)

draw(
  ht_mut %v%
    ht_cnv %v%
    ht_rna,
  merge_legends = TRUE
)
#ust hata
# Mutation
mut_integrated <- mut_mat[, patient_order]

# CNV
cnv_integrated <- cnv_top[, patient_order]

# RNA
rna_integrated <- rna_top_z[, patient_order]

# CNV renk skalası
cnv_col <- circlize::colorRamp2(
  c(0, 2, 4, 6, 10),
  c("navy", "blue", "white", "orange", "red")
)

# RNA renk skalası
rna_col <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

# Mutation
ht_mut <- oncoPrint(
  mut_integrated,
  alter_fun = alter_fun,
  col = c(Mutation = "#D73027"),
  show_column_names = FALSE,
  show_pct = FALSE,
  top_annotation = clinical_anno,
  column_order = patient_order,
  row_title = "Mutation"
)

# CNV
ht_cnv <- Heatmap(
  cnv_integrated,
  name = "CNV",
  col = cnv_col,
  show_column_names = FALSE,
  column_order = patient_order,
  row_title = "CNV"
)

# RNA
ht_rna <- Heatmap(
  rna_integrated,
  name = "RNA Z-score",
  col = rna_col,
  show_column_names = FALSE,
  column_order = patient_order,
  row_title = "RNA expression"
)

# Hepsini birleştir
draw(
  ht_mut %v%
    ht_cnv %v%
    ht_rna,
  merge_legends = TRUE
)
pdf(
  "TCGA_OV_277_MultiOmics_Integrated_Heatmap.pdf",
  width = 16,
  height = 18,
  onefile = TRUE
)

draw(
  ht_mut %v%
    ht_cnv %v%
    ht_rna,
  merge_legends = TRUE
)

dev.off()
png(
  "TCGA_OV_277_MultiOmics_Integrated_Heatmap.png",
  width = 4000,
  height = 5000,
  res = 300,
  type = "quartz"
)

draw(
  ht_mut %v%
    ht_cnv %v%
    ht_rna,
  merge_legends = TRUE
)

dev.off()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")

library(MOFA2)
packageVersion("MOFA2")
# RNA: sadece protein-coding
rna_mofa <- rna_pc[, patient_order]

# CNV: autosomal + protein-coding
cnv_mofa <- cnv_pc_auto[, patient_order]

# Mutation: gene x patient binary matrix
mut_genes <- unique(mutation_final$Hugo_Symbol)

mut_mofa <- matrix(
  0,
  nrow = length(mut_genes),
  ncol = length(patient_order),
  dimnames = list(mut_genes, patient_order)
)

for (g in mut_genes) {
  ids <- unique(
    mutation_final$patient_id[
      mutation_final$Hugo_Symbol == g
    ]
  )
  
  mut_mofa[g, patient_order %in% ids] <- 1
}

# Kontrol
dim(rna_mofa)
dim(cnv_mofa)
dim(mut_mofa)

identical(colnames(rna_mofa), patient_order)
identical(colnames(cnv_mofa), patient_order)
identical(colnames(mut_mofa), patient_order)
# RNA - top %10 variable genes
rna_var <- apply(rna_mofa, 1, var, na.rm = TRUE)
rna_keep <- rna_var >= quantile(rna_var, 0.90, na.rm = TRUE)
rna_mofa_f <- rna_mofa[rna_keep, ]

# CNV - top %10 variable genes
cnv_var <- apply(cnv_mofa, 1, var, na.rm = TRUE)
cnv_keep <- cnv_var >= quantile(cnv_var, 0.90, na.rm = TRUE)
cnv_mofa_f <- cnv_mofa[cnv_keep, ]

# Mutation - en az %2 hastada bulunan
mut_freq <- rowMeans(mut_mofa)

mut_mofa_f <- mut_mofa[
  mut_freq >= 0.02,
]

# Kontrol
dim(rna_mofa_f)
dim(cnv_mofa_f)
dim(mut_mofa_f)
#MOFA
mofa_data <- list(
  RNA = rna_mofa_f,
  CNV = cnv_mofa_f,
  Mutation = mut_mofa_f
)

mofa_model <- create_mofa(mofa_data)

mofa_model
library(MOFA2)

"create_mofa" %in% getNamespaceExports("MOFA2")
rownames(rna_mofa_f) <- paste0("RNA_", rownames(rna_mofa_f))
rownames(cnv_mofa_f) <- paste0("CNV_", rownames(cnv_mofa_f))
rownames(mut_mofa_f) <- paste0("MUT_", rownames(mut_mofa_f))

mofa_data <- list(
  RNA = rna_mofa_f,
  CNV = cnv_mofa_f,
  Mutation = mut_mofa_f
)

mofa_model <- MOFA2::create_mofa(mofa_data)

mofa_model
#error veriyor
sum(duplicated(rownames(rna_mofa_f)))
sum(duplicated(rownames(cnv_mofa_f)))
sum(duplicated(rownames(mut_mofa_f)))
unique(rownames(rna_mofa_f)[duplicated(rownames(rna_mofa_f))])
unique(rownames(cnv_mofa_f)[duplicated(rownames(cnv_mofa_f))])
unique(rownames(mut_mofa_f)[duplicated(rownames(mut_mofa_f))])
#
# NA / geçersiz gene isimlerini çıkar
cnv_mofa_f <- cnv_mofa_f[
  !is.na(rownames(cnv_mofa_f)) &
    !grepl("_NA$", rownames(cnv_mofa_f)),
]

# Aynı gene sahip CNV satırlarını ortalama ile birleştir
cnv_mofa_f <- rowsum(
  cnv_mofa_f,
  group = rownames(cnv_mofa_f),
  reorder = FALSE
) / as.vector(table(factor(
  rownames(cnv_mofa_f),
  levels = unique(rownames(cnv_mofa_f))
)))

# Kontrol
dim(cnv_mofa_f)
sum(duplicated(rownames(cnv_mofa_f)))
mofa_data <- list(
  RNA = rna_mofa_f,
  CNV = cnv_mofa_f,
  Mutation = mut_mofa_f
)

mofa_model <- MOFA2::create_mofa(mofa_data)

mofa_model
model_opts <- get_default_model_options(mofa_model)

# Başlangıç latent factor sayısı
model_opts$num_factors <- 20

# Her omics katmanı için likelihood
model_opts$likelihoods["RNA"] <- "gaussian"
model_opts$likelihoods["CNV"] <- "gaussian"
model_opts$likelihoods["Mutation"] <- "bernoulli"

model_opts
train_opts <- get_default_training_options(mofa_model)

train_opts$convergence_mode <- "medium"
train_opts$seed <- 123

train_opts
data_opts <- get_default_data_options(mofa_model)

mofa_model <- prepare_mofa(
  object = mofa_model,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

mofa_model
#egitiyoruz
mofa_model <- run_mofa(
  mofa_model,
  outfile = "TCGA_OV_MOFA_model.hdf5"
)
#
mofa_model <- run_mofa(
  mofa_model,
  outfile = "TCGA_OV_MOFA_model.hdf5",
  use_basilisk = TRUE
)
# factor 1 uyari veriyormus
plot_factor_cor(
  mofa_model,
  factors = 1
)
plot_variance_explained(
  mofa_model,
  plot_total = FALSE
)
#factor 1 CNV tarafindan domine ediliyor
r2 <- get_variance_explained(mofa_model)

r2$r2_per_factor
#DESEQ2 gerekiyor
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(rna_final),
  colData = data.frame(
    dummy = rep(1, ncol(rna_final)),
    row.names = colnames(rna_final)
  ),
  design = ~1
)

# Çok düşük ifade edilen genleri çıkar
keep <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(dds))
dds <- dds[keep, ]

vsd <- vst(dds, blind = TRUE)

rna_vst <- assay(vsd)

dim(rna_vst)
# Protein-coding genleri belirle
pc_ids <- rowData(rna_tcga)$gene_id[
  rowData(rna_tcga)$gene_type == "protein_coding"
]

pc_names <- rowData(rna_tcga)$gene_name[
  rowData(rna_tcga)$gene_type == "protein_coding"
]

# VST matrisinde bulunanları eşleştir
idx <- match(rownames(rna_vst), pc_ids)
keep_pc <- !is.na(idx)

rna_vst_pc <- rna_vst[keep_pc, ]
rownames(rna_vst_pc) <- pc_names[idx[keep_pc]]

# Aynı hasta sırası
rna_vst_pc <- rna_vst_pc[, patient_order]

# En değişken %10
rna_var <- apply(rna_vst_pc, 1, var, na.rm = TRUE)
cutoff <- quantile(rna_var, 0.90, na.rm = TRUE)

rna_mofa_new <- rna_vst_pc[rna_var >= cutoff, ]

dim(rna_mofa_new)
cnv_sample_stats <- data.frame(
  patient = colnames(cnv_mofa_f),
  median_CN = apply(cnv_mofa_f, 2, median, na.rm = TRUE),
  mean_CN   = colMeans(cnv_mofa_f, na.rm = TRUE)
)

summary(cnv_sample_stats[, c("median_CN", "mean_CN")])
# Her hastanın median CN değerini çıkar
cnv_centered <- sweep(
  cnv_mofa_f,
  2,
  apply(cnv_mofa_f, 2, median, na.rm = TRUE),
  FUN = "-"
)

# Kontrol
summary(apply(cnv_centered, 2, median, na.rm = TRUE))
cnv_var_new <- apply(cnv_centered, 1, var, na.rm = TRUE)

cutoff_cnv <- quantile(
  cnv_var_new,
  0.90,
  na.rm = TRUE
)

cnv_mofa_new <- cnv_centered[
  cnv_var_new >= cutoff_cnv,
]

dim(cnv_mofa_new)
#yanlis filtrelenmis
# Hasta sırasını eşitle
cnv_base <- cnv_pc_auto[, patient_order]

# Ploidy/general CN seviyesine göre center
cnv_centered <- sweep(
  cnv_base,
  2,
  apply(cnv_base, 2, median, na.rm = TRUE),
  FUN = "-"
)

# Tek kez top %10 variance
cnv_var_new <- apply(cnv_centered, 1, var, na.rm = TRUE)
cutoff_cnv <- quantile(cnv_var_new, 0.90, na.rm = TRUE)

cnv_mofa_new <- cnv_centered[cnv_var_new >= cutoff_cnv, ]

dim(cnv_mofa_new)

#simdi mutasyno katmani duzeltiliyor niye katkisi dusuk
mut_counts <- rowSums(mut_mofa > 0)

summary(mut_counts)

quantile(
  mut_counts,
  probs = c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1)
)
#mutasyonlar asiri sparse, o yuzden tekrarlayan mutasyonlara odaklanmak gerekiyor, yani >5 gibi 
mut_mofa_new <- mut_mofa[
  rowSums(mut_mofa) >= 5,
]

dim(mut_mofa_new)

sort(
  rowSums(mut_mofa_new),
  decreasing = TRUE
)
#yine cok geldi
mut_mofa_new <- mut_mofa[
  rowSums(mut_mofa) >= 10,
]

dim(mut_mofa_new)

sort(rowSums(mut_mofa_new), decreasing = TRUE)
# Feature isimlerini omics'e özgü yap
rownames(rna_mofa_new) <- paste0("RNA_", rownames(rna_mofa_new))
rownames(cnv_mofa_new) <- paste0("CNV_", rownames(cnv_mofa_new))
rownames(mut_mofa_new) <- paste0("MUT_", rownames(mut_mofa_new))

# Hasta sırası kontrolü
identical(colnames(rna_mofa_new), patient_order)
identical(colnames(cnv_mofa_new), patient_order)
identical(colnames(mut_mofa_new), patient_order)

# Duplicate kontrolü
sum(duplicated(rownames(rna_mofa_new)))
sum(duplicated(rownames(cnv_mofa_new)))
sum(duplicated(rownames(mut_mofa_new)))
###
# NA isimleri çıkar
cnv_mofa_new <- cnv_mofa_new[
  !is.na(rownames(cnv_mofa_new)) &
    !grepl("_NA$", rownames(cnv_mofa_new)),
]

# Duplicate genleri ortalama ile birleştir
cnv_mofa_new <- rowsum(
  cnv_mofa_new,
  group = rownames(cnv_mofa_new),
  reorder = FALSE
) / as.vector(table(factor(
  rownames(cnv_mofa_new),
  levels = unique(rownames(cnv_mofa_new))
)))

# Kontrol
dim(cnv_mofa_new)
sum(duplicated(rownames(cnv_mofa_new)))
#
mofa_data_v2 <- list(
  RNA = rna_mofa_new,
  CNV = cnv_mofa_new,
  Mutation = mut_mofa_new
)

mofa_v2 <- MOFA2::create_mofa(mofa_data_v2)

mofa_v2
###
model_opts_v2 <- get_default_model_options(mofa_v2)

model_opts_v2$num_factors <- 20
model_opts_v2$likelihoods["RNA"] <- "gaussian"
model_opts_v2$likelihoods["CNV"] <- "gaussian"
model_opts_v2$likelihoods["Mutation"] <- "bernoulli"

train_opts_v2 <- get_default_training_options(mofa_v2)
train_opts_v2$convergence_mode <- "medium"
train_opts_v2$seed <- 123

data_opts_v2 <- get_default_data_options(mofa_v2)

mofa_v2 <- prepare_mofa(
  mofa_v2,
  data_options = data_opts_v2,
  model_options = model_opts_v2,
  training_options = train_opts_v2
)

mofa_v2
#model egitiyoruz
mofa_v2 <- run_mofa(
  mofa_v2,
  outfile = "TCGA_OV_MOFA_v2.hdf5",
  use_basilisk = TRUE
)
#converge olmus ama hala 1 faktor uyarisi var
r2_v2 <- get_variance_explained(mofa_v2)

r2_v2$r2_per_factor
#%28.95 → %12.90 e dustu bu modelde ilkine gore
#mutation hala 0 - TCGA_OV'de recurrent mutation yapisi cok zayif + TP53 de ayirici degil
# o yuzden artik mutationu latent yapisina sokmuyorum sadece RNA+CNV ile latent molecular stati ogrenip sonra mutasyonlari ekleriz
#faktor 1 uyarisi?"
plot_factor_cor(
  mofa_v2,
  factors = 1
)
factor_scores <- get_factors(mofa_v2, factors = "all")$group1

dim(factor_scores)
head(factor_scores[, 1:3])
##
# Her hastada ifade edilen RNA geni sayısı
rna_detected <- colSums(rna_final[, patient_order] > 0)

# Hastanın genel CNV seviyesi
cnv_level <- apply(cnv_final[, patient_order], 2, median, na.rm = TRUE)

# Hastanın toplam somatik mutation burden'ı
mutation_burden <- sapply(patient_order, function(id) {
  sum(mutation_final$patient_id == id)
})

# Factor 1 ile korelasyon
c(
  RNA_detected = cor(factor_scores[, "Factor1"], rna_detected,
                     method = "spearman"),
  CNV_level = cor(factor_scores[, "Factor1"], cnv_level,
                  method = "spearman"),
  Mutation_burden = cor(factor_scores[, "Factor1"], mutation_burden,
                        method = "spearman")
)
r2_total_v2 <- get_variance_explained(mofa_v2)$r2_total

r2_total_v2

r2_pf <- r2_v2$r2_per_factor$group1

factor_summary <- data.frame(
  Factor = rownames(r2_pf),
  RNA = r2_pf[, "RNA"],
  CNV = r2_pf[, "CNV"],
  Total_RNA_CNV = r2_pf[, "RNA"] + r2_pf[, "CNV"]
)

factor_summary <- factor_summary[
  order(-factor_summary$Total_RNA_CNV),
]

factor_summary
# Tabloyu kaydet
write.csv(
  factor_summary,
  "MOFA_v2_factor_variance_summary.csv",
  row.names = FALSE
)

# Görselleştir
library(ggplot2)

ggplot(factor_summary,
       aes(x = reorder(Factor, Total_RNA_CNV),
           y = Total_RNA_CNV)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "MOFA Factor",
    y = "Variance explained (%)",
    title = "Variance explained by MOFA factors",
    subtitle = "Combined RNA + CNV variance"
  ) +
  theme_classic()
library(tidyr)
library(ggplot2)

factor_long <- factor_summary |>
  select(Factor, RNA, CNV) |>
  pivot_longer(
    cols = c(RNA, CNV),
    names_to = "Omics",
    values_to = "Variance"
  )

factor_long$Factor <- factor(
  factor_long$Factor,
  levels = rev(factor_summary$Factor)
)

p <- ggplot(factor_long,
            aes(x = Factor, y = Variance, fill = Omics)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "MOFA Factor",
    y = "Variance explained (%)",
    title = "MOFA factor contributions by omics layer"
  ) +
  theme_classic()

p

ggsave(
  "MOFA_v2_factor_variance_RNA_CNV.pdf",
  p,
  width = 7,
  height = 6
)
library(dplyr)
library(tidyr)
library(ggplot2)

factor_long <- factor_summary |>
  dplyr::select(Factor, RNA, CNV) |>
  tidyr::pivot_longer(
    cols = c(RNA, CNV),
    names_to = "Omics",
    values_to = "Variance"
  )
head(factor_long)
factor_long$Factor <- factor(
  factor_long$Factor,
  levels = rev(factor_summary$Factor)
)

p <- ggplot(
  factor_long,
  aes(x = Factor, y = Variance, fill = Omics)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "MOFA Factor",
    y = "Variance explained (%)",
    title = "MOFA factor contributions by omics layer"
  ) +
  theme_classic()

ggsave(
  "MOFA_v2_factor_variance_RNA_CNV.pdf",
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  "MOFA_v2_factor_variance_RNA_CNV.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
weights_v2 <- get_weights(
  mofa_v2,
  views = "all",
  factors = "all"
)

names(weights_v2)

dim(weights_v2$RNA)
dim(weights_v2$CNV)
dim(weights_v2$Mutation)
# Factor1 RNA weights
rna_f1 <- sort(weights_v2$RNA[, "Factor1"], decreasing = TRUE)

cat("RNA positive:\n")
print(head(rna_f1, 10))

cat("\nRNA negative:\n")
print(tail(rna_f1, 10))


# Factor1 CNV weights
cnv_f1 <- sort(weights_v2$CNV[, "Factor1"], decreasing = TRUE)

cat("\nCNV positive:\n")
print(head(cnv_f1, 10))

cat("\nCNV negative:\n")
print(tail(cnv_f1, 10))
#top features
top_features <- function(weight_matrix, n = 10) {
  
  do.call(rbind, lapply(colnames(weight_matrix), function(f) {
    
    w <- weight_matrix[, f]
    idx <- order(abs(w), decreasing = TRUE)[1:n]
    
    data.frame(
      Factor = f,
      Feature = rownames(weight_matrix)[idx],
      Weight = w[idx]
    )
    
  }))
}

top_rna <- top_features(weights_v2$RNA, 10)
top_cnv <- top_features(weights_v2$CNV, 10)

head(top_rna, 20)
write.csv(
  top_rna,
  "MOFA_v2_top10_RNA_features_per_factor.csv",
  row.names = FALSE
)

write.csv(
  top_cnv,
  "MOFA_v2_top10_CNV_features_per_factor.csv",
  row.names = FALSE
)
plot_weights_heatmap(
  mofa_v2,
  view = "RNA",
  factors = 1:10,
  nfeatures = 5,
  scale = TRUE
)
plot_weights_heatmap(
  mofa_v2,
  view = "RNA",
  factors = 1:10,
  nfeatures = 5,
  scale = "row"
)
top10_factors <- paste0("Factor", 1:10)

genes_top <- unique(
  top_rna$Feature[top_rna$Factor %in% top10_factors]
)

length(genes_top)

library(ComplexHeatmap)

# Prefix zaten mevcut, doğrudan eşleştiriyoruz
rna_weight_heatmap <- weights_v2$RNA[
  genes_top,
  top10_factors,
  drop = FALSE
]

Heatmap(
  rna_weight_heatmap,
  name = "Weight",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = "Top RNA features across MOFA factors",
  row_title = "Genes"
)
rownames(rna_weight_heatmap) <- sub("^RNA_", "", rownames(rna_weight_heatmap))

Heatmap(
  rna_weight_heatmap,
  name = "Weight",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = "Top RNA features across MOFA factors",
  row_title = "Genes"
)
#CNV heatmap
cnv_genes_top <- unique(
  top_cnv$Feature[top_cnv$Factor %in% top10_factors]
)

length(cnv_genes_top)
cnv_weight_heatmap <- weights_v2$CNV[
  cnv_genes_top,
  top10_factors,
  drop = FALSE
]

# CNV_ prefixini kaldır
rownames(cnv_weight_heatmap) <- sub(
  "^CNV_", "",
  rownames(cnv_weight_heatmap)
)

Heatmap(
  cnv_weight_heatmap,
  name = "Weight",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = "Top CNV features across MOFA factors",
  row_title = "Genes"
)
summary(abs(weights_v2$RNA[, "Factor1"]))
quantile(
  abs(weights_v2$RNA[, "Factor1"]),
  probs = c(0.90, 0.95, 0.975, 0.99)
)
factor_gene_sets <- lapply(top10_factors, function(f) {
  
  w <- weights_v2$RNA[, f]
  cutoff <- quantile(abs(w), 0.95)
  
  data.frame(
    Factor = f,
    Gene = sub("^RNA_", "", names(w)[abs(w) >= cutoff]),
    Weight = w[abs(w) >= cutoff],
    Direction = ifelse(w[abs(w) >= cutoff] > 0, "Positive", "Negative")
  )
})

factor_gene_sets <- do.call(rbind, factor_gene_sets)

table(
  factor_gene_sets$Factor,
  factor_gene_sets$Direction
)
library(clusterProfiler)
library(org.Hs.eg.db)
go_results <- compareCluster(
  Gene ~ Factor + Direction,
  data = factor_gene_sets,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_results
dotplot(
  go_results,
  showCategory = 3,
  includeAll = FALSE
) +
  ggtitle("GO Biological Process enrichment of MOFA factors")
p_go <- dotplot(
  go_results,
  showCategory = 3,
  includeAll = FALSE
) +
  ggtitle("GO Biological Process enrichment of MOFA factors")

ggsave(
  "MOFA_v2_GO_BP_enrichment.pdf",
  plot = p_go,
  width = 14,
  height = 9
)

ggsave(
  "MOFA_v2_GO_BP_enrichment.png",
  plot = p_go,
  width = 14,
  height = 9,
  dpi = 300
)
write.csv(
  as.data.frame(go_results),
  "MOFA_v2_GO_BP_enrichment_results.csv",
  row.names = FALSE
)



go_df <- as.data.frame(go_results)

go_top5 <- go_df |>
  dplyr::group_by(Factor, Direction) |>
  dplyr::arrange(p.adjust, .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::ungroup() |>
  dplyr::select(
    Factor,
    Direction,
    Description,
    GeneRatio,
    FoldEnrichment,
    p.adjust,
    Count
  )

go_top5
go_top5$Factor_num <- as.numeric(
  sub("Factor", "", go_top5$Factor)
)

go_top5 <- go_top5 |>
  dplyr::arrange(Factor_num, Direction, p.adjust)

write.csv(
  go_top5,
  "MOFA_v2_top5_GO_per_factor_direction.csv",
  row.names = FALSE
)

go_top5[, c(
  "Factor", "Direction", "Description",
  "FoldEnrichment", "p.adjust"
)]
print(
  go_top5[, c(
    "Factor", "Direction", "Description",
    "FoldEnrichment", "p.adjust"
  )],
  n = Inf
)
#molekuler veriyi hastayla eslestirme
factor_clinical <- data.frame(
  submitter_id = rownames(factor_scores),
  factor_scores,
  row.names = NULL
)

factor_clinical <- merge(
  factor_clinical,
  clinical_final,
  by = "submitter_id",
  all.x = TRUE,
  sort = FALSE
)

# MOFA hasta sırasını geri koru
factor_clinical <- factor_clinical[
  match(rownames(factor_scores), factor_clinical$submitter_id),
]

dim(factor_clinical)

# Kontrol
factor_clinical[
  1:5,
  c("submitter_id", "Factor1", "Factor2",
    "age_years", "figo_stage", "tumor_grade",
    "BRCA_any_mut", "vital_status")
]
age_cor <- sapply(
  paste0("Factor", 1:10),
  function(f) {
    cor(
      factor_clinical[[f]],
      factor_clinical$age_years,
      method = "spearman",
      use = "complete.obs"
    )
  }
)

sort(age_cor, decreasing = TRUE)
age_results <- do.call(rbind, lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    test <- cor.test(
      factor_clinical[[f]],
      factor_clinical$age_years,
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Factor = f,
      Rho = unname(test$estimate),
      P_value = test$p.value
    )
  }
))

age_results$FDR <- p.adjust(
  age_results$P_value,
  method = "BH"
)

age_results <- age_results[
  order(age_results$FDR),
]

age_results
factor_clinical$FIGO_main <- dplyr::case_when(
  grepl("Stage I[^I]", factor_clinical$figo_stage) ~ "I",
  grepl("Stage II[^I]", factor_clinical$figo_stage) ~ "II",
  grepl("Stage III", factor_clinical$figo_stage) ~ "III",
  grepl("Stage IV", factor_clinical$figo_stage) ~ "IV",
  TRUE ~ NA_character_
)

table(
  factor_clinical$figo_stage,
  factor_clinical$FIGO_main,
  useNA = "ifany"
)
factor_clinical$FIGO_main <- dplyr::case_when(
  grepl("^Stage IV",  factor_clinical$figo_stage) ~ "IV",
  grepl("^Stage III", factor_clinical$figo_stage) ~ "III",
  grepl("^Stage II",  factor_clinical$figo_stage) ~ "II",
  grepl("^Stage I",   factor_clinical$figo_stage) ~ "I",
  TRUE ~ NA_character_
)

table(
  factor_clinical$figo_stage,
  factor_clinical$FIGO_main,
  useNA = "ifany"
)
figo_34 <- factor_clinical[
  factor_clinical$FIGO_main %in% c("III", "IV"),
]

figo_results <- do.call(rbind, lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    test <- wilcox.test(
      figo_34[[f]] ~ figo_34$FIGO_main,
      exact = FALSE
    )
    
    data.frame(
      Factor = f,
      Median_III = median(
        figo_34[[f]][figo_34$FIGO_main == "III"],
        na.rm = TRUE
      ),
      Median_IV = median(
        figo_34[[f]][figo_34$FIGO_main == "IV"],
        na.rm = TRUE
      ),
      P_value = test$p.value
    )
  }
))

figo_results$FDR <- p.adjust(
  figo_results$P_value,
  method = "BH"
)

figo_results <- figo_results[
  order(figo_results$FDR),
]

figo_results
grade_23 <- factor_clinical[
  factor_clinical$tumor_grade %in% c("G2", "G3"),
]

grade_results <- do.call(rbind, lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    test <- wilcox.test(
      grade_23[[f]] ~ grade_23$tumor_grade,
      exact = FALSE
    )
    
    data.frame(
      Factor = f,
      Median_G2 = median(
        grade_23[[f]][grade_23$tumor_grade == "G2"],
        na.rm = TRUE
      ),
      Median_G3 = median(
        grade_23[[f]][grade_23$tumor_grade == "G3"],
        na.rm = TRUE
      ),
      P_value = test$p.value
    )
  }
))

grade_results$FDR <- p.adjust(
  grade_results$P_value,
  method = "BH"
)

grade_results <- grade_results[
  order(grade_results$FDR),
]

grade_results
brca_results <- do.call(rbind, lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    test <- wilcox.test(
      factor_clinical[[f]] ~ factor_clinical$BRCA_any_mut,
      exact = FALSE
    )
    
    data.frame(
      Factor = f,
      Median_BRCA_WT = median(
        factor_clinical[[f]][!factor_clinical$BRCA_any_mut],
        na.rm = TRUE
      ),
      Median_BRCA_Mut = median(
        factor_clinical[[f]][factor_clinical$BRCA_any_mut],
        na.rm = TRUE
      ),
      P_value = test$p.value
    )
  }
))

brca_results$FDR <- p.adjust(
  brca_results$P_value,
  method = "BH"
)

brca_results <- brca_results[
  order(brca_results$FDR),
]

brca_results
grep(
  "survival|death|follow|progress|recurr|days_to",
  colnames(factor_clinical),
  value = TRUE,
  ignore.case = TRUE
)
factor_clinical$OS_event <- ifelse(
  factor_clinical$vital_status == "Dead", 1, 0
)

factor_clinical$OS_days <- ifelse(
  factor_clinical$OS_event == 1,
  factor_clinical$days_to_death,
  factor_clinical$days_to_last_follow_up
)

summary(factor_clinical$OS_days)
table(factor_clinical$OS_event, useNA = "ifany")
sum(is.na(factor_clinical$OS_days))
#cox to ask which molecular factors related w patient survival
library(survival)

cox_results <- do.call(rbind, lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    fit <- coxph(
      Surv(OS_days, OS_event) ~ factor_clinical[[f]],
      data = factor_clinical
    )
    
    s <- summary(fit)
    
    data.frame(
      Factor = f,
      HR = s$coefficients[1, "exp(coef)"],
      CI_lower = s$conf.int[1, "lower .95"],
      CI_upper = s$conf.int[1, "upper .95"],
      P_value = s$coefficients[1, "Pr(>|z|)"]
    )
  }
))

cox_results$FDR <- p.adjust(
  cox_results$P_value,
  method = "BH"
)

cox_results <- cox_results[
  order(cox_results$FDR),
]

cox_results
association_matrix <- matrix(
  NA,
  nrow = 10,
  ncol = 5,
  dimnames = list(
    paste0("Factor", 1:10),
    c("Age", "FIGO III-IV", "Grade G2-G3", "BRCA", "Overall Survival")
  )
)

# Age = Spearman rho
association_matrix[age_results$Factor, "Age"] <-
  age_results$Rho

# FIGO = Stage IV median - Stage III median
association_matrix[figo_results$Factor, "FIGO III-IV"] <-
  figo_results$Median_IV - figo_results$Median_III

# Grade = G3 median - G2 median
association_matrix[grade_results$Factor, "Grade G2-G3"] <-
  grade_results$Median_G3 - grade_results$Median_G2

# BRCA = mutant median - WT median
association_matrix[brca_results$Factor, "BRCA"] <-
  brca_results$Median_BRCA_Mut - brca_results$Median_BRCA_WT

# Survival = log(HR)
association_matrix[cox_results$Factor, "Overall Survival"] <-
  log(cox_results$HR)

round(association_matrix, 2)
# FDR matrisi
fdr_matrix <- matrix(
  NA,
  nrow = 10,
  ncol = 5,
  dimnames = dimnames(association_matrix)
)

fdr_matrix[age_results$Factor, "Age"] <- age_results$FDR
fdr_matrix[figo_results$Factor, "FIGO III-IV"] <- figo_results$FDR
fdr_matrix[grade_results$Factor, "Grade G2-G3"] <- grade_results$FDR
fdr_matrix[brca_results$Factor, "BRCA"] <- brca_results$FDR
fdr_matrix[cox_results$Factor, "Overall Survival"] <- cox_results$FDR

# Her klinik değişkeni kendi içinde standardize et
association_scaled <- scale(association_matrix)

# FDR < 0.05 için yıldız
sig_labels <- ifelse(fdr_matrix < 0.05, "*", "")

library(ComplexHeatmap)
library(circlize)

Heatmap(
  association_scaled,
  name = "Association\n(Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      sig_labels[i, j],
      x, y,
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  },
  column_title = "Clinical associations of MOFA factors",
  row_title = "MOFA Factors"
)
library(ComplexHeatmap)
library(circlize)

# Hücrelerde gösterilecek gerçek effect size
effect_labels <- matrix(
  sprintf("%.2f", association_matrix),
  nrow = nrow(association_matrix),
  dimnames = dimnames(association_matrix)
)

# Anlamlı olanlara yıldız ekle
effect_labels[fdr_matrix < 0.05] <-
  paste0(effect_labels[fdr_matrix < 0.05], "*")

# Heatmap
ht_clinical <- Heatmap(
  association_scaled,
  name = "Relative\nassociation",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#3B4CC0", "white", "#B40426")
  ),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_title = "Clinical associations of MOFA factors",
  row_title = "MOFA factors",
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      effect_labels[i, j],
      x, y,
      gp = grid::gpar(
        fontsize = 9,
        fontface = ifelse(fdr_matrix[i, j] < 0.05, "bold", "plain")
      )
    )
  },
  
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9)
)

draw(ht_clinical)
#####
#artik hastalara state atiyoruz

library(mclust)

state_data <- factor_scores[, paste0("Factor", 1:10)]

# Her faktörü aynı ölçeğe getir
state_data_scaled <- scale(state_data)

dim(state_data_scaled)
summary(state_data_scaled[, 1:3])
mclust_fit <- Mclust(
  state_data_scaled,
  G = 2:8
)

summary(mclust_fit)

plot(
  mclust_fit,
  what = "BIC"
)
summary(mclust_fit)
state_probabilities <- mclust_fit$z

colnames(state_probabilities) <- paste0("State", 1:5)
rownames(state_probabilities) <- rownames(state_data_scaled)

head(round(state_probabilities, 3))

# Her hastanın en yüksek state olasılığı
summary(apply(state_probabilities, 1, max))
patient_state <- mclust_fit$classification

state_factor_profile <- aggregate(
  state_data_scaled,
  by = list(State = patient_state),
  FUN = mean
)

round(state_factor_profile, 2)
library(ComplexHeatmap)
library(circlize)

# Hastaları atandıkları state'e göre sırala
patient_state <- mclust_fit$classification
ord <- order(patient_state)

Heatmap(
  state_probabilities[ord, ],
  name = "Probability",
  col = colorRamp2(
    c(0, 0.5, 1),
    c("white", "orange", "red")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_split = factor(
    patient_state[ord],
    levels = 1:5,
    labels = paste0("State ", 1:5)
  ),
  column_title = "Probabilistic Patient-State Membership",
  row_title = "Patients",
  column_names_gp = grid::gpar(fontsize = 11)
)
library(ComplexHeatmap)
library(circlize)

state_mat <- as.matrix(
  state_factor_profile[, paste0("Factor", 1:10)]
)

rownames(state_mat) <- paste0("State ", state_factor_profile$State)

Heatmap(
  state_mat,
  name = "Mean\nFactor score",
  col = colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "Molecular Profiles of Patient States",
  row_title = "Patient States",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      sprintf("%.2f", state_mat[i, j]),
      x, y,
      gp = grid::gpar(fontsize = 9)
    )
  }
)
#klinik varyasyonlarin state based dagilimi
factor_clinical$State <- factor(
  patient_state[
    match(factor_clinical$submitter_id, names(patient_state))
  ],
  levels = 1:5
)

# Hasta sayıları
table(factor_clinical$State)

# Yaş
aggregate(age_years ~ State, factor_clinical, 
          function(x) c(n = length(x),
                        median = median(x, na.rm = TRUE),
                        mean = mean(x, na.rm = TRUE)))

# BRCA
table(factor_clinical$State,
      factor_clinical$BRCA_any_mut)

# FIGO
table(factor_clinical$State,
      factor_clinical$FIGO_main)

# Grade
table(factor_clinical$State,
      factor_clinical$tumor_grade)
library(ggplot2)

ggplot(
  factor_clinical,
  aes(x = State, y = age_years, fill = State)
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.4,
    alpha = 0.45
  ) +
  labs(
    title = "Age Distribution Across Patient States",
    x = "Patient State",
    y = "Age (years)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none"
  )
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# -------------------------
# A) AGE
# -------------------------
p_age <- ggplot(
  factor_clinical,
  aes(x = State, y = age_years, fill = State)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
  labs(
    title = "A. Age",
    x = "Patient State",
    y = "Age (years)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")


# -------------------------
# B) BRCA
# -------------------------
brca_plot_df <- factor_clinical %>%
  filter(!is.na(BRCA_any_mut)) %>%
  count(State, BRCA_any_mut) %>%
  group_by(State) %>%
  mutate(prop = n / sum(n))

p_brca <- ggplot(
  brca_plot_df,
  aes(x = State, y = prop, fill = BRCA_any_mut)
) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "B. BRCA mutation",
    x = "Patient State",
    y = "Patients (%)",
    fill = "BRCA mutated"
  ) +
  theme_classic(base_size = 12)


# -------------------------
# C) FIGO
# -------------------------
figo_plot_df <- factor_clinical %>%
  filter(!is.na(FIGO_main)) %>%
  count(State, FIGO_main) %>%
  group_by(State) %>%
  mutate(prop = n / sum(n))

p_figo <- ggplot(
  figo_plot_df,
  aes(x = State, y = prop, fill = FIGO_main)
) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "C. FIGO stage",
    x = "Patient State",
    y = "Patients (%)",
    fill = "FIGO"
  ) +
  theme_classic(base_size = 12)


# -------------------------
# D) GRADE
# -------------------------
grade_plot_df <- factor_clinical %>%
  filter(!is.na(tumor_grade)) %>%
  count(State, tumor_grade) %>%
  group_by(State) %>%
  mutate(prop = n / sum(n))

p_grade <- ggplot(
  grade_plot_df,
  aes(x = State, y = prop, fill = tumor_grade)
) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "D. Tumor grade",
    x = "Patient State",
    y = "Patients (%)",
    fill = "Grade"
  ) +
  theme_classic(base_size = 12)


# -------------------------
# COMBINE
# -------------------------
(p_age | p_brca) /
  (p_figo | p_grade) +
  plot_annotation(
    title = "Clinical Characteristics of Probabilistic Patient States"
  )
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ============================================================
# A) AGE
# ============================================================

p_age <- ggplot(
  factor_clinical,
  aes(x = State, y = age_years, fill = State)
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,
    alpha = 0.35
  ) +
  labs(
    title = "A. Age",
    x = "Patient State",
    y = "Age (years)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )


# ============================================================
# B) BRCA MUTATION
# ============================================================

brca_plot_df <- factor_clinical %>%
  dplyr::filter(!is.na(BRCA_any_mut)) %>%
  dplyr::count(State, BRCA_any_mut) %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p_brca <- ggplot(
  brca_plot_df,
  aes(
    x = State,
    y = prop,
    fill = factor(BRCA_any_mut)
  )
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "B. BRCA mutation",
    x = "Patient State",
    y = "Patients (%)",
    fill = "BRCA mutated"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )


# ============================================================
# C) FIGO STAGE
# ============================================================

figo_plot_df <- factor_clinical %>%
  dplyr::filter(!is.na(FIGO_main)) %>%
  dplyr::count(State, FIGO_main) %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p_figo <- ggplot(
  figo_plot_df,
  aes(
    x = State,
    y = prop,
    fill = FIGO_main
  )
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "C. FIGO stage",
    x = "Patient State",
    y = "Patients (%)",
    fill = "FIGO stage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )


# ============================================================
# D) TUMOR GRADE
# ============================================================

grade_plot_df <- factor_clinical %>%
  dplyr::filter(!is.na(tumor_grade)) %>%
  dplyr::count(State, tumor_grade) %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p_grade <- ggplot(
  grade_plot_df,
  aes(
    x = State,
    y = prop,
    fill = tumor_grade
  )
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = "D. Tumor grade",
    x = "Patient State",
    y = "Patients (%)",
    fill = "Grade"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )


# ============================================================
# COMBINE ALL PANELS
# ============================================================

clinical_state_panel <-
  (p_age | p_brca) /
  (p_figo | p_grade) +
  patchwork::plot_annotation(
    title = "Clinical Characteristics of Probabilistic Patient States",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      )
    )
  )

clinical_state_panel
#clinical state ssignmentleri significant mi bakiyrozu
# ============================================================
# CLINICAL ASSOCIATION WITH 5 PATIENT STATES
# ============================================================

library(dplyr)
library(ggplot2)

# 1. AGE — Kruskal-Wallis
age_test <- kruskal.test(
  age_years ~ State,
  data = factor_clinical
)

# 2. BRCA — Fisher exact
brca_tab <- table(
  factor_clinical$State,
  factor_clinical$BRCA_any_mut
)

brca_test <- fisher.test(brca_tab)

# 3. FIGO — Fisher exact with simulation
figo_tab <- table(
  factor_clinical$State,
  factor_clinical$FIGO_main
)

figo_test <- fisher.test(
  figo_tab,
  simulate.p.value = TRUE,
  B = 100000
)

# 4. GRADE — Fisher exact with simulation
grade_tab <- table(
  factor_clinical$State,
  factor_clinical$tumor_grade
)

grade_test <- fisher.test(
  grade_tab,
  simulate.p.value = TRUE,
  B = 100000
)


# ============================================================
# COLLECT RESULTS
# ============================================================

clinical_tests <- data.frame(
  Clinical_feature = c(
    "Age",
    "BRCA mutation",
    "FIGO stage",
    "Tumor grade"
  ),
  Test = c(
    "Kruskal-Wallis",
    "Fisher exact",
    "Fisher exact",
    "Fisher exact"
  ),
  P_value = c(
    age_test$p.value,
    brca_test$p.value,
    figo_test$p.value,
    grade_test$p.value
  )
)

# Multiple-testing correction
clinical_tests$FDR <- p.adjust(
  clinical_tests$P_value,
  method = "BH"
)

clinical_tests$minus_log10_FDR <- -log10(clinical_tests$FDR)

clinical_tests


# ============================================================
# VISUALIZATION
# ============================================================

ggplot(
  clinical_tests,
  aes(
    x = reorder(Clinical_feature, minus_log10_FDR),
    y = minus_log10_FDR
  )
) +
  geom_col(width = 0.65) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text(
    aes(
      label = paste0(
        "FDR = ",
        formatC(FDR, format = "g", digits = 3)
      )
    ),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  expand_limits(
    y = max(clinical_tests$minus_log10_FDR) * 1.25
  ) +
  labs(
    title = "Clinical Associations with Probabilistic Patient States",
    x = NULL,
    y = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 13)
#simdi statelere biyolojik anlamlilik veriyoruz
state_factor_long <- state_factor_profile %>%
  tidyr::pivot_longer(
    cols = starts_with("Factor"),
    names_to = "Factor",
    values_to = "MeanScore"
  ) %>%
  group_by(State) %>%
  arrange(desc(abs(MeanScore)), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup()

state_factor_long
library(dplyr)

# ============================================================
# PATIENT - STATE MATCHING
# ============================================================

state_assignment <- data.frame(
  submitter_id = names(patient_state),
  State = patient_state
)

# RNA hasta-state eşleşmesi
state_rna <- state_assignment$State[
  match(colnames(rna_mofa_new),
        state_assignment$submitter_id)
]

# CNV hasta-state eşleşmesi
state_cnv <- state_assignment$State[
  match(colnames(cnv_mofa_new),
        state_assignment$submitter_id)
]


# ============================================================
# FUNCTION: ONE STATE vs ALL OTHER STATES
# ============================================================

get_state_signature <- function(mat, states, target_state) {
  
  in_state <- states == target_state
  
  mean_state <- rowMeans(
    mat[, in_state, drop = FALSE],
    na.rm = TRUE
  )
  
  mean_others <- rowMeans(
    mat[, !in_state, drop = FALSE],
    na.rm = TRUE
  )
  
  data.frame(
    Feature = rownames(mat),
    Mean_State = mean_state,
    Mean_Others = mean_others,
    Difference = mean_state - mean_others
  ) %>%
    arrange(desc(abs(Difference)))
}


# ============================================================
# RNA SIGNATURES — STATES 1-5
# ============================================================

rna_state_signatures <- lapply(
  1:5,
  function(s) {
    
    result <- get_state_signature(
      rna_mofa_new,
      state_rna,
      s
    )
    
    result$State <- paste0("State", s)
    result
  }
)

names(rna_state_signatures) <- paste0("State", 1:5)


# ============================================================
# CNV SIGNATURES — STATES 1-5
# ============================================================

cnv_state_signatures <- lapply(
  1:5,
  function(s) {
    
    result <- get_state_signature(
      cnv_mofa_new,
      state_cnv,
      s
    )
    
    result$State <- paste0("State", s)
    result
  }
)

names(cnv_state_signatures) <- paste0("State", 1:5)


# ============================================================
# TOP 10 FEATURES FOR EACH STATE
# ============================================================

for (s in 1:5) {
  
  cat("\n============================\n")
  cat("STATE", s, "- RNA\n")
  cat("============================\n")
  
  print(
    head(rna_state_signatures[[s]], 10)
  )
  
  cat("\nSTATE", s, "- CNV\n")
  
  print(
    head(cnv_state_signatures[[s]], 10)
  )
}
#gorsel
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# TOP 10 RNA + CNV FEATURES FROM EACH STATE
# ============================================================

rna_top <- bind_rows(rna_state_signatures) %>%
  group_by(State) %>%
  slice_max(order_by = abs(Difference), n = 10, with_ties = FALSE) %>%
  ungroup()

cnv_top <- bind_rows(cnv_state_signatures) %>%
  group_by(State) %>%
  slice_max(order_by = abs(Difference), n = 10, with_ties = FALSE) %>%
  ungroup()


# ============================================================
# CREATE STATE × FEATURE MATRICES
# For each selected feature calculate Difference in ALL states
# ============================================================

rna_features <- unique(rna_top$Feature)
cnv_features <- unique(cnv_top$Feature)

rna_sig_matrix <- sapply(1:5, function(s) {
  
  sig <- rna_state_signatures[[s]]
  
  sig$Difference[
    match(rna_features, sig$Feature)
  ]
})

rownames(rna_sig_matrix) <- sub("^RNA_", "", rna_features)
colnames(rna_sig_matrix) <- paste0("State ", 1:5)


cnv_sig_matrix <- sapply(1:5, function(s) {
  
  sig <- cnv_state_signatures[[s]]
  
  sig$Difference[
    match(cnv_features, sig$Feature)
  ]
})

rownames(cnv_sig_matrix) <- sub("^CNV_", "", cnv_features)
colnames(cnv_sig_matrix) <- paste0("State ", 1:5)


# ============================================================
# HEATMAP COLOR SCALES
# ============================================================

rna_lim <- max(abs(rna_sig_matrix), na.rm = TRUE)

rna_col <- colorRamp2(
  c(-rna_lim, 0, rna_lim),
  c("#2166AC", "white", "#B2182B")
)

cnv_lim <- max(abs(cnv_sig_matrix), na.rm = TRUE)

cnv_col <- colorRamp2(
  c(-cnv_lim, 0, cnv_lim),
  c("#2166AC", "white", "#B2182B")
)


# ============================================================
# RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  rna_sig_matrix,
  name = "RNA\nDifference",
  col = rna_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "RNA State Signatures",
  
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  width = unit(5, "cm")
)


# ============================================================
# CNV HEATMAP
# ============================================================

ht_cnv <- Heatmap(
  cnv_sig_matrix,
  name = "CNV\nDifference",
  col = cnv_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "CNV State Signatures",
  
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  width = unit(5, "cm")
)


# ============================================================
# DRAW TOGETHER
# ============================================================

draw(
  ht_rna + ht_cnv,
  column_title = "Molecular Signatures of the Five Patient States"
)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

# ============================================================
# TOP 5 RNA FEATURE / STATE
# ============================================================

rna_top <- bind_rows(rna_state_signatures) %>%
  group_by(State) %>%
  slice_max(abs(Difference), n = 5, with_ties = FALSE) %>%
  ungroup()

rna_features <- unique(rna_top$Feature)

rna_state_mat <- sapply(1:5, function(s) {
  
  x <- rna_state_signatures[[s]]
  
  x$Difference[
    match(rna_features, x$Feature)
  ]
})

rownames(rna_state_mat) <- sub("^RNA_", "", rna_features)
colnames(rna_state_mat) <- paste0("State ", 1:5)


# ============================================================
# TOP 5 CNV FEATURE / STATE
# ============================================================

cnv_top <- bind_rows(cnv_state_signatures) %>%
  group_by(State) %>%
  slice_max(abs(Difference), n = 5, with_ties = FALSE) %>%
  ungroup()

cnv_features <- unique(cnv_top$Feature)

cnv_state_mat <- sapply(1:5, function(s) {
  
  x <- cnv_state_signatures[[s]]
  
  x$Difference[
    match(cnv_features, x$Feature)
  ]
})

rownames(cnv_state_mat) <- sub("^CNV_", "", cnv_features)
colnames(cnv_state_mat) <- paste0("State ", 1:5)


# ============================================================
# RNA HEATMAP
# ============================================================

rna_lim <- max(abs(rna_state_mat))

ht_rna <- Heatmap(
  rna_state_mat,
  name = "RNA\nΔ",
  col = colorRamp2(
    c(-rna_lim, 0, rna_lim),
    c("blue", "white", "red")
  ),
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "RNA signatures",
  row_title = "Genes",
  
  column_names_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 9),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.1f", rna_state_mat[i,j]),
      x, y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# CNV HEATMAP
# ============================================================

cnv_lim <- max(abs(cnv_state_mat))

ht_cnv <- Heatmap(
  cnv_state_mat,
  name = "CNV\nΔ",
  col = colorRamp2(
    c(-cnv_lim, 0, cnv_lim),
    c("blue", "white", "red")
  ),
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "CNV signatures",
  row_title = "Genes",
  
  column_names_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 9),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.1f", cnv_state_mat[i,j]),
      x, y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# DRAW
# ============================================================
# RNA üstte, CNV altta
ht_list <- ht_rna %v% ht_cnv

draw(
  ht_list,
  column_title = "Molecular Signatures of the Five Patient States"
)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. TOP 5 RNA FEATURES PER STATE
# ============================================================

rna_top <- bind_rows(rna_state_signatures) %>%
  group_by(State) %>%
  slice_max(
    order_by = abs(Difference),
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()

rna_features <- unique(rna_top$Feature)

# Her seçilmiş RNA geni için State 1-5 Difference değerleri
rna_state_mat <- sapply(1:5, function(s) {
  
  x <- rna_state_signatures[[s]]
  
  x$Difference[
    match(rna_features, x$Feature)
  ]
})

rownames(rna_state_mat) <- sub("^RNA_", "", rna_features)
colnames(rna_state_mat) <- paste0("State ", 1:5)


# ============================================================
# 2. TOP 5 CNV FEATURES PER STATE
# ============================================================

cnv_top <- bind_rows(cnv_state_signatures) %>%
  group_by(State) %>%
  slice_max(
    order_by = abs(Difference),
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()

cnv_features <- unique(cnv_top$Feature)

# Her seçilmiş CNV geni için State 1-5 Difference değerleri
cnv_state_mat <- sapply(1:5, function(s) {
  
  x <- cnv_state_signatures[[s]]
  
  x$Difference[
    match(cnv_features, x$Feature)
  ]
})

rownames(cnv_state_mat) <- sub("^CNV_", "", cnv_features)
colnames(cnv_state_mat) <- paste0("State ", 1:5)


# ============================================================
# 3. RNA COLOR SCALE
# ============================================================

rna_lim <- max(abs(rna_state_mat), na.rm = TRUE)

rna_col <- colorRamp2(
  c(-rna_lim, 0, rna_lim),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 4. CNV COLOR SCALE
# ============================================================

cnv_lim <- max(abs(cnv_state_mat), na.rm = TRUE)

cnv_col <- colorRamp2(
  c(-cnv_lim, 0, cnv_lim),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 5. RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  
  rna_state_mat,
  
  name = "RNA\nΔ",
  col = rna_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # ----- CLEAR PANEL TITLE -----
  column_title = "RNA EXPRESSION SIGNATURES",
  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),
  
  row_title = "RNA genes",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  # State labels
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  # Numbers inside cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      sprintf("%.1f", rna_state_mat[i, j]),
      x,
      y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# 6. CNV HEATMAP
# ============================================================

ht_cnv <- Heatmap(
  
  cnv_state_mat,
  
  name = "CNV\nΔ",
  col = cnv_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # ----- CLEAR PANEL TITLE -----
  column_title = "COPY NUMBER (CNV) SIGNATURES",
  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),
  
  row_title = "CNV genes",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  # State labels
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  # Numbers inside cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      sprintf("%.1f", cnv_state_mat[i, j]),
      x,
      y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# 7. COMBINE VERTICALLY
# ============================================================

ht_list <- ht_rna %v% ht_cnv


# ============================================================
# 8. DRAW
# ============================================================

draw(
  ht_list,
  
  column_title =
    "Molecular Signatures of the Five Patient States",
  
  column_title_gp = gpar(
    fontsize = 17,
    fontface = "bold"
  ),
  
  # Clear visual separation between RNA and CNV
  ht_gap = unit(12, "mm"),
  
  merge_legends = FALSE
)
##
# ============================================================
# RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  rna_state_mat,
  
  name = "RNA\nΔ",
  col = rna_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "RNA EXPRESSION SIGNATURES",
  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),
  
  row_title = "RNA genes",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 9),
  
  heatmap_legend_param = list(
    title = "RNA\nΔ",
    legend_height = unit(22, "mm")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.1f", rna_state_mat[i, j]),
      x, y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# CNV HEATMAP
# ============================================================

ht_cnv <- Heatmap(
  cnv_state_mat,
  
  name = "CNV\nΔ",
  col = cnv_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "COPY NUMBER (CNV) SIGNATURES",
  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),
  
  row_title = "CNV genes",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 9),
  
  heatmap_legend_param = list(
    title = "CNV\nΔ",
    legend_height = unit(22, "mm")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.1f", cnv_state_mat[i, j]),
      x, y,
      gp = gpar(fontsize = 7)
    )
  }
)


# ============================================================
# COMBINE
# ============================================================

ht_list <- ht_rna %v% ht_cnv


# ============================================================
# DRAW
# ============================================================

draw(
  ht_list,
  
  column_title =
    "Molecular Signatures of the Five Patient States",
  
  column_title_gp = gpar(
    fontsize = 17,
    fontface = "bold"
  ),
  
  ht_gap = unit(14, "mm"),
  
  merge_legends = FALSE,
  
  # RNA legend upper, CNV legend lower
  heatmap_legend_side = "right"
)


# ============================================================
# DASHED SEPARATOR BETWEEN RNA AND CNV
# ============================================================

decorate_heatmap_body("RNA\nΔ", {
  
  # line slightly below RNA heatmap
  grid.lines(
    x = unit(c(0, 1), "npc"),
    y = unit(-7, "mm"),
    gp = gpar(
      lty = 2,
      lwd = 1,
      col = "grey40"
    )
  )
})
# ============================================================
# SAVE FIGURE
# ============================================================

plot_state_heatmap <- function() {
  
  draw(
    ht_list,
    column_title = "Molecular Signatures of the Five Patient States",
    column_title_gp = gpar(
      fontsize = 17,
      fontface = "bold"
    ),
    ht_gap = unit(14, "mm"),
    merge_legends = FALSE,
    heatmap_legend_side = "right"
  )
  
  # Dashed separator
  decorate_heatmap_body("RNA\nΔ", {
    grid.lines(
      x = unit(c(0, 1), "npc"),
      y = unit(-7, "mm"),
      gp = gpar(
        lty = 2,
        lwd = 1,
        col = "grey40"
      )
    )
  })
}


# PDF
pdf(
  "Molecular_Signatures_Patient_States.pdf",
  width = 14,
  height = 10
)

plot_state_heatmap()

dev.off()


# High-resolution PNG
png(
  "Molecular_Signatures_Patient_States.png",
  width = 4200,
  height = 3000,
  res = 300
)

plot_state_heatmap()

dev.off()
#pathway
rna_universe <- sub("^RNA_", "", rownames(weights_v2$RNA))

length(rna_universe)
head(rna_universe)
library(clusterProfiler)
library(org.Hs.eg.db)

go_results_corrected <- compareCluster(
  Gene ~ Factor + Direction,
  data = factor_gene_sets,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  universe = rna_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_corrected_df <- as.data.frame(go_results_corrected)

dim(go_corrected_df)
head(go_corrected_df)
library(dplyr)

go_top5_corrected <- go_corrected_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Factor, Direction) %>%
  arrange(p.adjust, desc(FoldEnrichment)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(
    Factor,
    Direction,
    Description,
    FoldEnrichment,
    p.adjust,
    geneID
  )

go_top5_corrected
go_top5_corrected <- go_corrected_df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(Factor, Direction) %>%
  dplyr::arrange(p.adjust, dplyr::desc(FoldEnrichment)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Factor,
    Direction,
    Description,
    FoldEnrichment,
    p.adjust,
    geneID
  )

go_top5_corrected
library(dplyr)
library(ggplot2)

# ============================================================
# 1. SAVE FULL CORRECTED GO RESULTS
# ============================================================

write.csv(
  go_corrected_df,
  "GO_enrichment_corrected_full.csv",
  row.names = FALSE
)

write.csv(
  go_top5_corrected,
  "GO_enrichment_corrected_top5.csv",
  row.names = FALSE
)


# ============================================================
# 2. PREPARE TOP GO TERMS FOR VISUALIZATION
# ============================================================

go_plot <- go_top5_corrected %>%
  mutate(
    Factor = factor(
      Factor,
      levels = paste0("Factor", 1:10)
    ),
    
    Direction = factor(
      Direction,
      levels = c("Positive", "Negative")
    ),
    
    minus_log10_FDR = -log10(p.adjust),
    
    Label = paste(Factor, Direction, sep = " | ")
  )


# ============================================================
# 3. GO DOTPLOT
# ============================================================

p_go <- ggplot(
  go_plot,
  aes(
    x = FoldEnrichment,
    y = reorder(Description, FoldEnrichment)
  )
) +
  
  geom_point(
    aes(
      size = minus_log10_FDR,
      fill = Direction
    ),
    shape = 21
  ) +
  
  facet_wrap(
    ~ Factor,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    title = "Biological Programs Associated with MOFA Factors",
    subtitle = "GO Biological Process enrichment using the MOFA RNA feature universe",
    x = "Fold enrichment",
    y = NULL,
    size = "-log10(FDR)",
    fill = "Weight direction"
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text.y = element_text(
      size = 8
    ),
    
    legend.position = "right"
  )

p_go
ggsave(
  "MOFA_Factor_GO_Enrichment.pdf",
  p_go,
  width = 14,
  height = 18
)

ggsave(
  "MOFA_Factor_GO_Enrichment.png",
  p_go,
  width = 14,
  height = 18,
  dpi = 300
)
library(dplyr)
library(tidyr)

# ============================================================
# 1. STATE FACTOR PROFILES -> LONG FORMAT
# ============================================================

state_factor_all <- state_factor_profile %>%
  pivot_longer(
    cols = starts_with("Factor"),
    names_to = "Factor",
    values_to = "MeanScore"
  ) %>%
  mutate(
    Direction = ifelse(
      MeanScore >= 0,
      "Positive",
      "Negative"
    ),
    AbsScore = abs(MeanScore)
  )


# ============================================================
# 2. KEEP 4 STRONGEST FACTORS PER STATE
# ============================================================

state_factor_top <- state_factor_all %>%
  group_by(State) %>%
  arrange(desc(AbsScore), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup()


# ============================================================
# 3. JOIN WITH CORRECTED GO RESULTS
# ============================================================

state_go <- state_factor_top %>%
  left_join(
    go_top5_corrected,
    by = c(
      "Factor" = "Factor",
      "Direction" = "Direction"
    )
  )


# ============================================================
# 4. KEEP INFORMATIVE COLUMNS
# ============================================================

state_go_table <- state_go %>%
  dplyr::select(
    State,
    Factor,
    MeanScore,
    Direction,
    Description,
    FoldEnrichment,
    p.adjust
  ) %>%
  arrange(
    State,
    desc(abs(MeanScore)),
    p.adjust
  )


# ============================================================
# 5. VIEW
# ============================================================

print(
  state_go_table,
  n = Inf
)


# ============================================================
# 6. SAVE
# ============================================================

write.csv(
  state_go_table,
  "Patient_State_Biological_Programs.csv",
  row.names = FALSE
)
##
library(dplyr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. ASSIGN GO TERMS TO BROADER BIOLOGICAL PROGRAMS
# ============================================================

state_go_programs <- state_go_table %>%
  filter(!is.na(Description)) %>%
  
  mutate(
    Program = case_when(
      
      str_detect(
        Description,
        regex(
          "extracellular matrix|extracellular structure|encapsulating structure|collagen|connective tissue|fiber organization",
          ignore_case = TRUE
        )
      ) ~ "ECM / Collagen",
      
      str_detect(
        Description,
        regex(
          "immune|defense|cytokine|cell killing",
          ignore_case = TRUE
        )
      ) ~ "Immune / Defense",
      
      str_detect(
        Description,
        regex(
          "granulocyte|neutrophil|leukocyte|chemotaxis|migration",
          ignore_case = TRUE
        )
      ) ~ "Myeloid / Chemotaxis",
      
      str_detect(
        Description,
        regex(
          "adhesion",
          ignore_case = TRUE
        )
      ) ~ "Cell Adhesion",
      
      str_detect(
        Description,
        regex(
          "vasculature|circulatory|angiogenesis",
          ignore_case = TRUE
        )
      ) ~ "Vasculature",
      
      str_detect(
        Description,
        regex(
          "development|morphogenesis|neurogenesis",
          ignore_case = TRUE
        )
      ) ~ "Developmental",
      
      str_detect(
        Description,
        regex(
          "signaling|receptor",
          ignore_case = TRUE
        )
      ) ~ "Signaling",
      
      TRUE ~ "Other"
    )
  )


# ============================================================
# 2. CALCULATE STATE-PROGRAM SCORE
#
# MeanScore:
# direction + magnitude of the MOFA factor in that state
#
# -log10(FDR):
# strength of pathway enrichment
# ============================================================

state_go_programs <- state_go_programs %>%
  mutate(
    GO_strength = -log10(p.adjust),
    ProgramScore = MeanScore * GO_strength
  )


# ============================================================
# 3. COLLAPSE REDUNDANT GO TERMS
# ============================================================

state_program_score <- state_go_programs %>%
  filter(Program != "Other") %>%
  group_by(State, Program) %>%
  summarise(
    Score = mean(ProgramScore, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 4. CREATE STATE × PROGRAM MATRIX
# ============================================================

state_program_matrix <- state_program_score %>%
  pivot_wider(
    names_from = State,
    values_from = Score,
    values_fill = 0
  )

program_names <- state_program_matrix$Program

state_program_matrix <- as.data.frame(state_program_matrix)

rownames(state_program_matrix) <- program_names

state_program_matrix$Program <- NULL

state_program_matrix <- as.matrix(state_program_matrix)

colnames(state_program_matrix) <- paste0(
  "State ",
  colnames(state_program_matrix)
)


# ============================================================
# 5. HEATMAP COLOR SCALE
# ============================================================

lim <- max(abs(state_program_matrix), na.rm = TRUE)

program_col <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 6. HEATMAP
# ============================================================

ht_program <- Heatmap(
  
  state_program_matrix,
  
  name = "Program\nscore",
  
  col = program_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title =
    "Biological Programs Across Probabilistic Patient States",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = "Biological program",
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 10
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      sprintf(
        "%.1f",
        state_program_matrix[i, j]
      ),
      x,
      y,
      gp = gpar(fontsize = 8)
    )
  }
)


# ============================================================
# 7. DRAW
# ============================================================

draw(ht_program)


# ============================================================
# 8. SAVE MATRIX
# ============================================================

write.csv(
  state_program_matrix,
  "Patient_State_Biological_Program_Matrix.csv"
)


# ============================================================
# 9. SAVE FIGURE
# ============================================================

pdf(
  "Patient_State_Biological_Programs.pdf",
  width = 9,
  height = 6
)

draw(ht_program)

dev.off()


png(
  "Patient_State_Biological_Programs.png",
  width = 2700,
  height = 1800,
  res = 300
)

draw(ht_program)

dev.off()
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. STATE + EXACT GO TERMS
# ============================================================

state_go_exact <- state_go_table %>%
  filter(
    !is.na(Description),
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  mutate(
    GO_score = MeanScore * (-log10(p.adjust))
  )


# ============================================================
# 2. SELECT MOST INFORMATIVE GO TERMS
# Top 5 strongest GO associations per state
# ============================================================

state_go_top <- state_go_exact %>%
  group_by(State) %>%
  arrange(desc(abs(GO_score)), .by_group = TRUE) %>%
  distinct(Description, .keep_all = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

selected_GO <- unique(state_go_top$Description)


# ============================================================
# 3. CALCULATE SCORE FOR SELECTED GO TERMS IN EACH STATE
# ============================================================

go_state_scores <- state_go_exact %>%
  filter(Description %in% selected_GO) %>%
  group_by(State, Description) %>%
  summarise(
    Score = mean(GO_score),
    .groups = "drop"
  )


# ============================================================
# 4. MATRIX
# ============================================================

go_state_matrix <- go_state_scores %>%
  pivot_wider(
    names_from = State,
    values_from = Score,
    values_fill = 0
  )

go_names <- go_state_matrix$Description

go_state_matrix <- as.data.frame(go_state_matrix)

rownames(go_state_matrix) <- go_names
go_state_matrix$Description <- NULL

go_state_matrix <- as.matrix(go_state_matrix)

# Force all five states to appear
for(s in 1:5){
  if(!as.character(s) %in% colnames(go_state_matrix)){
    go_state_matrix <- cbind(
      go_state_matrix,
      setNames(
        data.frame(rep(0, nrow(go_state_matrix))),
        as.character(s)
      )
    )
  }
}

go_state_matrix <- go_state_matrix[, as.character(1:5)]

colnames(go_state_matrix) <- paste0("State ", 1:5)


# ============================================================
# 5. COLOR
# ============================================================

lim <- max(abs(go_state_matrix), na.rm = TRUE)

go_col <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 6. HEATMAP — EXACT GO NAMES
# ============================================================

ht_go_state <- Heatmap(
  
  go_state_matrix,
  
  name = "GO\nscore",
  col = go_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title =
    "GO Biological Processes Across Patient States",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = "GO Biological Process",
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_max_width = unit(10, "cm"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    if(go_state_matrix[i,j] != 0){
      grid.text(
        sprintf("%.1f", go_state_matrix[i,j]),
        x, y,
        gp = gpar(fontsize = 7)
      )
    }
  }
)


# ============================================================
# 7. DRAW
# ============================================================

draw(ht_go_state)


# ============================================================
# 8. SAVE
# ============================================================

write.csv(
  go_state_matrix,
  "Patient_State_Exact_GO_Matrix.csv"
)

pdf(
  "Patient_State_Exact_GO_Heatmap.pdf",
  width = 13,
  height = 9
)

draw(ht_go_state)

dev.off()
library(dplyr)
library(tidyr)

# ============================================================
# 1. ALL 10 FACTORS FOR EVERY STATE
# ============================================================

state_factor_all <- state_factor_profile %>%
  pivot_longer(
    cols = starts_with("Factor"),
    names_to = "Factor",
    values_to = "MeanScore"
  ) %>%
  mutate(
    Direction = ifelse(
      MeanScore >= 0,
      "Positive",
      "Negative"
    )
  )


# ============================================================
# 2. JOIN ALL FACTORS WITH CORRECTED GO RESULTS
# ============================================================

state_go_all <- state_factor_all %>%
  left_join(
    go_corrected_df %>%
      filter(p.adjust < 0.05),
    by = c(
      "Factor" = "Factor",
      "Direction" = "Direction"
    )
  )


# ============================================================
# 3. CALCULATE STATE-GO SCORE
# ============================================================

state_go_all <- state_go_all %>%
  filter(
    !is.na(Description),
    !is.na(p.adjust)
  ) %>%
  mutate(
    GO_strength = -log10(p.adjust),
    GO_score = MeanScore * GO_strength
  )


# ============================================================
# 4. TOP 10 EXACT GO TERMS PER STATE
# ============================================================

state_go_top10 <- state_go_all %>%
  group_by(State) %>%
  arrange(desc(abs(GO_score)), .by_group = TRUE) %>%
  distinct(Description, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()


# ============================================================
# 5. VIEW RESULTS
# ============================================================

state_go_top10 %>%
  dplyr::select(
    State,
    Factor,
    Direction,
    MeanScore,
    Description,
    p.adjust,
    GO_score
  ) %>%
  print(n = 50)


# ============================================================
# 6. SAVE
# ============================================================

write.csv(
  state_go_top10,
  "Patient_States_Top10_GO_AllFactors.csv",
  row.names = FALSE
)
# Tüm anlamlı GO'ları kullan
all_go_terms <- unique(state_go_all$Description)

go_matrix_all <- state_go_all %>%
  dplyr::group_by(State, Description) %>%
  dplyr::summarise(
    Score = mean(GO_score),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = Score,
    values_fill = 0
  )

go_names <- go_matrix_all$Description

go_matrix_all <- as.data.frame(go_matrix_all)
rownames(go_matrix_all) <- go_names
go_matrix_all$Description <- NULL
go_matrix_all <- as.matrix(go_matrix_all)

# State sırası
go_matrix_all <- go_matrix_all[, as.character(1:5)]
colnames(go_matrix_all) <- paste0("State ", 1:5)

dim(go_matrix_all)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# ALL 312 SIGNIFICANT GO TERMS
# ============================================================

lim <- max(abs(go_matrix_all), na.rm = TRUE)

go_col_all <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B4CC0", "white", "#B40426")
)

ht_go_all <- Heatmap(
  go_matrix_all,
  
  name = "GO\nscore",
  col = go_col_all,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "All Significant GO Biological Processes Across Patient States",
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = "GO Biological Process",
  
  row_names_gp = gpar(fontsize = 5),
  row_names_max_width = unit(12, "cm"),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  show_row_dend = TRUE,
  border = TRUE
)

# ============================================================
# SHOW
# ============================================================

draw(ht_go_all)

# ============================================================
# SAVE LARGE PDF
# ============================================================

pdf(
  "Patient_States_ALL_Significant_GO.pdf",
  width = 15,
  height = 35
)

draw(ht_go_all)

dev.off()

# Save matrix
write.csv(
  go_matrix_all,
  "Patient_States_ALL_Significant_GO_Matrix.csv"
)
if (!requireNamespace("rrvgo", quietly = TRUE)) {
  BiocManager::install("rrvgo")
}

library(rrvgo)
library(org.Hs.eg.db)
# ============================================================
# 1. INSTALL / LOAD PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("rrvgo", quietly = TRUE))
  BiocManager::install("rrvgo")

library(rrvgo)
library(org.Hs.eg.db)
library(dplyr)


# ============================================================
# 2. GET GO IDs USED IN OUR SIGNIFICANT RESULTS
# ============================================================

# GO descriptions -> corresponding GO IDs
go_lookup <- go_corrected_df %>%
  dplyr::filter(
    p.adjust < 0.05,
    Description %in% rownames(go_matrix_all)
  ) %>%
  dplyr::select(
    ID,
    Description,
    p.adjust
  ) %>%
  dplyr::distinct(ID, .keep_all = TRUE)

length(unique(go_lookup$ID))
head(go_lookup)


# ============================================================
# 3. SCORE EACH GO TERM
# Higher score = stronger enrichment
# ============================================================

go_scores <- -log10(go_lookup$p.adjust)

names(go_scores) <- go_lookup$ID


# ============================================================
# 4. CALCULATE GO SEMANTIC SIMILARITY
# ============================================================

simMatrix <- calculateSimMatrix(
  unique(go_lookup$ID),
  orgdb = "org.Hs.eg.db",
  ont = "BP",
  method = "Rel"
)


# ============================================================
# 5. AUTOMATICALLY REDUCE / CLUSTER GO TERMS
# ============================================================

reducedTerms <- reduceSimMatrix(
  simMatrix,
  scores = go_scores,
  threshold = 0.7,
  orgdb = "org.Hs.eg.db"
)


# ============================================================
# 6. CHECK RESULTS
# ============================================================

dim(reducedTerms)

head(
  reducedTerms[
    order(reducedTerms$parentTerm),
  ],
  20
)

length(unique(reducedTerms$parentTerm))


# ============================================================
# 7. SAVE
# ============================================================

write.csv(
  reducedTerms,
  "GO_semantic_clusters_rrvgo.csv",
  row.names = FALSE
)

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. MAP GO IDs -> rrvgo SEMANTIC CLUSTERS
# ============================================================

go_cluster_map <- reducedTerms %>%
  dplyr::select(
    ID = go,
    SemanticCluster = parentTerm,
    Representative = parent
  )


# ============================================================
# 2. ADD GO IDs TO STATE RESULTS
# ============================================================

state_go_clustered <- state_go_all %>%
  dplyr::left_join(
    go_corrected_df %>%
      dplyr::select(ID, Description) %>%
      dplyr::distinct(),
    by = "Description"
  ) %>%
  dplyr::left_join(
    go_cluster_map,
    by = "ID"
  ) %>%
  dplyr::filter(!is.na(SemanticCluster))


# ============================================================
# 3. COLLAPSE GO TERMS WITHIN EACH SEMANTIC CLUSTER
# ============================================================

state_cluster_score <- state_go_clustered %>%
  dplyr::group_by(
    State,
    SemanticCluster,
    Representative
  ) %>%
  dplyr::summarise(
    Score = mean(GO_score, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 4. CREATE 35 × 5 MATRIX
# ============================================================

cluster_matrix <- state_cluster_score %>%
  dplyr::select(
    State,
    Representative,
    Score
  ) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = Score,
    values_fill = 0
  )

cluster_names <- cluster_matrix$Representative

cluster_matrix <- as.data.frame(cluster_matrix)

rownames(cluster_matrix) <- cluster_names
cluster_matrix$Representative <- NULL

cluster_matrix <- as.matrix(cluster_matrix)

# State order
cluster_matrix <- cluster_matrix[, as.character(1:5)]

colnames(cluster_matrix) <- paste0("State ", 1:5)

dim(cluster_matrix)


# ============================================================
# 5. HEATMAP
# ============================================================

lim <- max(abs(cluster_matrix), na.rm = TRUE)

cluster_col <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B4CC0", "white", "#B40426")
)

ht_cluster <- Heatmap(
  cluster_matrix,
  
  name = "GO cluster\nscore",
  col = cluster_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title =
    "Semantic GO Programs Across Patient States",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = "GO semantic cluster",
  
  row_names_gp = gpar(fontsize = 8),
  
  row_names_max_width = unit(10, "cm"),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  border = TRUE
)

draw(ht_cluster)


# ============================================================
# 6. SAVE
# ============================================================

write.csv(
  cluster_matrix,
  "Patient_State_GO_Semantic_Clusters.csv"
)

pdf(
  "Patient_State_GO_Semantic_Clusters.pdf",
  width = 13,
  height = 12
)

draw(ht_cluster)

dev.off()
library(dplyr)
library(tidyr)

# ============================================================
# 1. CHECK rrvgo COLUMN NAMES
# ============================================================

colnames(reducedTerms)
# GO ID -> semantic cluster mapping

go_cluster_map <- reducedTerms %>%
  dplyr::transmute(
    ID = go,
    SemanticCluster = cluster,
    Representative = parentTerm
  )

# State GO sonuçlarıyla birleştir
state_go_clustered <- state_go_all %>%
  dplyr::left_join(
    go_cluster_map,
    by = "ID"
  ) %>%
  dplyr::filter(!is.na(SemanticCluster))

# Kontrol
dim(state_go_clustered)

length(unique(state_go_clustered$SemanticCluster))

state_go_clustered %>%
  dplyr::select(
    State,
    Factor,
    Direction,
    ID,
    Description,
    SemanticCluster,
    Representative,
    GO_score
  ) %>%
  head(20)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. COLLAPSE GO TERMS INTO SEMANTIC CLUSTERS
# ============================================================

state_semantic_score <- state_go_clustered %>%
  dplyr::group_by(
    State,
    SemanticCluster,
    Representative
  ) %>%
  dplyr::summarise(
    Score = mean(GO_score, na.rm = TRUE),
    N_GO = dplyr::n_distinct(ID),
    .groups = "drop"
  )


# ============================================================
# 2. CREATE CLUSTER × STATE MATRIX
# ============================================================

semantic_matrix_df <- state_semantic_score %>%
  dplyr::select(
    State,
    SemanticCluster,
    Representative,
    Score
  ) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = Score,
    values_fill = 0
  )

# Representative GO names
semantic_names <- semantic_matrix_df$Representative

semantic_matrix <- semantic_matrix_df %>%
  dplyr::select(-SemanticCluster, -Representative) %>%
  as.data.frame()

rownames(semantic_matrix) <- make.unique(semantic_names)

semantic_matrix <- as.matrix(semantic_matrix)

# State order
semantic_matrix <- semantic_matrix[, as.character(1:5)]

colnames(semantic_matrix) <- paste0("State ", 1:5)

dim(semantic_matrix)


# ============================================================
# 3. COLOR SCALE
# ============================================================

lim <- max(abs(semantic_matrix), na.rm = TRUE)

semantic_col <- colorRamp2(
  c(-lim, 0, lim),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 4. HEATMAP
# ============================================================

ht_semantic <- Heatmap(
  semantic_matrix,
  
  name = "Semantic\nGO score",
  col = semantic_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title =
    "Semantic GO Programs Across Patient States",
  
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_title = "Representative GO Biological Process",
  
  row_names_gp = gpar(fontsize = 8),
  
  row_names_max_width = unit(11, "cm"),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  border = TRUE
)

draw(ht_semantic)


# ============================================================
# 5. SAVE
# ============================================================

write.csv(
  semantic_matrix,
  "Patient_State_Semantic_GO_Matrix.csv"
)

write.csv(
  state_semantic_score,
  "Patient_State_Semantic_GO_Clusters.csv",
  row.names = FALSE
)

pdf(
  "Patient_State_Semantic_GO_Heatmap.pdf",
  width = 13,
  height = 12
)

draw(ht_semantic)

dev.off()
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. TOP RNA FEATURES ACROSS STATES
# ============================================================

top_rna_state <- lapply(1:5, function(s) {
  
  rna_state_signatures[[s]] %>%
    slice_max(
      order_by = abs(Difference),
      n = 5
    ) %>%
    pull(Feature)
  
}) %>%
  unlist() %>%
  unique()


# RNA state × feature matrix
rna_summary <- sapply(1:5, function(s) {
  
  x <- rna_state_signatures[[s]]
  
  x$Difference[
    match(top_rna_state, x$Feature)
  ]
  
})

rownames(rna_summary) <- top_rna_state
colnames(rna_summary) <- paste0("State ", 1:5)


# ============================================================
# 2. TOP CNV FEATURES ACROSS STATES
# ============================================================

top_cnv_state <- lapply(1:5, function(s) {
  
  cnv_state_signatures[[s]] %>%
    slice_max(
      order_by = abs(Difference),
      n = 5
    ) %>%
    pull(Feature)
  
}) %>%
  unlist() %>%
  unique()


cnv_summary <- sapply(1:5, function(s) {
  
  x <- cnv_state_signatures[[s]]
  
  x$Difference[
    match(top_cnv_state, x$Feature)
  ]
  
})

rownames(cnv_summary) <- top_cnv_state
colnames(cnv_summary) <- paste0("State ", 1:5)


# ============================================================
# 3. SEMANTIC GO MATRIX
# already created previously
# ============================================================

go_summary <- semantic_matrix


# ============================================================
# 4. AGE PER STATE
# ============================================================

age_summary <- factor_clinical %>%
  filter(!is.na(State)) %>%
  group_by(State) %>%
  summarise(
    Median_Age = median(age_years, na.rm = TRUE),
    .groups = "drop"
  )

age_matrix <- matrix(
  age_summary$Median_Age,
  nrow = 1
)

rownames(age_matrix) <- "Median age"
colnames(age_matrix) <- paste0(
  "State ",
  age_summary$State
)


# ============================================================
# 5. BRCA MUTATION FREQUENCY
# ============================================================

brca_summary <- factor_clinical %>%
  filter(!is.na(State)) %>%
  group_by(State) %>%
  summarise(
    BRCA_mut_percent =
      mean(BRCA_any_mut, na.rm = TRUE) * 100,
    .groups = "drop"
  )

brca_matrix <- matrix(
  brca_summary$BRCA_mut_percent,
  nrow = 1
)

rownames(brca_matrix) <- "BRCA1/2 mut (%)"
colnames(brca_matrix) <- paste0(
  "State ",
  brca_summary$State
)


# ============================================================
# 6. CHECK
# ============================================================

dim(rna_summary)
dim(cnv_summary)
dim(go_summary)

age_matrix
brca_matrix
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. ROW-WISE Z-SCORE FUNCTION
# ============================================================

row_zscore <- function(x) {
  t(scale(t(x)))
}

rna_z <- row_zscore(rna_summary)
cnv_z <- row_zscore(cnv_summary)
go_z  <- row_zscore(go_summary)

# remove possible NA rows caused by zero variance
rna_z[is.na(rna_z)] <- 0
cnv_z[is.na(cnv_z)] <- 0
go_z[is.na(go_z)]   <- 0


# ============================================================
# 2. CLINICAL VALUES -> Z-SCORE FOR COLOR ONLY
# ============================================================

age_z  <- row_zscore(age_matrix)
brca_z <- row_zscore(brca_matrix)


# ============================================================
# 3. COMMON COLOR SCALE
# ============================================================

z_col <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "white", "#B40426")
)


# ============================================================
# 4. RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  rna_z,
  name = "Relative\nstate profile",
  col = z_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "MOLECULAR AND BIOLOGICAL CHARACTERIZATION OF PATIENT STATES",
  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),
  
  row_title = "RNA expression",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 8),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  )
)


# ============================================================
# 5. CNV HEATMAP
# ============================================================

ht_cnv <- Heatmap(
  cnv_z,
  name = "CNV Z-score",
  col = z_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  row_title = "Copy number",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 8),
  
  show_column_names = FALSE,
  
  show_heatmap_legend = FALSE
)


# ============================================================
# 6. SEMANTIC GO HEATMAP
# ============================================================

ht_go <- Heatmap(
  go_z,
  name = "GO Z-score",
  col = z_col,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  row_title = "Semantic GO programs",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(fontsize = 7),
  
  row_names_max_width = unit(10, "cm"),
  
  show_column_names = FALSE,
  
  show_heatmap_legend = FALSE
)


# ============================================================
# 7. AGE
# ============================================================

ht_age <- Heatmap(
  age_z,
  name = "Age",
  col = z_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title = "Clinical",
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  show_column_names = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      sprintf("%.1f", age_matrix[i, j]),
      x,
      y,
      gp = gpar(
        fontsize = 9,
        fontface = "bold"
      )
    )
  },
  
  show_heatmap_legend = FALSE
)


# ============================================================
# 8. BRCA
# ============================================================

ht_brca <- Heatmap(
  brca_z,
  name = "BRCA",
  col = z_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      sprintf("%.1f%%", brca_matrix[i, j]),
      x,
      y,
      gp = gpar(
        fontsize = 9,
        fontface = "bold"
      )
    )
  },
  
  show_heatmap_legend = FALSE
)


# ============================================================
# 9. COMBINE VERTICALLY
# ============================================================

ht_final <-
  ht_rna %v%
  ht_cnv %v%
  ht_go %v%
  ht_age %v%
  ht_brca


draw(
  ht_final,
  heatmap_legend_side = "right"
)


# ============================================================
# 10. SAVE
# ============================================================

pdf(
  "Patient_State_Integrated_Characterization.pdf",
  width = 14,
  height = 20
)

draw(
  ht_final,
  heatmap_legend_side = "right"
)

dev.off()
library(dplyr)

rna_state_stats <- lapply(1:5, function(s) {
  
  in_state <- patient_state == s
  
  res <- apply(rna_mofa_new, 1, function(x) {
    
    test <- wilcox.test(
      x[in_state],
      x[!in_state],
      exact = FALSE
    )
    
    c(
      Median_State  = median(x[in_state], na.rm = TRUE),
      Median_Others = median(x[!in_state], na.rm = TRUE),
      P_value       = test$p.value
    )
  })
  
  res <- as.data.frame(t(res))
  
  res$Feature <- rownames(res)
  
  res$Difference <-
    res$Median_State - res$Median_Others
  
  res$FDR <- p.adjust(
    res$P_value,
    method = "BH"
  )
  
  res$State <- s
  
  res %>%
    arrange(FDR, desc(abs(Difference)))
  
})

rna_state_stats_all <- bind_rows(rna_state_stats)

# Kaç anlamlı RNA feature var?
rna_state_stats_all %>%
  filter(FDR < 0.05) %>%
  count(State)

# Her state için en güçlü 10
rna_state_stats_all %>%
  filter(FDR < 0.05) %>%
  group_by(State) %>%
  slice_max(
    order_by = abs(Difference),
    n = 10
  ) %>%
  arrange(State, desc(abs(Difference))) %>%
  select(State, Feature, Difference, FDR) %>%
  print(n = 50)
##
rna_state_stats <- lapply(1:5, function(s) {
  
  in_state <- patient_state == s
  
  res <- lapply(rownames(rna_mofa_new), function(gene) {
    
    x <- as.numeric(rna_mofa_new[gene, ])
    
    test <- wilcox.test(
      x[in_state],
      x[!in_state],
      exact = FALSE
    )
    
    data.frame(
      Feature = gene,
      Median_State = median(x[in_state], na.rm = TRUE),
      Median_Others = median(x[!in_state], na.rm = TRUE),
      Difference =
        median(x[in_state], na.rm = TRUE) -
        median(x[!in_state], na.rm = TRUE),
      P_value = test$p.value
    )
  })
  
  res <- dplyr::bind_rows(res)
  
  res$FDR <- p.adjust(
    res$P_value,
    method = "BH"
  )
  
  res$State <- s
  
  res
})


rna_state_stats_all <- dplyr::bind_rows(rna_state_stats)


# kontrol
str(rna_state_stats_all$FDR)


# FDR < 0.05 feature sayısı
rna_state_stats_all %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::count(State)
top_rna_markers <- rna_state_stats_all %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::group_by(State) %>%
  dplyr::slice_max(
    order_by = abs(Difference),
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(
    State,
    dplyr::desc(abs(Difference))
  ) %>%
  dplyr::select(
    State,
    Feature,
    Median_State,
    Median_Others,
    Difference,
    FDR
  )

print(top_rna_markers, n = 50)
cnv_state_stats <- lapply(1:5, function(s) {
  
  in_state <- patient_state == s
  
  res <- lapply(rownames(cnv_mofa_new), function(gene) {
    
    x <- as.numeric(cnv_mofa_new[gene, ])
    
    test <- wilcox.test(
      x[in_state],
      x[!in_state],
      exact = FALSE
    )
    
    data.frame(
      Feature = gene,
      Median_State = median(x[in_state], na.rm = TRUE),
      Median_Others = median(x[!in_state], na.rm = TRUE),
      Difference =
        median(x[in_state], na.rm = TRUE) -
        median(x[!in_state], na.rm = TRUE),
      P_value = test$p.value
    )
  })
  
  res <- dplyr::bind_rows(res)
  
  res$FDR <- p.adjust(
    res$P_value,
    method = "BH"
  )
  
  res$State <- s
  
  res
})

cnv_state_stats_all <- dplyr::bind_rows(cnv_state_stats)


# Number of significant CNV features per state
cnv_state_stats_all %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::count(State)


# Top 10 significant CNV features per state
top_cnv_markers <- cnv_state_stats_all %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::group_by(State) %>%
  dplyr::slice_max(
    order_by = abs(Difference),
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(
    State,
    dplyr::desc(abs(Difference))
  ) %>%
  dplyr::select(
    State,
    Feature,
    Median_State,
    Median_Others,
    Difference,
    FDR
  )

print(top_cnv_markers, n = 50)
cnv_state_stats_all %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::count(State)
# CNV annotation kontrolü

cnv_annotation <- cnv_test %>%
  dplyr::select(
    gene_name,
    chromosome,
    start,
    end
  ) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

# Prefix'i kaldır
cnv_state_stats_annotated <- cnv_state_stats_all %>%
  dplyr::mutate(
    gene_name = sub("^CNV_", "", Feature)
  ) %>%
  dplyr::left_join(
    cnv_annotation,
    by = "gene_name"
  )

# Her state için anlamlı CNV'lerin kromozom dağılımı
cnv_chr_summary <- cnv_state_stats_annotated %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::count(
    State,
    chromosome,
    sort = TRUE
  )

cnv_chr_summary %>%
  dplyr::group_by(State) %>%
  dplyr::slice_max(
    order_by = n,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(State, dplyr::desc(n)) %>%
  print(n = 25)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# STATE × CHROMOSOME CNV MATRIX
# ============================================================

cnv_chr_matrix <- cnv_chr_summary %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(cnv_chr_matrix) <- cnv_chr_matrix$chromosome
cnv_chr_matrix$chromosome <- NULL

cnv_chr_matrix <- as.matrix(cnv_chr_matrix)

# chromosome order
chr_order <- paste0("chr", 1:22)

cnv_chr_matrix <- cnv_chr_matrix[
  intersect(chr_order, rownames(cnv_chr_matrix)),
  ,
  drop = FALSE
]

colnames(cnv_chr_matrix) <- paste0("State ", 1:5)


# ============================================================
# COLOR
# ============================================================

max_n <- max(cnv_chr_matrix)

chr_col <- colorRamp2(
  c(0, max_n/2, max_n),
  c("white", "#FDB863", "#B2182B")
)


# ============================================================
# HEATMAP
# ============================================================

ht_chr <- Heatmap(
  cnv_chr_matrix,
  
  name = "Significant\nCNV genes",
  
  col = chr_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_title =
    "Chromosomal Distribution of State-Specific CNV Features",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  row_title = "Chromosome",
  
  row_names_gp = gpar(fontsize = 9),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      cnv_chr_matrix[i, j],
      x, y,
      gp = gpar(fontsize = 8)
    )
  },
  
  border = TRUE
)

draw(ht_chr)
# ============================================================
# 1. TOTAL TESTED CNV FEATURES PER CHROMOSOME
# ============================================================

cnv_tested_chr <- cnv_state_stats_annotated %>%
  dplyr::select(Feature, chromosome) %>%
  dplyr::distinct() %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22)
  ) %>%
  dplyr::count(
    chromosome,
    name = "Total_tested"
  )


# ============================================================
# 2. SIGNIFICANT / TOTAL TESTED × 100
# ============================================================

cnv_chr_percent <- cnv_chr_summary %>%
  dplyr::left_join(
    cnv_tested_chr,
    by = "chromosome"
  ) %>%
  dplyr::mutate(
    Percent_significant =
      100 * n / Total_tested
  )


# ============================================================
# 3. MATRIX
# ============================================================

cnv_chr_percent_matrix <- cnv_chr_percent %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22)
  ) %>%
  dplyr::select(
    State,
    chromosome,
    Percent_significant
  ) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = Percent_significant,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(cnv_chr_percent_matrix) <-
  cnv_chr_percent_matrix$chromosome

cnv_chr_percent_matrix$chromosome <- NULL

cnv_chr_percent_matrix <-
  as.matrix(cnv_chr_percent_matrix)

# chromosome order
chr_order <- paste0("chr", 1:22)

cnv_chr_percent_matrix <-
  cnv_chr_percent_matrix[
    intersect(
      chr_order,
      rownames(cnv_chr_percent_matrix)
    ),
    ,
    drop = FALSE
  ]

colnames(cnv_chr_percent_matrix) <-
  paste0("State ", 1:5)


# ============================================================
# 4. HEATMAP
# ============================================================

max_pct <- max(
  cnv_chr_percent_matrix,
  na.rm = TRUE
)

pct_col <- circlize::colorRamp2(
  c(0, max_pct / 2, max_pct),
  c("white", "#FDB863", "#B2182B")
)

ht_chr_normalized <- ComplexHeatmap::Heatmap(
  
  cnv_chr_percent_matrix,
  
  name = "% significant\nCNV features",
  
  col = pct_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_title =
    "Chromosome-Level Enrichment of State-Specific CNV Features",
  
  column_title_gp = grid::gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  row_title = "Chromosome",
  
  cell_fun = function(
    j, i, x, y,
    width, height, fill
  ) {
    
    grid::grid.text(
      sprintf(
        "%.1f%%",
        cnv_chr_percent_matrix[i, j]
      ),
      x, y,
      gp = grid::gpar(fontsize = 8)
    )
  },
  
  border = TRUE
)

ComplexHeatmap::draw(ht_chr_normalized)
cnv_chr_percent %>%
  dplyr::filter(chromosome %in% paste0("chr", 1:22)) %>%
  dplyr::select(
    State,
    chromosome,
    n,
    Total_tested,
    Percent_significant
  ) %>%
  dplyr::arrange(
    State,
    factor(chromosome, levels = paste0("chr", 1:22))
  ) %>%
  print(n = Inf)
cnv_chr_check <- cnv_chr_percent %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22)
  ) %>%
  dplyr::select(
    State,
    chromosome,
    n,
    Total_tested,
    Percent_significant
  ) %>%
  dplyr::arrange(
    State,
    factor(
      chromosome,
      levels = paste0("chr", 1:22)
    )
  )

cnv_chr_check
# ============================================================
# PREPARE GENOMIC CNV EFFECT DATA
# ============================================================

cnv_genomic <- cnv_state_stats_annotated %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22),
    !is.na(start),
    !is.na(end)
  ) %>%
  dplyr::mutate(
    chr_num = as.integer(
      sub("chr", "", chromosome)
    ),
    position = (start + end) / 2,
    Significant = FDR < 0.05
  ) %>%
  dplyr::arrange(
    State,
    chr_num,
    position
  )

# kontrol
cnv_genomic %>%
  dplyr::select(
    State,
    Feature,
    chromosome,
    position,
    Difference,
    FDR,
    Significant
  ) %>%
  head(20)
library(ggplot2)
library(dplyr)

# Significant CNV features
cnv_plot <- cnv_genomic %>%
  filter(Significant)

# Chromosome lengths from our observed coordinates
chr_info <- cnv_genomic %>%
  group_by(chr_num) %>%
  summarise(chr_length = max(end, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr_num) %>%
  mutate(
    offset = cumsum(lag(chr_length, default = 0)),
    center = offset + chr_length / 2
  )

# Add cumulative genomic position
cnv_plot <- cnv_plot %>%
  left_join(
    chr_info %>% select(chr_num, offset),
    by = "chr_num"
  ) %>%
  mutate(
    genome_position = position + offset
  )

# Chromosome boundaries
chr_boundaries <- chr_info$offset[-1]

# ============================================================
# PLOT
# ============================================================

p_cnv_landscape <- ggplot(
  cnv_plot,
  aes(
    x = genome_position,
    y = Difference
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  
  geom_vline(
    xintercept = chr_boundaries,
    linewidth = 0.2,
    linetype = "dashed"
  ) +
  
  geom_point(
    aes(color = Difference),
    size = 0.7,
    alpha = 0.7
  ) +
  
  facet_grid(
    State ~ .,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = chr_info$center,
    labels = paste0("chr", chr_info$chr_num),
    expand = c(0.01, 0.01)
  ) +
  
  scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  
  labs(
    title = "State-Specific CNV Landscape Across the Genome",
    x = "Genomic position",
    y = "CNV difference\n(State vs Others)",
    color = "CNV\ndifference"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    
    strip.text.y = element_text(
      face = "bold",
      size = 10
    ),
    
    panel.spacing.y = unit(0.3, "cm")
  )

p_cnv_landscape
# ============================================================
# STATE-SPECIFIC CNV GENOMIC LANDSCAPE
# ============================================================

library(dplyr)
library(ggplot2)
library(grid)


# ============================================================
# 1. PREPARE GENOMIC CNV DATA
# ============================================================

cnv_genomic <- cnv_state_stats_annotated %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22),
    !is.na(start),
    !is.na(end)
  ) %>%
  dplyr::mutate(
    chr_num = as.integer(
      sub("chr", "", chromosome)
    ),
    position = (start + end) / 2,
    Significant = FDR < 0.05
  ) %>%
  dplyr::arrange(
    State,
    chr_num,
    position
  )


# ============================================================
# 2. KEEP SIGNIFICANT CNV FEATURES
# ============================================================

cnv_plot <- cnv_genomic %>%
  dplyr::filter(Significant)


# ============================================================
# 3. CALCULATE CHROMOSOME OFFSETS
# ============================================================

chr_info <- cnv_genomic %>%
  dplyr::group_by(chr_num) %>%
  dplyr::summarise(
    chr_length = max(end, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(chr_num) %>%
  dplyr::mutate(
    offset = cumsum(
      dplyr::lag(
        chr_length,
        default = 0
      )
    ),
    center = offset + chr_length / 2
  )


# ============================================================
# 4. ADD CUMULATIVE GENOMIC POSITION
# ============================================================

cnv_plot <- cnv_plot %>%
  dplyr::left_join(
    chr_info %>%
      dplyr::select(
        chr_num,
        offset
      ),
    by = "chr_num"
  ) %>%
  dplyr::mutate(
    genome_position = position + offset
  )


# ============================================================
# 5. CHROMOSOME BOUNDARIES
# ============================================================

chr_boundaries <- chr_info$offset[-1]


# ============================================================
# 6. CHECK
# ============================================================

head(
  cnv_plot[
    ,
    c(
      "State",
      "Feature",
      "chromosome",
      "Difference",
      "FDR",
      "genome_position"
    )
  ]
)


# ============================================================
# 7. GENOMIC CNV LANDSCAPE
# ============================================================

p_cnv_landscape <- ggplot(
  cnv_plot,
  aes(
    x = genome_position,
    y = Difference
  )
) +
  
  # zero line
  geom_hline(
    yintercept = 0,
    linewidth = 0.3
  ) +
  
  # chromosome boundaries
  geom_vline(
    xintercept = chr_boundaries,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  
  # CNV features
  geom_point(
    aes(color = Difference),
    size = 0.8,
    alpha = 0.75
  ) +
  
  # one panel per state
  facet_grid(
    State ~ .,
    scales = "free_y",
    labeller = labeller(
      State = function(x) paste0("State ", x)
    )
  ) +
  
  # chromosome labels
  scale_x_continuous(
    breaks = chr_info$center,
    labels = paste0(
      "chr",
      chr_info$chr_num
    ),
    expand = c(0.01, 0.01)
  ) +
  
  # relative CNV difference
  scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  
  labs(
    title = "State-Specific CNV Landscape Across the Genome",
    subtitle = "Significant CNV features (FDR < 0.05)",
    x = "Genomic position",
    y = "Relative CNV difference\n(State vs Others)",
    color = "CNV\ndifference"
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    strip.text.y = element_text(
      face = "bold",
      size = 10
    ),
    
    panel.spacing.y = unit(
      0.25,
      "cm"
    ),
    
    legend.position = "right"
  )


# ============================================================
# 8. SHOW
# ============================================================

p_cnv_landscape


# ============================================================
# 9. SAVE
# ============================================================

ggsave(
  "Patient_State_CNV_Genomic_Landscape.pdf",
  p_cnv_landscape,
  width = 14,
  height = 10
)

ggsave(
  "Patient_State_CNV_Genomic_Landscape.png",
  p_cnv_landscape,
  width = 14,
  height = 10,
  dpi = 300
)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. PREPARE SEGMENTS
# ============================================================

cnv_track <- cnv_segments %>%
  dplyr::filter(
    n_genes >= 2,
    chromosome %in% paste0("chr", 1:22)
  ) %>%
  dplyr::left_join(
    chr_info %>%
      dplyr::select(chr_num, offset),
    by = "chr_num"
  ) %>%
  dplyr::mutate(
    genome_start = start + offset,
    genome_end   = end + offset,
    State_label  = paste0("State ", State)
  )

# State order
cnv_track$State_label <- factor(
  cnv_track$State_label,
  levels = paste0("State ", 5:1)
)


# ============================================================
# 2. GENOME TRACK PLOT
# ============================================================

p_cnv_track <- ggplot(cnv_track) +
  
  # CNV genomic segments
  geom_segment(
    aes(
      x = genome_start,
      xend = genome_end,
      y = State_label,
      yend = State_label,
      color = Difference
    ),
    linewidth = 8,
    lineend = "butt"
  ) +
  
  # chromosome boundaries
  geom_vline(
    xintercept = chr_boundaries,
    linetype = "dashed",
    linewidth = 0.25
  ) +
  
  # chromosome labels
  scale_x_continuous(
    breaks = chr_info$center,
    labels = paste0("chr", chr_info$chr_num),
    expand = c(0.005, 0.005)
  ) +
  
  # CNV difference
  scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = "Relative CNV\ndifference"
  ) +
  
  labs(
    title = "State-Specific CNV Genomic Landscape",
    subtitle = "Significant CNV regions (FDR < 0.05)",
    x = "Genomic position",
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    
    axis.text.y = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.ticks.y = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.position = "right"
  )


p_cnv_track


# ============================================================
# SAVE
# ============================================================

ggsave(
  "Patient_State_CNV_Genome_Tracks.pdf",
  p_cnv_track,
  width = 15,
  height = 5
)

ggsave(
  "Patient_State_CNV_Genome_Tracks.png",
  p_cnv_track,
  width = 15,
  height = 5,
  dpi = 300
)
library(mclust)

factor_numbers <- c(8, 10, 12, 15)

stability_results <- lapply(factor_numbers, function(k) {
  
  # Select MOFA factors
  x <- factor_scores[, paste0("Factor", 1:k)]
  
  # Standardize
  x_scaled <- scale(x)
  
  # Fit GMM
  fit <- Mclust(
    x_scaled,
    G = 2:8
  )
  
  # Compare with original 10-factor state solution
  ari <- adjustedRandIndex(
    patient_state,
    fit$classification
  )
  
  data.frame(
    N_Factors = k,
    Selected_States = fit$G,
    Model = fit$modelName,
    BIC = fit$bic,
    ARI_vs_10Factor = ari
  )
})

factor_stability <- dplyr::bind_rows(stability_results)

factor_stability
get_variance_explained(mofa_v2)$r2_per_factor
library(mclust)
library(dplyr)

factor_numbers <- c(8, 10, 12, 15)

unscaled_results <- lapply(factor_numbers, function(k) {
  
  # IMPORTANT: no scale()
  x <- factor_scores[, paste0("Factor", 1:k)]
  
  fit <- Mclust(
    x,
    G = 2:8
  )
  
  data.frame(
    N_Factors = k,
    Selected_States = fit$G,
    Model = fit$modelName,
    BIC = fit$bic
  )
})

unscaled_stability <- dplyr::bind_rows(unscaled_results)

unscaled_stability
library(mclust)

# Unscaled GMM models
fit_8 <- Mclust(
  factor_scores[, paste0("Factor", 1:8)],
  G = 2:8
)

fit_10 <- Mclust(
  factor_scores[, paste0("Factor", 1:10)],
  G = 2:8
)

fit_12 <- Mclust(
  factor_scores[, paste0("Factor", 1:12)],
  G = 2:8
)

# Pairwise ARI
ari_8_10 <- adjustedRandIndex(
  fit_8$classification,
  fit_10$classification
)

ari_8_12 <- adjustedRandIndex(
  fit_8$classification,
  fit_12$classification
)

ari_10_12 <- adjustedRandIndex(
  fit_10$classification,
  fit_12$classification
)

# Results
ari_results <- data.frame(
  Comparison = c("8 vs 10", "8 vs 12", "10 vs 12"),
  ARI = c(ari_8_10, ari_8_12, ari_10_12)
)

ari_results
library(mclust)
library(dplyr)

set.seed(123)

# ============================================================
# ORIGINAL 10-FACTOR UNSCALED MODEL
# ============================================================

X <- factor_scores[, paste0("Factor", 1:10)]

original_fit <- Mclust(
  X,
  G = 5,
  modelNames = "VVE"
)

original_cluster <- original_fit$classification


# ============================================================
# BOOTSTRAP
# ============================================================

n_boot <- 100

bootstrap_results <- lapply(1:n_boot, function(i) {
  
  # Sample patients with replacement
  boot_idx <- sample(
    1:nrow(X),
    size = nrow(X),
    replace = TRUE
  )
  
  # Unique patients present in bootstrap sample
  unique_idx <- unique(boot_idx)
  
  X_boot <- X[boot_idx, , drop = FALSE]
  
  # Fit same 5-state model
  fit_boot <- tryCatch(
    Mclust(
      X_boot,
      G = 5,
      modelNames = "VVE",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_boot)) {
    return(NA_real_)
  }
  
  # Predict states for all original patients
  pred <- predict(
    fit_boot,
    newdata = X
  )$classification
  
  # Compare with original clustering
  adjustedRandIndex(
    original_cluster,
    pred
  )
})


# ============================================================
# RESULTS
# ============================================================

bootstrap_ari <- unlist(bootstrap_results)

summary(bootstrap_ari)

mean(bootstrap_ari, na.rm = TRUE)
median(bootstrap_ari, na.rm = TRUE)
sd(bootstrap_ari, na.rm = TRUE)

quantile(
  bootstrap_ari,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)
library(mclust)

set.seed(123)

X <- as.matrix(
  factor_scores[, paste0("Factor", 1:10)]
)

# Reference model
original_fit <- Mclust(
  X,
  G = 5,
  modelNames = "VVE",
  verbose = FALSE
)

original_cluster <- original_fit$classification

# Bootstrap
n_boot <- 100
bootstrap_ari <- rep(NA_real_, n_boot)

for (i in 1:n_boot) {
  
  boot_idx <- sample(
    seq_len(nrow(X)),
    size = nrow(X),
    replace = TRUE
  )
  
  X_boot <- X[boot_idx, , drop = FALSE]
  
  fit_boot <- tryCatch(
    Mclust(
      X_boot,
      G = 5,
      modelNames = "VVE",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_boot)) {
    
    pred <- mclust::predict.Mclust(
      fit_boot,
      X
    )$classification
    
    bootstrap_ari[i] <- adjustedRandIndex(
      original_cluster,
      pred
    )
  }
}

# Results
summary(bootstrap_ari)

mean(bootstrap_ari, na.rm = TRUE)
median(bootstrap_ari, na.rm = TRUE)
sd(bootstrap_ari, na.rm = TRUE)

quantile(
  bootstrap_ari,
  c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)
library(mclust)
library(dplyr)

set.seed(123)

# ============================================================
# DATA
# ============================================================

X <- as.matrix(
  factor_scores[, paste0("Factor", 1:10)]
)

K_values <- 2:8
n_boot <- 100


# ============================================================
# TEST EACH K
# ============================================================

stability_results <- lapply(K_values, function(K) {
  
  cat("Testing K =", K, "\n")
  
  # ----------------------------------------------------------
  # Original model
  # ----------------------------------------------------------
  
  original_fit <- tryCatch(
    Mclust(
      X,
      G = K,
      modelNames = "VVE",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(original_fit)) {
    return(NULL)
  }
  
  original_cluster <- original_fit$classification
  
  # ----------------------------------------------------------
  # Bootstrap
  # ----------------------------------------------------------
  
  boot_ari <- rep(NA_real_, n_boot)
  
  for (i in 1:n_boot) {
    
    boot_idx <- sample(
      seq_len(nrow(X)),
      size = nrow(X),
      replace = TRUE
    )
    
    X_boot <- X[boot_idx, , drop = FALSE]
    
    fit_boot <- tryCatch(
      Mclust(
        X_boot,
        G = K,
        modelNames = "VVE",
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_boot)) {
      
      pred <- tryCatch(
        mclust::predict.Mclust(
          fit_boot,
          X
        )$classification,
        error = function(e) NULL
      )
      
      if (!is.null(pred)) {
        
        boot_ari[i] <- adjustedRandIndex(
          original_cluster,
          pred
        )
      }
    }
  }
  
  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------
  
  data.frame(
    K = K,
    
    BIC = as.numeric(original_fit$bic),
    
    Mean_ARI = mean(
      boot_ari,
      na.rm = TRUE
    ),
    
    Median_ARI = median(
      boot_ari,
      na.rm = TRUE
    ),
    
    SD_ARI = sd(
      boot_ari,
      na.rm = TRUE
    ),
    
    Q25_ARI = quantile(
      boot_ari,
      0.25,
      na.rm = TRUE
    ),
    
    Q75_ARI = quantile(
      boot_ari,
      0.75,
      na.rm = TRUE
    ),
    
    Successful_Bootstraps = sum(
      !is.na(boot_ari)
    )
  )
})


# ============================================================
# FINAL TABLE
# ============================================================

stability_table <- dplyr::bind_rows(
  stability_results
)

stability_table
library(mclust)
library(dplyr)

# ============================================================
# 1. FINAL ROBUST K=2 MODEL
# ============================================================

X <- as.matrix(
  factor_scores[, paste0("Factor", 1:10)]
)

final_state_fit <- Mclust(
  X,
  G = 2,
  modelNames = "VVE",
  verbose = FALSE
)

# State assignment
robust_state <- final_state_fit$classification

table(robust_state)


# ============================================================
# 2. POSTERIOR PROBABILITIES
# ============================================================

robust_prob <- final_state_fit$z

colnames(robust_prob) <- c(
  "P_State1",
  "P_State2"
)

# Maximum assignment probability
max_probability <- apply(
  robust_prob,
  1,
  max
)

summary(max_probability)


# ============================================================
# 3. CREATE PATIENT STATE TABLE
# ============================================================

robust_state_df <- data.frame(
  submitter_id = rownames(factor_scores),
  Robust_State = robust_state,
  P_State1 = robust_prob[, 1],
  P_State2 = robust_prob[, 2],
  Max_Probability = max_probability
)

head(robust_state_df)


# ============================================================
# 4. ADD CLINICAL DATA
# ============================================================

robust_state_clinical <- clinical_final %>%
  dplyr::left_join(
    robust_state_df,
    by = "submitter_id"
  )


# ============================================================
# 5. STATE SIZE
# ============================================================

table(
  robust_state_clinical$Robust_State
)


# ============================================================
# 6. AGE
# ============================================================

robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Median_Age = median(age_years, na.rm = TRUE),
    Mean_Age = mean(age_years, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 7. BRCA
# ============================================================

robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    BRCA_mut_N = sum(BRCA_any_mut, na.rm = TRUE),
    BRCA_mut_percent =
      100 * mean(BRCA_any_mut, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 8. FIGO
# ============================================================

table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$FIGO_main,
  useNA = "ifany"
)


# ============================================================
# 9. GRADE
# ============================================================

table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$tumor_grade,
  useNA = "ifany"
)
# FIGO
table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$figo_stage,
  useNA = "ifany"
)

# Grade
table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$tumor_grade,
  useNA = "ifany"
)

# Age
robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Median_Age = median(age_years, na.rm = TRUE),
    Mean_Age = mean(age_years, na.rm = TRUE),
    .groups = "drop"
  )

# BRCA
robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    BRCA_mut_N = sum(BRCA_any_mut, na.rm = TRUE),
    BRCA_mut_percent = 100 * mean(BRCA_any_mut, na.rm = TRUE),
    .groups = "drop"
  )
library(ggplot2)
library(dplyr)

brca_plot <- robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    BRCA_mut_percent = 100 * mean(BRCA_any_mut, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  brca_plot,
  aes(
    x = factor(Robust_State),
    y = BRCA_mut_percent
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(round(BRCA_mut_percent, 1), "%")),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = "Robust Molecular State",
    y = "BRCA1/2-mutated patients (%)",
    title = "Somatic BRCA Mutation Frequency Across Molecular States"
  ) +
  ylim(0, 10) +
  theme_classic(base_size = 13)
# FIGO distribution (% within each state)
figo_dist <- robust_state_clinical %>%
  dplyr::count(Robust_State, figo_stage) %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::mutate(
    Percent = 100 * n / sum(n)
  )

figo_dist


# Grade distribution (% within each state)
grade_dist <- robust_state_clinical %>%
  dplyr::count(Robust_State, tumor_grade) %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::mutate(
    Percent = 100 * n / sum(n)
  )

grade_dist
grade_test <- robust_state_clinical %>%
  dplyr::filter(tumor_grade %in% c("G2", "G3"))

grade_table <- table(
  grade_test$Robust_State,
  grade_test$tumor_grade
)

grade_table

fisher.test(grade_table)
library(dplyr)

# ============================================================
# 1. AGE — Wilcoxon
# ============================================================

age_test <- wilcox.test(
  age_years ~ Robust_State,
  data = robust_state_clinical
)

age_p <- age_test$p.value


# ============================================================
# 2. BRCA — Fisher
# ============================================================

brca_table <- table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$BRCA_any_mut
)

brca_test <- fisher.test(brca_table)
brca_p <- brca_test$p.value


# ============================================================
# 3. FIGO — I/II vs III/IV
# ============================================================

clinical_test <- robust_state_clinical %>%
  dplyr::mutate(
    FIGO_group = dplyr::case_when(
      grepl("^Stage I($|[^I])|^I[A-C]?$", figo_stage, ignore.case = TRUE) ~ "I",
      grepl("^Stage II($|[^I])|^II[A-C]?$", figo_stage, ignore.case = TRUE) ~ "II",
      grepl("^Stage III|^III[A-C]?$", figo_stage, ignore.case = TRUE) ~ "III",
      grepl("^Stage IV|^IV$", figo_stage, ignore.case = TRUE) ~ "IV",
      TRUE ~ NA_character_
    ),
    FIGO_binary = dplyr::case_when(
      FIGO_group %in% c("I", "II") ~ "I/II",
      FIGO_group %in% c("III", "IV") ~ "III/IV",
      TRUE ~ NA_character_
    )
  )

figo_table <- table(
  clinical_test$Robust_State,
  clinical_test$FIGO_binary
)

figo_test <- fisher.test(figo_table)
figo_p <- figo_test$p.value


# ============================================================
# 4. GRADE — G2 vs G3
# ============================================================

grade_test_data <- robust_state_clinical %>%
  dplyr::filter(tumor_grade %in% c("G2", "G3"))

grade_table <- table(
  grade_test_data$Robust_State,
  grade_test_data$tumor_grade
)

grade_test <- fisher.test(grade_table)
grade_p <- grade_test$p.value


# ============================================================
# 5. VITAL STATUS — Alive vs Dead
# ============================================================

vital_table <- table(
  robust_state_clinical$Robust_State,
  robust_state_clinical$vital_status
)

vital_test <- fisher.test(vital_table)
vital_p <- vital_test$p.value


# ============================================================
# 6. SUMMARY + BH FDR
# ============================================================

clinical_associations <- data.frame(
  Variable = c(
    "Age",
    "BRCA mutation",
    "FIGO stage (I/II vs III/IV)",
    "Grade (G2 vs G3)",
    "Vital status"
  ),
  P_value = c(
    age_p,
    brca_p,
    figo_p,
    grade_p,
    vital_p
  )
) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

clinical_associations
#two state analizleri
library(dplyr)

# State bilgisi RNA sample sırasına getir
state_vector <- robust_state[
  match(colnames(rna_mofa_new), rownames(factor_scores))
]

table(state_vector)


# ============================================================
# RNA: STATE 1 vs STATE 2
# ============================================================

rna_robust_stats <- lapply(
  seq_len(nrow(rna_mofa_new)),
  function(i) {

    x <- rna_mofa_new[i, ]

    group1 <- x[state_vector == 1]
    group2 <- x[state_vector == 2]

    test <- wilcox.test(
      group1,
      group2,
      exact = FALSE
    )

    data.frame(
      Feature = rownames(rna_mofa_new)[i],
      Median_State1 = median(group1, na.rm = TRUE),
      Median_State2 = median(group2, na.rm = TRUE),
      Difference = median(group2, na.rm = TRUE) -
                   median(group1, na.rm = TRUE),
      P_value = test$p.value
    )
  }
)

rna_robust_stats <- dplyr::bind_rows(rna_robust_stats) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)


# ============================================================
# RESULTS
# ============================================================

# Number significant
sum(rna_robust_stats$FDR < 0.05)

# Top significant genes
rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(20)
sum(rna_robust_stats$FDR < 0.05)
# State vector
state_vector_cnv <- robust_state[
  match(colnames(cnv_mofa_new), rownames(factor_scores))
]

table(state_vector_cnv)


# ============================================================
# CNV: STATE 1 vs STATE 2
# ============================================================

cnv_robust_stats <- lapply(
  seq_len(nrow(cnv_mofa_new)),
  function(i) {
    
    x <- cnv_mofa_new[i, ]
    
    group1 <- x[state_vector_cnv == 1]
    group2 <- x[state_vector_cnv == 2]
    
    test <- wilcox.test(
      group1,
      group2,
      exact = FALSE
    )
    
    data.frame(
      Feature = rownames(cnv_mofa_new)[i],
      Median_State1 = median(group1, na.rm = TRUE),
      Median_State2 = median(group2, na.rm = TRUE),
      Difference = median(group2, na.rm = TRUE) -
        median(group1, na.rm = TRUE),
      P_value = test$p.value
    )
  }
)

cnv_robust_stats <- dplyr::bind_rows(cnv_robust_stats) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)


# Number significant
sum(cnv_robust_stats$FDR < 0.05)


# Top 20 CNV features
cnv_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(20)
sum(cnv_robust_stats$FDR < 0.05)
# Gene annotation ekle
cnv_robust_annotated <- cnv_robust_stats %>%
  dplyr::mutate(
    gene_name = sub("^CNV_", "", Feature)
  ) %>%
  dplyr::left_join(
    cnv_annotation,
    by = "gene_name"
  )

# Sadece anlamlı CNV'ler
cnv_sig_robust <- cnv_robust_annotated %>%
  dplyr::filter(
    FDR < 0.05,
    chromosome %in% paste0("chr", 1:22)
  )

# Her kromozomda:
# kaç feature test edildi / kaçı anlamlı / yüzde kaç?
cnv_chr_robust <- cnv_robust_annotated %>%
  dplyr::filter(
    chromosome %in% paste0("chr", 1:22)
  ) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(
    Total_tested = dplyr::n(),
    Significant = sum(FDR < 0.05),
    Percent_significant = 100 * Significant / Total_tested,
    Median_difference = median(
      Difference[FDR < 0.05],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    chr_num = as.integer(sub("chr", "", chromosome))
  ) %>%
  dplyr::arrange(chr_num)

cnv_chr_robust
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# ============================================================
# 1. SIGNIFICANT RNA GENES
# Difference > 0 = State 2 high
# Difference < 0 = State 1 high
# ============================================================

rna_sig <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::mutate(
    Gene = sub("^RNA_", "", Feature),
    Direction = ifelse(
      Difference > 0,
      "State2_high",
      "State1_high"
    )
  )

table(rna_sig$Direction)


# ============================================================
# 2. CORRECT BACKGROUND
# All RNA features entering MOFA
# ============================================================

rna_background <- sub(
  "^RNA_",
  "",
  rownames(rna_mofa_new)
)


# ============================================================
# 3. GO ENRICHMENT — STATE 2 HIGH
# ============================================================

genes_state2 <- rna_sig$Gene[
  rna_sig$Direction == "State2_high"
]

go_state2 <- enrichGO(
  gene          = genes_state2,
  universe      = rna_background,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)


# ============================================================
# 4. GO ENRICHMENT — STATE 1 HIGH
# ============================================================

genes_state1 <- rna_sig$Gene[
  rna_sig$Direction == "State1_high"
]

go_state1 <- enrichGO(
  gene          = genes_state1,
  universe      = rna_background,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)


# ============================================================
# 5. RESULTS
# ============================================================

go_state2_df <- as.data.frame(go_state2)
go_state1_df <- as.data.frame(go_state1)

cat("State 2-high significant GO terms:",
    sum(go_state2_df$p.adjust < 0.05), "\n")

cat("State 1-high significant GO terms:",
    sum(go_state1_df$p.adjust < 0.05), "\n")


go_state2_df %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, GeneRatio, FoldEnrichment, p.adjust) %>%
  head(15)

go_state1_df %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, GeneRatio, FoldEnrichment, p.adjust) %>%
  head(15)
#GO sonucu alamadik
# Factor scores + robust state
factor_robust <- as.data.frame(factor_scores[, paste0("Factor", 1:10)])

factor_robust$Robust_State <- robust_state


# State mean factor scores
factor_state_means <- factor_robust %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    dplyr::across(
      starts_with("Factor"),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

factor_state_means


# Statistical comparison of each factor
factor_tests <- lapply(
  paste0("Factor", 1:10),
  function(f) {
    
    test <- wilcox.test(
      factor_robust[[f]] ~ factor_robust$Robust_State,
      exact = FALSE
    )
    
    data.frame(
      Factor = f,
      Mean_State1 = mean(
        factor_robust[[f]][factor_robust$Robust_State == 1]
      ),
      Mean_State2 = mean(
        factor_robust[[f]][factor_robust$Robust_State == 2]
      ),
      Difference =
        mean(factor_robust[[f]][factor_robust$Robust_State == 2]) -
        mean(factor_robust[[f]][factor_robust$Robust_State == 1]),
      P_value = test$p.value
    )
  }
)

factor_tests <- dplyr::bind_rows(factor_tests) %>%
  dplyr::mutate(
    FDR = p.adjust(P_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

factor_tests
library(dplyr)

# ============================================================
# 1. SAMPLE-LEVEL RAW CNV / PLOIDY PROXIES
# ============================================================

cnv_raw_test <- cnv_pc_auto[, rownames(factor_scores)]

cnv_level_df <- data.frame(
  submitter_id = colnames(cnv_raw_test),
  
  Median_CNV = apply(
    cnv_raw_test,
    2,
    median,
    na.rm = TRUE
  ),
  
  Mean_CNV = colMeans(
    cnv_raw_test,
    na.rm = TRUE
  ),
  
  Factor1 = factor_scores[, "Factor1"],
  
  Robust_State = robust_state
)


# ============================================================
# 2. FACTOR 1 vs GLOBAL CNV LEVEL
# ============================================================

cor.test(
  cnv_level_df$Factor1,
  cnv_level_df$Median_CNV,
  method = "spearman"
)

cor.test(
  cnv_level_df$Factor1,
  cnv_level_df$Mean_CNV,
  method = "spearman"
)


# ============================================================
# 3. GLOBAL CNV LEVEL: STATE 1 vs STATE 2
# ============================================================

cnv_level_df %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    Median_of_Median_CNV =
      median(Median_CNV, na.rm = TRUE),
    
    Median_of_Mean_CNV =
      median(Mean_CNV, na.rm = TRUE),
    
    .groups = "drop"
  )


# ============================================================
# 4. STATISTICAL TEST
# ============================================================

wilcox.test(
  Median_CNV ~ Robust_State,
  data = cnv_level_df,
  exact = FALSE
)

wilcox.test(
  Mean_CNV ~ Robust_State,
  data = cnv_level_df,
  exact = FALSE
)
library(mclust)

# Factor 1'i çıkar
X_noF1 <- as.matrix(
  factor_scores[, paste0("Factor", 2:10)]
)

# Aynı model: K = 2, VVE
fit_noF1 <- Mclust(
  X_noF1,
  G = 2,
  modelNames = "VVE",
  verbose = FALSE
)

state_noF1 <- fit_noF1$classification

# Yeni state boyutları
table(state_noF1)

# Mevcut state vs F1 çıkarılmış state
table(
  Original_State = robust_state,
  NoF1_State = state_noF1
)

# Adjusted Rand Index
ARI_noF1 <- adjustedRandIndex(
  robust_state,
  state_noF1
)

ARI_noF1
library(mclust)

set.seed(123)

# Factor1 excluded
X_noF1 <- as.matrix(
  factor_scores[, paste0("Factor", 2:10)]
)

# Reference K=2 model
reference_fit_noF1 <- Mclust(
  X_noF1,
  G = 2,
  modelNames = "VVE",
  verbose = FALSE
)

reference_state_noF1 <- reference_fit_noF1$classification


# ============================================================
# BOOTSTRAP STABILITY
# ============================================================

n_boot <- 100

bootstrap_ARI_noF1 <- rep(NA_real_, n_boot)

for (b in seq_len(n_boot)) {
  
  # bootstrap samples
  idx <- sample(
    seq_len(nrow(X_noF1)),
    size = nrow(X_noF1),
    replace = TRUE
  )
  
  X_boot <- X_noF1[idx, , drop = FALSE]
  
  # Fit model
  fit_boot <- tryCatch(
    Mclust(
      X_boot,
      G = 2,
      modelNames = "VVE",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_boot)) {
    
    # Predict ALL original patients
    pred <- tryCatch(
      mclust::predict.Mclust(
        fit_boot,
        X_noF1
      )$classification,
      error = function(e) NULL
    )
    
    if (!is.null(pred)) {
      
      bootstrap_ARI_noF1[b] <- adjustedRandIndex(
        reference_state_noF1,
        pred
      )
    }
  }
}


# ============================================================
# RESULTS
# ============================================================

summary(bootstrap_ARI_noF1)

sd(
  bootstrap_ARI_noF1,
  na.rm = TRUE
)

quantile(
  bootstrap_ARI_noF1,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

sum(!is.na(bootstrap_ARI_noF1))
# ============================================================
# ROBUST STATE CHARACTERIZATION
# ============================================================

# Factor means
factor_summary <- data.frame(
  Robust_State = robust_state,
  factor_scores[, paste0("Factor", 1:10)]
) %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    dplyr::across(
      starts_with("Factor"),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

factor_summary


# ============================================================
# CLINICAL SUMMARY
# ============================================================

clinical_state_summary <- robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    Median_age = median(age_years, na.rm = TRUE),
    
    BRCA_mutated_n = sum(BRCA_any_mut == TRUE, na.rm = TRUE),
    
    BRCA_mutated_percent =
      100 * mean(BRCA_any_mut == TRUE, na.rm = TRUE),
    
    Dead_n = sum(vital_status == "Dead", na.rm = TRUE),
    
    Dead_percent =
      100 * mean(vital_status == "Dead", na.rm = TRUE),
    
    .groups = "drop"
  )

clinical_state_summary


# ============================================================
# TOP RNA FEATURES FOR EACH STATE
# ============================================================

state2_rna <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05, Difference > 0) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(10)

state1_rna <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05, Difference < 0) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(10)

state1_rna
state2_rna


# ============================================================
# TOP CNV FEATURES FOR EACH STATE
# ============================================================

state2_cnv <- cnv_robust_stats %>%
  dplyr::filter(FDR < 0.05, Difference > 0) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(10)

state1_cnv <- cnv_robust_stats %>%
  dplyr::filter(FDR < 0.05, Difference < 0) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(10)

state1_cnv
state2_cnv
state1_rna
state2_rna
library(survival)
library(dplyr)

# ============================================================
# 1. SURVIVAL DATA
# ============================================================

surv_df <- robust_state_clinical %>%
  dplyr::mutate(
    OS_time = days_to_death,
    OS_event = ifelse(vital_status == "Dead", 1, 0),
    
    # yaşayan hastalarda follow-up kullan
    OS_time = ifelse(
      OS_event == 0,
      days_to_last_follow_up,
      OS_time
    ),
    
    Robust_State = factor(
      Robust_State,
      levels = c(1, 2)
    )
  )

# Kontrol
surv_df %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events = sum(OS_event, na.rm = TRUE),
    Median_followup = median(OS_time, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 2. UNADJUSTED COX
# State 1 = reference
# ============================================================

cox_unadjusted <- coxph(
  Surv(OS_time, OS_event) ~ Robust_State,
  data = surv_df
)

summary(cox_unadjusted)


# ============================================================
# 3. FIGO'YU SADELEŞTİR
# I/II vs III vs IV
# ============================================================

surv_df <- surv_df %>%
  dplyr::mutate(
    FIGO_group = dplyr::case_when(
      grepl("^Stage I[^V]", figo_stage, ignore.case = TRUE) ~ "I_II",
      grepl("^Stage II", figo_stage, ignore.case = TRUE) ~ "I_II",
      grepl("^Stage III", figo_stage, ignore.case = TRUE) ~ "III",
      grepl("^Stage IV", figo_stage, ignore.case = TRUE) ~ "IV",
      TRUE ~ NA_character_
    ),
    
    FIGO_group = factor(
      FIGO_group,
      levels = c("I_II", "III", "IV")
    )
  )

table(
  surv_df$FIGO_group,
  useNA = "ifany"
)


# ============================================================
# 4. AGE + FIGO ADJUSTED COX
# ============================================================

cox_adjusted <- coxph(
  Surv(OS_time, OS_event) ~
    Robust_State +
    age_years +
    FIGO_group,
  data = surv_df
)

summary(cox_adjusted)


# ============================================================
# 5. PROPORTIONAL HAZARDS ASSUMPTION
# ============================================================

cox.zph(cox_adjusted)
table(
  surv_df$figo_stage,
  useNA = "ifany"
)
surv_df <- surv_df %>%
  dplyr::mutate(
    FIGO_group = dplyr::case_when(
      figo_stage %in% c(
        "Stage IC",
        "Stage IIA",
        "Stage IIB",
        "Stage IIC"
      ) ~ "I_II",
      
      figo_stage %in% c(
        "Stage IIIA",
        "Stage IIIB",
        "Stage IIIC"
      ) ~ "III",
      
      figo_stage == "Stage IV" ~ "IV",
      
      TRUE ~ NA_character_
    ),
    
    FIGO_group = factor(
      FIGO_group,
      levels = c("I_II", "III", "IV")
    )
  )

table(surv_df$FIGO_group, useNA = "ifany")
cox_adjusted <- coxph(
  Surv(OS_time, OS_event) ~
    Robust_State +
    age_years +
    FIGO_group,
  data = surv_df
)

summary(cox_adjusted)

cox.zph(cox_adjusted)
cox_timevarying <- coxph(
  Surv(OS_time, OS_event) ~
    Robust_State +
    FIGO_group +
    age_years +
    tt(age_years),
  data = surv_df,
  tt = function(x, t, ...) {
    x * log(t + 1)
  }
)

summary(cox_timevarying)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

# ============================================================
# 1. FACTOR MATRIX
# ============================================================

factor_mat <- factor_state_means %>%
  tibble::column_to_rownames("Robust_State") %>%
  as.matrix()

factor_mat <- t(factor_mat)

colnames(factor_mat) <- paste0("State ", colnames(factor_mat))

# sadece anlamlı factorlar
factor_keep <- c(
  "Factor1",
  "Factor7",
  "Factor4",
  "Factor6",
  "Factor3"
)

factor_mat <- factor_mat[factor_keep, , drop = FALSE]


# ============================================================
# 2. RNA MATRIX
# Top discriminatory RNA genes
# ============================================================

rna_top <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(20)

rna_mat <- as.matrix(
  rna_top[, c("Median_State1", "Median_State2")]
)

rownames(rna_mat) <- sub("^RNA_", "", rna_top$Feature)

# row-wise Z score
rna_mat <- t(scale(t(rna_mat)))

colnames(rna_mat) <- c("State 1", "State 2")


# ============================================================
# 3. CNV CHROMOSOME MATRIX
# use % significant selected CNV features
# ============================================================

cnv_chr_plot <- cnv_chr_robust %>%
  dplyr::filter(Significant > 0) %>%
  dplyr::arrange(chr_num)

# Directional representation:
# State1 = negative, State2 = positive
cnv_chr_mat <- cbind(
  State1 = -cnv_chr_plot$Percent_significant,
  State2 =  cnv_chr_plot$Percent_significant
)

rownames(cnv_chr_mat) <- cnv_chr_plot$chromosome

colnames(cnv_chr_mat) <- c("State 1", "State 2")


# ============================================================
# 4. CLINICAL MATRIX
# actual values for annotation
# ============================================================

clinical_plot <- robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    
    Age = median(age_years, na.rm = TRUE),
    
    BRCA_percent =
      100 * mean(BRCA_any_mut == TRUE, na.rm = TRUE),
    
    FIGO_IV_percent =
      100 * mean(figo_stage == "Stage IV", na.rm = TRUE),
    
    G3_percent =
      100 * mean(tumor_grade == "G3", na.rm = TRUE),
    
    .groups = "drop"
  )

clinical_mat <- t(
  as.matrix(
    clinical_plot[, -1]
  )
)

colnames(clinical_mat) <- paste0(
  "State ",
  clinical_plot$Robust_State
)

rownames(clinical_mat) <- c(
  "Median age",
  "BRCA mutation (%)",
  "FIGO IV (%)",
  "Grade G3 (%)"
)

# Z-score ONLY for heatmap colors
clinical_z <- t(scale(t(clinical_mat)))


# ============================================================
# 5. COLOR FUNCTIONS
# ============================================================

factor_lim <- max(abs(factor_mat), na.rm = TRUE)

factor_col <- circlize::colorRamp2(
  c(-factor_lim, 0, factor_lim),
  c("navy", "white", "firebrick")
)

rna_col <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick")
)

cnv_col <- circlize::colorRamp2(
  c(-100, 0, 100),
  c("navy", "white", "firebrick")
)

clinical_col <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("navy", "white", "firebrick")
)


# ============================================================
# 6. HEATMAPS
# ============================================================

ht_factor <- Heatmap(
  factor_mat,
  name = "Mean factor\nscore",
  col = factor_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "Robust Molecular States",
  row_title = "MOFA factors",
  rect_gp = gpar(col = "white"),
  column_names_gp = gpar(fontface = "bold"),
  row_names_gp = gpar(fontsize = 9)
)


ht_rna <- Heatmap(
  rna_mat,
  name = "RNA\nZ-score",
  col = rna_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "RNA signature",
  rect_gp = gpar(col = "white"),
  row_names_gp = gpar(fontsize = 8)
)


ht_cnv <- Heatmap(
  cnv_chr_mat,
  name = "CNV\nsignal",
  col = cnv_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "CNV chromosome pattern",
  rect_gp = gpar(col = "white"),
  row_names_gp = gpar(fontsize = 8),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.text(
      paste0(
        round(abs(cnv_chr_mat[i, j]), 1),
        "%"
      ),
      x,
      y,
      gp = gpar(fontsize = 7)
    )
  }
)


ht_clinical <- Heatmap(
  clinical_z,
  name = "Clinical\nZ-score",
  col = clinical_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Clinical characteristics",
  rect_gp = gpar(col = "white"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    value <- clinical_mat[i, j]
    
    label <- if (
      rownames(clinical_mat)[i] == "Median age"
    ) {
      sprintf("%.1f", value)
    } else {
      sprintf("%.1f%%", value)
    }
    
    grid.text(
      label,
      x,
      y,
      gp = gpar(fontsize = 8)
    )
  }
)


# ============================================================
# 7. DRAW FINAL FIGURE
# ============================================================

final_ht <-
  ht_factor %v%
  ht_rna %v%
  ht_cnv %v%
  ht_clinical

draw(
  final_ht,
  heatmap_legend_side = "right"
)
# ============================================================
# CORRECT CNV PANEL
# ============================================================

cnv_chr_plot <- cnv_chr_robust %>%
  dplyr::filter(Significant > 0) %>%
  dplyr::arrange(chr_num)

# State-level relative CNV direction
cnv_chr_mat <- cbind(
  "State 1" = rep(0, nrow(cnv_chr_plot)),
  "State 2" = cnv_chr_plot$Median_difference
)

rownames(cnv_chr_mat) <- cnv_chr_plot$chromosome

cnv_lim <- max(abs(cnv_chr_mat), na.rm = TRUE)

cnv_col <- circlize::colorRamp2(
  c(-cnv_lim, 0, cnv_lim),
  c("navy", "white", "firebrick")
)

ht_cnv <- Heatmap(
  cnv_chr_mat,
  name = "Relative CNV\ndifference",
  col = cnv_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "CNV chromosome pattern",
  rect_gp = gpar(col = "white"),
  row_names_gp = gpar(fontsize = 8),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    # State 2 hücresinde % significant feature göster
    if (j == 2) {
      grid.text(
        paste0(
          round(cnv_chr_plot$Percent_significant[i], 1),
          "%"
        ),
        x, y,
        gp = gpar(fontsize = 7)
      )
    }
  }
)
###daha duzgun 
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

# ============================================================
# 1. MOFA FACTORS
# ============================================================

factor_keep <- c(
  "Factor1",
  "Factor7",
  "Factor4",
  "Factor6",
  "Factor3"
)

factor_mat <- factor_state_means %>%
  tibble::column_to_rownames("Robust_State") %>%
  as.matrix() %>%
  t()

factor_mat <- factor_mat[factor_keep, , drop = FALSE]

colnames(factor_mat) <- c("State 1", "State 2")


# ============================================================
# 2. RNA — TOP 20 DISCRIMINATORY GENES
# ============================================================

rna_top <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(dplyr::desc(abs(Difference))) %>%
  head(20)

rna_mat_raw <- as.matrix(
  rna_top[, c("Median_State1", "Median_State2")]
)

rownames(rna_mat_raw) <- sub("^RNA_", "", rna_top$Feature)

rna_mat <- t(scale(t(rna_mat_raw)))

colnames(rna_mat) <- c("State 1", "State 2")


# ============================================================
# 3. CNV — CHROMOSOME LEVEL
# ============================================================

cnv_chr_plot <- cnv_chr_robust %>%
  dplyr::filter(Significant > 0) %>%
  dplyr::arrange(chr_num)

# Difference relative to State 1
cnv_chr_mat <- cbind(
  "State 1" = 0,
  "State 2" = cnv_chr_plot$Median_difference
)

rownames(cnv_chr_mat) <- cnv_chr_plot$chromosome


# ============================================================
# 4. CLINICAL CHARACTERISTICS
# ============================================================

clinical_plot <- robust_state_clinical %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    
    Age = median(age_years, na.rm = TRUE),
    
    BRCA_percent =
      100 * mean(BRCA_any_mut == TRUE, na.rm = TRUE),
    
    FIGO_IV_percent =
      100 * mean(figo_stage == "Stage IV", na.rm = TRUE),
    
    G3_percent =
      100 * mean(tumor_grade == "G3", na.rm = TRUE),
    
    .groups = "drop"
  )

clinical_mat <- t(
  as.matrix(clinical_plot[, -1])
)

colnames(clinical_mat) <- c("State 1", "State 2")

rownames(clinical_mat) <- c(
  "Median age",
  "BRCA mutation (%)",
  "FIGO IV (%)",
  "Grade G3 (%)"
)

clinical_z <- t(scale(t(clinical_mat)))


# ============================================================
# 5. NEW COLOR PALETTE
# ============================================================

# Blue → very light grey → orange/red
factor_lim <- max(abs(factor_mat), na.rm = TRUE)

factor_col <- circlize::colorRamp2(
  c(-factor_lim, 0, factor_lim),
  c("#3B6FB6", "#F7F7F7", "#D95F45")
)

rna_col <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("#4575B4", "#F7F7F7", "#D6604D")
)

cnv_lim <- max(abs(cnv_chr_mat), na.rm = TRUE)

cnv_col <- circlize::colorRamp2(
  c(-cnv_lim, 0, cnv_lim),
  c("#4C78A8", "#F7F7F7", "#E07B39")
)

clinical_col <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("#5B8DB8", "#F7F7F7", "#D9825B")
)


# ============================================================
# Helper: dashed line between State 1 and State 2
# ============================================================

add_state_separator <- function() {
  
  grid.lines(
    x = unit(c(0.5, 0.5), "npc"),
    y = unit(c(0, 1), "npc"),
    gp = gpar(
      col = "#555555",
      lty = 2,
      lwd = 1.2
    )
  )
}


# ============================================================
# 6. MOFA FACTOR HEATMAP
# ============================================================

ht_factor <- Heatmap(
  
  factor_mat,
  
  name = "Mean factor\nscore",
  
  col = factor_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_title = "Robust Molecular States",
  
  row_title = "MOFA factors",
  
  rect_gp = gpar(
    col = "white",
    lwd = 1
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9,
    fontface = "bold"
  ),
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    value <- factor_mat[i, j]
    
    # white text for strong/dark cells
    text_col <- ifelse(
      abs(value) > factor_lim * 0.45,
      "white",
      "#222222"
    )
    
    grid.text(
      sprintf("%.2f", value),
      x, y,
      gp = gpar(
        fontsize = 8,
        col = text_col,
        fontface = "bold"
      )
    )
  },
  
  layer_fun = function(j, i, x, y, width, height, fill) {
    add_state_separator()
  }
)


# ============================================================
# 7. RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  
  rna_mat,
  
  name = "RNA\nZ-score",
  
  col = rna_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  row_title = "RNA signature",
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  row_names_gp = gpar(
    fontsize = 8
  ),
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  layer_fun = function(j, i, x, y, width, height, fill) {
    add_state_separator()
  }
)


# ============================================================
# 8. CNV HEATMAP
# ============================================================

ht_cnv <- Heatmap(
  
  cnv_chr_mat,
  
  name = "Relative CNV\ndifference",
  
  col = cnv_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  row_title = "CNV chromosome pattern",
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  row_names_gp = gpar(
    fontsize = 8
  ),
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    # Percentage only shown for State 2
    if (j == 2) {
      
      value <- cnv_chr_plot$Percent_significant[i]
      
      # automatic text colour
      text_col <- ifelse(
        cnv_chr_mat[i, j] > cnv_lim * 0.45,
        "white",
        "#222222"
      )
      
      grid.text(
        sprintf("%.1f%%", value),
        x, y,
        gp = gpar(
          fontsize = 7.5,
          col = text_col,
          fontface = "bold"
        )
      )
    }
  },
  
  layer_fun = function(j, i, x, y, width, height, fill) {
    add_state_separator()
  }
)


# ============================================================
# 9. CLINICAL HEATMAP
# ============================================================

ht_clinical <- Heatmap(
  
  clinical_z,
  
  name = "Clinical\nZ-score",
  
  col = clinical_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title = "Clinical characteristics",
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    value <- clinical_mat[i, j]
    
    label <- if (
      rownames(clinical_mat)[i] == "Median age"
    ) {
      sprintf("%.1f", value)
    } else {
      sprintf("%.1f%%", value)
    }
    
    z <- clinical_z[i, j]
    
    text_col <- ifelse(
      abs(z) > 0.7,
      "white",
      "#222222"
    )
    
    grid.text(
      label,
      x, y,
      gp = gpar(
        fontsize = 8.5,
        col = text_col,
        fontface = "bold"
      )
    )
  },
  
  layer_fun = function(j, i, x, y, width, height, fill) {
    add_state_separator()
  }
)


# ============================================================
# 10. COMBINE
# ============================================================

final_ht <-
  
  ht_factor %v%
  
  ht_rna %v%
  
  ht_cnv %v%
  
  ht_clinical


# ============================================================
# 11. DRAW
# ============================================================

draw(
  
  final_ht,
  
  heatmap_legend_side = "right",
  
  padding = unit(
    c(5, 12, 5, 5),
    "mm"
  )
)

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(tibble)


# ============================================================
# 1. MOFA FACTOR PANEL
# ============================================================

# State'ler arasında anlamlı farklılık gösteren factorlar
factor_keep <- c(
  "Factor1",
  "Factor7",
  "Factor4",
  "Factor6",
  "Factor3"
)

factor_mat <- factor_state_means %>%
  tibble::column_to_rownames("Robust_State") %>%
  as.matrix() %>%
  t()

factor_mat <- factor_mat[
  factor_keep,
  ,
  drop = FALSE
]

colnames(factor_mat) <- c(
  "State 1",
  "State 2"
)


# ============================================================
# 2. RNA PANEL
# Top 20 discriminatory genes
# ============================================================

rna_top <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(
    dplyr::desc(abs(Difference))
  ) %>%
  head(20)

rna_mat_raw <- as.matrix(
  rna_top[
    ,
    c(
      "Median_State1",
      "Median_State2"
    )
  ]
)

rownames(rna_mat_raw) <- sub(
  "^RNA_",
  "",
  rna_top$Feature
)

# Row-wise Z-score for visualization
rna_mat <- t(
  scale(
    t(rna_mat_raw)
  )
)

colnames(rna_mat) <- c(
  "State 1",
  "State 2"
)


# ============================================================
# 3. CNV PANEL
#
# IMPORTANT:
# This panel shows State2 - State1 relative CNV difference.
#
# It does NOT mean State1 has no CNV.
#
# Percentage shown inside each cell =
# percentage of tested MOFA CNV features on that chromosome
# significantly different between the two states.
# ============================================================

cnv_chr_plot <- cnv_chr_robust %>%
  dplyr::filter(
    Significant > 0
  ) %>%
  dplyr::arrange(
    chr_num
  )

cnv_chr_mat <- matrix(
  cnv_chr_plot$Median_difference,
  ncol = 1
)

rownames(cnv_chr_mat) <-
  cnv_chr_plot$chromosome

colnames(cnv_chr_mat) <-
  "State 2 - State 1"


# ============================================================
# 4. CLINICAL PANEL
# ============================================================

clinical_plot <- robust_state_clinical %>%
  dplyr::group_by(
    Robust_State
  ) %>%
  dplyr::summarise(
    
    Age =
      median(
        age_years,
        na.rm = TRUE
      ),
    
    BRCA_percent =
      100 *
      mean(
        BRCA_any_mut == TRUE,
        na.rm = TRUE
      ),
    
    FIGO_IV_percent =
      100 *
      mean(
        figo_stage == "Stage IV",
        na.rm = TRUE
      ),
    
    G3_percent =
      100 *
      mean(
        tumor_grade == "G3",
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


clinical_mat <- t(
  as.matrix(
    clinical_plot[, -1]
  )
)

colnames(clinical_mat) <- c(
  "State 1",
  "State 2"
)

rownames(clinical_mat) <- c(
  "Median age",
  "BRCA mutation (%)",
  "FIGO IV (%)",
  "Grade G3 (%)"
)

# Z-score ONLY for heatmap colors
clinical_z <- t(
  scale(
    t(clinical_mat)
  )
)


# ============================================================
# 5. COLOR PALETTES
# ============================================================

# ---------- MOFA ----------

factor_lim <- max(
  abs(factor_mat),
  na.rm = TRUE
)

factor_col <- circlize::colorRamp2(
  c(
    -factor_lim,
    0,
    factor_lim
  ),
  c(
    "#4C78A8",
    "#F7F7F7",
    "#D65F4A"
  )
)


# ---------- RNA ----------

rna_col <- circlize::colorRamp2(
  c(
    -1,
    0,
    1
  ),
  c(
    "#5B7DB1",
    "#F7F7F7",
    "#DE7868"
  )
)


# ---------- CNV ----------

cnv_lim <- max(
  abs(cnv_chr_mat),
  na.rm = TRUE
)

if (cnv_lim == 0) {
  cnv_lim <- 1
}

cnv_col <- circlize::colorRamp2(
  c(
    -cnv_lim,
    0,
    cnv_lim
  ),
  c(
    "#4C78A8",
    "#F7F7F7",
    "#E17C35"
  )
)


# ---------- Clinical ----------

clinical_col <- circlize::colorRamp2(
  c(
    -1,
    0,
    1
  ),
  c(
    "#6F94B8",
    "#F7F7F7",
    "#DF8A68"
  )
)


# ============================================================
# 6. HELPER:
# DASHED LINE BETWEEN STATE 1 AND STATE 2
# ============================================================

add_state_separator <- function() {
  
  grid.lines(
    x = unit(
      c(0.5, 0.5),
      "npc"
    ),
    
    y = unit(
      c(0, 1),
      "npc"
    ),
    
    gp = gpar(
      col = "#666666",
      lty = 2,
      lwd = 1.2
    )
  )
}


# ============================================================
# 7. MOFA FACTOR HEATMAP
# ============================================================

ht_factor <- Heatmap(
  
  factor_mat,
  
  name = "Mean factor\nscore",
  
  col = factor_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  column_title =
    "Robust Molecular States",
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  row_title =
    "MOFA factors",
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  rect_gp = gpar(
    col = "white",
    lwd = 1
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9,
    fontface = "bold"
  ),
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    value <- factor_mat[i, j]
    
    text_col <- ifelse(
      abs(value) >
        factor_lim * 0.45,
      "white",
      "#222222"
    )
    
    grid.text(
      sprintf(
        "%.2f",
        value
      ),
      
      x,
      y,
      
      gp = gpar(
        fontsize = 8,
        col = text_col,
        fontface = "bold"
      )
    )
  },
  
  layer_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    add_state_separator()
  }
)


# ============================================================
# 8. RNA HEATMAP
# ============================================================

ht_rna <- Heatmap(
  
  rna_mat,
  
  name = "RNA\nZ-score",
  
  col = rna_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  row_title =
    "RNA signature",
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  row_names_gp = gpar(
    fontsize = 8
  ),
  
  layer_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    add_state_separator()
  }
)


# ============================================================
# 9. CNV HEATMAP
#
# ONE COLUMN ONLY:
# State 2 - State 1
# ============================================================

ht_cnv <- Heatmap(
  
  cnv_chr_mat,
  
  name =
    "Relative CNV\nDifference",
  
  col = cnv_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title =
    "CNV chromosome pattern",
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  column_names_gp = gpar(
    fontsize = 9,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 8
  ),
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    diff_value <-
      cnv_chr_mat[i, j]
    
    percentage <-
      cnv_chr_plot$Percent_significant[i]
    
    text_col <- ifelse(
      abs(diff_value) >
        cnv_lim * 0.45,
      "white",
      "#222222"
    )
    
    grid.text(
      
      sprintf(
        "%.1f%%",
        percentage
      ),
      
      x,
      y,
      
      gp = gpar(
        fontsize = 8,
        col = text_col,
        fontface = "bold"
      )
    )
  }
)


# ============================================================
# 10. CLINICAL HEATMAP
# ============================================================

ht_clinical <- Heatmap(
  
  clinical_z,
  
  name =
    "Clinical\nZ-score",
  
  col = clinical_col,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title =
    "Clinical characteristics",
  
  row_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.8
  ),
  
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "bold"
  ),
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    value <-
      clinical_mat[i, j]
    
    label <- if (
      rownames(clinical_mat)[i] ==
      "Median age"
    ) {
      
      sprintf(
        "%.1f",
        value
      )
      
    } else {
      
      sprintf(
        "%.1f%%",
        value
      )
    }
    
    z <- clinical_z[i, j]
    
    text_col <- ifelse(
      abs(z) > 0.7,
      "white",
      "#222222"
    )
    
    grid.text(
      
      label,
      
      x,
      y,
      
      gp = gpar(
        fontsize = 8.5,
        col = text_col,
        fontface = "bold"
      )
    )
  },
  
  layer_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    add_state_separator()
  }
)


# ============================================================
# 11. COMBINE ALL PANELS
# ============================================================

final_ht <-
  
  ht_factor %v%
  
  ht_rna %v%
  
  ht_cnv %v%
  
  ht_clinical


# ============================================================
# 12. DRAW FINAL FIGURE
# ============================================================

draw(
  
  final_ht,
  
  heatmap_legend_side = "right",
  
  padding = unit(
    c(
      5,
      15,
      5,
      5
    ),
    "mm"
  )
)

# 110 significant RNA genes
rna_110 <- rna_robust_stats %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      Difference > 0 ~ "State2_high",
      Difference < 0 ~ "State1_high",
      TRUE ~ "No_difference"
    )
  ) %>%
  dplyr::arrange(FDR)

# Gör
View(rna_110)

# Kontrol
dim(rna_110)
table(rna_110$Direction)

# Excel olarak kaydet
writexl::write_xlsx(
  rna_110,
  "Robust_State1_vs_State2_110_RNA_genes.xlsx"
)
library(writexl)

write_xlsx(
  rna_110,
  "Robust_State1_vs_State2_110_RNA_genes.xlsx"
)


workflow_plot <- grViz("
digraph molecular_state_workflow {

# SENİN AYNI KODUN BURADA

}
")

library(DiagrammeRsvg)
library(rsvg)

# SVG'ye çevir
svg_code <- DiagrammeRsvg::export_svg(workflow_plot)

# SVG kaydet
writeLines(
  svg_code,
  "Robust_Molecular_State_Workflow.svg"
)

# Yüksek çözünürlüklü PNG kaydet
rsvg::rsvg_png(
  charToRaw(svg_code),
  file = "Robust_Molecular_State_Workflow.png",
  width = 3000
)

# PDF kaydet
rsvg::rsvg_pdf(
  charToRaw(svg_code),
  file = "Robust_Molecular_State_Workflow.pdf"
)
# ============================================================
# MOLECULAR STATE DISCOVERY WORKFLOW - GGPLOT VERSION
# ============================================================

library(ggplot2)
library(grid)


# ============================================================
# 1. BOX POSITIONS
# ============================================================

nodes <- data.frame(
  
  id = c(
    "A",
    "B1", "B2", "B3",
    "C1", "C2", "C3",
    "D", "E", "F",
    "G", "H",
    "I",
    "J1", "J2", "J3",
    "K1", "K2", "K3", "K4",
    "L", "M"
  ),
  
  x = c(
    6,
    2, 6, 10,
    2, 6, 10,
    6, 6, 6,
    3.5, 8.5,
    6,
    3.5, 8.5, 8.5,
    1.5, 4.5, 7.5, 10.5,
    6, 6
  ),
  
  y = c(
    16,
    14.5, 14.5, 14.5,
    13, 13, 13,
    11.5, 10, 8.5,
    7, 7,
    5.5,
    3.8, 3.8, 2.4,
    0.5, 0.5, 0.5, 0.5,
    -1.5, -3
  ),
  
  label = c(
    
    "TCGA-OV Primary Tumors\n277 patients",
    
    "RNA-seq\nVST normalization\nProtein-coding genes",
    
    "Copy Number Variation\nAutosomal protein-coding genes\nSample-median centering",
    
    "Somatic Mutations\nBinary mutation matrix",
    
    
    "RNA\nTop 10% variable\n1,722 features",
    
    "CNV\nTop 10% variable\n1,902 features",
    
    "Mutation\n>=10 patients\n107 features",
    
    
    "MOFA2\nBayesian multi-omics integration\n20 latent factors",
    
    "Factor Evaluation\nVariance explained + biology\nFactors 1-10 retained",
    
    "Gaussian Mixture Modeling\nCandidate K = 2-8",
    
    
    "Initial 5-state solution\nLow bootstrap stability\nMedian ARI ~ 0.35",
    
    "Stability-based model selection\nBootstrap + factor sensitivity",
    
    
    "ROBUST K = 2 SOLUTION\nState 1: n = 226 | State 2: n = 51",
    
    
    "K=2 Bootstrap Stability\nMean ARI ~ 0.78\nMedian ARI ~ 0.80",
    
    "Factor 1 Sensitivity\nRemove Factor 1\nOriginal vs no-F1 ARI = 0.79",
    
    "No-F1 Bootstrap\nMean ARI = 0.83\nMedian ARI = 0.84\n99/100 successful",
    
    
    "RNA Characterization\n110 genes, FDR < 0.05\n34 State1-high | 76 State2-high\nState2: ZIC1 / ZIC2 / ZIC4",
    
    "CNV Characterization\n1,152 / 1,902 features significant\nBroad relative CNV differences",
    
    "Clinical Association\nAge: significant\nBRCA / FIGO / Grade: NS",
    
    "GO Biological Process\nState1 vs State2 RNA\nNo significant GO terms",
    
    
    "Overall Survival Analysis\nAdjusted Cox model\nState2 HR = 1.20\n95% CI 0.79-1.81; p = 0.396",
    
    
    "ROBUST MOLECULAR PATIENT STATES\nReproducible multi-omic configurations\nNot OS-defined prognostic subgroups"
  ),
  
  group = c(
    "cohort",
    
    rep("input", 3),
    rep("feature", 3),
    
    "mofa",
    "mofa",
    "cluster",
    
    "failed",
    "cluster",
    
    "robust",
    
    rep("validation", 3),
    
    rep("characterization", 4),
    
    "survival",
    "final"
  )
)


# ============================================================
# 2. COLORS
# ============================================================

node_colors <- c(
  
  cohort = "#DCEAF7",
  input = "#E8F1FA",
  feature = "#EDF5FB",
  
  mofa = "#D9EAD3",
  
  cluster = "#FFF0D6",
  
  failed = "#FBE3E3",
  
  robust = "#CDECCF",
  
  validation = "#E2F3E3",
  
  characterization = "#E6E0F3",
  
  survival = "#F4E4F2",
  
  final = "#BFE3D0"
)


# ============================================================
# 3. CONNECTIONS
# ============================================================

edges <- data.frame(
  
  from = c(
    
    "A", "A", "A",
    
    "B1", "B2", "B3",
    
    "C1", "C2", "C3",
    
    "D",
    "E",
    
    "F", "F",
    
    "G",
    
    "H",
    
    "I", "I",
    
    "J2",
    
    "J1", "J3",
    
    "I", "I", "I",
    
    "K1", "K2", "K3", "K4",
    
    "L"
  ),
  
  to = c(
    
    "B1", "B2", "B3",
    
    "C1", "C2", "C3",
    
    "D", "D", "D",
    
    "E",
    "F",
    
    "G", "H",
    
    "H",
    
    "I",
    
    "J1", "J2",
    
    "J3",
    
    "K1", "K1",
    
    "K2", "K3", "K4",
    
    "L", "L", "L", "L",
    
    "M"
  )
)


# ============================================================
# 4. ADD COORDINATES TO EDGES
# ============================================================

edges <- merge(
  edges,
  nodes[, c("id", "x", "y")],
  by.x = "from",
  by.y = "id"
)

colnames(edges)[
  colnames(edges) %in% c("x", "y")
] <- c("x_from", "y_from")


edges <- merge(
  edges,
  nodes[, c("id", "x", "y")],
  by.x = "to",
  by.y = "id"
)

colnames(edges)[
  colnames(edges) %in% c("x", "y")
] <- c("x_to", "y_to")


# ============================================================
# 5. PLOT
# ============================================================

workflow_plot <- ggplot() +
  
  # ----------------------------------------------------------
# ARROWS
# ----------------------------------------------------------

geom_segment(
  
  data = edges,
  
  aes(
    x = x_from,
    y = y_from - 0.38,
    xend = x_to,
    yend = y_to + 0.38
  ),
  
  linewidth = 0.45,
  colour = "#718096",
  
  arrow = arrow(
    length = unit(0.13, "cm"),
    type = "closed"
  )
) +
  
  
  # ----------------------------------------------------------
# BOXES
# ----------------------------------------------------------

geom_label(
  
  data = nodes,
  
  aes(
    x = x,
    y = y,
    label = label,
    fill = group
  ),
  
  size = 3.0,
  
  label.size = 0.35,
  
  label.padding = unit(
    0.25,
    "lines"
  ),
  
  label.r = unit(
    0.15,
    "lines"
  ),
  
  colour = "#263238",
  
  lineheight = 0.95,
  
  fontface = "plain"
) +
  
  
  # ----------------------------------------------------------
# COLORS
# ----------------------------------------------------------

scale_fill_manual(
  values = node_colors
) +
  
  
  # ----------------------------------------------------------
# LIMITS
# ----------------------------------------------------------

coord_cartesian(
  xlim = c(0, 12),
  ylim = c(-4, 17)
) +
  
  
  # ----------------------------------------------------------
# CLEAN THEME
# ----------------------------------------------------------

theme_void() +
  
  theme(
    
    legend.position = "none",
    
    plot.margin = margin(
      20,
      20,
      20,
      20
    )
  )


# ============================================================
# 6. SHOW
# ============================================================

workflow_plot


# ============================================================
# 7. SAVE PDF
# ============================================================

ggsave(
  
  filename =
    "Robust_Molecular_State_Workflow.pdf",
  
  plot =
    workflow_plot,
  
  width = 14,
  
  height = 18,
  
  units = "in",
  
  device = cairo_pdf
)


# ============================================================
# 8. ALSO SAVE HIGH-RESOLUTION PNG
# ============================================================

ggsave(
  
  filename =
    "Robust_Molecular_State_Workflow.png",
  
  plot =
    workflow_plot,
  
  width = 14,
  
  height = 18,
  
  units = "in",
  
  dpi = 400,
  
  bg = "white"
)


# ============================================================
# 9. CHECK LOCATION
# ============================================================

getwd()

# ============================================================
# ACTIONABLE / TREATMENT-RELEVANT MUTATION PANEL
# ============================================================

actionable_genes <- c(
  "BRCA1",
  "BRCA2",
  "CDK12",
  "NF1",
  "RB1",
  "PTEN",
  "PIK3CA",
  "KRAS"
)

# Patient x gene binary mutation matrix
actionable_mut <- sapply(
  actionable_genes,
  function(gene) {
    
    rownames(factor_scores) %in%
      unique(
        mutation_final$patient_id[
          mutation_final$Hugo_Symbol == gene
        ]
      )
  }
)

actionable_mut <- as.data.frame(actionable_mut)

rownames(actionable_mut) <- rownames(factor_scores)

# Add patient ID
actionable_mut$submitter_id <- rownames(actionable_mut)

# Add robust state
actionable_mut$Robust_State <- robust_state[
  match(
    actionable_mut$submitter_id,
    rownames(factor_scores)
  )
]

# Put ID + state first
actionable_mut <- actionable_mut %>%
  dplyr::select(
    submitter_id,
    Robust_State,
    dplyr::everything()
  )

View(actionable_mut)

# Number of mutated patients for each gene
colSums(
  actionable_mut[, actionable_genes]
)

# Mutation frequency by molecular state
actionable_mut_summary <- actionable_mut %>%
  dplyr::group_by(Robust_State) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(actionable_genes),
      ~ 100 * mean(.x, na.rm = TRUE)
    )
  )

View(actionable_mut_summary)

actionable_mut_summary

# ============================================================
# MUTATION ASSOCIATION WITH ROBUST STATE
# ============================================================

mutation_tests <- lapply(
  actionable_genes,
  function(gene) {
    
    tab <- table(
      State = actionable_mut$Robust_State,
      Mutated = actionable_mut[[gene]]
    )
    
    test <- fisher.test(tab)
    
    data.frame(
      Gene = gene,
      
      State1_mutated = sum(
        actionable_mut$Robust_State == 1 &
          actionable_mut[[gene]]
      ),
      
      State1_percent = 100 * mean(
        actionable_mut[[gene]][
          actionable_mut$Robust_State == 1
        ]
      ),
      
      State2_mutated = sum(
        actionable_mut$Robust_State == 2 &
          actionable_mut[[gene]]
      ),
      
      State2_percent = 100 * mean(
        actionable_mut[[gene]][
          actionable_mut$Robust_State == 2
        ]
      ),
      
      P_value = test$p.value
    )
  }
)

mutation_tests <- dplyr::bind_rows(
  mutation_tests
)

mutation_tests$FDR <- p.adjust(
  mutation_tests$P_value,
  method = "BH"
)

mutation_tests <- mutation_tests %>%
  dplyr::arrange(FDR)

mutation_tests

# ============================================================
# STRING-CONNECTED STATE SIGNATURE
# ============================================================

state1_string <- c(
  "GABRG3",
  "GABRE",
  "GABBR2",
  "DEFB1",
  "MMP11"
)

state2_string <- c(
  "ZIC4",
  "ZIC1",
  "FOXG1",
  "SIX3",
  "ZIC2",
  "ZIC5",
  "CD177",
  "LYPD3",
  "LIN28B",
  "MKRN3",
  "POU6F2",
  "PALM3",
  "LRRN1",
  "ADGRL3",
  "NLGN1",
  "ATP1A2",
  "MAGEA11",
  "CTAG2",
  "NFE2",
  "GATA2",
  "INA",
  "NEFH",
  "UCHL1",
  "SULT1C2",
  "SULT1C4",
  "UGT2B17",
  "MROH2A",
  "SLCO1A2",
  "SLC15A1",
  "TDRD12",
  "PIWIL1",
  "STK31"
)

string_signature <- unique(c(
  state1_string,
  state2_string
))

length(state1_string)
length(state2_string)
length(string_signature)
# ============================================================
# KEEP ONLY OUR 110 SIGNIFICANT GENES
# ============================================================

sig110 <- sub(
  "^RNA_",
  "",
  rna_110$Feature
)

state1_string_sig <- intersect(
  state1_string,
  sig110
)

state2_string_sig <- intersect(
  state2_string,
  sig110
)

string_signature_sig <- unique(c(
  state1_string_sig,
  state2_string_sig
))

cat(
  "State 1 connected + significant:",
  length(state1_string_sig),
  "\n"
)

print(state1_string_sig)

cat(
  "\nState 2 connected + significant:",
  length(state2_string_sig),
  "\n"
)

print(state2_string_sig)

cat(
  "\nTotal STRING-connected significant signature:",
  length(string_signature_sig),
  "\n"
)

print(string_signature_sig)
# STRING-connected significant genes
string_rna <- rna_mofa_new[
  rownames(rna_mofa_new) %in% string_signature_sig,
  ,
  drop = FALSE
]

dim(string_rna)
rownames(string_rna)

# Patient x gene
X_string <- t(string_rna)

# Patient states
y_string <- robust_state[
  match(
    rownames(X_string),
    rownames(factor_scores)
  )
]

table(y_string)

# Standardize genes
X_string_scaled <- scale(X_string)

# PCA
pca_string <- prcomp(
  X_string_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_df <- data.frame(
  PC1 = pca_string$x[, 1],
  PC2 = pca_string$x[, 2],
  State = factor(y_string)
)

library(ggplot2)

ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = State
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.75
  ) +
  stat_ellipse(
    aes(group = State),
    linewidth = 0.8
  ) +
  theme_classic() +
  labs(
    title = "STRING-connected RNA signature",
    subtitle = paste0(
      ncol(X_string),
      " network-connected state-associated genes"
    ),
    x = "PC1",
    y = "PC2",
    color = "Molecular State"
  )

library(glmnet)

# ============================================================
# DATA
# ============================================================

X <- as.matrix(X_string_scaled)

# State 1 = 0
# State 2 = 1
y <- ifelse(y_string == 2, 1, 0)

table(y)


# ============================================================
# CROSS-VALIDATED LASSO LOGISTIC REGRESSION
# ============================================================

set.seed(123)

cv_lasso <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "auc"
)

plot(cv_lasso)

cv_lasso$lambda.min
cv_lasso$lambda.1se
coef_1se <- coef(
  cv_lasso,
  s = "lambda.1se"
)

selected_genes_1se <- data.frame(
  Gene = rownames(coef_1se),
  Coefficient = as.numeric(coef_1se)
) %>%
  dplyr::filter(
    Coefficient != 0,
    Gene != "(Intercept)"
  ) %>%
  dplyr::arrange(
    dplyr::desc(abs(Coefficient))
  )

selected_genes_1se

nrow(selected_genes_1se)
# ============================================================
# lambda.min SELECTED GENES
# ============================================================

coef_min <- coef(
  cv_lasso,
  s = "lambda.min"
)

selected_genes <- data.frame(
  Gene = rownames(coef_min),
  Coefficient = as.numeric(coef_min)
) %>%
  dplyr::filter(
    Coefficient != 0,
    Gene != "(Intercept)"
  ) %>%
  dplyr::arrange(
    dplyr::desc(abs(Coefficient))
  )

selected_genes

nrow(selected_genes)

# ============================================================
# STRING-CONNECTED 37-GENE SIGNATURE
# OUT-OF-FOLD LASSO VALIDATION
# ============================================================

library(glmnet)
library(pROC)
library(dplyr)
library(ggplot2)


# ============================================================
# 1. PREPARE DATA
# ============================================================

# X_string_scaled:
# rows    = patients
# columns = 37 STRING-connected genes

X <- as.matrix(X_string_scaled)

# State 1 = 0
# State 2 = 1
y <- ifelse(y_string == 2, 1, 0)

# Checks
dim(X)
table(y)
colnames(X)


# ============================================================
# 2. CREATE STRATIFIED 10-FOLD SPLITS
# ============================================================

set.seed(123)

n <- nrow(X)

foldid <- rep(
  NA_integer_,
  n
)

# State 1
idx0 <- which(y == 0)

# State 2
idx1 <- which(y == 1)


# Distribute each state across 10 folds
foldid[idx0] <- sample(
  rep(
    1:10,
    length.out = length(idx0)
  )
)

foldid[idx1] <- sample(
  rep(
    1:10,
    length.out = length(idx1)
  )
)


# Check class distribution
table(
  Fold = foldid,
  State = y
)


# ============================================================
# 3. OUT-OF-FOLD PREDICTIONS
# ============================================================

oof_prob <- rep(
  NA_real_,
  n
)

# Store lambda selected in each fold
lambda_per_fold <- rep(
  NA_real_,
  10
)


set.seed(123)

for (k in 1:10) {
  
  cat(
    "Running fold",
    k,
    "/ 10\n"
  )
  
  # ----------------------------------------------------------
  # Training / testing
  # ----------------------------------------------------------
  
  train_idx <- which(
    foldid != k
  )
  
  test_idx <- which(
    foldid == k
  )
  
  
  # ----------------------------------------------------------
  # INNER CV
  # Lambda selection using TRAINING patients only
  # ----------------------------------------------------------
  
  cv_fit <- glmnet::cv.glmnet(
    
    x = X[
      train_idx,
      ,
      drop = FALSE
    ],
    
    y = y[
      train_idx
    ],
    
    family = "binomial",
    
    alpha = 1,
    
    nfolds = 9,
    
    type.measure = "auc"
  )
  
  
  # More conservative lambda
  lambda_use <- cv_fit$lambda.1se
  
  lambda_per_fold[k] <- lambda_use
  
  
  # ----------------------------------------------------------
  # FIT LASSO ON TRAINING DATA
  # ----------------------------------------------------------
  
  fit <- glmnet::glmnet(
    
    x = X[
      train_idx,
      ,
      drop = FALSE
    ],
    
    y = y[
      train_idx
    ],
    
    family = "binomial",
    
    alpha = 1,
    
    lambda = lambda_use
  )
  
  
  # ----------------------------------------------------------
  # PREDICT HELD-OUT PATIENTS
  #
  # Explicit glmnet method used to avoid predict() conflicts
  # ----------------------------------------------------------
  
  oof_prob[test_idx] <- as.numeric(
    
    glmnet:::predict.glmnet(
      
      object = fit,
      
      newx = X[
        test_idx,
        ,
        drop = FALSE
      ],
      
      s = lambda_use,
      
      type = "response"
    )
  )
}


# ============================================================
# 4. CHECK OOF PREDICTIONS
# ============================================================

summary(oof_prob)

sum(
  is.na(oof_prob)
)

lambda_per_fold


# IMPORTANT:
# Number of NA predictions should be 0.


# ============================================================
# 5. ROC ANALYSIS
# ============================================================

roc_string <- pROC::roc(
  
  response = y,
  
  predictor = oof_prob,
  
  levels = c(
    0,
    1
  ),
  
  direction = "<",
  
  quiet = TRUE
)


# ============================================================
# 6. AUC
# ============================================================

auc_string <- as.numeric(
  pROC::auc(
    roc_string
  )
)

auc_string


# ============================================================
# 7. 95% CI FOR AUC
# ============================================================

ci_auc <- pROC::ci.auc(
  roc_string
)

ci_auc


# ============================================================
# 8. FIND OPTIMAL THRESHOLD
# Youden index
# ============================================================

best_coords <- pROC::coords(
  
  roc_string,
  
  x = "best",
  
  best.method = "youden",
  
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  
  transpose = FALSE
)

best_coords


# ============================================================
# 9. EXTRACT THRESHOLD
# ============================================================

threshold <- as.numeric(
  best_coords["threshold", ]
)

threshold


# ============================================================
# 10. CLASSIFY PATIENTS
# ============================================================

pred_class <- ifelse(
  oof_prob >= threshold,
  1,
  0
)


# ============================================================
# 11. CONFUSION MATRIX
# ============================================================

confusion <- table(
  
  Actual = y,
  
  Predicted = pred_class
)

confusion


# ============================================================
# 12. PERFORMANCE METRICS
# ============================================================

TP <- sum(
  pred_class == 1 &
    y == 1
)

TN <- sum(
  pred_class == 0 &
    y == 0
)

FP <- sum(
  pred_class == 1 &
    y == 0
)

FN <- sum(
  pred_class == 0 &
    y == 1
)


# Sensitivity
sensitivity <- TP / (
  TP + FN
)


# Specificity
specificity <- TN / (
  TN + FP
)


# Accuracy
accuracy <- (
  TP + TN
) / length(y)


# Balanced accuracy
balanced_accuracy <- (
  sensitivity +
    specificity
) / 2


# PPV / Precision
precision <- TP / (
  TP + FP
)


# NPV
npv <- TN / (
  TN + FN
)


# F1 score
f1_score <- 2 * (
  precision *
    sensitivity
) / (
  precision +
    sensitivity
)


# ============================================================
# 13. FINAL PERFORMANCE TABLE
# ============================================================

performance <- data.frame(
  
  Model =
    "37-gene STRING-LASSO",
  
  N_patients =
    length(y),
  
  N_State1 =
    sum(y == 0),
  
  N_State2 =
    sum(y == 1),
  
  N_genes =
    ncol(X),
  
  AUC =
    auc_string,
  
  AUC_low =
    as.numeric(
      ci_auc[1]
    ),
  
  AUC_high =
    as.numeric(
      ci_auc[3]
    ),
  
  Sensitivity =
    sensitivity,
  
  Specificity =
    specificity,
  
  Balanced_Accuracy =
    balanced_accuracy,
  
  Accuracy =
    accuracy,
  
  Precision =
    precision,
  
  NPV =
    npv,
  
  F1 =
    f1_score,
  
  Threshold =
    threshold
)

performance


# ============================================================
# 14. ROC CURVE
# ============================================================

plot(
  
  roc_string,
  
  legacy.axes = TRUE,
  
  print.auc = TRUE,
  
  print.auc.cex = 1.2,
  
  lwd = 2,
  
  main =
    "OOF ROC: STRING-connected RNA signature"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)


# ============================================================
# 15. SAVE PATIENT-LEVEL OOF RESULTS
# ============================================================

oof_results <- data.frame(
  
  Patient =
    rownames(X),
  
  Actual_State =
    ifelse(
      y == 1,
      "State2",
      "State1"
    ),
  
  P_State2 =
    oof_prob,
  
  Predicted_State =
    ifelse(
      pred_class == 1,
      "State2",
      "State1"
    ),
  
  Correct =
    pred_class == y,
  
  Fold =
    foldid
)

View(oof_results)


# ============================================================
# 16. OPTIONAL: SAVE RESULTS
# ============================================================

write.csv(
  
  performance,
  
  "STRING_37gene_LASSO_performance.csv",
  
  row.names = FALSE
)

write.csv(
  
  oof_results,
  
  "STRING_37gene_LASSO_OOF_predictions.csv",
  
  row.names = FALSE
)


# ============================================================
# 17. PRINT MAIN RESULTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "37-GENE STRING-LASSO RESULTS\n"
)

cat(
  "========================================\n"
)

cat(
  "AUC:",
  round(
    auc_string,
    3
  ),
  "\n"
)

cat(
  "95% CI:",
  round(
    ci_auc[1],
    3
  ),
  "-",
  round(
    ci_auc[3],
    3
  ),
  "\n"
)

cat(
  "Sensitivity:",
  round(
    sensitivity,
    3
  ),
  "\n"
)

cat(
  "Specificity:",
  round(
    specificity,
    3
  ),
  "\n"
)

cat(
  "Balanced Accuracy:",
  round(
    balanced_accuracy,
    3
  ),
  "\n"
)

cat(
  "Accuracy:",
  round(
    accuracy,
    3
  ),
  "\n"
)

cat(
  "F1:",
  round(
    f1_score,
    3
  ),
  "\n"
)

cat(
  "Threshold:",
  round(
    threshold,
    3
  ),
  "\n"
)

cat(
  "========================================\n"
)

performance

confusion
# ============================================================
# FIX THRESHOLD + PERFORMANCE METRICS
# ============================================================

# Optimal threshold
best_coords <- pROC::coords(
  roc_string,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

best_coords


# ============================================================
# CORRECTLY EXTRACT VALUES
# ============================================================

threshold <- as.numeric(best_coords$threshold)

roc_sensitivity <- as.numeric(best_coords$sensitivity)
roc_specificity <- as.numeric(best_coords$specificity)

threshold
roc_sensitivity
roc_specificity


# ============================================================
# CLASSIFY OOF PREDICTIONS
# ============================================================

pred_class <- ifelse(
  oof_prob >= threshold,
  1,
  0
)


# ============================================================
# CONFUSION MATRIX
# ============================================================

confusion <- table(
  Actual = factor(y, levels = c(0, 1)),
  Predicted = factor(pred_class, levels = c(0, 1))
)

confusion


# ============================================================
# EXTRACT COUNTS
# ============================================================

TN <- confusion["0", "0"]
FP <- confusion["0", "1"]

FN <- confusion["1", "0"]
TP <- confusion["1", "1"]


# ============================================================
# PERFORMANCE METRICS
# ============================================================

sensitivity <- TP / (TP + FN)

specificity <- TN / (TN + FP)

accuracy <- (TP + TN) / sum(confusion)

balanced_accuracy <- (
  sensitivity + specificity
) / 2

precision <- ifelse(
  (TP + FP) > 0,
  TP / (TP + FP),
  NA
)

npv <- ifelse(
  (TN + FN) > 0,
  TN / (TN + FN),
  NA
)

f1_score <- ifelse(
  !is.na(precision) &&
    (precision + sensitivity) > 0,
  2 * precision * sensitivity /
    (precision + sensitivity),
  NA
)


# ============================================================
# FINAL TABLE
# ============================================================

performance_fixed <- data.frame(
  
  Model = "37-gene STRING-LASSO",
  
  N_patients = length(y),
  
  N_State1 = sum(y == 0),
  
  N_State2 = sum(y == 1),
  
  N_genes = ncol(X),
  
  AUC = auc_string,
  
  AUC_low = as.numeric(ci_auc[1]),
  
  AUC_high = as.numeric(ci_auc[3]),
  
  Threshold = threshold,
  
  Sensitivity = as.numeric(sensitivity),
  
  Specificity = as.numeric(specificity),
  
  Balanced_Accuracy =
    as.numeric(balanced_accuracy),
  
  Accuracy =
    as.numeric(accuracy),
  
  Precision =
    as.numeric(precision),
  
  NPV =
    as.numeric(npv),
  
  F1 =
    as.numeric(f1_score)
)


performance_fixed

confusion
# ============================================================
# STABILITY OF GENE SELECTION ACROSS THE 10 OUTER FOLDS
# ============================================================

gene_selection <- matrix(
  0,
  nrow = ncol(X),
  ncol = 10,
  dimnames = list(
    colnames(X),
    paste0("Fold", 1:10)
  )
)

gene_coef <- matrix(
  0,
  nrow = ncol(X),
  ncol = 10,
  dimnames = list(
    colnames(X),
    paste0("Fold", 1:10)
  )
)

set.seed(123)

for (k in 1:10) {
  
  train_idx <- which(foldid != k)
  
  # ------------------------------------------
  # Inner CV: choose lambda using training only
  # ------------------------------------------
  
  cv_fit <- glmnet::cv.glmnet(
    x = X[train_idx, , drop = FALSE],
    y = y[train_idx],
    family = "binomial",
    alpha = 1,
    nfolds = 9,
    type.measure = "auc"
  )
  
  lambda_use <- cv_fit$lambda.1se
  
  # ------------------------------------------
  # Fit LASSO
  # ------------------------------------------
  
  fit <- glmnet::glmnet(
    x = X[train_idx, , drop = FALSE],
    y = y[train_idx],
    family = "binomial",
    alpha = 1,
    lambda = lambda_use
  )
  
  # ------------------------------------------
  # Extract coefficients
  # ------------------------------------------
  
  cf <- as.matrix(
    coef(
      fit,
      s = lambda_use
    )
  )
  
  cf <- cf[
    rownames(cf) != "(Intercept)",
    ,
    drop = FALSE
  ]
  
  # Put coefficients in correct gene order
  cf_vec <- cf[, 1]
  
  gene_coef[
    names(cf_vec),
    k
  ] <- cf_vec
  
  gene_selection[
    names(cf_vec),
    k
  ] <- as.numeric(
    cf_vec != 0
  )
}


# ============================================================
# SUMMARY
# ============================================================

selection_frequency <- rowSums(
  gene_selection
)

selection_percent <- 100 *
  selection_frequency / 10

mean_coefficient <- rowMeans(
  gene_coef
)

mean_abs_coefficient <- rowMeans(
  abs(gene_coef)
)


stable_genes <- data.frame(
  
  Gene = rownames(gene_selection),
  
  Selected_Folds =
    selection_frequency,
  
  Selection_Percent =
    selection_percent,
  
  Mean_Coefficient =
    mean_coefficient,
  
  Mean_Abs_Coefficient =
    mean_abs_coefficient
  
) %>%
  
  dplyr::arrange(
    dplyr::desc(Selected_Folds),
    dplyr::desc(Mean_Abs_Coefficient)
  )


stable_genes
# Genes selected in >= 80% of outer folds
stable_core_80 <- stable_genes %>%
  dplyr::filter(
    Selection_Percent >= 80
  )

stable_core_80

nrow(stable_core_80)
# Genes selected in >= 90% of folds
stable_core_90 <- stable_genes %>%
  dplyr::filter(
    Selection_Percent >= 90
  )

stable_core_90

nrow(stable_core_90)
# ============================================================
# 10-GENE STABLE CORE SIGNATURE
# OOF LOGISTIC REGRESSION
# ============================================================

library(glmnet)
library(pROC)

# ------------------------------------------------------------
# 1. Stable 10 genes
# ------------------------------------------------------------

core10 <- stable_core_80$Gene

core10


# ------------------------------------------------------------
# 2. Extract 10 genes from existing standardized matrix
# ------------------------------------------------------------

X10 <- X[, core10, drop = FALSE]

dim(X10)
colnames(X10)


# ============================================================
# 3. OOF PREDICTIONS
#
# IMPORTANT:
# Use the SAME outer folds as the 37-gene analysis
# ============================================================

oof_prob_10 <- rep(NA_real_, nrow(X10))


for (k in 1:10) {
  
  train_idx <- which(foldid != k)
  test_idx  <- which(foldid == k)
  
  # ----------------------------------------------------------
  # Logistic regression with small ridge penalty
  # ----------------------------------------------------------
  
  cv_fit10 <- glmnet::cv.glmnet(
    x = X10[train_idx, , drop = FALSE],
    y = y[train_idx],
    family = "binomial",
    alpha = 0,                 # ridge
    nfolds = 9,
    type.measure = "auc"
  )
  
  lambda_use <- cv_fit10$lambda.1se
  
  fit10 <- glmnet::glmnet(
    x = X10[train_idx, , drop = FALSE],
    y = y[train_idx],
    family = "binomial",
    alpha = 0,
    lambda = lambda_use
  )
  
  oof_prob_10[test_idx] <- as.numeric(
    glmnet:::predict.glmnet(
      object = fit10,
      newx = X10[test_idx, , drop = FALSE],
      s = lambda_use,
      type = "response"
    )
  )
}


# ============================================================
# 4. CHECK
# ============================================================

summary(oof_prob_10)

sum(is.na(oof_prob_10))


# ============================================================
# 5. ROC
# ============================================================

roc10 <- pROC::roc(
  response = y,
  predictor = oof_prob_10,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc10 <- as.numeric(
  pROC::auc(roc10)
)

ci10 <- pROC::ci.auc(roc10)

auc10
ci10


# ============================================================
# 6. OPTIMAL THRESHOLD
# ============================================================

coords10 <- pROC::coords(
  roc10,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

threshold10 <- as.numeric(
  coords10$threshold
)

threshold10


# ============================================================
# 7. PREDICTED CLASSES
# ============================================================

pred10 <- ifelse(
  oof_prob_10 >= threshold10,
  1,
  0
)


confusion10 <- table(
  Actual = factor(y, levels = c(0, 1)),
  Predicted = factor(pred10, levels = c(0, 1))
)

confusion10


# ============================================================
# 8. PERFORMANCE
# ============================================================

TN10 <- confusion10["0", "0"]
FP10 <- confusion10["0", "1"]

FN10 <- confusion10["1", "0"]
TP10 <- confusion10["1", "1"]


sens10 <- TP10 / (TP10 + FN10)

spec10 <- TN10 / (TN10 + FP10)

bal10 <- (
  sens10 + spec10
) / 2

acc10 <- (
  TP10 + TN10
) / sum(confusion10)

precision10 <- TP10 / (
  TP10 + FP10
)

npv10 <- TN10 / (
  TN10 + FN10
)

f1_10 <- 2 *
  precision10 *
  sens10 /
  (
    precision10 +
      sens10
  )


performance10 <- data.frame(
  
  Model = "10-gene stable core",
  
  N_genes = 10,
  
  AUC = auc10,
  
  AUC_low = as.numeric(ci10[1]),
  
  AUC_high = as.numeric(ci10[3]),
  
  Sensitivity = as.numeric(sens10),
  
  Specificity = as.numeric(spec10),
  
  Balanced_Accuracy = as.numeric(bal10),
  
  Accuracy = as.numeric(acc10),
  
  Precision = as.numeric(precision10),
  
  NPV = as.numeric(npv10),
  
  F1 = as.numeric(f1_10)
)

performance10


# ============================================================
# 9. COMPARE 37 vs 10 GENES
# ============================================================

comparison <- data.frame(
  
  Model = c(
    "37-gene STRING",
    "10-gene stable core"
  ),
  
  N_genes = c(
    37,
    10
  ),
  
  AUC = c(
    auc_string,
    auc10
  ),
  
  Sensitivity = c(
    sensitivity,
    sens10
  ),
  
  Specificity = c(
    specificity,
    spec10
  ),
  
  Balanced_Accuracy = c(
    balanced_accuracy,
    bal10
  )
)

comparison


# ============================================================
# 10. STATISTICAL COMPARISON OF ROC CURVES
# ============================================================

pROC::roc.test(
  roc_string,
  roc10,
  paired = TRUE,
  method = "delong"
)


# ============================================================
# 11. ROC FIGURE
# ============================================================

plot(
  roc_string,
  legacy.axes = TRUE,
  lwd = 2,
  main = "Molecular State RNA Classifiers"
)

plot(
  roc10,
  add = TRUE,
  lwd = 2,
  lty = 2
)

abline(
  a = 0,
  b = 1,
  lty = 3
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "37 genes: AUC = ",
      round(auc_string, 3)
    ),
    paste0(
      "10-gene core: AUC = ",
      round(auc10, 3)
    )
  ),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)
# ============================================================
# FULL NESTED CROSS-VALIDATION
#
# Training fold only:
# RNA differential analysis
#        ↓
# FDR < 0.05 genes
#        ↓
# STRING-connected filter (predefined 37 genes)
#        ↓
# LASSO + inner CV
#        ↓
# Held-out patients
#
# Outcome:
# State 1 = 0
# State 2 = 1
# ============================================================


library(glmnet)
library(pROC)
library(dplyr)


# ============================================================
# 1. PREPARE FULL RNA DATA
# ============================================================

# rna_mofa_new:
# rows    = RNA features
# columns = patients

rna_nested <- rna_mofa_new

# Remove RNA_ prefix
rownames(rna_nested) <- sub(
  "^RNA_",
  "",
  rownames(rna_nested)
)


# Patient x gene matrix
X_full <- t(rna_nested)


# Align patients with factor scores
common_patients <- intersect(
  rownames(X_full),
  rownames(factor_scores)
)

X_full <- X_full[
  common_patients,
  ,
  drop = FALSE
]


# State labels
state_nested <- robust_state[
  match(
    common_patients,
    rownames(factor_scores)
  )
]


# Binary outcome
# State1 = 0
# State2 = 1

y_nested <- ifelse(
  state_nested == 2,
  1,
  0
)


# Checks
dim(X_full)

table(y_nested)

stopifnot(
  nrow(X_full) ==
    length(y_nested)
)


# ============================================================
# 2. PREDEFINED STRING NETWORK GENE SET
# ============================================================

string_network_genes <- c(
  
  "GABRG3",
  "GABRE",
  "GABBR2",
  "DEFB1",
  "MMP11",
  
  "ZIC4",
  "ZIC1",
  "FOXG1",
  "SIX3",
  "ZIC2",
  "ZIC5",
  
  "CD177",
  "LYPD3",
  
  "LIN28B",
  "MKRN3",
  
  "POU6F2",
  "PALM3",
  
  "LRRN1",
  "ADGRL3",
  "NLGN1",
  "ATP1A2",
  
  "MAGEA11",
  "CTAG2",
  
  "NFE2",
  "GATA2",
  
  "INA",
  "NEFH",
  "UCHL1",
  
  "SULT1C2",
  "SULT1C4",
  "UGT2B17",
  "MROH2A",
  "SLCO1A2",
  "SLC15A1",
  
  "TDRD12",
  "PIWIL1",
  "STK31"
)


# Keep only genes actually available
string_network_genes <- intersect(
  string_network_genes,
  colnames(X_full)
)

length(string_network_genes)

string_network_genes


# ============================================================
# 3. STRATIFIED OUTER 10-FOLD CV
# ============================================================

set.seed(123)

n <- nrow(X_full)

outer_fold <- rep(
  NA_integer_,
  n
)


idx_state1 <- which(
  y_nested == 0
)

idx_state2 <- which(
  y_nested == 1
)


outer_fold[idx_state1] <- sample(
  rep(
    1:10,
    length.out = length(idx_state1)
  )
)

outer_fold[idx_state2] <- sample(
  rep(
    1:10,
    length.out = length(idx_state2)
  )
)


# Check distribution
table(
  Fold = outer_fold,
  State = y_nested
)


# ============================================================
# 4. STORAGE OBJECTS
# ============================================================

nested_score <- rep(
  NA_real_,
  n
)

nested_pred_prob <- rep(
  NA_real_,
  n
)


# Store which genes were significant / selected
DE_gene_list <- vector(
  "list",
  10
)

STRING_gene_list <- vector(
  "list",
  10
)

LASSO_gene_list <- vector(
  "list",
  10
)

lambda_list <- rep(
  NA_real_,
  10
)


# ============================================================
# 5. OUTER NESTED CV LOOP
# ============================================================

set.seed(123)

for (k in 1:10) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "OUTER FOLD:",
    k,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  # ----------------------------------------------------------
  # TRAIN / TEST
  # ----------------------------------------------------------
  
  train_idx <- which(
    outer_fold != k
  )
  
  test_idx <- which(
    outer_fold == k
  )
  
  
  X_train_full <- X_full[
    train_idx,
    ,
    drop = FALSE
  ]
  
  X_test_full <- X_full[
    test_idx,
    ,
    drop = FALSE
  ]
  
  
  y_train <- y_nested[
    train_idx
  ]
  
  y_test <- y_nested[
    test_idx
  ]
  
  
  # ==========================================================
  # 5A. DIFFERENTIAL RNA ANALYSIS
  #
  # IMPORTANT:
  # TRAINING PATIENTS ONLY
  # ==========================================================
  
  p_values <- rep(
    NA_real_,
    ncol(X_train_full)
  )
  
  differences <- rep(
    NA_real_,
    ncol(X_train_full)
  )
  
  
  names(p_values) <- colnames(
    X_train_full
  )
  
  names(differences) <- colnames(
    X_train_full
  )
  
  
  for (j in seq_len(
    ncol(X_train_full)
  )) {
    
    x1 <- X_train_full[
      y_train == 0,
      j
    ]
    
    x2 <- X_train_full[
      y_train == 1,
      j
    ]
    
    
    differences[j] <-
      median(
        x2,
        na.rm = TRUE
      ) -
      median(
        x1,
        na.rm = TRUE
      )
    
    
    p_values[j] <- tryCatch(
      
      wilcox.test(
        x1,
        x2,
        exact = FALSE
      )$p.value,
      
      error = function(e) NA_real_
    )
  }
  
  
  # ----------------------------------------------------------
  # FDR
  # ----------------------------------------------------------
  
  fdr_values <- p.adjust(
    p_values,
    method = "BH"
  )
  
  
  de_table <- data.frame(
    
    Gene =
      names(p_values),
    
    Difference =
      differences,
    
    P_value =
      p_values,
    
    FDR =
      fdr_values
  )
  
  
  de_sig <- de_table %>%
    
    dplyr::filter(
      !is.na(FDR),
      FDR < 0.05
    )
  
  
  DE_gene_list[[k]] <-
    de_sig$Gene
  
  
  cat(
    "Training DE genes:",
    nrow(de_sig),
    "\n"
  )
  
  
  # ==========================================================
  # 5B. STRING NETWORK FILTER
  # ==========================================================
  
  candidate_genes <- intersect(
    
    de_sig$Gene,
    
    string_network_genes
  )
  
  
  STRING_gene_list[[k]] <-
    candidate_genes
  
  
  cat(
    "STRING-connected significant genes:",
    length(candidate_genes),
    "\n"
  )
  
  
  # Safety check
  if (
    length(candidate_genes) < 2
  ) {
    
    stop(
      paste(
        "Too few candidate genes in fold",
        k
      )
    )
  }
  
  
  # ==========================================================
  # 5C. EXTRACT TRAIN / TEST MATRICES
  # ==========================================================
  
  X_train <- X_train_full[
    ,
    candidate_genes,
    drop = FALSE
  ]
  
  X_test <- X_test_full[
    ,
    candidate_genes,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 5D. STANDARDIZATION
  #
  # CRITICAL:
  # Calculate mean/SD from TRAINING ONLY
  # ==========================================================
  
  train_mean <- apply(
    X_train,
    2,
    mean,
    na.rm = TRUE
  )
  
  train_sd <- apply(
    X_train,
    2,
    sd,
    na.rm = TRUE
  )
  
  
  # Remove zero-variance genes
  valid_genes <- names(
    train_sd[
      is.finite(train_sd) &
        train_sd > 0
    ]
  )
  
  
  X_train <- X_train[
    ,
    valid_genes,
    drop = FALSE
  ]
  
  X_test <- X_test[
    ,
    valid_genes,
    drop = FALSE
  ]
  
  
  train_mean <- train_mean[
    valid_genes
  ]
  
  train_sd <- train_sd[
    valid_genes
  ]
  
  
  # Standardize TRAIN
  X_train_scaled <- sweep(
    X_train,
    2,
    train_mean,
    "-"
  )
  
  X_train_scaled <- sweep(
    X_train_scaled,
    2,
    train_sd,
    "/"
  )
  
  
  # Standardize TEST using TRAIN parameters
  X_test_scaled <- sweep(
    X_test,
    2,
    train_mean,
    "-"
  )
  
  X_test_scaled <- sweep(
    X_test_scaled,
    2,
    train_sd,
    "/"
  )
  
  
  # ==========================================================
  # 5E. INNER CV FOR LAMBDA
  # ==========================================================
  
  inner_cv <- glmnet::cv.glmnet(
    
    x = as.matrix(
      X_train_scaled
    ),
    
    y = y_train,
    
    family = "binomial",
    
    alpha = 1,
    
    nfolds = 9,
    
    type.measure = "auc",
    
    standardize = FALSE
  )
  
  
  lambda_use <-
    inner_cv$lambda.1se
  
  
  lambda_list[k] <-
    lambda_use
  
  
  # ==========================================================
  # 5F. FINAL LASSO MODEL
  # ==========================================================
  
  fit <- glmnet::glmnet(
    
    x = as.matrix(
      X_train_scaled
    ),
    
    y = y_train,
    
    family = "binomial",
    
    alpha = 1,
    
    lambda =
      lambda_use,
    
    standardize = FALSE
  )
  
  
  # ==========================================================
  # 5G. STORE SELECTED GENES
  # ==========================================================
  
  cf <- as.matrix(
    coef(
      fit,
      s = lambda_use
    )
  )
  
  
  selected <- rownames(cf)[
    cf[, 1] != 0
  ]
  
  
  selected <- setdiff(
    selected,
    "(Intercept)"
  )
  
  
  LASSO_gene_list[[k]] <-
    selected
  
  
  cat(
    "LASSO-selected genes:",
    length(selected),
    "\n"
  )
  
  
  cat(
    paste(
      selected,
      collapse = ", "
    ),
    "\n"
  )
  
  
  # ==========================================================
  # 5H. PREDICT HELD-OUT PATIENTS
  #
  # We obtain LINK scores first.
  # Then manually transform to probabilities.
  # This avoids the earlier predict.glmnet issue.
  # ==========================================================
  
  test_link <- as.numeric(
    
    glmnet:::predict.glmnet(
      
      object = fit,
      
      newx = as.matrix(
        X_test_scaled
      ),
      
      s = lambda_use,
      
      type = "link"
    )
  )
  
  
  # Store linear predictor
  nested_score[
    test_idx
  ] <- test_link
  
  
  # Convert logit -> probability
  nested_pred_prob[
    test_idx
  ] <- plogis(
    test_link
  )
  
  
  cat(
    "Test probability range:",
    round(
      min(
        nested_pred_prob[test_idx]
      ),
      3
    ),
    "-",
    round(
      max(
        nested_pred_prob[test_idx]
      ),
      3
    ),
    "\n"
  )
}


# ============================================================
# 6. CHECK NESTED PREDICTIONS
# ============================================================

summary(
  nested_pred_prob
)

range(
  nested_pred_prob
)

sum(
  is.na(
    nested_pred_prob
  )
)


# Should be:
# probability range between 0 and 1
# NA = 0


# ============================================================
# 7. NESTED ROC / AUC
# ============================================================

nested_roc <- pROC::roc(
  
  response =
    y_nested,
  
  predictor =
    nested_pred_prob,
  
  levels = c(
    0,
    1
  ),
  
  direction = "<",
  
  quiet = TRUE
)


nested_auc <- as.numeric(
  pROC::auc(
    nested_roc
  )
)


nested_ci <- pROC::ci.auc(
  nested_roc
)


nested_auc

nested_ci


# ============================================================
# 8. OPTIMAL THRESHOLD
# ============================================================

nested_coords <- pROC::coords(
  
  nested_roc,
  
  x = "best",
  
  best.method = "youden",
  
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  
  transpose = FALSE
)


nested_coords


nested_threshold <-
  as.numeric(
    nested_coords$threshold
  )


# ============================================================
# 9. CLASSIFICATION
# ============================================================

nested_pred_class <- ifelse(
  
  nested_pred_prob >=
    nested_threshold,
  
  1,
  
  0
)


nested_confusion <- table(
  
  Actual = factor(
    y_nested,
    levels = c(0, 1)
  ),
  
  Predicted = factor(
    nested_pred_class,
    levels = c(0, 1)
  )
)


nested_confusion


# ============================================================
# 10. PERFORMANCE METRICS
# ============================================================

TN <- nested_confusion[
  "0",
  "0"
]

FP <- nested_confusion[
  "0",
  "1"
]

FN <- nested_confusion[
  "1",
  "0"
]

TP <- nested_confusion[
  "1",
  "1"
]


nested_sensitivity <-
  TP / (
    TP + FN
  )


nested_specificity <-
  TN / (
    TN + FP
  )


nested_balanced_accuracy <-
  (
    nested_sensitivity +
      nested_specificity
  ) / 2


nested_accuracy <-
  (
    TP + TN
  ) /
  sum(
    nested_confusion
  )


nested_precision <-
  TP / (
    TP + FP
  )


nested_npv <-
  TN / (
    TN + FN
  )


nested_f1 <-
  2 *
  nested_precision *
  nested_sensitivity /
  (
    nested_precision +
      nested_sensitivity
  )


# ============================================================
# 11. FINAL PERFORMANCE TABLE
# ============================================================

nested_performance <- data.frame(
  
  Model =
    "Nested CV STRING-LASSO",
  
  N_patients =
    length(
      y_nested
    ),
  
  N_State1 =
    sum(
      y_nested == 0
    ),
  
  N_State2 =
    sum(
      y_nested == 1
    ),
  
  AUC =
    nested_auc,
  
  AUC_low =
    as.numeric(
      nested_ci[1]
    ),
  
  AUC_high =
    as.numeric(
      nested_ci[3]
    ),
  
  Threshold =
    nested_threshold,
  
  Sensitivity =
    as.numeric(
      nested_sensitivity
    ),
  
  Specificity =
    as.numeric(
      nested_specificity
    ),
  
  Balanced_Accuracy =
    as.numeric(
      nested_balanced_accuracy
    ),
  
  Accuracy =
    as.numeric(
      nested_accuracy
    ),
  
  Precision =
    as.numeric(
      nested_precision
    ),
  
  NPV =
    as.numeric(
      nested_npv
    ),
  
  F1 =
    as.numeric(
      nested_f1
    )
)


nested_performance


# ============================================================
# 12. GENE SELECTION STABILITY ACROSS OUTER FOLDS
# ============================================================

all_nested_genes <- sort(
  unique(
    unlist(
      LASSO_gene_list
    )
  )
)


nested_gene_stability <- data.frame(
  
  Gene =
    all_nested_genes,
  
  Selected_Folds =
    sapply(
      all_nested_genes,
      function(g) {
        
        sum(
          sapply(
            LASSO_gene_list,
            function(x)
              g %in% x
          )
        )
      }
    )
)


nested_gene_stability <-
  nested_gene_stability %>%
  
  dplyr::mutate(
    
    Selection_Percent =
      100 *
      Selected_Folds /
      10
    
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Selected_Folds
    )
  )


nested_gene_stability


# ============================================================
# 13. GENES AVAILABLE AFTER DE + STRING FILTER
# ============================================================

all_string_genes <- sort(
  unique(
    unlist(
      STRING_gene_list
    )
  )
)


string_filter_stability <- data.frame(
  
  Gene =
    all_string_genes,
  
  Significant_STRING_Folds =
    sapply(
      all_string_genes,
      function(g) {
        
        sum(
          sapply(
            STRING_gene_list,
            function(x)
              g %in% x
          )
        )
      }
    )
)


string_filter_stability <-
  string_filter_stability %>%
  
  dplyr::mutate(
    
    Percent =
      100 *
      Significant_STRING_Folds /
      10
    
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Significant_STRING_Folds
    )
  )


string_filter_stability


# ============================================================
# 14. ROC FIGURE
# ============================================================

plot(
  
  nested_roc,
  
  legacy.axes = TRUE,
  
  lwd = 2,
  
  main =
    "Nested CV: Network-guided RNA classifier",
  
  print.auc = TRUE,
  
  print.auc.cex = 1.2
)


abline(
  
  a = 0,
  
  b = 1,
  
  lty = 2
)


# ============================================================
# 15. PATIENT-LEVEL RESULTS
# ============================================================

nested_patient_results <- data.frame(
  
  Patient =
    rownames(
      X_full
    ),
  
  Actual_State =
    ifelse(
      y_nested == 1,
      "State2",
      "State1"
    ),
  
  P_State2 =
    nested_pred_prob,
  
  Predicted_State =
    ifelse(
      nested_pred_class == 1,
      "State2",
      "State1"
    ),
  
  Correct =
    nested_pred_class ==
    y_nested,
  
  Outer_Fold =
    outer_fold
)


View(
  nested_patient_results
)


# ============================================================
# 16. SAVE RESULTS
# ============================================================

write.csv(
  
  nested_performance,
  
  "NestedCV_STRING_LASSO_performance.csv",
  
  row.names = FALSE
)


write.csv(
  
  nested_gene_stability,
  
  "NestedCV_gene_stability.csv",
  
  row.names = FALSE
)


write.csv(
  
  string_filter_stability,
  
  "NestedCV_STRING_filter_stability.csv",
  
  row.names = FALSE
)


write.csv(
  
  nested_patient_results,
  
  "NestedCV_patient_predictions.csv",
  
  row.names = FALSE
)


# ============================================================
# 17. PRINT IMPORTANT RESULTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "FULL NESTED CV RESULTS\n"
)

cat(
  "========================================\n"
)

print(
  nested_performance
)


cat(
  "\nCONFUSION MATRIX\n"
)

print(
  nested_confusion
)


cat(
  "\nLASSO GENE STABILITY\n"
)

print(
  nested_gene_stability
)


cat(
  "\nSTRING + DE STABILITY\n"
)

print(
  string_filter_stability
)

cat(
  "========================================\n"
)
# ============================================================
# ROBUST MOLECULAR STATE → RNA SIGNATURE WORKFLOW
# ============================================================

library(ggplot2)
library(grid)


# ============================================================
# 1. NODES
# ============================================================

nodes <- data.frame(
  
  id = c(
    "TCGA",
    "MOFA",
    "K5",
    "STABILITY",
    "K2",
    "RNA110",
    "STRING37",
    "NESTED",
    "GENE9",
    "CORE4"
  ),
  
  x = rep(5, 10),
  
  y = c(
    10,
    8.8,
    7.6,
    6.4,
    5.2,
    4.0,
    2.8,
    1.6,
    0.4,
    -0.8
  ),
  
  label = c(
    
    "TCGA-OV PRIMARY TUMORS\n277 patients\nMatched RNA + CNV + mutation data",
    
    "MOFA2 MULTI-OMICS INTEGRATION\n1,722 RNA + 1,902 CNV + 107 mutation features\n20 latent factors → Factors 1–10 retained",
    
    "EXPLORATORY 5-STATE SOLUTION\nState 1–5\nApparently distinct molecular substructure",
    
    "ROBUSTNESS TESTING\nFactor-number sensitivity + bootstrap validation\n5-state solution unstable\nMedian bootstrap ARI ≈ 0.35",
    
    "ROBUST 2-STATE SOLUTION\nState 1: n = 226   |   State 2: n = 51\nK=2 bootstrap median ARI ≈ 0.80",
    
    "STATE-ASSOCIATED RNA SIGNAL\n110 genes at FDR < 0.05\n34 State 1-high   |   76 State 2-high",
    
    "STRING NETWORK FILTERING\n37 state-associated genes connected\nwithin protein-interaction networks",
    
    "NETWORK-GUIDED NESTED CV LASSO\nOuter-fold DE selection + STRING filtering + inner LASSO\nAUC = 0.763   (95% CI: 0.688–0.838)\nBalanced accuracy = 0.729",
    
    "STABLE 9-GENE CORE\nSelected in ≥80% of outer folds\nADGRL3  •  PIWIL1  •  SLC15A1  •  ZIC5  •  GATA2\nFOXG1  •  MAGEA11  •  PALM3  •  SLCO1A2",
    
    "HIGHEST-STABILITY CORE\nSelected in 10/10 outer folds\nADGRL3  •  PIWIL1  •  SLC15A1  •  ZIC5"
  ),
  
  group = c(
    "input",
    "mofa",
    "exploratory",
    "validation",
    "robust",
    "rna",
    "string",
    "model",
    "signature",
    "core"
  )
)


# ============================================================
# 2. COLORS
# ============================================================

node_colors <- c(
  
  input       = "#DCEAF7",
  
  mofa        = "#D9EAD3",
  
  exploratory = "#F8DEDE",
  
  validation  = "#FFF0D6",
  
  robust      = "#CDECCF",
  
  rna         = "#E8E1F2",
  
  string      = "#E5DDF0",
  
  model       = "#DCE6F2",
  
  signature   = "#C9E6D4",
  
  core        = "#A9D5BC"
)


# ============================================================
# 3. EDGES
# ============================================================

edges <- data.frame(
  
  from = c(
    "TCGA",
    "MOFA",
    "K5",
    "STABILITY",
    "K2",
    "RNA110",
    "STRING37",
    "NESTED",
    "GENE9"
  ),
  
  to = c(
    "MOFA",
    "K5",
    "STABILITY",
    "K2",
    "RNA110",
    "STRING37",
    "NESTED",
    "GENE9",
    "CORE4"
  )
)


# ============================================================
# 4. ADD COORDINATES
# ============================================================

edges <- merge(
  edges,
  nodes[, c("id", "x", "y")],
  by.x = "from",
  by.y = "id"
)

names(edges)[
  names(edges) == "x"
] <- "x_from"

names(edges)[
  names(edges) == "y"
] <- "y_from"


edges <- merge(
  edges,
  nodes[, c("id", "x", "y")],
  by.x = "to",
  by.y = "id"
)

names(edges)[
  names(edges) == "x"
] <- "x_to"

names(edges)[
  names(edges) == "y"
] <- "y_to"


# ============================================================
# 5. DRAW WORKFLOW
# ============================================================

signature_workflow <- ggplot() +
  
  # ----------------------------------------------------------
# Arrows
# ----------------------------------------------------------

geom_segment(
  
  data = edges,
  
  aes(
    x = x_from,
    y = y_from - 0.38,
    xend = x_to,
    yend = y_to + 0.38
  ),
  
  linewidth = 0.55,
  
  colour = "#667785",
  
  arrow = arrow(
    length = unit(
      0.15,
      "cm"
    ),
    type = "closed"
  )
) +
  
  
  # ----------------------------------------------------------
# Boxes
# ----------------------------------------------------------

geom_label(
  
  data = nodes,
  
  aes(
    x = x,
    y = y,
    label = label,
    fill = group
  ),
  
  size = 3.5,
  
  colour = "#263238",
  
  label.size = 0.4,
  
  label.padding = unit(
    0.32,
    "lines"
  ),
  
  label.r = unit(
    0.15,
    "lines"
  ),
  
  lineheight = 1.05
) +
  
  
  # ----------------------------------------------------------
# Colors
# ----------------------------------------------------------

scale_fill_manual(
  values = node_colors
) +
  
  
  # ----------------------------------------------------------
# Layout
# ----------------------------------------------------------

coord_cartesian(
  
  xlim = c(
    1,
    9
  ),
  
  ylim = c(
    -1.6,
    10.8
  ),
  
  clip = "off"
) +
  
  
  theme_void() +
  
  theme(
    
    legend.position = "none",
    
    plot.margin = margin(
      25,
      30,
      25,
      30
    ),
    
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5,
      colour = "#263238"
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      colour = "#546E7A",
      margin = margin(
        b = 15
      )
    )
  ) +
  
  labs(
    
    title =
      "Discovery of a Network-Supported Molecular State RNA Signature",
    
    subtitle =
      "Multi-omics patient-state discovery, robustness testing and network-guided signature reduction"
  )


# ============================================================
# 6. SHOW FIGURE
# ============================================================

signature_workflow


# ============================================================
# 7. SAVE PDF
# ============================================================

ggsave(
  
  filename =
    "Molecular_State_RNA_Signature_Workflow.pdf",
  
  plot =
    signature_workflow,
  
  width = 11,
  
  height = 16,
  
  units = "in",
  
  device = cairo_pdf
)


# ============================================================
# 8. SAVE HIGH-RESOLUTION PNG
# ============================================================

ggsave(
  
  filename =
    "Molecular_State_RNA_Signature_Workflow.png",
  
  plot =
    signature_workflow,
  
  width = 11,
  
  height = 16,
  
  units = "in",
  
  dpi = 400,
  
  bg = "white"
)


# ============================================================
# 9. FILE LOCATION
# ============================================================

getwd()
# ============================================================
# GSEA CORE ENRICHMENT GENES
# State1 vs State2
# ============================================================

library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)


# ============================================================
# 1. CORE ENRICHMENT GENES -> LONG FORMAT
# ============================================================

core_gene_long <- gsea_go_sig %>%
  
  dplyr::mutate(
    
    State = ifelse(
      NES > 0,
      "State2",
      "State1"
    )
    
  ) %>%
  
  dplyr::select(
    ID,
    Description,
    NES,
    p.adjust,
    State,
    core_enrichment
  ) %>%
  
  tidyr::separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  
  dplyr::rename(
    ENTREZID = core_enrichment
  )


# ============================================================
# 2. ENTREZ -> GENE SYMBOL
# ============================================================

entrez_to_symbol <- AnnotationDbi::select(
  
  org.Hs.eg.db,
  
  keys = unique(
    core_gene_long$ENTREZID
  ),
  
  keytype = "ENTREZID",
  
  columns = "SYMBOL"
)


core_gene_long <- core_gene_long %>%
  
  dplyr::left_join(
    entrez_to_symbol,
    by = "ENTREZID"
  ) %>%
  
  dplyr::filter(
    !is.na(SYMBOL)
  )


# ============================================================
# 3. ADD ORIGINAL RNA DIFFERENCE
# ============================================================

core_gene_long <- core_gene_long %>%
  
  dplyr::left_join(
    
    gsea_stats %>%
      
      dplyr::select(
        Gene,
        Difference,
        Z,
        FDR
      ),
    
    by = c(
      "SYMBOL" = "Gene"
    )
  )


# ============================================================
# 4. STATE1 CORE GENES
# ============================================================

state1_core_genes <- core_gene_long %>%
  
  dplyr::filter(
    State == "State1"
  ) %>%
  
  dplyr::group_by(
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    N_pathways =
      dplyr::n_distinct(
        Description
      ),
    
    Mean_NES =
      mean(
        NES,
        na.rm = TRUE
      ),
    
    Best_FDR =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    RNA_Difference =
      dplyr::first(
        Difference
      ),
    
    RNA_Z =
      dplyr::first(
        Z
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      N_pathways
    ),
    Best_FDR
  )


# ============================================================
# 5. STATE2 CORE GENES
# ============================================================

state2_core_genes <- core_gene_long %>%
  
  dplyr::filter(
    State == "State2"
  ) %>%
  
  dplyr::group_by(
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    N_pathways =
      dplyr::n_distinct(
        Description
      ),
    
    Mean_NES =
      mean(
        NES,
        na.rm = TRUE
      ),
    
    Best_FDR =
      min(
        p.adjust,
        na.rm = TRUE
      ),
    
    RNA_Difference =
      dplyr::first(
        Difference
      ),
    
    RNA_Z =
      dplyr::first(
        Z
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      N_pathways
    ),
    Best_FDR
  )


# ============================================================
# 6. VIEW TOP CORE GENES
# ============================================================

cat("\nSTATE 1 CORE GENES\n")

print(
  head(
    state1_core_genes,
    30
  )
)


cat("\nSTATE 2 CORE GENES\n")

print(
  head(
    state2_core_genes,
    30
  )
)


# ============================================================
# 7. CHECK THERAPY-RELEVANT IMMUNE GENES
#
# This is annotation only, NOT yet evidence of drug response.
# ============================================================

immune_targets <- c(
  
  "PDCD1",
  "CD274",
  "PDCD1LG2",
  
  "CTLA4",
  
  "LAG3",
  "TIGIT",
  "HAVCR2",
  
  "CD27",
  "CD28",
  "ICOS",
  
  "TNFRSF4",
  "TNFRSF9",
  
  "JAK1",
  "JAK2",
  "STAT1",
  "STAT3",
  
  "IFNG",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  
  "CD3D",
  "CD3E",
  
  "CD8A",
  "CD8B",
  
  "GZMA",
  "GZMB",
  "PRF1"
)


state1_immune_targets <- state1_core_genes %>%
  
  dplyr::filter(
    SYMBOL %in%
      immune_targets
  )


state1_immune_targets


# ============================================================
# 8. CHECK CANCER / DRUG-RELATED TARGET CANDIDATES
#
# Preliminary annotation panel only.
# ============================================================

candidate_targets <- c(
  
  # DNA repair
  "BRCA1",
  "BRCA2",
  "ATM",
  "ATR",
  "CHEK1",
  "CHEK2",
  "WEE1",
  "PARP1",
  
  # PI3K / AKT / mTOR
  "PIK3CA",
  "PIK3CB",
  "AKT1",
  "AKT2",
  "MTOR",
  "PTEN",
  
  # MAPK
  "KRAS",
  "NRAS",
  "BRAF",
  "RAF1",
  "MAP2K1",
  "MAP2K2",
  
  # RTKs
  "ERBB2",
  "ERBB3",
  "EGFR",
  "MET",
  "FGFR1",
  "FGFR2",
  "FGFR3",
  
  # Cell cycle
  "CCNE1",
  "CDK2",
  "CDK4",
  "CDK6",
  "AURKA",
  "AURKB",
  
  # Transcription / chromatin
  "BET1",
  "BRD2",
  "BRD3",
  "BRD4",
  "EZH2",
  "HDAC1",
  "HDAC2",
  
  # Apoptosis
  "BCL2",
  "BCL2L1",
  "MCL1",
  
  # Immune
  immune_targets
)


# ============================================================
# 9. FIND CANDIDATE TARGETS IN BOTH STATES
# ============================================================

state1_candidate_targets <- state1_core_genes %>%
  
  dplyr::filter(
    SYMBOL %in%
      candidate_targets
  )


state2_candidate_targets <- state2_core_genes %>%
  
  dplyr::filter(
    SYMBOL %in%
      candidate_targets
  )


cat("\nSTATE 1 CANDIDATE TARGETS\n")

print(
  state1_candidate_targets
)


cat("\nSTATE 2 CANDIDATE TARGETS\n")

print(
  state2_candidate_targets
)


# ============================================================
# 10. MOST RECURRENT CORE GENES
# ============================================================

top_state1_genes <- state1_core_genes %>%
  
  dplyr::slice_head(
    n = 20
  ) %>%
  
  dplyr::mutate(
    State = "State 1"
  )


top_state2_genes <- state2_core_genes %>%
  
  dplyr::slice_head(
    n = 20
  ) %>%
  
  dplyr::mutate(
    State = "State 2"
  )


top_core_plot <- dplyr::bind_rows(
  
  top_state1_genes,
  top_state2_genes
  
) %>%
  
  dplyr::mutate(
    
    Gene_State =
      paste(
        SYMBOL,
        State,
        sep = "_"
      )
    
  )


# ============================================================
# 11. PLOT
# ============================================================

top_core_plot$Gene_State <- factor(
  
  top_core_plot$Gene_State,
  
  levels =
    top_core_plot$Gene_State[
      order(
        top_core_plot$N_pathways
      )
    ]
)


core_plot <- ggplot(
  
  top_core_plot,
  
  aes(
    x = N_pathways,
    y = Gene_State,
    shape = State,
    size = abs(RNA_Z)
  )
  
) +
  
  geom_point(
    alpha = 0.85
  ) +
  
  labs(
    
    title =
      "Core Genes Driving State-Associated Pathway Programs",
    
    subtitle =
      "Genes ranked by recurrence across significant GSEA pathways",
    
    x =
      "Number of significant pathways",
    
    y =
      NULL,
    
    size =
      "|RNA Z|"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


print(
  core_plot
)


# ============================================================
# 12. SAVE TABLES
# ============================================================

write.csv(
  
  state1_core_genes,
  
  "State1_GSEA_core_genes.csv",
  
  row.names = FALSE
)


write.csv(
  
  state2_core_genes,
  
  "State2_GSEA_core_genes.csv",
  
  row.names = FALSE
)


write.csv(
  
  state1_candidate_targets,
  
  "State1_candidate_drug_targets.csv",
  
  row.names = FALSE
)


write.csv(
  
  state2_candidate_targets,
  
  "State2_candidate_drug_targets.csv",
  
  row.names = FALSE
)


# ============================================================
# 13. IMPORTANT OUTPUT
# ============================================================

cat(
  "\n===== STATE 1 TOP CORE GENES =====\n"
)

print(
  head(
    state1_core_genes,
    30
  )
)


cat(
  "\n===== STATE 2 TOP CORE GENES =====\n"
)

print(
  head(
    state2_core_genes,
    30
  )
)


cat(
  "\n===== STATE 1 THERAPY-RELATED TARGETS =====\n"
)

print(
  state1_candidate_targets
)


cat(
  "\n===== STATE 2 THERAPY-RELATED TARGETS =====\n"
)

print(
  state2_candidate_targets
)
library(dplyr)
library(ggplot2)

# ============================================================
# TOP CORE GENES PER STATE
# ============================================================

top_n_core <- 25

plot_state1 <- state1_core_genes %>%
  arrange(desc(N_pathways), Best_FDR) %>%
  slice_head(n = top_n_core) %>%
  mutate(State = "State 1")

plot_state2 <- state2_core_genes %>%
  arrange(desc(N_pathways), Best_FDR) %>%
  slice_head(n = top_n_core) %>%
  mutate(State = "State 2")

core_plot_df <- bind_rows(
  plot_state1,
  plot_state2
)

# ============================================================
# PLOT
# ============================================================

p_core <- ggplot(
  core_plot_df,
  aes(
    x = reorder(SYMBOL, N_pathways),
    y = N_pathways,
    size = abs(RNA_Z),
    fill = RNA_Difference
  )
) +
  geom_point(
    shape = 21,
    alpha = 0.85
  ) +
  coord_flip() +
  facet_wrap(
    ~ State,
    scales = "free_y"
  ) +
  labs(
    title = "Core Genes Driving Molecular State Programs",
    subtitle = "Leading-edge genes from significant GO-GSEA pathways",
    x = NULL,
    y = "Number of significant pathways",
    size = "|RNA Z|",
    fill = "State2 - State1\nRNA difference"
  ) +
  theme_classic(base_size = 12)

print(p_core)


# ============================================================
# PATHWAY <-> CORE GENE NETWORK
# State 1 vs State 2
# ============================================================

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)


# ============================================================
# FUNCTION TO BUILD NETWORK FOR ONE STATE
# ============================================================

make_pathway_gene_network <- function(
    pathway_df,
    state_name,
    top_n_pathways = 8
) {
  
  # ----------------------------------------------------------
  # 1. Select most significant pathways
  # ----------------------------------------------------------
  
  top_pathways <- pathway_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n_pathways)
  
  
  # ----------------------------------------------------------
  # 2. Extract pathway-core gene relationships
  # ----------------------------------------------------------
  
  edges <- top_pathways %>%
    select(
      ID,
      Description,
      NES,
      p.adjust,
      core_enrichment
    ) %>%
    separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%
    rename(
      ENTREZID = core_enrichment
    )
  
  
  # ----------------------------------------------------------
  # 3. Convert ENTREZ IDs to gene symbols
  # ----------------------------------------------------------
  
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(edges$ENTREZID),
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  
  
  edges <- edges %>%
    left_join(
      gene_map,
      by = "ENTREZID"
    ) %>%
    filter(
      !is.na(SYMBOL)
    )
  
  
  # ----------------------------------------------------------
  # 4. Count how many selected pathways each gene connects to
  # ----------------------------------------------------------
  
  gene_degree <- edges %>%
    group_by(SYMBOL) %>%
    summarise(
      N_pathways = n_distinct(Description),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # 5. Create edge table
  # ----------------------------------------------------------
  
  network_edges <- edges %>%
    distinct(
      Description,
      SYMBOL
    ) %>%
    transmute(
      from = Description,
      to = SYMBOL
    )
  
  
  # ----------------------------------------------------------
  # 6. Pathway nodes
  # ----------------------------------------------------------
  
  pathway_nodes <- top_pathways %>%
    transmute(
      name = Description,
      Type = "Pathway",
      Score = -log10(p.adjust)
    )
  
  
  # ----------------------------------------------------------
  # 7. Gene nodes
  # ----------------------------------------------------------
  
  gene_nodes <- gene_degree %>%
    transmute(
      name = SYMBOL,
      Type = "Gene",
      Score = N_pathways
    )
  
  
  # ----------------------------------------------------------
  # 8. Combine nodes
  # ----------------------------------------------------------
  
  nodes <- bind_rows(
    pathway_nodes,
    gene_nodes
  ) %>%
    distinct(
      name,
      .keep_all = TRUE
    )
  
  
  # ----------------------------------------------------------
  # 9. Create graph
  # ----------------------------------------------------------
  
  graph <- graph_from_data_frame(
    network_edges,
    vertices = nodes,
    directed = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 10. Plot
  # ----------------------------------------------------------
  
  p <- ggraph(
    graph,
    layout = "fr"
  ) +
    
    geom_edge_link(
      alpha = 0.25,
      linewidth = 0.5
    ) +
    
    geom_node_point(
      aes(
        shape = Type,
        size = Score
      ),
      alpha = 0.9
    ) +
    
    geom_node_text(
      aes(
        label = name,
        filter = Type == "Pathway"
      ),
      repel = TRUE,
      size = 3.5,
      fontface = "bold"
    ) +
    
    geom_node_text(
      aes(
        label = name,
        filter = Type == "Gene"
      ),
      repel = TRUE,
      size = 2.7
    ) +
    
    scale_shape_manual(
      values = c(
        "Pathway" = 15,
        "Gene" = 16
      )
    ) +
    
    labs(
      title = paste0(
        state_name,
        ": Pathway-Core Gene Network"
      ),
      subtitle = paste0(
        "Top ",
        top_n_pathways,
        " significant GO-GSEA pathways and their leading-edge genes"
      ),
      size = "Connectivity /\nSignificance",
      shape = NULL
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold"
      ),
      
      plot.subtitle = element_text(
        size = 11
      ),
      
      legend.position = "right"
    )
  
  
  return(
    list(
      plot = p,
      edges = edges,
      graph = graph,
      top_pathways = top_pathways,
      gene_degree = gene_degree
    )
  )
}


# ============================================================
# STATE 1 NETWORK
# ============================================================

network_state1 <- make_pathway_gene_network(
  pathway_df = state1_pathways,
  state_name = "State 1",
  top_n_pathways = 8
)

print(
  network_state1$plot
)


# ============================================================
# STATE 2 NETWORK
# ============================================================

network_state2 <- make_pathway_gene_network(
  pathway_df = state2_pathways,
  state_name = "State 2",
  top_n_pathways = 8
)

print(
  network_state2$plot
)


# ============================================================
# SAVE FIGURES
# ============================================================

ggsave(
  "State1_Pathway_CoreGene_Network.pdf",
  network_state1$plot,
  width = 12,
  height = 10
)

ggsave(
  "State2_Pathway_CoreGene_Network.pdf",
  network_state2$plot,
  width = 12,
  height = 10
)


# ============================================================
# MOST CONNECTED GENES
# ============================================================

cat("\n===== STATE 1 MOST CONNECTED CORE GENES =====\n")

print(
  network_state1$gene_degree %>%
    arrange(desc(N_pathways)) %>%
    head(30)
)


cat("\n===== STATE 2 MOST CONNECTED CORE GENES =====\n")

print(
  network_state2$gene_degree %>%
    arrange(desc(N_pathways)) %>%
    head(30)
)

make_pathway_gene_network <- function(
    pathway_df,
    state_name,
    top_n_pathways = 8
) {
  
  # 1. Top pathways
  top_pathways <- pathway_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n_pathways)
  
  # 2. Core genes
  edges <- top_pathways %>%
    dplyr::select(
      ID,
      Description,
      NES,
      p.adjust,
      core_enrichment
    ) %>%
    tidyr::separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%
    dplyr::rename(
      ENTREZID = core_enrichment
    )
  
  # 3. ENTREZ -> SYMBOL
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(edges$ENTREZID),
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  
  edges <- edges %>%
    dplyr::left_join(
      gene_map,
      by = "ENTREZID"
    ) %>%
    dplyr::filter(!is.na(SYMBOL))
  
  # 4. Gene recurrence
  gene_degree <- edges %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      N_pathways = dplyr::n_distinct(Description),
      .groups = "drop"
    )
  
  # 5. Edges
  network_edges <- edges %>%
    dplyr::distinct(
      Description,
      SYMBOL
    ) %>%
    dplyr::transmute(
      from = Description,
      to = SYMBOL
    )
  
  # 6. Pathway nodes
  pathway_nodes <- top_pathways %>%
    dplyr::transmute(
      name = Description,
      Type = "Pathway",
      Score = -log10(p.adjust)
    )
  
  # 7. Gene nodes
  gene_nodes <- gene_degree %>%
    dplyr::transmute(
      name = SYMBOL,
      Type = "Gene",
      Score = N_pathways
    )
  
  # 8. Combine nodes
  nodes <- dplyr::bind_rows(
    pathway_nodes,
    gene_nodes
  ) %>%
    dplyr::distinct(
      name,
      .keep_all = TRUE
    )
  
  # 9. Graph
  graph <- igraph::graph_from_data_frame(
    network_edges,
    vertices = nodes,
    directed = FALSE
  )
  
  # 10. Plot
  p <- ggraph::ggraph(
    graph,
    layout = "fr"
  ) +
    
    ggraph::geom_edge_link(
      alpha = 0.25,
      linewidth = 0.5
    ) +
    
    ggraph::geom_node_point(
      ggplot2::aes(
        shape = Type,
        size = Score
      ),
      alpha = 0.9
    ) +
    
    ggraph::geom_node_text(
      ggplot2::aes(
        label = name,
        filter = Type == "Pathway"
      ),
      repel = TRUE,
      size = 3.5,
      fontface = "bold"
    ) +
    
    ggraph::geom_node_text(
      ggplot2::aes(
        label = name,
        filter = Type == "Gene"
      ),
      repel = TRUE,
      size = 2.7
    ) +
    
    ggplot2::scale_shape_manual(
      values = c(
        "Pathway" = 15,
        "Gene" = 16
      )
    ) +
    
    ggplot2::labs(
      title = paste0(
        state_name,
        ": Pathway-Core Gene Network"
      ),
      subtitle = paste0(
        "Top ",
        top_n_pathways,
        " significant GO-GSEA pathways and their leading-edge genes"
      ),
      size = "Connectivity /\nSignificance",
      shape = NULL
    ) +
    
    ggplot2::theme_void()
  
  return(
    list(
      plot = p,
      edges = edges,
      graph = graph,
      top_pathways = top_pathways,
      gene_degree = gene_degree
    )
  )
}
network_state1 <- make_pathway_gene_network(
  pathway_df = state1_pathways,
  state_name = "State 1",
  top_n_pathways = 8
)

network_state2 <- make_pathway_gene_network(
  pathway_df = state2_pathways,
  state_name = "State 2",
  top_n_pathways = 8
)

print(network_state1$plot)
print(network_state2$plot)

###immune tiplerine bakiyorum
# ============================================================
# IMMUNE HOT-LIKE vs COLD-LIKE ANALYSIS
# Robust State 1 vs State 2
#
# Full protein-coding VST RNA
# + MSigDB Hallmark immune programs
# + sample-level ssGSEA
# ============================================================


# ============================================================
# 0. INSTALL PACKAGES IF NEEDED
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GSVA", quietly = TRUE)) {
  BiocManager::install("GSVA")
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(GSVA)
library(msigdbr)
library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 2. USE FULL PROTEIN-CODING VST MATRIX
#
# NOT rna_mofa_new.
# We want the broader RNA matrix so immune genes are not lost
# because of MOFA variance filtering.
# ============================================================

immune_expr <- rna_vst_pc


# Remove RNA_ prefix if present
rownames(immune_expr) <- sub(
  "^RNA_",
  "",
  rownames(immune_expr)
)


# ============================================================
# 3. ALIGN PATIENTS WITH ROBUST STATES
# ============================================================

common_patients <- intersect(
  colnames(immune_expr),
  rownames(factor_scores)
)


immune_expr <- immune_expr[
  ,
  common_patients,
  drop = FALSE
]


state_vector_immune <- robust_state[
  match(
    common_patients,
    rownames(factor_scores)
  )
]


names(state_vector_immune) <- common_patients


cat(
  "\nNumber of patients:",
  length(state_vector_immune),
  "\n"
)


cat(
  "\nState distribution:\n"
)

print(
  table(state_vector_immune)
)


# ============================================================
# 4. REMOVE DUPLICATED GENE SYMBOLS
#
# If the same gene symbol occurs more than once,
# retain the row with the largest variance.
# ============================================================

gene_var <- apply(
  immune_expr,
  1,
  var,
  na.rm = TRUE
)


immune_expr <- immune_expr[
  order(
    gene_var,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


immune_expr <- immune_expr[
  !duplicated(
    rownames(immune_expr)
  ),
  ,
  drop = FALSE
]


cat(
  "\nGenes in expression matrix:",
  nrow(immune_expr),
  "\n"
)


# ============================================================
# 5. GET MSigDB HALLMARK GENE SETS
# ============================================================

hallmark_df <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)


cat(
  "\nMSigDB Hallmark rows:",
  nrow(hallmark_df),
  "\n"
)


# ============================================================
# 6. DEFINE IMMUNE-RELATED HALLMARK PROGRAMS
#
# IMPORTANT:
# We select established MSigDB gene SETS here.
# Individual genes are NOT manually selected.
# ============================================================

wanted_sets <- c(
  
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  
  "HALLMARK_COMPLEMENT",
  
  "HALLMARK_IL2_STAT5_SIGNALING",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)


# Keep only sets actually available in installed MSigDB version
wanted_sets_available <- intersect(
  wanted_sets,
  unique(hallmark_df$gs_name)
)


cat(
  "\nRequested immune Hallmark sets:",
  length(wanted_sets),
  "\n"
)


cat(
  "Available immune Hallmark sets:",
  length(wanted_sets_available),
  "\n\n"
)


print(
  wanted_sets_available
)


# ============================================================
# 7. STOP IF NOTHING WAS FOUND
# ============================================================

if (length(wanted_sets_available) == 0) {
  
  cat(
    "\nAvailable Hallmark names containing immune keywords:\n"
  )
  
  print(
    unique(hallmark_df$gs_name)[
      grepl(
        "INTERFERON|INFLAMMATORY|IL6|IL2|TNFA|COMPLEMENT",
        unique(hallmark_df$gs_name)
      )
    ]
  )
  
  stop(
    "None of the requested Hallmark sets were found."
  )
}


# ============================================================
# 8. BUILD GENE-SET TABLE
# ============================================================

immune_sets_df <- hallmark_df %>%
  
  dplyr::filter(
    gs_name %in% wanted_sets_available
  ) %>%
  
  dplyr::select(
    gs_name,
    gene_symbol
  ) %>%
  
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  
  dplyr::distinct()


cat(
  "\nRows in immune gene-set table:",
  nrow(immune_sets_df),
  "\n"
)


# ============================================================
# 9. CONVERT TO NAMED LIST
#
# This fixes the previous empty geneSets problem.
# ============================================================

immune_sets_tbl <- immune_sets_df %>%
  
  dplyr::group_by(
    gs_name
  ) %>%
  
  dplyr::summarise(
    
    genes = list(
      unique(gene_symbol)
    ),
    
    .groups = "drop"
  )


immune_sets_list <- immune_sets_tbl$genes


names(immune_sets_list) <-
  immune_sets_tbl$gs_name


cat(
  "\nGene sets successfully created:",
  length(immune_sets_list),
  "\n\n"
)


print(
  names(immune_sets_list)
)


# ============================================================
# 10. CHECK GENE OVERLAP WITH TCGA RNA
# ============================================================

overlap_table <- data.frame(
  
  Program =
    names(immune_sets_list),
  
  MSigDB_genes =
    sapply(
      immune_sets_list,
      length
    ),
  
  Genes_in_TCGA =
    sapply(
      
      immune_sets_list,
      
      function(x) {
        
        length(
          intersect(
            x,
            rownames(immune_expr)
          )
        )
        
      }
      
    ),
  
  stringsAsFactors = FALSE
)


overlap_table$Coverage_percent <-
  
  100 *
  
  overlap_table$Genes_in_TCGA /
  
  overlap_table$MSigDB_genes


cat(
  "\n============================================\n"
)

cat(
  "GENE-SET COVERAGE\n"
)

cat(
  "============================================\n\n"
)


print(
  overlap_table
)


# ============================================================
# 11. KEEP SETS WITH >= 10 OBSERVED GENES
# ============================================================

keep_sets <- overlap_table$Program[
  overlap_table$Genes_in_TCGA >= 10
]


immune_sets_use <- immune_sets_list[
  names(immune_sets_list) %in% keep_sets
]


cat(
  "\nGene sets retained for ssGSEA:",
  length(immune_sets_use),
  "\n\n"
)


print(
  names(immune_sets_use)
)


# ============================================================
# 12. SAFETY CHECK
# ============================================================

if (length(immune_sets_use) == 0) {
  
  stop(
    paste0(
      "No immune gene sets have >=10 genes overlapping ",
      "with the expression matrix. Inspect overlap_table."
    )
  )
  
}


# ============================================================
# 13. ssGSEA
#
# New GSVA API
# ============================================================

ssgsea_par <- GSVA::ssgseaParam(
  
  exprData =
    as.matrix(immune_expr),
  
  geneSets =
    immune_sets_use,
  
  minSize = 10,
  
  maxSize = 500,
  
  normalize = TRUE
)


ssgsea_scores <- GSVA::gsva(
  
  ssgsea_par,
  
  verbose = TRUE
  
)


# ============================================================
# 14. CHECK ssGSEA OUTPUT
# ============================================================

cat(
  "\nssGSEA matrix dimensions:\n"
)

print(
  dim(ssgsea_scores)
)


cat(
  "\nPrograms scored:\n"
)

print(
  rownames(ssgsea_scores)
)


# ============================================================
# 15. CONVERT TO SAMPLE × PROGRAM DATAFRAME
# ============================================================

score_df <- as.data.frame(
  t(ssgsea_scores)
)


score_df$submitter_id <-
  rownames(score_df)


score_df$Robust_State <-
  state_vector_immune[
    score_df$submitter_id
  ]


# ============================================================
# 16. LONG FORMAT
# ============================================================

score_long <- score_df %>%
  
  tidyr::pivot_longer(
    
    cols =
      -c(
        submitter_id,
        Robust_State
      ),
    
    names_to =
      "Program",
    
    values_to =
      "Score"
    
  ) %>%
  
  dplyr::mutate(
    
    State =
      ifelse(
        Robust_State == 1,
        "State 1",
        "State 2"
      )
    
  )


# ============================================================
# 17. STATE1 vs STATE2 STATISTICAL TEST
# ============================================================

immune_stats <- score_long %>%
  
  dplyr::group_by(
    Program
  ) %>%
  
  dplyr::summarise(
    
    N_State1 =
      sum(
        Robust_State == 1
      ),
    
    N_State2 =
      sum(
        Robust_State == 2
      ),
    
    Median_State1 =
      median(
        Score[
          Robust_State == 1
        ],
        na.rm = TRUE
      ),
    
    Median_State2 =
      median(
        Score[
          Robust_State == 2
        ],
        na.rm = TRUE
      ),
    
    Difference =
      Median_State2 -
      Median_State1,
    
    P_value =
      wilcox.test(
        Score[
          Robust_State == 1
        ],
        Score[
          Robust_State == 2
        ],
        exact = FALSE
      )$p.value,
    
    .groups = "drop"
  )


# ============================================================
# 18. FDR CORRECTION
# ============================================================

immune_stats <- immune_stats %>%
  
  dplyr::mutate(
    
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    Direction =
      dplyr::case_when(
        
        Difference < 0 &
          FDR < 0.05 ~
          "Higher in State 1",
        
        Difference > 0 &
          FDR < 0.05 ~
          "Higher in State 2",
        
        TRUE ~
          "Not significant"
        
      )
    
  ) %>%
  
  dplyr::arrange(
    FDR
  )


# ============================================================
# 19. PRINT RESULTS
# ============================================================

cat(
  "\n\n============================================\n"
)

cat(
  "HOT/COLD IMMUNE TEST\n"
)

cat(
  "============================================\n\n"
)


print(
  
  immune_stats %>%
    
    dplyr::select(
      
      Program,
      
      Median_State1,
      
      Median_State2,
      
      Difference,
      
      P_value,
      
      FDR,
      
      Direction
      
    )
  
)


# ============================================================
# 20. CLEAN PROGRAM NAMES
# ============================================================

score_long <- score_long %>%
  
  dplyr::mutate(
    
    Program_clean =
      gsub(
        "^HALLMARK_",
        "",
        Program
      ),
    
    Program_clean =
      gsub(
        "_",
        " ",
        Program_clean
      )
    
  )


immune_stats_plot <- immune_stats %>%
  
  dplyr::mutate(
    
    Program_clean =
      gsub(
        "^HALLMARK_",
        "",
        Program
      ),
    
    Program_clean =
      gsub(
        "_",
        " ",
        Program_clean
      ),
    
    Significance =
      -log10(FDR)
    
  )


# ============================================================
# 21. BOXPLOT
# ============================================================

p_hotcold <- ggplot2::ggplot(
  
  score_long,
  
  ggplot2::aes(
    x = State,
    y = Score,
    fill = State
  )
  
) +
  
  ggplot2::geom_boxplot(
    
    width = 0.60,
    
    outlier.shape = NA,
    
    alpha = 0.70
    
  ) +
  
  ggplot2::geom_jitter(
    
    width = 0.12,
    
    alpha = 0.25,
    
    size = 0.7
    
  ) +
  
  ggplot2::facet_wrap(
    
    ~ Program_clean,
    
    scales = "free_y",
    
    ncol = 3
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Immune Programs Across Robust Molecular States",
    
    subtitle =
      "Sample-level ssGSEA using curated MSigDB Hallmark gene sets",
    
    x = NULL,
    
    y =
      "ssGSEA score"
    
  ) +
  
  ggplot2::theme_classic(
    base_size = 11
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "none",
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 9
      ),
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 15
      )
    
  )


print(
  p_hotcold
)


# ============================================================
# 22. EFFECT PLOT
#
# Negative difference = State1 higher
# Positive difference = State2 higher
# ============================================================

p_effect <- ggplot2::ggplot(
  
  immune_stats_plot,
  
  ggplot2::aes(
    
    x = Difference,
    
    y = reorder(
      Program_clean,
      Difference
    ),
    
    size = Significance
    
  )
  
) +
  
  ggplot2::geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.4
    
  ) +
  
  ggplot2::geom_point(
    alpha = 0.85
  ) +
  
  ggplot2::labs(
    
    title =
      "State-Specific Immune Program Differences",
    
    subtitle =
      "Negative values = higher in State 1; positive values = higher in State 2",
    
    x =
      "Median ssGSEA difference (State 2 - State 1)",
    
    y = NULL,
    
    size =
      "-log10(FDR)"
    
  ) +
  
  ggplot2::theme_classic(
    base_size = 11
  )


print(
  p_effect
)


# ============================================================
# 23. SAVE FIGURES
# ============================================================

ggplot2::ggsave(
  
  filename =
    "RobustStates_Immune_ssGSEA_Boxplots.pdf",
  
  plot =
    p_hotcold,
  
  width = 11,
  
  height = 8,
  
  units = "in"
  
)


ggplot2::ggsave(
  
  filename =
    "RobustStates_Immune_ssGSEA_Effect.pdf",
  
  plot =
    p_effect,
  
  width = 9,
  
  height = 5.5,
  
  units = "in"
  
)


# ============================================================
# 24. SAVE RESULTS
# ============================================================

write.csv(
  
  immune_stats,
  
  "RobustStates_Immune_ssGSEA_Statistics.csv",
  
  row.names = FALSE
  
)


write.csv(
  
  overlap_table,
  
  "Immune_GeneSet_Coverage.csv",
  
  row.names = FALSE
  
)


write.csv(
  
  score_df,
  
  "Patient_Immune_ssGSEA_Scores.csv",
  
  row.names = FALSE
  
)


# ============================================================
# 25. FINAL OUTPUT
# ============================================================

cat(
  "\n\n============================================\n"
)

cat(
  "FINAL HOT/COLD RESULTS\n"
)

cat(
  "============================================\n\n"
)


print(
  
  immune_stats %>%
    
    dplyr::select(
      
      Program,
      
      Median_State1,
      
      Median_State2,
      
      Difference,
      
      FDR,
      
      Direction
      
    )
  
)
install.packages("remotes")

remotes::install_github(
  "omnideconv/immunedeconv",
  dependencies = TRUE
)

library(immunedeconv)
immunedeconv::deconvolution_methods

if (!requireNamespace("estimate", quietly = TRUE)) {
  install.packages(
    "estimate",
    repos = c(
      "https://cloud.r-project.org",
      "https://r-forge.r-project.org"
    )
  )
}

library(estimate)


url <- paste0(
  "https://api.gdc.cancer.gov/data/",
  "b3df502e-3594-46ef-9f94-d041a20a0b9a"
)

download.file(
  url,
  destfile = "TCGA_CIBERSORT_relative.tsv",
  mode = "wb",
  method = "libcurl"
)

cibersort_tcga <- read.delim(
  "TCGA_CIBERSORT_relative.tsv",
  check.names = FALSE
)

dim(cibersort_tcga)
colnames(cibersort_tcga)
head(cibersort_tcga[, 1:5])


# ============================================================
# TCGA PANIMMUNE CIBERSORT
# OV Robust State1 vs State2
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 1. FILTER OVARIAN CANCER
# ============================================================

ov_ciber <- cibersort_tcga %>%
  dplyr::filter(CancerType == "OV")


cat("OV CIBERSORT samples:", nrow(ov_ciber), "\n")


# ============================================================
# 2. CREATE TCGA PATIENT ID
#
# Example:
# TCGA.XX.XXXX.01A...  -> TCGA-XX-XXXX
# ============================================================

ov_ciber$submitter_id <- gsub(
  "\\.",
  "-",
  substr(ov_ciber$SampleID, 1, 12)
)


cat(
  "Unique OV patients:",
  length(unique(ov_ciber$submitter_id)),
  "\n"
)


# ============================================================
# 3. CHECK OVERLAP WITH OUR 277 PATIENTS
# ============================================================

state_df <- data.frame(
  
  submitter_id = rownames(factor_scores),
  
  Robust_State = robust_state,
  
  stringsAsFactors = FALSE
)


common_ids <- intersect(
  state_df$submitter_id,
  ov_ciber$submitter_id
)


cat(
  "\nOur molecular-state patients:",
  nrow(state_df),
  "\n"
)

cat(
  "CIBERSORT matched patients:",
  length(common_ids),
  "\n\n"
)


# ============================================================
# 4. KEEP MATCHED PATIENTS
# ============================================================

ov_ciber_match <- ov_ciber %>%
  
  dplyr::filter(
    submitter_id %in% common_ids
  ) %>%
  
  dplyr::left_join(
    state_df,
    by = "submitter_id"
  )


cat("\nState distribution:\n")

print(
  table(ov_ciber_match$Robust_State)
)


# ============================================================
# 5. CHECK WHETHER PATIENTS HAVE DUPLICATE SAMPLES
# ============================================================

dup_check <- ov_ciber_match %>%
  
  dplyr::count(
    submitter_id
  ) %>%
  
  dplyr::filter(
    n > 1
  )


cat(
  "\nPatients with >1 CIBERSORT sample:",
  nrow(dup_check),
  "\n"
)


# ============================================================
# 6. IF DUPLICATES EXIST:
# KEEP SAMPLE WITH BEST CIBERSORT P-VALUE
# ============================================================

ov_ciber_match <- ov_ciber_match %>%
  
  dplyr::arrange(
    submitter_id,
    P.value
  ) %>%
  
  dplyr::group_by(
    submitter_id
  ) %>%
  
  dplyr::slice(1) %>%
  
  dplyr::ungroup()


cat(
  "\nFinal matched patients:",
  nrow(ov_ciber_match),
  "\n"
)

print(
  table(ov_ciber_match$Robust_State)
)


# ============================================================
# 7. DEFINE THE 22 CIBERSORT CELL TYPES
# ============================================================

cell_types <- c(
  
  "B.cells.naive",
  "B.cells.memory",
  "Plasma.cells",
  
  "T.cells.CD8",
  "T.cells.CD4.naive",
  "T.cells.CD4.memory.resting",
  "T.cells.CD4.memory.activated",
  "T.cells.follicular.helper",
  "T.cells.regulatory..Tregs.",
  "T.cells.gamma.delta",
  
  "NK.cells.resting",
  "NK.cells.activated",
  
  "Monocytes",
  
  "Macrophages.M0",
  "Macrophages.M1",
  "Macrophages.M2",
  
  "Dendritic.cells.resting",
  "Dendritic.cells.activated",
  
  "Mast.cells.resting",
  "Mast.cells.activated",
  
  "Eosinophils",
  "Neutrophils"
)


# ============================================================
# 8. LONG FORMAT
# ============================================================

ciber_long <- ov_ciber_match %>%
  
  dplyr::select(
    submitter_id,
    Robust_State,
    dplyr::all_of(cell_types)
  ) %>%
  
  tidyr::pivot_longer(
    
    cols =
      dplyr::all_of(cell_types),
    
    names_to =
      "Cell_Type",
    
    values_to =
      "Fraction"
    
  ) %>%
  
  dplyr::mutate(
    
    State =
      ifelse(
        Robust_State == 1,
        "State 1",
        "State 2"
      )
    
  )


# ============================================================
# 9. STATE1 vs STATE2
# ============================================================

ciber_stats <- ciber_long %>%
  
  dplyr::group_by(
    Cell_Type
  ) %>%
  
  dplyr::summarise(
    
    Median_State1 =
      median(
        Fraction[Robust_State == 1],
        na.rm = TRUE
      ),
    
    Median_State2 =
      median(
        Fraction[Robust_State == 2],
        na.rm = TRUE
      ),
    
    Difference =
      Median_State2 -
      Median_State1,
    
    P_value =
      wilcox.test(
        Fraction[Robust_State == 1],
        Fraction[Robust_State == 2],
        exact = FALSE
      )$p.value,
    
    .groups = "drop"
    
  ) %>%
  
  dplyr::mutate(
    
    FDR =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    Direction =
      dplyr::case_when(
        
        FDR < 0.05 & Difference < 0 ~
          "Higher in State 1",
        
        FDR < 0.05 & Difference > 0 ~
          "Higher in State 2",
        
        TRUE ~
          "Not significant"
        
      )
    
  ) %>%
  
  dplyr::arrange(
    FDR
  )


# ============================================================
# 10. PRINT ALL RESULTS
# ============================================================

cat(
  "\n\n===========================================\n"
)

cat(
  "CIBERSORT STATE1 vs STATE2\n"
)

cat(
  "===========================================\n\n"
)


print(
  ciber_stats,
  n = 22
)


# ============================================================
# 11. HOT/COLD-RELEVANT CELLS
# ============================================================

hot_cold_cells <- c(
  
  "T.cells.CD8",
  
  "T.cells.CD4.memory.activated",
  
  "T.cells.follicular.helper",
  
  "NK.cells.activated",
  
  "Macrophages.M1",
  
  "Macrophages.M2",
  
  "Dendritic.cells.activated",
  
  "T.cells.regulatory..Tregs."
)


hotcold_stats <- ciber_stats %>%
  
  dplyr::filter(
    Cell_Type %in% hot_cold_cells
  )


cat(
  "\n\n===========================================\n"
)

cat(
  "HOT/COLD RELEVANT POPULATIONS\n"
)

cat(
  "===========================================\n\n"
)


print(
  hotcold_stats
)


# ============================================================
# 12. BOXPLOT — HOT/COLD RELEVANT CELLS
# ============================================================

hotcold_long <- ciber_long %>%
  
  dplyr::filter(
    Cell_Type %in% hot_cold_cells
  )


p_ciber <- ggplot2::ggplot(
  
  hotcold_long,
  
  ggplot2::aes(
    x = State,
    y = Fraction,
    fill = State
  )
  
) +
  
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  
  ggplot2::geom_jitter(
    width = 0.12,
    alpha = 0.25,
    size = 0.7
  ) +
  
  ggplot2::facet_wrap(
    ~ Cell_Type,
    scales = "free_y",
    ncol = 4
  ) +
  
  ggplot2::labs(
    
    title =
      "Immune Cell Composition Across Robust Molecular States",
    
    subtitle =
      "TCGA PanImmune CIBERSORT relative fractions",
    
    x = NULL,
    
    y =
      "CIBERSORT relative fraction"
    
  ) +
  
  ggplot2::theme_classic(
    base_size = 11
  ) +
  
  ggplot2::theme(
    legend.position = "none",
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 8
      )
  )


print(p_ciber)


# ============================================================
# 13. SAVE
# ============================================================

write.csv(
  ciber_stats,
  "TCGA_OV_CIBERSORT_State1_vs_State2.csv",
  row.names = FALSE
)


write.csv(
  ov_ciber_match,
  "TCGA_OV_CIBERSORT_MatchedPatients.csv",
  row.names = FALSE
)


ggplot2::ggsave(
  "TCGA_OV_CIBERSORT_HotCold.pdf",
  p_ciber,
  width = 12,
  height = 6
)


# ============================================================
# 14. FINAL OUTPUT
# ============================================================

cat(
  "\n\n===== FINAL HOT/COLD CIBERSORT RESULTS =====\n\n"
)

print(
  hotcold_stats %>%
    dplyr::select(
      Cell_Type,
      Median_State1,
      Median_State2,
      Difference,
      FDR,
      Direction
    )
)
# State 1 - all core genes from representative programs
state1_drug_genes <- reduced_state1$core_long %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# State 2
state2_drug_genes <- reduced_state2$core_long %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

cat("State 1 genes:", nrow(state1_drug_genes), "\n")
cat("State 2 genes:", nrow(state2_drug_genes), "\n")

head(state1_drug_genes)
head(state2_drug_genes)
###
#drug
# ============================================================
# OPEN TARGETS DRUG MAPPING
# State 1 vs State 2 core-enrichment genes
# ============================================================

# ------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------

packages <- c(
  "httr",
  "jsonlite",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


# ============================================================
# 1. PREPARE STATE GENE TABLES
# ============================================================

state1_input <- state1_drug_genes %>%
  dplyr::mutate(State = "State1")

state2_input <- state2_drug_genes %>%
  dplyr::mutate(State = "State2")


drug_input <- dplyr::bind_rows(
  state1_input,
  state2_input
) %>%
  dplyr::select(
    State,
    Program,
    SYMBOL,
    ENTREZID,
    Difference,
    Z,
    FDR
  ) %>%
  dplyr::distinct()


cat("Total input rows:", nrow(drug_input), "\n")
cat("Unique genes:", length(unique(drug_input$SYMBOL)), "\n")


# ============================================================
# 2. OPEN TARGETS API
# ============================================================

ot_url <- "https://api.platform.opentargets.org/api/v4/graphql"


# ============================================================
# 3. FUNCTION: SYMBOL -> OPEN TARGETS / ENSEMBL ID
# ============================================================

get_ot_target <- function(symbol) {
  
  query <- '
  query searchTarget($queryString: String!) {
    search(queryString: $queryString, entityNames: ["target"]) {
      hits {
        id
        name
        entity
      }
    }
  }
  '
  
  response <- httr::POST(
    url = ot_url,
    body = list(
      query = query,
      variables = list(
        queryString = symbol
      )
    ),
    encode = "json"
  )
  
  if (httr::status_code(response) != 200) {
    return(NULL)
  }
  
  result <- httr::content(
    response,
    as = "parsed",
    simplifyVector = FALSE
  )
  
  hits <- result$data$search$hits
  
  if (length(hits) == 0) {
    return(NULL)
  }
  
  # Prefer exact gene-symbol match
  hit_names <- sapply(
    hits,
    function(x) {
      if (is.null(x$name)) NA_character_ else x$name
    }
  )
  
  exact <- which(
    toupper(hit_names) == toupper(symbol)
  )
  
  if (length(exact) > 0) {
    chosen <- hits[[exact[1]]]
  } else {
    chosen <- hits[[1]]
  }
  
  tibble::tibble(
    SYMBOL = symbol,
    ensembl_id = chosen$id,
    OT_name = chosen$name
  )
}


# ============================================================
# 4. TEST API FIRST
# ============================================================

test_target <- get_ot_target("TNF")

print(test_target)

if (is.null(test_target)) {
  stop("Open Targets API test failed.")
}


# ============================================================
# 5. MAP ALL GENES TO OPEN TARGETS IDs
# ============================================================

all_symbols <- unique(drug_input$SYMBOL)

cat(
  "\nMapping",
  length(all_symbols),
  "genes to Open Targets...\n"
)


target_map_list <- vector(
  "list",
  length(all_symbols)
)


for (i in seq_along(all_symbols)) {
  
  gene <- all_symbols[i]
  
  cat(
    "[",
    i,
    "/",
    length(all_symbols),
    "] ",
    gene,
    "\n",
    sep = ""
  )
  
  target_map_list[[i]] <- tryCatch(
    
    get_ot_target(gene),
    
    error = function(e) {
      NULL
    }
    
  )
  
  Sys.sleep(0.05)
}


target_map <- dplyr::bind_rows(
  target_map_list
)


cat(
  "\nGenes mapped:",
  nrow(target_map),
  "\n"
)


# ============================================================
# 6. FUNCTION: GET KNOWN DRUGS FOR TARGET
# ============================================================

get_target_drugs <- function(ensembl_id, symbol) {
  
  query <- '
  query targetDrugs($ensemblId: String!) {

    target(ensemblId: $ensemblId) {

      id
      approvedSymbol

      knownDrugs {

        rows {

          drugId
          prefName
          drugType

          mechanismOfAction
          actionType

          phase

          status

        }

      }

    }

  }
  '
  
  response <- httr::POST(
    
    url = ot_url,
    
    body = list(
      
      query = query,
      
      variables = list(
        ensemblId = ensembl_id
      )
      
    ),
    
    encode = "json"
    
  )
  
  
  if (httr::status_code(response) != 200) {
    return(NULL)
  }
  
  
  result <- httr::content(
    response,
    as = "parsed",
    simplifyVector = FALSE
  )
  
  
  target <- result$data$target
  
  
  if (is.null(target)) {
    return(NULL)
  }
  
  
  rows <- target$knownDrugs$rows
  
  
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  
  
  out <- lapply(
    
    rows,
    
    function(x) {
      
      tibble::tibble(
        
        SYMBOL =
          symbol,
        
        ensembl_id =
          ensembl_id,
        
        Drug_ID =
          ifelse(
            is.null(x$drugId),
            NA,
            x$drugId
          ),
        
        Drug =
          ifelse(
            is.null(x$prefName),
            NA,
            x$prefName
          ),
        
        Drug_Type =
          ifelse(
            is.null(x$drugType),
            NA,
            x$drugType
          ),
        
        Mechanism =
          ifelse(
            is.null(x$mechanismOfAction),
            NA,
            x$mechanismOfAction
          ),
        
        Action_Type =
          ifelse(
            is.null(x$actionType),
            NA,
            paste(
              unlist(x$actionType),
              collapse = "; "
            )
          ),
        
        Max_Phase =
          ifelse(
            is.null(x$phase),
            NA,
            x$phase
          ),
        
        Status =
          ifelse(
            is.null(x$status),
            NA,
            x$status
          )
        
      )
      
    }
    
  )
  
  
  dplyr::bind_rows(out)
}


# ============================================================
# 7. DOWNLOAD DRUG-TARGET RELATIONSHIPS
# ============================================================

drug_results_list <- vector(
  "list",
  nrow(target_map)
)


for (i in seq_len(nrow(target_map))) {
  
  cat(
    "[",
    i,
    "/",
    nrow(target_map),
    "] Drugs for ",
    target_map$SYMBOL[i],
    "\n",
    sep = ""
  )
  
  
  drug_results_list[[i]] <- tryCatch(
    
    get_target_drugs(
      
      target_map$ensembl_id[i],
      
      target_map$SYMBOL[i]
      
    ),
    
    error = function(e) {
      
      cat(
        "   Error:",
        conditionMessage(e),
        "\n"
      )
      
      NULL
      
    }
    
  )
  
  
  Sys.sleep(0.05)
}


drug_target_raw <- dplyr::bind_rows(
  drug_results_list
)


cat(
  "\nDrug-target rows:",
  nrow(drug_target_raw),
  "\n"
)


cat(
  "Genes with at least one drug:",
  length(unique(drug_target_raw$SYMBOL)),
  "\n"
)


# ============================================================
# 8. JOIN BACK TO OUR BIOLOGICAL EVIDENCE
# ============================================================

state_drug_table <- drug_input %>%
  
  dplyr::inner_join(
    drug_target_raw,
    by = "SYMBOL"
  ) %>%
  
  dplyr::distinct()


# ============================================================
# 9. ADD GENE-LEVEL EVIDENCE CLASS
# ============================================================

state_drug_table <- state_drug_table %>%
  
  dplyr::mutate(
    
    Gene_Evidence =
      dplyr::case_when(
        
        FDR < 0.05 ~
          "Core + gene-level significant",
        
        TRUE ~
          "Pathway-core only"
        
      )
    
  )


# ============================================================
# 10. CLINICAL STAGE CLASSIFICATION
# ============================================================

state_drug_table <- state_drug_table %>%
  
  dplyr::mutate(
    
    Clinical_Status =
      dplyr::case_when(
        
        Max_Phase >= 4 ~
          "Approved / Phase 4",
        
        Max_Phase >= 3 ~
          "Phase 3",
        
        Max_Phase >= 2 ~
          "Phase 2",
        
        Max_Phase >= 1 ~
          "Phase 1",
        
        TRUE ~
          "Preclinical / unknown"
        
      )
    
  )


# ============================================================
# 11. SUMMARY BY STATE
# ============================================================

state_summary <- state_drug_table %>%
  
  dplyr::group_by(State) %>%
  
  dplyr::summarise(
    
    Core_genes =
      n_distinct(SYMBOL),
    
    Drugs =
      n_distinct(Drug),
    
    Approved_or_Phase4 =
      n_distinct(
        Drug[
          Max_Phase >= 4
        ]
      ),
    
    Gene_level_significant_targets =
      n_distinct(
        SYMBOL[
          FDR < 0.05
        ]
      ),
    
    .groups = "drop"
  )


cat(
  "\n============================\n"
)

cat(
  "DRUG MAPPING SUMMARY\n"
)

cat(
  "============================\n\n"
)


print(state_summary)


# ============================================================
# 12. HIGHER-CONFIDENCE TARGETS
#
# Core-enrichment gene
# +
# gene-level FDR < 0.05
# +
# at least one drug
# ============================================================

high_confidence <- state_drug_table %>%
  
  dplyr::filter(
    FDR < 0.05
  ) %>%
  
  dplyr::arrange(
    State,
    FDR,
    dplyr::desc(Max_Phase)
  )


cat(
  "\n\n====================================\n"
)

cat(
  "HIGH-CONFIDENCE DRUGGABLE GENES\n"
)

cat(
  "====================================\n\n"
)


print(
  
  high_confidence %>%
    
    dplyr::select(
      
      State,
      Program,
      SYMBOL,
      Difference,
      Z,
      FDR,
      Drug,
      Mechanism,
      Action_Type,
      Max_Phase,
      Clinical_Status
      
    ),
  
  n = 100
  
)


# ============================================================
# 13. APPROVED / PHASE 4 DRUGS
# ============================================================

approved_drugs <- state_drug_table %>%
  
  dplyr::filter(
    Max_Phase >= 4
  ) %>%
  
  dplyr::arrange(
    State,
    FDR
  )


cat(
  "\n\n====================================\n"
)

cat(
  "APPROVED / PHASE 4 DRUG CONNECTIONS\n"
)

cat(
  "====================================\n\n"
)


print(
  
  approved_drugs %>%
    
    dplyr::select(
      
      State,
      Program,
      SYMBOL,
      FDR,
      Gene_Evidence,
      Drug,
      Mechanism,
      Action_Type
      
    ),
  
  n = 100
  
)


# ============================================================
# 14. TARGET-LEVEL SUMMARY
# ============================================================

target_summary <- state_drug_table %>%
  
  dplyr::group_by(
    State,
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    Program =
      paste(
        unique(Program),
        collapse = "; "
      ),
    
    Gene_FDR =
      min(
        FDR,
        na.rm = TRUE
      ),
    
    Difference =
      first(Difference),
    
    Z =
      first(Z),
    
    N_drugs =
      n_distinct(Drug),
    
    Highest_phase =
      max(
        Max_Phase,
        na.rm = TRUE
      ),
    
    Drugs =
      paste(
        head(
          unique(Drug),
          10
        ),
        collapse = "; "
      ),
    
    .groups = "drop"
    
  ) %>%
  
  dplyr::arrange(
    State,
    Gene_FDR
  )


# ============================================================
# 15. SAVE EVERYTHING
# ============================================================

write.csv(
  target_map,
  "OpenTargets_Gene_ID_Mapping.csv",
  row.names = FALSE
)


write.csv(
  drug_target_raw,
  "OpenTargets_Raw_Drug_Targets.csv",
  row.names = FALSE
)


write.csv(
  state_drug_table,
  "State1_State2_OpenTargets_DrugMapping.csv",
  row.names = FALSE
)


write.csv(
  high_confidence,
  "State1_State2_HighConfidence_DrugTargets.csv",
  row.names = FALSE
)


write.csv(
  approved_drugs,
  "State1_State2_ApprovedDrug_Targets.csv",
  row.names = FALSE
)


write.csv(
  target_summary,
  "State1_State2_DruggableTarget_Summary.csv",
  row.names = FALSE
)


# ============================================================
# 16. FINAL QUICK VIEW
# ============================================================

cat(
  "\n\n===== FINAL TARGET SUMMARY =====\n\n"
)


print(
  target_summary,
  n = 50
)



# ============================================================
# CURRENT OPEN TARGETS DRUG MAPPING
# ============================================================

get_target_drugs_v2 <- function(ensembl_id, symbol) {
  
  query <- '
  query targetDrugs($ensemblId: String!) {
    target(ensemblId: $ensemblId) {
      id
      approvedSymbol

      drugAndClinicalCandidates {
        rows {
          maxClinicalStage

          drug {
            id
            name
            drugType
          }

          diseases {
            diseaseFromSource
            disease {
              id
              name
            }
          }
        }
      }
    }
  }
  '
  
  response <- httr::POST(
    ot_url,
    body = list(
      query = query,
      variables = list(ensemblId = ensembl_id)
    ),
    encode = "json"
  )
  
  if (httr::status_code(response) != 200) {
    cat("HTTP error:", symbol, httr::status_code(response), "\n")
    return(NULL)
  }
  
  result <- httr::content(
    response,
    as = "parsed",
    simplifyVector = FALSE
  )
  
  rows <- result$data$target$drugAndClinicalCandidates$rows
  
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  
  out <- lapply(rows, function(x) {
    
    diseases <- x$diseases
    
    disease_names <- character(0)
    disease_ids   <- character(0)
    
    if (!is.null(diseases) && length(diseases) > 0) {
      
      disease_names <- sapply(diseases, function(d) {
        if (!is.null(d$disease) && !is.null(d$disease$name)) {
          d$disease$name
        } else if (!is.null(d$diseaseFromSource)) {
          d$diseaseFromSource
        } else {
          NA_character_
        }
      })
      
      disease_ids <- sapply(diseases, function(d) {
        if (!is.null(d$disease) && !is.null(d$disease$id)) {
          d$disease$id
        } else {
          NA_character_
        }
      })
    }
    
    tibble::tibble(
      SYMBOL = symbol,
      ensembl_id = ensembl_id,
      
      Drug_ID = ifelse(
        is.null(x$drug$id),
        NA_character_,
        x$drug$id
      ),
      
      Drug = ifelse(
        is.null(x$drug$name),
        NA_character_,
        x$drug$name
      ),
      
      Drug_Type = ifelse(
        is.null(x$drug$drugType),
        NA_character_,
        x$drug$drugType
      ),
      
      Clinical_Stage = ifelse(
        is.null(x$maxClinicalStage),
        NA_character_,
        x$maxClinicalStage
      ),
      
      Diseases = paste(
        unique(na.omit(disease_names)),
        collapse = "; "
      ),
      
      Disease_IDs = paste(
        unique(na.omit(disease_ids)),
        collapse = "; "
      )
    )
  })
  
  dplyr::bind_rows(out)
}


# ============================================================
# QUERY ALL MAPPED GENES
# ============================================================

drug_results_v2 <- vector("list", nrow(target_map))

for (i in seq_len(nrow(target_map))) {
  
  cat(
    "[", i, "/", nrow(target_map), "] ",
    target_map$SYMBOL[i], "\n",
    sep = ""
  )
  
  drug_results_v2[[i]] <- tryCatch(
    
    get_target_drugs_v2(
      target_map$ensembl_id[i],
      target_map$SYMBOL[i]
    ),
    
    error = function(e) {
      cat(
        "ERROR:",
        target_map$SYMBOL[i],
        conditionMessage(e),
        "\n"
      )
      NULL
    }
  )
  
  Sys.sleep(0.05)
}


# ============================================================
# COMBINE RESULTS
# ============================================================

drug_target_raw_v2 <- dplyr::bind_rows(drug_results_v2)

cat("\nDrug-target rows:",
    nrow(drug_target_raw_v2), "\n")

cat("Genes with drugs:",
    dplyr::n_distinct(drug_target_raw_v2$SYMBOL), "\n")

cat("Unique drugs:",
    dplyr::n_distinct(drug_target_raw_v2$Drug), "\n")


# ============================================================
# JOIN WITH OUR STATE EVIDENCE
# ============================================================

state_drug_table <- drug_input %>%
  
  dplyr::inner_join(
    drug_target_raw_v2,
    by = "SYMBOL"
  ) %>%
  
  dplyr::distinct()


# ============================================================
# EVIDENCE TYPE
# ============================================================

state_drug_table <- state_drug_table %>%
  
  dplyr::mutate(
    
    Gene_Evidence = dplyr::case_when(
      
      FDR < 0.05 ~
        "Core + gene-level significant",
      
      TRUE ~
        "Pathway-core only"
    )
  )


# ============================================================
# OVARIAN-CANCER RELATED RECORDS
# ============================================================

ovarian_drugs <- state_drug_table %>%
  
  dplyr::filter(
    grepl(
      "ovarian|ovary",
      Diseases,
      ignore.case = TRUE
    )
  )


# ============================================================
# HIGHER-CONFIDENCE MOLECULAR TARGETS
# ============================================================

high_confidence <- state_drug_table %>%
  
  dplyr::filter(
    FDR < 0.05
  ) %>%
  
  dplyr::arrange(
    State,
    FDR
  )


# ============================================================
# SUMMARY
# ============================================================

state_summary <- state_drug_table %>%
  
  dplyr::group_by(State) %>%
  
  dplyr::summarise(
    
    Druggable_genes =
      dplyr::n_distinct(SYMBOL),
    
    Unique_drugs =
      dplyr::n_distinct(Drug),
    
    Gene_level_significant_targets =
      dplyr::n_distinct(SYMBOL[FDR < 0.05]),
    
    Ovarian_related_drugs =
      dplyr::n_distinct(
        Drug[
          grepl(
            "ovarian|ovary",
            Diseases,
            ignore.case = TRUE
          )
        ]
      ),
    
    .groups = "drop"
  )


cat("\n\n===== STATE DRUG SUMMARY =====\n\n")
print(state_summary)


cat("\n\n===== HIGH-CONFIDENCE TARGETS =====\n\n")

print(
  high_confidence %>%
    dplyr::select(
      State,
      Program,
      SYMBOL,
      Difference,
      Z,
      FDR,
      Drug,
      Drug_Type,
      Clinical_Stage,
      Diseases
    ),
  n = 100
)


cat("\n\n===== OVARIAN-CANCER RELATED =====\n\n")

print(
  ovarian_drugs %>%
    dplyr::select(
      State,
      SYMBOL,
      FDR,
      Gene_Evidence,
      Drug,
      Clinical_Stage,
      Diseases
    ),
  n = 100
)


# ============================================================
# SAVE
# ============================================================

write.csv(
  drug_target_raw_v2,
  "OpenTargets_DrugTargets_CurrentAPI.csv",
  row.names = FALSE
)

write.csv(
  state_drug_table,
  "State1_State2_DrugMapping.csv",
  row.names = FALSE
)

write.csv(
  high_confidence,
  "State1_State2_HighConfidence_DrugTargets.csv",
  row.names = FALSE
)

write.csv(
  ovarian_drugs,
  "State1_State2_OvarianCancer_DrugTargets.csv",
  row.names = FALSE
)

write.csv(
  state_summary,
  "State1_State2_DrugSummary.csv",
  row.names = FALSE
)

library(dplyr)
library(ggplot2)
library(forcats)

# ============================================================
# 1. TARGET-LEVEL SUMMARY
# ============================================================

drug_plot_df <- state_drug_table %>%
  dplyr::group_by(State, SYMBOL) %>%
  dplyr::summarise(
    Z = dplyr::first(Z),
    Difference = dplyr::first(Difference),
    FDR = dplyr::first(FDR),
    N_Drugs = dplyr::n_distinct(Drug),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Evidence = ifelse(
      FDR < 0.05,
      "Gene-level significant",
      "Pathway-core only"
    )
  )


# ============================================================
# 2. ORDER GENES BY Z
# ============================================================

drug_plot_df <- drug_plot_df %>%
  dplyr::arrange(Z) %>%
  dplyr::mutate(
    Gene = factor(SYMBOL, levels = unique(SYMBOL))
  )


# ============================================================
# 3. PLOT
# ============================================================

p_drug_targets <- ggplot(
  drug_plot_df,
  aes(
    x = Z,
    y = Gene,
    size = N_Drugs,
    shape = Evidence
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_point(
    aes(fill = State),
    alpha = 0.8
  ) +
  facet_wrap(
    ~State,
    scales = "free_y",
    ncol = 2
  ) +
  scale_shape_manual(
    values = c(
      "Gene-level significant" = 21,
      "Pathway-core only" = 1
    )
  ) +
  labs(
    title = "Drug-Target Landscape of Robust Molecular States",
    subtitle = "Data-derived GSEA core genes mapped to Open Targets",
    x = "RNA state association (Z-score)",
    y = NULL,
    size = "Number of drugs",
    shape = "Evidence",
    fill = "Molecular state"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

print(p_drug_targets)


# ============================================================
# 4. SAVE
# ============================================================

ggsave(
  "State1_State2_Druggable_Targets.pdf",
  p_drug_targets,
  width = 11,
  height = 9
)

ggsave(
  "State1_State2_Druggable_Targets.png",
  p_drug_targets,
  width = 11,
  height = 9,
  dpi = 300
)

########state 1 state 2 arasindaki farklar

