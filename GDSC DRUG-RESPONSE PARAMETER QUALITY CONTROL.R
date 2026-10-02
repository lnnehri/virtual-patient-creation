# ============================================================
# GDSC DRUG-RESPONSE PARAMETER QUALITY CONTROL
# Parameters:
# AAC / IC50 / EC50 / E_inf / HS
# ============================================================

library(PharmacoGx)
library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 1. GET ORIGINAL GDSC SENSITIVITY DATA
# ============================================================

gdsc_sens <- sensitivityProfiles(gdsc2)

cat("\n============================================\n")
cat("AVAILABLE GDSC PARAMETERS\n")
cat("============================================\n")

print(colnames(gdsc_sens))


# Parameters we want to investigate
param_names <- c(
  "aac_recomputed",
  "ic50_recomputed",
  "EC50",
  "E_inf",
  "HS"
)


# Check that all parameters exist
missing_params <- setdiff(
  param_names,
  colnames(gdsc_sens)
)

if(length(missing_params) > 0) {
  
  stop(
    paste(
      "Missing parameters:",
      paste(
        missing_params,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 2. CREATE PARAMETER DATA FRAME
# ============================================================

gdsc_params <- as.data.frame(
  gdsc_sens[
    ,
    param_names,
    drop = FALSE
  ]
)


cat("\nGDSC sensitivity dimensions:\n")
print(dim(gdsc_params))


# ============================================================
# 3. QC FUNCTION
# ============================================================

parameter_qc <- function(x, parameter_name) {
  
  x <- as.numeric(x)
  
  finite_x <- x[
    is.finite(x)
  ]
  
  if(length(finite_x) == 0) {
    
    return(
      data.frame(
        Parameter = parameter_name,
        N_total = length(x),
        N_available = 0,
        N_NA = sum(is.na(x)),
        Percent_NA = 100 * mean(is.na(x)),
        N_Inf = sum(is.infinite(x)),
        N_finite = 0,
        Min = NA_real_,
        Q01 = NA_real_,
        Q05 = NA_real_,
        Q25 = NA_real_,
        Median = NA_real_,
        Mean = NA_real_,
        Q75 = NA_real_,
        Q95 = NA_real_,
        Q99 = NA_real_,
        Max = NA_real_
      )
    )
  }
  
  data.frame(
    
    Parameter = parameter_name,
    
    N_total =
      length(x),
    
    N_available =
      sum(!is.na(x)),
    
    N_NA =
      sum(is.na(x)),
    
    Percent_NA =
      100 * mean(is.na(x)),
    
    N_Inf =
      sum(is.infinite(x)),
    
    N_finite =
      length(finite_x),
    
    Min =
      min(
        finite_x
      ),
    
    Q01 =
      as.numeric(
        quantile(
          finite_x,
          0.01,
          na.rm = TRUE
        )
      ),
    
    Q05 =
      as.numeric(
        quantile(
          finite_x,
          0.05,
          na.rm = TRUE
        )
      ),
    
    Q25 =
      as.numeric(
        quantile(
          finite_x,
          0.25,
          na.rm = TRUE
        )
      ),
    
    Median =
      median(
        finite_x,
        na.rm = TRUE
      ),
    
    Mean =
      mean(
        finite_x,
        na.rm = TRUE
      ),
    
    Q75 =
      as.numeric(
        quantile(
          finite_x,
          0.75,
          na.rm = TRUE
        )
      ),
    
    Q95 =
      as.numeric(
        quantile(
          finite_x,
          0.95,
          na.rm = TRUE
        )
      ),
    
    Q99 =
      as.numeric(
        quantile(
          finite_x,
          0.99,
          na.rm = TRUE
        )
      ),
    
    Max =
      max(
        finite_x
      )
  )
}


# ============================================================
# 4. QC — ALL GDSC
# ============================================================

qc_all <- dplyr::bind_rows(
  
  lapply(
    
    param_names,
    
    function(p) {
      
      parameter_qc(
        gdsc_params[[p]],
        p
      )
      
    }
  )
)


cat("\n")
cat("============================================\n")
cat("ALL GDSC PARAMETER QC\n")
cat("============================================\n")

print(
  as.data.frame(qc_all)
)


# ============================================================
# 5. CHECK HGSOC DATA EXISTS
# ============================================================

if(!exists("gdsc_hgsoc_all")) {
  
  stop(
    "Object 'gdsc_hgsoc_all' does not exist. Run the HGSOC matching step first."
  )
}


# ============================================================
# 6. QC — HGSOC SUBSET
# ============================================================

qc_hgsoc <- dplyr::bind_rows(
  
  lapply(
    
    param_names,
    
    function(p) {
      
      parameter_qc(
        gdsc_hgsoc_all[[p]],
        p
      )
      
    }
  )
)


cat("\n")
cat("============================================\n")
cat("HGSOC PARAMETER QC\n")
cat("============================================\n")

print(
  as.data.frame(qc_hgsoc)
)


# ============================================================
# 7. BASE-R SUMMARY — ALL GDSC
# ============================================================

cat("\n")
cat("============================================\n")
cat("ALL GDSC BASE-R SUMMARIES\n")
cat("============================================\n")


for(p in param_names) {
  
  cat("\n")
  cat("--------------------------------------------\n")
  cat(p, "\n")
  cat("--------------------------------------------\n")
  
  print(
    summary(
      gdsc_params[[p]]
    )
  )
}


# ============================================================
# 8. BASE-R SUMMARY — HGSOC
# ============================================================

cat("\n")
cat("============================================\n")
cat("HGSOC BASE-R SUMMARIES\n")
cat("============================================\n")


for(p in param_names) {
  
  cat("\n")
  cat("--------------------------------------------\n")
  cat(p, "\n")
  cat("--------------------------------------------\n")
  
  print(
    summary(
      gdsc_hgsoc_all[[p]]
    )
  )
}


# ============================================================
# 9. LONG FORMAT — ALL GDSC
# ============================================================

gdsc_long <- gdsc_params

gdsc_long$Experiment <-
  rownames(gdsc_params)


gdsc_long <- gdsc_long %>%
  
  tidyr::pivot_longer(
    
    cols =
      dplyr::all_of(
        param_names
      ),
    
    names_to =
      "Parameter",
    
    values_to =
      "Value"
  ) %>%
  
  dplyr::filter(
    is.finite(Value)
  )


# ============================================================
# 10. RAW DISTRIBUTIONS
# ============================================================

p_raw <- ggplot(
  gdsc_long,
  aes(
    x = Value
  )
) +
  
  geom_histogram(
    bins = 60
  ) +
  
  facet_wrap(
    ~ Parameter,
    scales = "free"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    
    title =
      "GDSC drug-response parameter distributions",
    
    subtitle =
      "Raw finite values",
    
    x =
      "Parameter value",
    
    y =
      "Number of experiments"
  )


print(p_raw)


# ============================================================
# 11. CENTRAL 98% DISTRIBUTIONS
#
# Visualization only.
# Data are NOT modified for analysis.
# ============================================================

gdsc_long_central <- gdsc_long %>%
  
  dplyr::group_by(
    Parameter
  ) %>%
  
  dplyr::filter(
    
    Value >=
      quantile(
        Value,
        0.01,
        na.rm = TRUE
      ),
    
    Value <=
      quantile(
        Value,
        0.99,
        na.rm = TRUE
      )
  ) %>%
  
  dplyr::ungroup()


p_central <- ggplot(
  gdsc_long_central,
  aes(
    x = Value
  )
) +
  
  geom_histogram(
    bins = 60
  ) +
  
  facet_wrap(
    ~ Parameter,
    scales = "free"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    
    title =
      "GDSC parameter distributions — central 98%",
    
    subtitle =
      "1st–99th percentile shown for visualization only",
    
    x =
      "Parameter value",
    
    y =
      "Number of experiments"
  )


print(p_central)


# ============================================================
# 12. EXTREME-TAIL SUMMARY
# ============================================================

extreme_summary <- gdsc_long %>%
  
  dplyr::group_by(
    Parameter
  ) %>%
  
  dplyr::summarise(
    
    N =
      dplyr::n(),
    
    Q01 =
      quantile(
        Value,
        0.01,
        na.rm = TRUE
      ),
    
    Q99 =
      quantile(
        Value,
        0.99,
        na.rm = TRUE
      ),
    
    N_below_Q01 =
      sum(
        Value <
          quantile(
            Value,
            0.01,
            na.rm = TRUE
          )
      ),
    
    N_above_Q99 =
      sum(
        Value >
          quantile(
            Value,
            0.99,
            na.rm = TRUE
          )
      ),
    
    .groups =
      "drop"
  )


cat("\n")
cat("============================================\n")
cat("EXTREME-TAIL SUMMARY\n")
cat("============================================\n")

print(
  as.data.frame(
    extreme_summary
  )
)


# ============================================================
# 13. ADD SIMPLE MAX / Q99 RATIO
#
# Very large ratios can flag extreme upper tails.
# Descriptive QC only.
# ============================================================

qc_all$Max_Q99_ratio <-
  qc_all$Max /
  qc_all$Q99


qc_hgsoc$Max_Q99_ratio <-
  qc_hgsoc$Max /
  qc_hgsoc$Q99


cat("\n")
cat("============================================\n")
cat("ALL GDSC — MAX / Q99\n")
cat("============================================\n")

print(
  qc_all[
    ,
    c(
      "Parameter",
      "Percent_NA",
      "Median",
      "Q95",
      "Q99",
      "Max",
      "Max_Q99_ratio"
    )
  ]
)


cat("\n")
cat("============================================\n")
cat("HGSOC — MAX / Q99\n")
cat("============================================\n")

print(
  qc_hgsoc[
    ,
    c(
      "Parameter",
      "Percent_NA",
      "Median",
      "Q95",
      "Q99",
      "Max",
      "Max_Q99_ratio"
    )
  ]
)


# ============================================================
# 14. SAVE RESULTS
# ============================================================

output_dir <-
  "~/Desktop/GSK_PK_State_Model"


if(!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


write.csv(
  
  qc_all,
  
  file.path(
    output_dir,
    "GDSC_parameter_QC_ALL.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  qc_hgsoc,
  
  file.path(
    output_dir,
    "GDSC_parameter_QC_HGSOC.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  extreme_summary,
  
  file.path(
    output_dir,
    "GDSC_parameter_extreme_tail_summary.csv"
  ),
  
  row.names = FALSE
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_parameter_distributions_RAW.pdf"
  ),
  
  p_raw,
  
  width = 11,
  height = 7
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_parameter_distributions_RAW.png"
  ),
  
  p_raw,
  
  width = 11,
  height = 7,
  dpi = 600
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_parameter_distributions_central98.pdf"
  ),
  
  p_central,
  
  width = 11,
  height = 7
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_parameter_distributions_central98.png"
  ),
  
  p_central,
  
  width = 11,
  height = 7,
  dpi = 600
)


# ============================================================
# 15. FINISHED
# ============================================================

cat("\n")
cat("============================================\n")
cat("QC COMPLETE\n")
cat("============================================\n")

cat(
  "Results saved to:\n",
  output_dir,
  "\n"
)

# ============================================================
# MULTI-PARAMETER GDSC DRUG SCREEN
#
# State2-likeness vs:
# 1. AAC
# 2. E_inf
# 3. EC50
# 4. HS
# 5. IC50_recomputed
#
# Requires:
#   gdsc_hgsoc_all
#   gdsc_hgsoc_match
#
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 1. ADD STATE INFORMATION IF NOT ALREADY PRESENT
# ============================================================

if(!"State2_Likeness" %in% colnames(gdsc_hgsoc_all)) {
  
  state_table <- gdsc_hgsoc_match %>%
    dplyr::select(
      Sample.Name,
      ModelID,
      State2_Likeness,
      State_Like
    )
  
  gdsc_hgsoc_all <- gdsc_hgsoc_all %>%
    dplyr::left_join(
      state_table,
      by = "Sample.Name"
    )
}


# ============================================================
# 2. EXTRACT DRUG NAME IF NEEDED
# ============================================================

if(!"Drug" %in% colnames(gdsc_hgsoc_all)) {
  
  extract_drug <- function(experiment, sample_name) {
    
    x <- sub(
      paste0("^", sample_name, "_"),
      "",
      experiment,
      ignore.case = TRUE
    )
    
    x <- sub(
      "_[0-9]+_.*$",
      "",
      x
    )
    
    return(x)
  }
  
  gdsc_hgsoc_all$Drug <- mapply(
    extract_drug,
    gdsc_hgsoc_all$Experiment,
    gdsc_hgsoc_all$Sample.Name,
    USE.NAMES = FALSE
  )
}


# ============================================================
# 3. COLLAPSE REPLICATES
#
# One value per:
# Drug × cell line
# ============================================================

drug_cell_5p <- gdsc_hgsoc_all %>%
  
  dplyr::group_by(
    Drug,
    ModelID,
    Sample.Name,
    State2_Likeness,
    State_Like
  ) %>%
  
  dplyr::summarise(
    
    AAC =
      if(all(is.na(aac_recomputed))) {
        NA_real_
      } else {
        median(
          aac_recomputed,
          na.rm = TRUE
        )
      },
    
    E_inf =
      if(all(is.na(E_inf))) {
        NA_real_
      } else {
        median(
          E_inf,
          na.rm = TRUE
        )
      },
    
    EC50 =
      if(all(is.na(EC50))) {
        NA_real_
      } else {
        median(
          EC50,
          na.rm = TRUE
        )
      },
    
    HS =
      if(all(is.na(HS))) {
        NA_real_
      } else {
        median(
          HS,
          na.rm = TRUE
        )
      },
    
    IC50 =
      if(all(is.na(ic50_recomputed))) {
        NA_real_
      } else {
        median(
          ic50_recomputed,
          na.rm = TRUE
        )
      },
    
    .groups = "drop"
  )


# ============================================================
# 4. REMOVE NON-FINITE VALUES
# ============================================================

for(v in c(
  "AAC",
  "E_inf",
  "EC50",
  "HS",
  "IC50"
)) {
  
  drug_cell_5p[[v]][
    !is.finite(
      drug_cell_5p[[v]]
    )
  ] <- NA
}


# ============================================================
# 5. SAFE SPEARMAN FUNCTION
# ============================================================

safe_cor <- function(x, y) {
  
  ok <-
    is.finite(x) &
    is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if(
    n < 6 ||
    length(unique(x)) < 3 ||
    length(unique(y)) < 3
  ) {
    
    return(
      data.frame(
        N = n,
        rho = NA_real_,
        P = NA_real_
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    N = n,
    rho = unname(test$estimate),
    P = test$p.value
  )
}


# ============================================================
# 6. GENERAL SCREEN FUNCTION
# ============================================================

screen_parameter <- function(data, parameter) {
  
  data %>%
    
    dplyr::group_by(Drug) %>%
    
    dplyr::group_modify(
      
      ~ safe_cor(
        .x$State2_Likeness,
        .x[[parameter]]
      )
    ) %>%
    
    dplyr::ungroup() %>%
    
    dplyr::mutate(
      
      FDR = p.adjust(
        P,
        method = "BH"
      ),
      
      Parameter = parameter
    )
}


# ============================================================
# 7. RUN ALL FIVE PARAMETERS
# ============================================================

res_AAC <-
  screen_parameter(
    drug_cell_5p,
    "AAC"
  )

res_Einf <-
  screen_parameter(
    drug_cell_5p,
    "E_inf"
  )

res_EC50 <-
  screen_parameter(
    drug_cell_5p,
    "EC50"
  )

res_HS <-
  screen_parameter(
    drug_cell_5p,
    "HS"
  )

res_IC50 <-
  screen_parameter(
    drug_cell_5p,
    "IC50"
  )


# ============================================================
# 8. COMBINE LONG-FORM RESULTS
# ============================================================

all_parameter_results <- dplyr::bind_rows(
  res_AAC,
  res_Einf,
  res_EC50,
  res_HS,
  res_IC50
)


# ============================================================
# 9. BIOLOGICAL DIRECTION
#
# AAC:
#   higher = greater sensitivity
#   rho > 0 => State2-sensitive
#
# E_inf:
#   lower residual viability = greater effect
#   rho < 0 => State2-sensitive
#
# EC50:
#   lower = greater sensitivity
#   rho < 0 => State2-sensitive
#
# IC50:
#   lower = greater sensitivity
#   rho < 0 => State2-sensitive
#
# HS:
#   curve slope/shape parameter
#   DO NOT classify directly as sensitivity.
# ============================================================

all_parameter_results <-
  all_parameter_results %>%
  
  dplyr::mutate(
    
    Direction =
      dplyr::case_when(
        
        Parameter == "AAC" &
          rho > 0 ~
          "State2-sensitive",
        
        Parameter == "AAC" &
          rho < 0 ~
          "State1-sensitive",
        
        
        Parameter == "E_inf" &
          rho < 0 ~
          "State2-sensitive",
        
        Parameter == "E_inf" &
          rho > 0 ~
          "State1-sensitive",
        
        
        Parameter == "EC50" &
          rho < 0 ~
          "State2-sensitive",
        
        Parameter == "EC50" &
          rho > 0 ~
          "State1-sensitive",
        
        
        Parameter == "IC50" &
          rho < 0 ~
          "State2-sensitive",
        
        Parameter == "IC50" &
          rho > 0 ~
          "State1-sensitive",
        
        
        Parameter == "HS" &
          rho > 0 ~
          "Positive slope association",
        
        Parameter == "HS" &
          rho < 0 ~
          "Negative slope association",
        
        TRUE ~
          NA_character_
      )
  )


# ============================================================
# 10. CREATE WIDE DRUG TABLE
# ============================================================

wide_results <-
  all_parameter_results %>%
  
  dplyr::select(
    Drug,
    Parameter,
    N,
    rho,
    P,
    FDR,
    Direction
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from =
      Parameter,
    
    values_from =
      c(
        N,
        rho,
        P,
        FDR,
        Direction
      ),
    
    names_sep = "_"
  )


# ============================================================
# 11. PRIMARY DIRECTIONAL SUPPORT
#
# Main biological endpoints:
# AAC
# E_inf
# EC50
#
# IC50 is retained separately because of QC problems.
# HS describes curve shape.
# ============================================================

wide_results <-
  wide_results %>%
  
  dplyr::mutate(
    
    AAC_State2 =
      rho_AAC > 0,
    
    Einf_State2 =
      rho_E_inf < 0,
    
    EC50_State2 =
      rho_EC50 < 0,
    
    IC50_State2 =
      rho_IC50 < 0,
    
    
    # Number of primary endpoints
    # supporting State2 sensitivity
    State2_support_primary =
      rowSums(
        cbind(
          AAC_State2,
          Einf_State2,
          EC50_State2
        ),
        na.rm = TRUE
      ),
    
    
    # Number supporting State1 sensitivity
    State1_support_primary =
      rowSums(
        cbind(
          rho_AAC < 0,
          rho_E_inf > 0,
          rho_EC50 > 0
        ),
        na.rm = TRUE
      )
  )


# ============================================================
# 12. NUMBER OF AVAILABLE PRIMARY ENDPOINTS
# ============================================================

wide_results <-
  wide_results %>%
  
  dplyr::mutate(
    
    Primary_available =
      rowSums(
        cbind(
          !is.na(rho_AAC),
          !is.na(rho_E_inf),
          !is.na(rho_EC50)
        )
      )
  )


# ============================================================
# 13. PRIMARY CONCORDANCE CLASS
# ============================================================

wide_results <-
  wide_results %>%
  
  dplyr::mutate(
    
    Primary_concordance =
      dplyr::case_when(
        
        Primary_available == 3 &
          State2_support_primary == 3 ~
          
          "3/3 State2-sensitive",
        
        
        Primary_available == 3 &
          State1_support_primary == 3 ~
          
          "3/3 State1-sensitive",
        
        
        Primary_available >= 2 &
          State2_support_primary ==
          Primary_available ~
          
          "Concordant State2-sensitive",
        
        
        Primary_available >= 2 &
          State1_support_primary ==
          Primary_available ~
          
          "Concordant State1-sensitive",
        
        
        Primary_available >= 2 ~
          
          "Discordant",
        
        
        TRUE ~
          
          "Insufficient"
      )
  )


# ============================================================
# 14. EFFECT-SIZE SCORE
#
# AAC:   rho
# E_inf: -rho
# EC50:  -rho
#
# Positive score =
# overall State2-sensitivity direction
#
# IC50 NOT included in primary score.
# HS NOT included in primary score.
# ============================================================

wide_results <-
  wide_results %>%
  
  dplyr::rowwise() %>%
  
  dplyr::mutate(
    
    State2_effect_score =
      
      mean(
        c(
          rho_AAC,
          -rho_E_inf,
          -rho_EC50
        ),
        na.rm = TRUE
      )
  ) %>%
  
  dplyr::ungroup()


# Replace NaN with NA
wide_results$State2_effect_score[
  !is.finite(
    wide_results$State2_effect_score
  )
] <- NA


# ============================================================
# 15. MINIMUM SAMPLE SIZE
#
# Require N >= 6 for primary parameters used.
# ============================================================

wide_results <-
  wide_results %>%
  
  dplyr::mutate(
    
    N_OK_AAC =
      !is.na(N_AAC) &
      N_AAC >= 6,
    
    N_OK_Einf =
      !is.na(N_E_inf) &
      N_E_inf >= 6,
    
    N_OK_EC50 =
      !is.na(N_EC50) &
      N_EC50 >= 6
  )


# ============================================================
# 16. HIGH-CONFIDENCE DIRECTIONAL CANDIDATES
#
# This is exploratory prioritization,
# NOT statistical significance.
# ============================================================

state2_candidates_5p <-
  wide_results %>%
  
  dplyr::filter(
    
    Primary_concordance %in%
      c(
        "3/3 State2-sensitive",
        "Concordant State2-sensitive"
      ),
    
    State2_effect_score > 0
    
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      State2_effect_score
    )
  )


state1_candidates_5p <-
  wide_results %>%
  
  dplyr::filter(
    
    Primary_concordance %in%
      c(
        "3/3 State1-sensitive",
        "Concordant State1-sensitive"
      ),
    
    State2_effect_score < 0
    
  ) %>%
  
  dplyr::arrange(
    State2_effect_score
  )


# ============================================================
# 17. SHOW TOP STATE2 CANDIDATES
# ============================================================

cat("\n")
cat("============================================\n")
cat("TOP MULTI-PARAMETER STATE2 CANDIDATES\n")
cat("============================================\n")

print(
  
  head(
    
    state2_candidates_5p %>%
      
      dplyr::select(
        
        Drug,
        
        Primary_concordance,
        
        State2_effect_score,
        
        N_AAC,
        rho_AAC,
        P_AAC,
        FDR_AAC,
        
        N_E_inf,
        rho_E_inf,
        P_E_inf,
        FDR_E_inf,
        
        N_EC50,
        rho_EC50,
        P_EC50,
        FDR_EC50,
        
        N_HS,
        rho_HS,
        
        N_IC50,
        rho_IC50
      ),
    
    30
  )
)


# ============================================================
# 18. SHOW TOP STATE1 CANDIDATES
# ============================================================

cat("\n")
cat("============================================\n")
cat("TOP MULTI-PARAMETER STATE1 CANDIDATES\n")
cat("============================================\n")

print(
  
  head(
    
    state1_candidates_5p %>%
      
      dplyr::select(
        
        Drug,
        
        Primary_concordance,
        
        State2_effect_score,
        
        rho_AAC,
        rho_E_inf,
        rho_EC50,
        rho_HS,
        rho_IC50
      ),
    
    30
  )
)


# ============================================================
# 19. FDR SIGNIFICANCE SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("FDR < 0.05 BY PARAMETER\n")
cat("============================================\n")

fdr_summary <-
  all_parameter_results %>%
  
  dplyr::group_by(
    Parameter
  ) %>%
  
  dplyr::summarise(
    
    Drugs_tested =
      sum(!is.na(P)),
    
    FDR_significant =
      sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(fdr_summary)


# ============================================================
# 20. SHOW ALL FDR-SIGNIFICANT ASSOCIATIONS
# ============================================================

significant_results <-
  all_parameter_results %>%
  
  dplyr::filter(
    FDR < 0.05
  ) %>%
  
  dplyr::arrange(
    Parameter,
    FDR
  )


cat("\n")
cat("============================================\n")
cat("ALL FDR-SIGNIFICANT ASSOCIATIONS\n")
cat("============================================\n")

print(significant_results)


# ============================================================
# 21. HEATMAP-LIKE DATA FRAME
#
# Direction normalized:
#
# Positive = State2-sensitive direction
# Negative = State1-sensitive direction
#
# HS retained as raw rho because it is not
# directly sensitivity direction.
# ============================================================

drug_fingerprint <-
  wide_results %>%
  
  dplyr::transmute(
    
    Drug,
    
    AAC =
      rho_AAC,
    
    E_inf =
      -rho_E_inf,
    
    EC50 =
      -rho_EC50,
    
    HS =
      rho_HS,
    
    IC50 =
      -rho_IC50,
    
    State2_effect_score
  )


# ============================================================
# 22. LONG FORMAT FOR PLOT
# ============================================================

fingerprint_long <-
  drug_fingerprint %>%
  
  dplyr::select(
    Drug,
    AAC,
    E_inf,
    EC50,
    HS,
    IC50
  ) %>%
  
  tidyr::pivot_longer(
    
    cols =
      c(
        AAC,
        E_inf,
        EC50,
        HS,
        IC50
      ),
    
    names_to =
      "Parameter",
    
    values_to =
      "Directional_rho"
  )


# ============================================================
# 23. TOP 20 STATE2 CANDIDATE FINGERPRINT PLOT
# ============================================================

top20_drugs <-
  state2_candidates_5p$Drug[
    seq_len(
      min(
        20,
        nrow(
          state2_candidates_5p
        )
      )
    )
  ]


plot_df <-
  fingerprint_long %>%
  
  dplyr::filter(
    Drug %in%
      top20_drugs
  )


plot_df$Drug <-
  factor(
    plot_df$Drug,
    levels =
      rev(
        top20_drugs
      )
  )


p_fingerprint <-
  ggplot(
    plot_df,
    aes(
      x = Parameter,
      y = Drug,
      fill = Directional_rho
    )
  ) +
  
  geom_tile() +
  
  geom_text(
    aes(
      label =
        ifelse(
          is.na(Directional_rho),
          "",
          sprintf(
            "%.2f",
            Directional_rho
          )
        )
    ),
    size = 3
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    
    title =
      "Multi-parameter drug-response fingerprints",
    
    subtitle =
      "Positive values indicate State2-sensitivity direction for AAC, E_inf, EC50 and IC50; HS represents slope association",
    
    x =
      "Drug-response parameter",
    
    y =
      NULL,
    
    fill =
      "Directional\nSpearman rho"
  )


print(p_fingerprint)


# ============================================================
# 24. SAVE RESULTS
# ============================================================

output_dir <-
  "~/Desktop/GSK_PK_State_Model"


if(!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


write.csv(
  
  all_parameter_results,
  
  file.path(
    output_dir,
    "GDSC_5parameter_all_associations.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  wide_results,
  
  file.path(
    output_dir,
    "GDSC_5parameter_drug_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  state2_candidates_5p,
  
  file.path(
    output_dir,
    "GDSC_5parameter_State2_candidates.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  state1_candidates_5p,
  
  file.path(
    output_dir,
    "GDSC_5parameter_State1_candidates.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  significant_results,
  
  file.path(
    output_dir,
    "GDSC_5parameter_FDR_significant.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  drug_fingerprint,
  
  file.path(
    output_dir,
    "GDSC_5parameter_fingerprint.csv"
  ),
  
  row.names = FALSE
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_5parameter_State2_fingerprint.pdf"
  ),
  
  p_fingerprint,
  
  width = 9,
  height = 9
)


ggsave(
  
  file.path(
    output_dir,
    "GDSC_5parameter_State2_fingerprint.png"
  ),
  
  p_fingerprint,
  
  width = 9,
  height = 9,
  dpi = 600
)


# ============================================================
# 25. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("FINAL SUMMARY\n")
cat("============================================\n")

cat(
  "Total drugs:",
  nrow(
    wide_results
  ),
  "\n"
)

cat(
  "3/3 State2-sensitive:",
  sum(
    wide_results$Primary_concordance ==
      "3/3 State2-sensitive",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "3/3 State1-sensitive:",
  sum(
    wide_results$Primary_concordance ==
      "3/3 State1-sensitive",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Any FDR < 0.05 associations:",
  nrow(
    significant_results
  ),
  "\n"
)

cat("\nAnalysis complete.\n")
state2_candidates_5p %>%
  dplyr::select(
    Drug,
    Primary_concordance,
    State2_effect_score,
    N_AAC, rho_AAC, P_AAC, FDR_AAC,
    N_E_inf, rho_E_inf, P_E_inf, FDR_E_inf,
    N_EC50, rho_EC50, P_EC50, FDR_EC50,
    N_HS, rho_HS, P_HS, FDR_HS,
    N_IC50, rho_IC50, P_IC50, FDR_IC50
  ) %>%
  head(20) %>%
  print()
gdc0810 <- state2_candidates_5p %>%
  dplyr::filter(Drug == "GDC0810")

data.frame(
  Parameter = c("AAC", "E_inf", "EC50", "HS", "IC50"),
  N = c(
    gdc0810$N_AAC,
    gdc0810$N_E_inf,
    gdc0810$N_EC50,
    gdc0810$N_HS,
    gdc0810$N_IC50
  ),
  rho = c(
    gdc0810$rho_AAC,
    gdc0810$rho_E_inf,
    gdc0810$rho_EC50,
    gdc0810$rho_HS,
    gdc0810$rho_IC50
  ),
  P = c(
    gdc0810$P_AAC,
    gdc0810$P_E_inf,
    gdc0810$P_EC50,
    gdc0810$P_HS,
    gdc0810$P_IC50
  ),
  FDR = c(
    gdc0810$FDR_AAC,
    gdc0810$FDR_E_inf,
    gdc0810$FDR_EC50,
    gdc0810$FDR_HS,
    gdc0810$FDR_IC50
  )
)
# ============================================================
# FULL RESULTS FOR ALL DRUGS SHOWN IN THE HEATMAP
# ============================================================

# Heatmapte kullanılan ilaçlar
heatmap_drugs <- unique(as.character(plot_df$Drug))

# Uzun format:
# Her drug için 5 parameter ayrı satır
heatmap_full_results <- all_parameter_results %>%
  dplyr::filter(
    Drug %in% heatmap_drugs
  ) %>%
  dplyr::select(
    Drug,
    Parameter,
    N,
    rho,
    P,
    FDR,
    Direction
  ) %>%
  dplyr::arrange(
    match(Drug, heatmap_drugs),
    match(
      Parameter,
      c(
        "AAC",
        "E_inf",
        "EC50",
        "HS",
        "IC50"
      )
    )
  )

cat("\n============================================\n")
cat("ALL HEATMAP DRUGS — ALL 5 PARAMETERS\n")
cat("============================================\n")

print(
  as.data.frame(heatmap_full_results),
  row.names = FALSE
)


# ============================================================
# CREATE COMPACT SUMMARY
# One row per drug
# ============================================================

heatmap_summary <- wide_results %>%
  dplyr::filter(
    Drug %in% heatmap_drugs
  ) %>%
  dplyr::arrange(
    match(Drug, heatmap_drugs)
  ) %>%
  dplyr::select(
    Drug,
    Primary_concordance,
    State2_effect_score,
    
    N_AAC,
    rho_AAC,
    P_AAC,
    FDR_AAC,
    
    N_E_inf,
    rho_E_inf,
    P_E_inf,
    FDR_E_inf,
    
    N_EC50,
    rho_EC50,
    P_EC50,
    FDR_EC50,
    
    N_HS,
    rho_HS,
    P_HS,
    FDR_HS,
    
    N_IC50,
    rho_IC50,
    P_IC50,
    FDR_IC50
  )

cat("\n============================================\n")
cat("COMPACT DRUG SUMMARY\n")
cat("============================================\n")

print(
  as.data.frame(heatmap_summary),
  row.names = FALSE
)


# ============================================================
# COUNT DIRECTIONAL SUPPORT
#
# Main sensitivity endpoints:
# AAC   > 0 = State2 sensitive
# E_inf < 0 = State2 sensitive
# EC50  < 0 = State2 sensitive
#
# HS is NOT counted as sensitivity.
# IC50 shown separately because of QC problems.
# ============================================================

comparison_table <- heatmap_summary %>%
  dplyr::mutate(
    
    AAC_support =
      !is.na(rho_AAC) &
      rho_AAC > 0,
    
    E_inf_support =
      !is.na(rho_E_inf) &
      rho_E_inf < 0,
    
    EC50_support =
      !is.na(rho_EC50) &
      rho_EC50 < 0,
    
    IC50_support =
      !is.na(rho_IC50) &
      rho_IC50 < 0,
    
    Main_support_N =
      rowSums(
        cbind(
          AAC_support,
          E_inf_support,
          EC50_support
        ),
        na.rm = TRUE
      ),
    
    Main_available_N =
      rowSums(
        cbind(
          !is.na(rho_AAC),
          !is.na(rho_E_inf),
          !is.na(rho_EC50)
        )
      ),
    
    Min_P_primary =
      apply(
        cbind(
          P_AAC,
          P_E_inf,
          P_EC50
        ),
        1,
        function(x) {
          if(all(is.na(x))) {
            NA_real_
          } else {
            min(x, na.rm = TRUE)
          }
        }
      ),
    
    Min_FDR_primary =
      apply(
        cbind(
          FDR_AAC,
          FDR_E_inf,
          FDR_EC50
        ),
        1,
        function(x) {
          if(all(is.na(x))) {
            NA_real_
          } else {
            min(x, na.rm = TRUE)
          }
        }
      )
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(Main_support_N),
    dplyr::desc(State2_effect_score)
  )


cat("\n============================================\n")
cat("FINAL COMPARISON OF HEATMAP DRUGS\n")
cat("============================================\n")

print(
  as.data.frame(
    comparison_table %>%
      dplyr::select(
        Drug,
        Main_support_N,
        Main_available_N,
        State2_effect_score,
        rho_AAC,
        rho_E_inf,
        rho_EC50,
        rho_HS,
        rho_IC50,
        Min_P_primary,
        Min_FDR_primary
      )
  ),
  row.names = FALSE
)


# ============================================================
# SAVE ALL TABLES
# ============================================================

output_dir <- "~/Desktop/GSK_PK_State_Model"

write.csv(
  heatmap_full_results,
  file.path(
    output_dir,
    "Heatmap_drugs_all_5_parameters_LONG.csv"
  ),
  row.names = FALSE
)

write.csv(
  heatmap_summary,
  file.path(
    output_dir,
    "Heatmap_drugs_all_5_parameters_WIDE.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparison_table,
  file.path(
    output_dir,
    "Heatmap_drugs_final_comparison.csv"
  ),
  row.names = FALSE
)

cat("\nDone.\n")

# ============================================================
# VISUAL SUMMARY OF ALL HEATMAP DRUGS
# Effect size + direction + statistical evidence
# ============================================================

library(dplyr)
library(ggplot2)

# Heatmapte bulunan ilaçlar
heatmap_drugs <- unique(as.character(plot_df$Drug))


# ============================================================
# 1. Prepare plotting table
# ============================================================

visual_df <- all_parameter_results %>%
  
  dplyr::filter(
    Drug %in% heatmap_drugs
  ) %>%
  
  dplyr::mutate(
    
    # Convert correlations into common biological direction
    #
    # Positive = State2-sensitive direction
    # Negative = State1-sensitive direction
    #
    # HS is kept as raw association because slope is not
    # directly a sensitivity endpoint.
    
    Directional_rho =
      dplyr::case_when(
        
        Parameter == "AAC"  ~ rho,
        
        Parameter == "E_inf" ~ -rho,
        
        Parameter == "EC50" ~ -rho,
        
        Parameter == "IC50" ~ -rho,
        
        Parameter == "HS" ~ rho,
        
        TRUE ~ NA_real_
      ),
    
    # P-value evidence for point size
    logP =
      ifelse(
        !is.na(P) & P > 0,
        -log10(P),
        NA_real_
      )
  )


# ============================================================
# 2. Order parameters
# ============================================================

visual_df$Parameter <- factor(
  visual_df$Parameter,
  levels = c(
    "AAC",
    "E_inf",
    "EC50",
    "HS",
    "IC50"
  )
)


# ============================================================
# 3. Order drugs by overall State2 effect score
# ============================================================

drug_order <- wide_results %>%
  
  dplyr::filter(
    Drug %in% heatmap_drugs
  ) %>%
  
  dplyr::arrange(
    State2_effect_score
  ) %>%
  
  dplyr::pull(Drug)


visual_df$Drug <- factor(
  visual_df$Drug,
  levels = drug_order
)


# ============================================================
# 4. Main visualization
# ============================================================

p_multi <- ggplot(
  visual_df,
  aes(
    x = Parameter,
    y = Drug
  )
) +
  
  # Background tile = direction/effect size
  geom_tile(
    aes(
      fill = Directional_rho
    ),
    linewidth = 0.5
  ) +
  
  # Circle = strength of nominal statistical evidence
  geom_point(
    aes(
      size = logP
    ),
    shape = 21,
    fill = "white",
    alpha = 0.8
  ) +
  
  # Show directional rho
  geom_text(
    aes(
      label =
        ifelse(
          is.na(Directional_rho),
          "",
          sprintf(
            "%.2f",
            Directional_rho
          )
        )
    ),
    size = 3
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    name = "Directional\nSpearman rho"
  ) +
  
  scale_size_continuous(
    range = c(1, 8),
    name = "-log10(P)"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    axis.text.y =
      element_text(
        size = 10
      ),
    
    legend.position =
      "right"
  ) +
  
  labs(
    
    title =
      "Multi-parameter drug sensitivity associated with State2-likeness",
    
    subtitle =
      "Positive directional rho = State2-sensitive direction; circle size = nominal statistical evidence",
    
    x =
      "Drug-response parameter",
    
    y =
      NULL
  )


print(p_multi)


# ============================================================
# 5. Save
# ============================================================

output_dir <- "~/Desktop/GSK_PK_State_Model"

ggsave(
  file.path(
    output_dir,
    "GDSC_State2_multparameter_summary.pdf"
  ),
  p_multi,
  width = 10,
  height = 8
)

ggsave(
  file.path(
    output_dir,
    "GDSC_State2_multparameter_summary.png"
  ),
  p_multi,
  width = 10,
  height = 8,
  dpi = 600
)
##Bayesian
library(brms)

test_df <- data.frame(
  y = rnorm(10),
  x = rnorm(10)
)

test_fit <- brm(
  y ~ x,
  data = test_df,
  chains = 2,
  iter = 1000,
  warmup = 500,
  cores = 2,
  seed = 123
)

summary(test_fit)
