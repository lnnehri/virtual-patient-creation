# ============================================================
# GSK-1904529A — EXPLORATORY PK/PD MODEL
# ============================================================

library(deSolve)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Published preclinical information
# ------------------------------------------------------------

Dose <- 30          # mg/kg
F    <- 0.45        # oral bioavailability
t_half <- 5.2       # hours

# Elimination rate
ke <- log(2) / t_half

# Ka and V are not identified from the information currently
# available to us.
# Use explicit exploratory assumptions for simulation.

Ka <- 1             # 1/h, assumed
V  <- 1             # L/kg, assumed apparent volume

cat("Elimination rate ke =", ke, "1/h\n")


# ------------------------------------------------------------
# 2. One-compartment oral PK model
# ------------------------------------------------------------

pk_model <- function(time, state, parameters) {
  
  with(
    as.list(c(state, parameters)),
    {
      
      dAgut <- -Ka * Agut
      
      dAcentral <- F * Ka * Agut -
        ke * Acentral
      
      list(c(dAgut, dAcentral))
    }
  )
}


# ------------------------------------------------------------
# 3. Initial conditions
# ------------------------------------------------------------

state <- c(
  Agut = Dose,
  Acentral = 0
)

parameters <- c(
  Ka = Ka,
  ke = ke,
  F = F
)

times <- seq(
  0,
  48,
  by = 0.1
)


# ------------------------------------------------------------
# 4. Simulate concentration-time profile
# ------------------------------------------------------------

pk_sim <- as.data.frame(
  ode(
    y = state,
    times = times,
    func = pk_model,
    parms = parameters
  )
)

pk_sim$Concentration <- pk_sim$Acentral / V


# ------------------------------------------------------------
# 5. Basic PK metrics
# ------------------------------------------------------------

Cmax <- max(pk_sim$Concentration)

Tmax <- pk_sim$time[
  which.max(pk_sim$Concentration)
]

AUC <- sum(
  diff(pk_sim$time) *
    (
      head(pk_sim$Concentration, -1) +
        tail(pk_sim$Concentration, -1)
    ) / 2
)

cat("\n============================\n")
cat("EXPLORATORY PK METRICS\n")
cat("============================\n")

cat("Cmax =", Cmax, "\n")
cat("Tmax =", Tmax, "hours\n")
cat("AUC0-48 =", AUC, "\n")
cat("Half-life =", t_half, "hours\n")


# ------------------------------------------------------------
# 6. Exploratory PK -> PD model
#
# Published preclinical work suggests strong IGF1R
# phosphorylation inhibition around ~3 uM exposure.
#
# This is NOT a fitted IC50.
# Use 3 uM only as an exploratory PD reference.
# ------------------------------------------------------------

EC50_PD <- 3

pk_sim$Target_Inhibition <-
  pk_sim$Concentration /
  (EC50_PD + pk_sim$Concentration)

pk_sim$Target_Inhibition_percent <-
  100 * pk_sim$Target_Inhibition


# ------------------------------------------------------------
# 7. Concentration-time plot
# ------------------------------------------------------------

p_pk <- ggplot(
  pk_sim,
  aes(
    x = time,
    y = Concentration
  )
) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 13) +
  labs(
    title = "Exploratory GSK-1904529A PK simulation",
    subtitle = "One-compartment oral model",
    x = "Time (hours)",
    y = "Simulated concentration"
  )

print(p_pk)


# ------------------------------------------------------------
# 8. PK-PD plot
# ------------------------------------------------------------

p_pd <- ggplot(
  pk_sim,
  aes(
    x = time,
    y = Target_Inhibition_percent
  )
) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 13) +
  labs(
    title = "Exploratory GSK-1904529A PK-PD simulation",
    x = "Time (hours)",
    y = "Predicted target inhibition (%)"
  )

print(p_pd)
# ============================================================
# GSK-1904529A
# EXPLORATORY PK × MOLECULAR-STATE RESPONSE MODEL
# ============================================================

library(deSolve)
library(ggplot2)
library(dplyr)


# ============================================================
# 1. PRECLINICAL PK PARAMETERS
# ============================================================

Dose   <- 30       # mg/kg
F      <- 0.45     # oral bioavailability
t_half <- 5.2      # hours

# Elimination rate constant
ke <- log(2) / t_half

# Exploratory assumptions
Ka <- 1            # absorption rate constant, 1/h
V  <- 1            # apparent volume, L/kg

cat(
  "Elimination rate ke =",
  ke,
  "1/h\n"
)


# ============================================================
# 2. ONE-COMPARTMENT ORAL PK MODEL
# ============================================================

pk_model <- function(time, state, parameters) {
  
  with(
    as.list(c(state, parameters)),
    {
      
      dAgut <- -Ka * Agut
      
      dAcentral <-
        F * Ka * Agut -
        ke * Acentral
      
      list(
        c(
          dAgut,
          dAcentral
        )
      )
    }
  )
}


# ============================================================
# 3. INITIAL CONDITIONS
# ============================================================

state <- c(
  Agut = Dose,
  Acentral = 0
)

parameters <- c(
  Ka = Ka,
  ke = ke,
  F = F
)

times <- seq(
  0,
  48,
  by = 0.1
)


# ============================================================
# 4. SIMULATE PK PROFILE
# ============================================================

pk_sim <- as.data.frame(
  deSolve::ode(
    y = state,
    times = times,
    func = pk_model,
    parms = parameters
  )
)

pk_sim$Concentration <-
  pk_sim$Acentral / V


# ============================================================
# 5. BASIC PK METRICS
# ============================================================

Cmax <- max(
  pk_sim$Concentration,
  na.rm = TRUE
)

Tmax <- pk_sim$time[
  which.max(pk_sim$Concentration)
]

AUC <- sum(
  diff(pk_sim$time) *
    (
      head(pk_sim$Concentration, -1) +
        tail(pk_sim$Concentration, -1)
    ) / 2
)

cat("\n============================\n")
cat("EXPLORATORY PK METRICS\n")
cat("============================\n")

cat(
  "Cmax =",
  round(Cmax, 4),
  "\n"
)

cat(
  "Tmax =",
  round(Tmax, 2),
  "hours\n"
)

cat(
  "AUC0-48 =",
  round(AUC, 4),
  "\n"
)

cat(
  "Half-life =",
  t_half,
  "hours\n"
)


# ============================================================
# 6. NORMALIZED EXPOSURE
# ============================================================

# IMPORTANT:
# Concentration is generated from an exploratory PK model.
# Therefore we normalize exposure rather than directly mixing
# these concentrations with experimental pharmacodynamic units.

pk_sim$Exposure_norm <-
  pk_sim$Concentration /
  max(
    pk_sim$Concentration,
    na.rm = TRUE
  )

summary(pk_sim$Exposure_norm)


# ============================================================
# 7. PLOT PK PROFILE
# ============================================================

p_pk <- ggplot(
  pk_sim,
  aes(
    x = time,
    y = Concentration
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Exploratory GSK-1904529A PK simulation",
    subtitle =
      "One-compartment oral model",
    x =
      "Time after dose (hours)",
    y =
      "Simulated concentration"
  )

print(p_pk)


# ============================================================
# 8. PLOT NORMALIZED EXPOSURE
# ============================================================

p_exposure <- ggplot(
  pk_sim,
  aes(
    x = time,
    y = Exposure_norm
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Normalized GSK-1904529A exposure",
    x =
      "Time after dose (hours)",
    y =
      "Normalized exposure"
  )

print(p_exposure)


# ============================================================
# 9. GET OBSERVED STATE2-LIKENESS SCORES
# ============================================================

# Uses the GSK-tested HGSOC cell lines from the previous
# GDSC and PRISM analyses.

all_state_scores <- c(
  gsk_gdsc_unique$State2_Likeness,
  gsk_prism_unique$State2_Likeness
)

all_state_scores <-
  all_state_scores[
    is.finite(all_state_scores)
  ]

state_min <- min(
  all_state_scores,
  na.rm = TRUE
)

state_max <- max(
  all_state_scores,
  na.rm = TRUE
)

state_mean <- mean(
  all_state_scores,
  na.rm = TRUE
)

state_sd <- sd(
  all_state_scores,
  na.rm = TRUE
)

cat("\n============================\n")
cat("STATE2-LIKENESS\n")
cat("============================\n")

cat(
  "Minimum =",
  state_min,
  "\n"
)

cat(
  "Maximum =",
  state_max,
  "\n"
)

cat(
  "Mean =",
  state_mean,
  "\n"
)

cat(
  "SD =",
  state_sd,
  "\n"
)


# ============================================================
# 10. CREATE CONTINUOUS STATE2-LIKENESS GRID
# ============================================================

state_grid <- seq(
  state_min,
  state_max,
  length.out = 100
)


# ============================================================
# 11. PK × MOLECULAR-STATE GRID
# ============================================================

pk_state_grid <- expand.grid(
  time = pk_sim$time,
  State2_Likeness = state_grid
)

pk_state_grid <- pk_state_grid %>%
  dplyr::left_join(
    pk_sim %>%
      dplyr::select(
        time,
        Concentration,
        Exposure_norm
      ),
    by = "time"
  )


# ============================================================
# 12. STANDARDIZE STATE2-LIKENESS
# ============================================================

pk_state_grid$State2_z <-
  (
    pk_state_grid$State2_Likeness -
      state_mean
  ) /
  state_sd


# ============================================================
# 13. EXPLORATORY EXPOSURE × STATE MODEL
# ============================================================

# IMPORTANT:
#
# These coefficients are NOT fitted clinical parameters.
#
# They are used only to explore the conceptual interaction
# between:
#
# drug exposure
# +
# molecular State2-likeness
# +
# exposure × molecular-state interaction


beta_exposure <- 3

beta_state <- 0.5

beta_interaction <- 1


# Linear predictor

pk_state_grid$Response_linear <-
  beta_exposure *
  pk_state_grid$Exposure_norm +
  
  beta_state *
  pk_state_grid$State2_z +
  
  beta_interaction *
  pk_state_grid$Exposure_norm *
  pk_state_grid$State2_z


# Logistic transformation

pk_state_grid$Predicted_Response <-
  plogis(
    pk_state_grid$Response_linear
  )


pk_state_grid$Predicted_Response_percent <-
  100 *
  pk_state_grid$Predicted_Response


# ============================================================
# 14. DEFINE REPRESENTATIVE STATE1-LIKE / STATE2-LIKE PROFILES
# ============================================================

state_low <- as.numeric(
  quantile(
    all_state_scores,
    0.25,
    na.rm = TRUE
  )
)

state_high <- as.numeric(
  quantile(
    all_state_scores,
    0.75,
    na.rm = TRUE
  )
)

cat("\n============================\n")
cat("REPRESENTATIVE PROFILES\n")
cat("============================\n")

cat(
  "State1-like reference =",
  state_low,
  "\n"
)

cat(
  "State2-like reference =",
  state_high,
  "\n"
)


representative_states <- data.frame(
  Group = c(
    "State1-like",
    "State2-like"
  ),
  State2_Likeness = c(
    state_low,
    state_high
  )
)


# ============================================================
# 15. CREATE RESPONSE CURVES
# ============================================================

response_curves <- expand.grid(
  time = pk_sim$time,
  Group = representative_states$Group
)

response_curves <- response_curves %>%
  dplyr::left_join(
    representative_states,
    by = "Group"
  ) %>%
  dplyr::left_join(
    pk_sim %>%
      dplyr::select(
        time,
        Concentration,
        Exposure_norm
      ),
    by = "time"
  )


# Standardize State2-likeness

response_curves$State2_z <-
  (
    response_curves$State2_Likeness -
      state_mean
  ) /
  state_sd


# Response model

response_curves$Response_linear <-
  beta_exposure *
  response_curves$Exposure_norm +
  
  beta_state *
  response_curves$State2_z +
  
  beta_interaction *
  response_curves$Exposure_norm *
  response_curves$State2_z


response_curves$Predicted_Response <-
  plogis(
    response_curves$Response_linear
  )


response_curves$Predicted_Response_percent <-
  100 *
  response_curves$Predicted_Response


# ============================================================
# 16. STATE-DEPENDENT RESPONSE OVER TIME
# ============================================================

p_state_response <- ggplot(
  response_curves,
  aes(
    x = time,
    y = Predicted_Response_percent,
    linetype = Group
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Exploratory state-dependent GSK-1904529A response",
    subtitle =
      "PK exposure integrated with molecular State2-likeness",
    x =
      "Time after dose (hours)",
    y =
      "Simulated response (%)",
    linetype =
      "Molecular profile"
  )

print(p_state_response)


# ============================================================
# 17. EXPOSURE × MOLECULAR-STATE RESPONSE LANDSCAPE
# ============================================================

p_response_map <- ggplot(
  pk_state_grid,
  aes(
    x = State2_Likeness,
    y = time,
    fill = Predicted_Response_percent
  )
) +
  geom_raster() +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "GSK-1904529A exposure × molecular-state response landscape",
    x =
      "State2-likeness",
    y =
      "Time after dose (hours)",
    fill =
      "Simulated\nresponse (%)"
  )

print(p_response_map)


# ============================================================
# 18. RESPONSE AT MAXIMUM EXPOSURE
# ============================================================

tmax_value <- pk_sim$time[
  which.max(
    pk_sim$Concentration
  )
]

response_at_tmax <- pk_state_grid[
  abs(
    pk_state_grid$time -
      tmax_value
  ) < 1e-8,
]


p_tmax <- ggplot(
  response_at_tmax,
  aes(
    x = State2_Likeness,
    y = Predicted_Response_percent
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Molecular-state effect at maximum GSK exposure",
    subtitle =
      paste0(
        "Tmax ≈ ",
        round(
          tmax_value,
          2
        ),
        " h"
      ),
    x =
      "State2-likeness",
    y =
      "Simulated response (%)"
  )

print(p_tmax)


# ============================================================
# 19. RESPONSE ACROSS EXPOSURE LEVELS
# ============================================================

exposure_grid <- seq(
  0,
  1,
  length.out = 100
)

exposure_state_plot <- expand.grid(
  Exposure_norm = exposure_grid,
  Group = representative_states$Group
)

exposure_state_plot <- exposure_state_plot %>%
  dplyr::left_join(
    representative_states,
    by = "Group"
  )

exposure_state_plot$State2_z <-
  (
    exposure_state_plot$State2_Likeness -
      state_mean
  ) /
  state_sd


exposure_state_plot$Response_linear <-
  beta_exposure *
  exposure_state_plot$Exposure_norm +
  
  beta_state *
  exposure_state_plot$State2_z +
  
  beta_interaction *
  exposure_state_plot$Exposure_norm *
  exposure_state_plot$State2_z


exposure_state_plot$Predicted_Response_percent <-
  100 *
  plogis(
    exposure_state_plot$Response_linear
  )


p_exposure_response <- ggplot(
  exposure_state_plot,
  aes(
    x = Exposure_norm,
    y = Predicted_Response_percent,
    linetype = Group
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title =
      "Exposure-response relationship by molecular profile",
    x =
      "Normalized GSK-1904529A exposure",
    y =
      "Simulated response (%)",
    linetype =
      "Molecular profile"
  )

print(p_exposure_response)


# ============================================================
# 20. SUMMARY TABLE
# ============================================================

simulation_summary <- data.frame(
  Parameter = c(
    "Dose_mg_kg",
    "Oral_bioavailability",
    "Half_life_hours",
    "Elimination_rate",
    "Assumed_Ka",
    "Assumed_V",
    "Simulated_Cmax",
    "Simulated_Tmax_hours",
    "Simulated_AUC_0_48",
    "State2_score_min",
    "State2_score_max",
    "State1_like_reference",
    "State2_like_reference"
  ),
  
  Value = c(
    Dose,
    F,
    t_half,
    ke,
    Ka,
    V,
    Cmax,
    Tmax,
    AUC,
    state_min,
    state_max,
    state_low,
    state_high
  )
)

print(simulation_summary)


# ============================================================
# 21. FINAL CHECK
# ============================================================

cat("\n========================================\n")
cat("PK × MOLECULAR STATE MODEL COMPLETE\n")
cat("========================================\n")

cat(
  "Dose:",
  Dose,
  "mg/kg\n"
)

cat(
  "Half-life:",
  t_half,
  "hours\n"
)

cat(
  "Simulated Tmax:",
  round(
    Tmax,
    2
  ),
  "hours\n"
)

cat(
  "State2-likeness range:",
  round(
    state_min,
    3
  ),
  "to",
  round(
    state_max,
    3
  ),
  "\n"
)

cat(
  "State1-like reference:",
  round(
    state_low,
    3
  ),
  "\n"
)

cat(
  "State2-like reference:",
  round(
    state_high,
    3
  ),
  "\n"
)

cat(
  "\nIMPORTANT:\n",
  "Response percentages are exploratory simulations,\n",
  "not experimentally fitted or clinically predicted responses.\n"
)
# ============================================================
# 22. SAVE ALL FIGURES TO DESKTOP
# ============================================================

# macOS Desktop
output_dir <- "~/Desktop/GSK_PK_State_Model"

# Create folder if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. PK concentration-time profile
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "01_GSK_PK_Profile.png"
  ),
  plot = p_pk,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "01_GSK_PK_Profile.pdf"
  ),
  plot = p_pk,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 2. Normalized exposure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "02_GSK_Normalized_Exposure.png"
  ),
  plot = p_exposure,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "02_GSK_Normalized_Exposure.pdf"
  ),
  plot = p_exposure,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 3. State-dependent response over time
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "03_State_Dependent_Response_Time.png"
  ),
  plot = p_state_response,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "03_State_Dependent_Response_Time.pdf"
  ),
  plot = p_state_response,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 4. Exposure × molecular-state response landscape
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "04_Exposure_State_Response_Landscape.png"
  ),
  plot = p_response_map,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "04_Exposure_State_Response_Landscape.pdf"
  ),
  plot = p_response_map,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# 5. State effect at Tmax
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "05_State_Effect_at_Tmax.png"
  ),
  plot = p_tmax,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "05_State_Effect_at_Tmax.pdf"
  ),
  plot = p_tmax,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 6. Exposure-response curves
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "06_Exposure_Response_by_State.png"
  ),
  plot = p_exposure_response,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "06_Exposure_Response_by_State.pdf"
  ),
  plot = p_exposure_response,
  width = 7,
  height = 5
)


# ============================================================
# SAVE SUMMARY TABLE AS CSV
# ============================================================

write.csv(
  simulation_summary,
  file = file.path(
    output_dir,
    "GSK_PK_State_Model_Summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# CONFIRM OUTPUT LOCATION
# ============================================================

cat(
  "\nAll figures saved to:\n",
  normalizePath(output_dir),
  "\n"
)

list.files(output_dir)
# ============================================================
# CHECK GDSC GSK-1904529A IC50 SCALE
# ============================================================

cat("\nGDSC GSK data:\n")
print(gsk_gdsc_unique)

cat("\nIC50 summary:\n")
print(summary(gsk_gdsc_unique$IC50))

cat("\nSorted IC50 values:\n")
print(
  gsk_gdsc_unique[
    order(gsk_gdsc_unique$IC50),
    c("ModelID", "State2_Likeness", "IC50")
  ]
)

cat("\nlog10(IC50) values:\n")

gsk_gdsc_unique$log10_IC50 <-
  log10(gsk_gdsc_unique$IC50)

print(
  gsk_gdsc_unique[
    order(gsk_gdsc_unique$log10_IC50),
    c(
      "ModelID",
      "State2_Likeness",
      "IC50",
      "log10_IC50"
    )
  ]
)

cat("\nlog10(IC50) summary:\n")
print(summary(gsk_gdsc_unique$log10_IC50))


# ============================================================
# CHECK ORIGINAL GDSC OBJECT / COLUMN NAMES
# ============================================================

cat("\nObjects containing 'gdsc':\n")
print(
  ls(
    pattern = "gdsc|GDSC"
  )
)

cat("\nStructure of gsk_gdsc_unique:\n")
str(gsk_gdsc_unique)


# ============================================================
# SPEARMAN CHECK
# ============================================================

test_raw <- cor.test(
  gsk_gdsc_unique$State2_Likeness,
  gsk_gdsc_unique$IC50,
  method = "spearman",
  exact = FALSE
)

test_log <- cor.test(
  gsk_gdsc_unique$State2_Likeness,
  gsk_gdsc_unique$log10_IC50,
  method = "spearman",
  exact = FALSE
)

cat(
  "\nRaw IC50 rho =",
  unname(test_raw$estimate),
  "P =",
  test_raw$p.value,
  "\n"
)

cat(
  "log10 IC50 rho =",
  unname(test_log$estimate),
  "P =",
  test_log$p.value,
  "\n"
)
ls(pattern = "gdsc|GDSC")
print(gsk_gdsc_unique)
# ============================================================
# TRACE ORIGINAL GDSC IC50 VALUES
# ============================================================

cat("\n--- gsk_gdsc2 ---\n")
print(gsk_gdsc2)

cat("\n--- Structure ---\n")
str(gsk_gdsc2)

cat("\n--- Column names ---\n")
print(colnames(gsk_gdsc2))

cat("\n--- gsk_gdsc ---\n")
print(gsk_gdsc)

cat("\n--- gdsc_score structure ---\n")
str(gdsc_score)

cat("\n--- gdsc_score columns ---\n")
print(colnames(gdsc_score))

cat("\n--- GDSC summary object ---\n")
print(gdsc2_summary)
# ============================================================
# FIND AVAILABLE PharmacoGx ACCESSORS
# ============================================================

library(PharmacoGx)

cat("\n--- Class ---\n")
print(class(gdsc2))

cat("\n--- Available sensitivity measures ---\n")
print(sensitivityMeasures(gdsc2))

cat("\n--- Methods containing sensitivity ---\n")
print(
  grep(
    "sensitiv|Sensitivity",
    getNamespaceExports("PharmacoGx"),
    value = TRUE
  )
)

cat("\n--- Slot names in gdsc2 ---\n")
print(slotNames(gdsc2))
# Original GDSC sensitivity profiles
gdsc_sens <- sensitivityProfiles(gdsc2)

class(gdsc_sens)
dim(gdsc_sens)
colnames(gdsc_sens)

# İlk satırlar
head(gdsc_sens)
# ============================================================
# EXTRACT ORIGINAL GSK-1904529A GDSC2 SENSITIVITY DATA
# ============================================================

# Find all rows containing GSK1904529A
gsk_rows <- grep(
  "GSK[-_ ]?1904529A",
  rownames(gdsc_sens),
  ignore.case = TRUE,
  value = TRUE
)

cat("\nNumber of GSK-1904529A experiments:", length(gsk_rows), "\n")

# Show first rows
head(gsk_rows, 20)


# ============================================================
# Extract sensitivity measurements
# ============================================================

gsk_original <- gdsc_sens[gsk_rows, , drop = FALSE]

cat("\nGSK sensitivity dimensions:\n")
print(dim(gsk_original))

cat("\nGSK sensitivity measurements:\n")
print(
  head(
    gsk_original[
      ,
      c(
        "aac_recomputed",
        "ic50_recomputed",
        "EC50",
        "E_inf",
        "HS"
      )
    ],
    20
  )
)


# ============================================================
# Extract cell-line names from row names
# ============================================================

gsk_original_df <- data.frame(
  Experiment = rownames(gsk_original),
  gsk_original,
  row.names = NULL
)

# Cell-line name is before "_GSK..."
gsk_original_df$CellLine <- sub(
  "_GSK.*$",
  "",
  gsk_original_df$Experiment,
  ignore.case = TRUE
)


# ============================================================
# Inspect clean table
# ============================================================

gsk_original_df <- gsk_original_df[
  ,
  c(
    "CellLine",
    "aac_recomputed",
    "ic50_recomputed",
    "EC50",
    "E_inf",
    "HS",
    "Experiment"
  )
]

cat("\n====================================\n")
cat("ORIGINAL GDSC GSK-1904529A DATA\n")
cat("====================================\n")

print(gsk_original_df)

cat("\nIC50 summary:\n")
print(summary(gsk_original_df$ic50_recomputed))

cat("\nEC50 summary:\n")
print(summary(gsk_original_df$EC50))

cat("\nAAC summary:\n")
print(summary(gsk_original_df$aac_recomputed))
# ============================================================
# MATCH ORIGINAL GDSC GSK DATA TO HGSOC CELL LINES
# ============================================================

# First inspect available HGSOC names
print(gdsc_hgsoc_names)

# Standardize names
clean_name <- function(x) {
  toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      x
    )
  )
}

gsk_original_df$CellLine_clean <-
  clean_name(gsk_original_df$CellLine)

# Inspect current HGSOC matching table
print(
  gdsc_hgsoc_match
)


# See column names
colnames(gdsc_hgsoc_match)
print(gdsc_hgsoc_match)
colnames(gdsc_hgsoc_match)

# ============================================================
# CLEAN GSK-1904529A GDSC ANALYSIS USING ORIGINAL DATA
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Standardize cell-line names
# ------------------------------------------------------------

clean_name <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", x))
}

gsk_original_df$CellLine_clean <-
  clean_name(gsk_original_df$CellLine)

gdsc_hgsoc_match$CellLine_clean <-
  clean_name(gdsc_hgsoc_match$Sample.Name)


# ------------------------------------------------------------
# 2. Match HGSOC lines to original GSK data
# ------------------------------------------------------------

gsk_hgsoc_clean <- gdsc_hgsoc_match %>%
  dplyr::left_join(
    gsk_original_df %>%
      dplyr::select(
        CellLine_clean,
        aac_recomputed,
        ic50_recomputed,
        EC50,
        E_inf,
        HS
      ),
    by = "CellLine_clean"
  )


# ------------------------------------------------------------
# 3. Inspect matched data
# ------------------------------------------------------------

cat("\n============================\n")
cat("GSK-1904529A HGSOC DATA\n")
cat("============================\n")

print(
  gsk_hgsoc_clean %>%
    dplyr::select(
      ModelID,
      Sample.Name,
      State2_Likeness,
      State_Like,
      aac_recomputed,
      EC50,
      E_inf,
      HS
    )
)

cat(
  "\nAAC available:",
  sum(!is.na(gsk_hgsoc_clean$aac_recomputed)),
  "/",
  nrow(gsk_hgsoc_clean),
  "\n"
)

cat(
  "EC50 available:",
  sum(!is.na(gsk_hgsoc_clean$EC50)),
  "/",
  nrow(gsk_hgsoc_clean),
  "\n"
)


# ------------------------------------------------------------
# 4. State2-likeness vs AAC
#
# IMPORTANT:
# Higher AAC = greater drug sensitivity
# ------------------------------------------------------------

aac_df <- gsk_hgsoc_clean %>%
  dplyr::filter(
    is.finite(State2_Likeness),
    is.finite(aac_recomputed)
  )

aac_test <- cor.test(
  aac_df$State2_Likeness,
  aac_df$aac_recomputed,
  method = "spearman",
  exact = FALSE
)

cat("\n============================\n")
cat("STATE2-LIKENESS vs AAC\n")
cat("============================\n")

cat(
  "N =",
  nrow(aac_df),
  "\nSpearman rho =",
  unname(aac_test$estimate),
  "\nP =",
  aac_test$p.value,
  "\n"
)


# ------------------------------------------------------------
# 5. State2-likeness vs EC50
#
# Lower EC50 = greater sensitivity
# ------------------------------------------------------------

ec50_df <- gsk_hgsoc_clean %>%
  dplyr::filter(
    is.finite(State2_Likeness),
    is.finite(EC50),
    EC50 > 0
  )

ec50_test <- cor.test(
  ec50_df$State2_Likeness,
  ec50_df$EC50,
  method = "spearman",
  exact = FALSE
)

cat("\n============================\n")
cat("STATE2-LIKENESS vs EC50\n")
cat("============================\n")

cat(
  "N =",
  nrow(ec50_df),
  "\nSpearman rho =",
  unname(ec50_test$estimate),
  "\nP =",
  ec50_test$p.value,
  "\n"
)


# ------------------------------------------------------------
# 6. Plot AAC
# ------------------------------------------------------------

p_gsk_aac <- ggplot(
  aac_df,
  aes(
    x = State2_Likeness,
    y = aac_recomputed
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "GSK-1904529A sensitivity in HGSOC models",
    subtitle = paste0(
      "AAC: Spearman rho = ",
      round(unname(aac_test$estimate), 2),
      ", P = ",
      signif(aac_test$p.value, 3)
    ),
    x = "State2-likeness",
    y = "AAC (higher = greater sensitivity)"
  )

print(p_gsk_aac)


# ------------------------------------------------------------
# 7. Plot EC50
# ------------------------------------------------------------

p_gsk_ec50 <- ggplot(
  ec50_df,
  aes(
    x = State2_Likeness,
    y = EC50
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  scale_y_log10() +
  theme_classic(base_size = 13) +
  labs(
    title = "GSK-1904529A EC50 in HGSOC models",
    subtitle = paste0(
      "Spearman rho = ",
      round(unname(ec50_test$estimate), 2),
      ", P = ",
      signif(ec50_test$p.value, 3)
    ),
    x = "State2-likeness",
    y = "EC50 (log10 scale)"
  )

print(p_gsk_ec50)


# ------------------------------------------------------------
# 8. Save corrected results
# ------------------------------------------------------------

output_dir <- "~/Desktop/GSK_PK_State_Model"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(
  gsk_hgsoc_clean,
  file.path(
    output_dir,
    "GSK_HGSOC_original_GDSC_sensitivity.csv"
  ),
  row.names = FALSE
)

ggsave(
  file.path(
    output_dir,
    "07_GSK_State2_Likeness_AAC.png"
  ),
  p_gsk_aac,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "07_GSK_State2_Likeness_AAC.pdf"
  ),
  p_gsk_aac,
  width = 7,
  height = 5
)

ggsave(
  file.path(
    output_dir,
    "08_GSK_State2_Likeness_EC50.png"
  ),
  p_gsk_ec50,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "08_GSK_State2_Likeness_EC50.pdf"
  ),
  p_gsk_ec50,
  width = 7,
  height = 5
)

#######ilaclari tekrar kontrol ediyoruz 
# ============================================================
# PAN-DRUG SCREEN:
# State2-likeness vs GDSC AAC + EC50
# ============================================================

library(PharmacoGx)
library(dplyr)
library(stringr)
library(ggplot2)

# ------------------------------------------------------------
# 1. GDSC sensitivity profiles
# ------------------------------------------------------------

gdsc_sens <- sensitivityProfiles(gdsc2)

cat("Sensitivity matrix:\n")
print(dim(gdsc_sens))

# ------------------------------------------------------------
# 2. Convert row names into experiment information
# ------------------------------------------------------------

gdsc_all <- data.frame(
  Experiment = rownames(gdsc_sens),
  gdsc_sens,
  row.names = NULL,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. Match experiments to our 13 HGSOC cell lines
#
# Use Sample.Name from gdsc_hgsoc_match.
# We match using the beginning of the experiment name.
# ------------------------------------------------------------

clean_name <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", x))
}

# Function to identify which HGSOC cell line each experiment belongs to
extract_hgsoc_line <- function(experiment, sample_names) {
  
  hits <- sample_names[
    sapply(
      sample_names,
      function(x) {
        startsWith(
          clean_name(experiment),
          clean_name(x)
        )
      }
    )
  ]
  
  if (length(hits) == 0) return(NA_character_)
  
  # In case of multiple matches use longest name
  hits[which.max(nchar(hits))]
}

gdsc_all$Sample.Name <- vapply(
  gdsc_all$Experiment,
  extract_hgsoc_line,
  FUN.VALUE = character(1),
  sample_names = gdsc_hgsoc_match$Sample.Name
)

# Keep HGSOC experiments
gdsc_hgsoc_all <- gdsc_all %>%
  filter(!is.na(Sample.Name))

cat("\nHGSOC experiments:", nrow(gdsc_hgsoc_all), "\n")


# ------------------------------------------------------------
# 4. Extract drug name from experiment name
#
# Structure:
# CellLine_Drug_DrugID_...
# ------------------------------------------------------------

extract_drug <- function(experiment, sample_name) {
  
  x <- sub(
    paste0("^", sample_name, "_"),
    "",
    experiment,
    ignore.case = TRUE
  )
  
  # Drug name = everything before first numeric drug ID
  sub("_[0-9]+_.*$", "", x)
}

gdsc_hgsoc_all$Drug <- mapply(
  extract_drug,
  gdsc_hgsoc_all$Experiment,
  gdsc_hgsoc_all$Sample.Name
)


# ------------------------------------------------------------
# 5. Add State2-likeness
# ------------------------------------------------------------

gdsc_hgsoc_all <- gdsc_hgsoc_all %>%
  left_join(
    gdsc_hgsoc_match %>%
      select(
        Sample.Name,
        ModelID,
        State2_Likeness,
        State_Like
      ),
    by = "Sample.Name"
  )


# ------------------------------------------------------------
# 6. Collapse technical/duplicate experiments
#
# One value per Drug × ModelID
# ------------------------------------------------------------

drug_cell <- gdsc_hgsoc_all %>%
  group_by(
    Drug,
    ModelID,
    Sample.Name,
    State2_Likeness,
    State_Like
  ) %>%
  summarise(
    AAC = median(aac_recomputed, na.rm = TRUE),
    EC50 = median(EC50, na.rm = TRUE),
    .groups = "drop"
  )

# Convert NaN generated by all-NA groups back to NA
drug_cell$AAC[!is.finite(drug_cell$AAC)] <- NA
drug_cell$EC50[!is.finite(drug_cell$EC50)] <- NA


# ------------------------------------------------------------
# 7. Function: Spearman test
# ------------------------------------------------------------

safe_spearman <- function(x, y) {
  
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if (
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


# ------------------------------------------------------------
# 8. AAC screen
#
# Higher AAC = greater sensitivity
#
# Therefore:
# rho > 0 = State2-like more sensitive
# rho < 0 = State1-like more sensitive
# ------------------------------------------------------------

aac_results <- drug_cell %>%
  group_by(Drug) %>%
  group_modify(
    ~ safe_spearman(
      .x$State2_Likeness,
      .x$AAC
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    AAC_direction = case_when(
      rho > 0 ~ "State2-sensitive",
      rho < 0 ~ "State1-sensitive",
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 9. EC50 screen
#
# Lower EC50 = greater sensitivity
#
# Therefore:
# rho < 0 = State2-like more sensitive
# rho > 0 = State1-like more sensitive
# ------------------------------------------------------------

ec50_results <- drug_cell %>%
  filter(is.na(EC50) | EC50 > 0) %>%
  group_by(Drug) %>%
  group_modify(
    ~ safe_spearman(
      .x$State2_Likeness,
      .x$EC50
    )
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    EC50_direction = case_when(
      rho < 0 ~ "State2-sensitive",
      rho > 0 ~ "State1-sensitive",
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 10. Combine AAC + EC50
# ------------------------------------------------------------

drug_screen <- aac_results %>%
  rename(
    N_AAC = N,
    rho_AAC = rho,
    P_AAC = P,
    FDR_AAC = FDR
  ) %>%
  full_join(
    ec50_results %>%
      rename(
        N_EC50 = N,
        rho_EC50 = rho,
        P_EC50 = P,
        FDR_EC50 = FDR
      ),
    by = "Drug"
  )


# ------------------------------------------------------------
# 11. Determine concordance
# ------------------------------------------------------------

drug_screen <- drug_screen %>%
  mutate(
    
    Concordance = case_when(
      
      AAC_direction == "State2-sensitive" &
        EC50_direction == "State2-sensitive" ~
        "Concordant State2-sensitive",
      
      AAC_direction == "State1-sensitive" &
        EC50_direction == "State1-sensitive" ~
        "Concordant State1-sensitive",
      
      !is.na(AAC_direction) &
        !is.na(EC50_direction) ~
        "Discordant",
      
      TRUE ~ "Insufficient"
    ),
    
    Mean_abs_rho =
      rowMeans(
        cbind(
          abs(rho_AAC),
          abs(rho_EC50)
        ),
        na.rm = TRUE
      )
  )


# ------------------------------------------------------------
# 12. Rank candidates
#
# Require at least 6 models for BOTH endpoints
# ------------------------------------------------------------

drug_ranked <- drug_screen %>%
  filter(
    N_AAC >= 6,
    N_EC50 >= 6
  ) %>%
  arrange(
    FDR_AAC,
    P_AAC,
    desc(Mean_abs_rho)
  )


# ------------------------------------------------------------
# 13. Show top concordant State2 candidates
# ------------------------------------------------------------

state2_candidates <- drug_ranked %>%
  filter(
    Concordance == "Concordant State2-sensitive"
  ) %>%
  arrange(
    FDR_AAC,
    P_AAC,
    desc(Mean_abs_rho)
  )

cat("\n========================================\n")
cat("TOP CONCORDANT STATE2-SENSITIVE DRUGS\n")
cat("========================================\n")

print(
  state2_candidates %>%
    select(
      Drug,
      N_AAC,
      rho_AAC,
      P_AAC,
      FDR_AAC,
      N_EC50,
      rho_EC50,
      P_EC50,
      FDR_EC50,
      Mean_abs_rho
    ) %>%
    head(30)
)


# ------------------------------------------------------------
# 14. FDR-significant results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("AAC FDR < 0.05\n")
cat("========================================\n")

print(
  drug_ranked %>%
    filter(FDR_AAC < 0.05) %>%
    select(
      Drug,
      N_AAC,
      rho_AAC,
      P_AAC,
      FDR_AAC,
      Concordance
    )
)


cat("\n========================================\n")
cat("EC50 FDR < 0.05\n")
cat("========================================\n")

print(
  drug_ranked %>%
    filter(FDR_EC50 < 0.05) %>%
    select(
      Drug,
      N_EC50,
      rho_EC50,
      P_EC50,
      FDR_EC50,
      Concordance
    )
)


# ------------------------------------------------------------
# 15. Strong concordant exploratory candidates
#
# This is NOT statistical significance.
# Just useful for candidate prioritization.
# ------------------------------------------------------------

strong_concordant <- drug_ranked %>%
  filter(
    Concordance != "Discordant",
    abs(rho_AAC) >= 0.5,
    abs(rho_EC50) >= 0.5
  ) %>%
  arrange(desc(Mean_abs_rho))

cat("\n========================================\n")
cat("STRONG CONCORDANT CANDIDATES\n")
cat("|rho AAC| >= 0.5 AND |rho EC50| >= 0.5\n")
cat("========================================\n")

print(
  strong_concordant %>%
    select(
      Drug,
      Concordance,
      N_AAC,
      rho_AAC,
      P_AAC,
      FDR_AAC,
      N_EC50,
      rho_EC50,
      P_EC50,
      FDR_EC50,
      Mean_abs_rho
    ) %>%
    head(30)
)


# ------------------------------------------------------------
# 16. Save everything
# ------------------------------------------------------------

output_dir <- "~/Desktop/GSK_PK_State_Model"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(
  drug_screen,
  file.path(
    output_dir,
    "GDSC_all_drugs_State2_AAC_EC50.csv"
  ),
  row.names = FALSE
)

write.csv(
  state2_candidates,
  file.path(
    output_dir,
    "GDSC_concordant_State2_candidates.csv"
  ),
  row.names = FALSE
)

write.csv(
  strong_concordant,
  file.path(
    output_dir,
    "GDSC_strong_concordant_candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")

cat(
  "Drugs tested:",
  nrow(drug_ranked),
  "\n"
)

cat(
  "Concordant State2-sensitive:",
  sum(
    drug_ranked$Concordance ==
      "Concordant State2-sensitive"
  ),
  "\n"
)

cat(
  "Concordant State1-sensitive:",
  sum(
    drug_ranked$Concordance ==
      "Concordant State1-sensitive"
  ),
  "\n"
)

cat(
  "Discordant:",
  sum(
    drug_ranked$Concordance ==
      "Discordant"
  ),
  "\n"
)

cat(
  "AAC FDR < 0.05:",
  sum(
    drug_ranked$FDR_AAC < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "EC50 FDR < 0.05:",
  sum(
    drug_ranked$FDR_EC50 < 0.05,
    na.rm = TRUE
  ),
  "\n"
)
# ============================================================
# PRISM VALIDATION OF TOP GDSC CANDIDATES
# Sepantronium / AZD4547 / Sorafenib
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Candidate drugs
# ------------------------------------------------------------

candidate_drugs <- c(
  "sepantronium",
  "azd4547",
  "sorafenib"
)

# ------------------------------------------------------------
# 2. Inspect PRISM drug names first
# Assumes prism data object used previously = prism
# ------------------------------------------------------------

cat("\nMatching PRISM drug names:\n")

for (drug in candidate_drugs) {
  
  cat("\n---", drug, "---\n")
  
  print(
    unique(
      prism$name[
        grepl(
          drug,
          prism$name,
          ignore.case = TRUE
        )
      ]
    )
  )
}


# ------------------------------------------------------------
# 3. Extract candidate drug measurements
# ------------------------------------------------------------

prism_candidates <- prism %>%
  dplyr::filter(
    grepl(
      paste(candidate_drugs, collapse = "|"),
      name,
      ignore.case = TRUE
    )
  )

cat("\nCandidate PRISM rows:\n")
print(dim(prism_candidates))

cat("\nMatched drug names:\n")
print(unique(prism_candidates$name))


# ------------------------------------------------------------
# 4. Add State2-likeness
#
# PRISM uses depmap_id
# gdsc_hgsoc_match uses ModelID
# ------------------------------------------------------------

state_lookup <- gdsc_hgsoc_match %>%
  dplyr::select(
    ModelID,
    CellLineName,
    State2_Likeness,
    State_Like
  )

prism_candidates <- prism_candidates %>%
  dplyr::inner_join(
    state_lookup,
    by = c(
      "depmap_id" = "ModelID"
    )
  )

cat(
  "\nRows after matching HGSOC models:",
  nrow(prism_candidates),
  "\n"
)


# ------------------------------------------------------------
# 5. Collapse replicate measurements
#
# One AUC per:
# Drug × cell line
# ------------------------------------------------------------

prism_candidate_cell <- prism_candidates %>%
  dplyr::group_by(
    name,
    depmap_id,
    CellLineName,
    State2_Likeness,
    State_Like
  ) %>%
  dplyr::summarise(
    
    AUC = median(
      auc,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(AUC),
    is.finite(State2_Likeness)
  )


cat("\nPRISM candidate-cell table:\n")

print(
  prism_candidate_cell %>%
    dplyr::arrange(
      name,
      State2_Likeness
    )
)


# ------------------------------------------------------------
# 6. Safe Spearman function
#
# Lower PRISM AUC = greater sensitivity
#
# rho < 0:
# State2-likeness ↑ -> AUC ↓
# = State2-sensitive
# ------------------------------------------------------------

safe_prism_cor <- function(df) {
  
  if (
    nrow(df) < 6 ||
    length(unique(df$State2_Likeness)) < 3 ||
    length(unique(df$AUC)) < 3
  ) {
    
    return(
      data.frame(
        N = nrow(df),
        rho = NA_real_,
        P = NA_real_
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      df$State2_Likeness,
      df$AUC,
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.frame(
    N = nrow(df),
    rho = unname(test$estimate),
    P = test$p.value
  )
}


# ------------------------------------------------------------
# 7. Test each drug
# ------------------------------------------------------------

prism_results <- prism_candidate_cell %>%
  
  dplyr::group_by(name) %>%
  
  dplyr::group_modify(
    ~ safe_prism_cor(.x)
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    
    FDR = p.adjust(
      P,
      method = "BH"
    ),
    
    Direction =
      dplyr::case_when(
        
        rho < 0 ~
          "State2-sensitive",
        
        rho > 0 ~
          "State1-sensitive",
        
        TRUE ~
          NA_character_
      )
  )


cat("\n")
cat("============================================\n")
cat("PRISM RESULTS\n")
cat("============================================\n")

print(prism_results)


# ------------------------------------------------------------
# 8. Plot each candidate
# ------------------------------------------------------------

plot_list <- list()

for (drug in unique(prism_candidate_cell$name)) {
  
  tmp <- prism_candidate_cell %>%
    dplyr::filter(
      name == drug
    )
  
  stat <- prism_results %>%
    dplyr::filter(
      name == drug
    )
  
  p <- ggplot(
    tmp,
    aes(
      x = State2_Likeness,
      y = AUC
    )
  ) +
    
    geom_point(
      size = 3
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    labs(
      
      title = paste0(
        drug,
        " sensitivity in HGSOC models"
      ),
      
      subtitle = paste0(
        "PRISM AUC: N = ",
        stat$N,
        ", rho = ",
        round(stat$rho, 2),
        ", P = ",
        signif(stat$P, 3)
      ),
      
      x = "State2-likeness",
      
      y = "PRISM AUC (lower = greater sensitivity)"
    )
  
  print(p)
  
  plot_list[[drug]] <- p
}


# ------------------------------------------------------------
# 9. Compare with GDSC results
# ------------------------------------------------------------

gdsc_top <- drug_ranked %>%
  
  dplyr::filter(
    grepl(
      paste(candidate_drugs, collapse = "|"),
      Drug,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::select(
    Drug,
    N_AAC,
    rho_AAC,
    P_AAC,
    FDR_AAC,
    N_EC50,
    rho_EC50,
    P_EC50,
    FDR_EC50,
    Concordance
  )


cat("\n")
cat("============================================\n")
cat("GDSC RESULTS\n")
cat("============================================\n")

print(gdsc_top)


cat("\n")
cat("============================================\n")
cat("PRISM RESULTS\n")
cat("============================================\n")

print(prism_results)


# ------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------

output_dir <- "~/Desktop/GSK_PK_State_Model"

if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


write.csv(
  prism_results,
  file.path(
    output_dir,
    "PRISM_top_GDSC_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  prism_candidate_cell,
  file.path(
    output_dir,
    "PRISM_top_GDSC_candidates_celllines.csv"
  ),
  row.names = FALSE
)


# Save plots
for (drug in names(plot_list)) {
  
  safe_name <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    drug
  )
  
  ggsave(
    file.path(
      output_dir,
      paste0(
        "PRISM_",
        safe_name,
        "_State2_AUC.png"
      )
    ),
    plot_list[[drug]],
    width = 7,
    height = 5,
    dpi = 600
  )
  
  ggsave(
    file.path(
      output_dir,
      paste0(
        "PRISM_",
        safe_name,
        "_State2_AUC.pdf"
      )
    ),
    plot_list[[drug]],
    width = 7,
    height = 5
  )
}


cat("\nAnalysis complete.\n")



