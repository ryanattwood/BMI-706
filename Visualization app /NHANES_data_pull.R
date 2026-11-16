# =========================================================================
# 01_fetch_and_clean_nhanes.R
# Standalone offline script to pull NHANES 2017-2020 (P_) & 2021-2023 (_L)
# directly from CDC/NCHS, compute clinical biomarkers (BMI, 2021 CKD-EPI eGFR,
# Non-HDL-C, Mean SBP, uACR, Glycemic/Smoking status), and save RDS files.
# =========================================================================

library(haven)
library(dplyr)
library(tidyr)
library(purrr)
library(survey)

# Prevent lonely PSU errors in stratified survey quantile calculations
options(survey.lonely.psu = "adjust")
options(timeout = 300) # Give CDC downloads plenty of time

# -------------------------------------------------------------------------
# 1. HELPER FUNCTION TO DOWNLOAD & MERGE ONE NHANES WAVE DIRECTLY FROM CDC
# -------------------------------------------------------------------------
fetch_nhanes_wave <- function(cycle_years, modules, is_prepandemic = TRUE) {
  base_url <- paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", cycle_years, "/DataFiles/")
  
  module_list <- imap(modules, function(file_stem, mod_name) {
    file_name <- if (is_prepandemic) {
      paste0("P_", file_stem, ".xpt")
    } else {
      paste0(file_stem, "_L.xpt")
    }
    full_url <- paste0(base_url, file_name)
    message(sprintf("[%s] Downloading %s from %s ...", mod_name, file_name, full_url))
    
    df <- read_xpt(full_url)
    return(df)
  })
  
  # Left-join all modules 1-to-1 onto DEMO via SEQN
  merged_df <- reduce(module_list, left_join, by = "SEQN")
  return(merged_df)
}

# Core modules needed across both waves
# Note: In 2021-2023, cholesterol/biochem/urine stems match standard NHANES names
modules_to_pull <- list(
  DEMO   = "DEMO",    # Demographics & Survey Weights
  BMX    = "BMX",     # Body Measures (BMI, Height, Weight)
  BPXO   = "BPXO",    # Oscillometric Blood Pressure (3 readings)
  TCHOL  = "TCHOL",   # Total Cholesterol
  HDL    = "HDL",     # HDL Cholesterol
  GHB    = "GHB",     # Glycohemoglobin (HbA1c)
  BIOPRO = "BIOPRO",  # Standard Biochemistry (Serum Creatinine LBXSCR)
  ALB_CR = "ALB_CR",  # Urine Albumin & Creatinine (URDACT)
  DIQ    = "DIQ",     # Diabetes Questionnaire
  BPQ    = "BPQ",     # Blood Pressure & Cholesterol Questionnaire
  SMQ    = "SMQ"      # Smoking - Cigarette Use Questionnaire
)

# -------------------------------------------------------------------------
# 2. DOWNLOAD BOTH WAVES INTO MEMORY
# -------------------------------------------------------------------------
message("\n--- Pulling 2017-March 2020 Pre-Pandemic Wave (P_) ---")
raw_2017_2020 <- fetch_nhanes_wave("2017", modules_to_pull, is_prepandemic = TRUE) %>%
  mutate(
    wave    = "2017-March 2020",
    # Rescale 3.2-year MEC exam weight to 5.2-year pooled denominator
    WT_POOL = (3.2 / 5.2) * WTMECPRP,
    # Offset stratum IDs by wave so strata never collide across cycles
    STRA_POOL = SDMVSTRA + 10000
  )

message("\n--- Pulling August 2021-August 2023 Wave (_L) ---")
raw_2021_2023 <- fetch_nhanes_wave("2021", modules_to_pull, is_prepandemic = FALSE) %>%
  mutate(
    wave    = "2021-2023",
    # Rescale 2.0-year MEC exam weight to 5.2-year pooled denominator
    WT_POOL = (2.0 / 5.2) * WTMEC2YR,
    STRA_POOL = SDMVSTRA + 20000
  )

# -------------------------------------------------------------------------
# 3. POOL WAVES & DERIVE ALL STATIC CLINICAL VARIABLES
# -------------------------------------------------------------------------
message("\n--- Cleaning & Calculating Derived Biomarkers (BMI, 2021 CKD-EPI eGFR, Non-HDL, SBP) ---")

nhanes_analytic_cohort <- bind_rows(raw_2017_2020, raw_2021_2023) %>%
  # Retain MEC-examined adults aged 30-79 (matching AHA PREVENT & KFRE domain)
  filter(!is.na(WT_POOL), WT_POOL > 0, RIDAGEYR >= 30, RIDAGEYR <= 79) %>%
  mutate(
    # --- A. Demographics & Stratification Bins ---
    age = as.numeric(RIDAGEYR),
    age_grp = cut(
      age,
      breaks = c(29, 39, 49, 59, 69, 79),
      labels = c("30-39", "40-49", "50-59", "60-69", "70-79")
    ),
    sex = ifelse(RIAGENDR == 1, "Male", "Female"),
    
    # --- B. Anthropometrics: BMI (Direct + Height/Weight Fallback) ---
    bmi_calc = ifelse(!is.na(BMXWT) & !is.na(BMXHT) & BMXHT > 0,
                      BMXWT / ((BMXHT / 100)^2), NA_real_),
    bmi = coalesce(BMXBMI, round(bmi_calc, 1)),
    bmi_tier = case_when(
      is.na(bmi) ~ NA_character_,
      bmi < 25   ~ "<25 (Normal/Under)",
      bmi < 30   ~ "25-29.9 (Overweight)",
      TRUE       ~ ">=30 (Obese)"
    ),
    
    # --- C. Hemodynamics: Mean Oscillometric Systolic & Diastolic BP ---
    sbp = rowMeans(pick(BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
    sbp = ifelse(is.nan(sbp), NA_real_, sbp),
    dbp = rowMeans(pick(BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE),
    dbp = ifelse(is.nan(dbp), NA_real_, dbp),
    
    # --- D. Lipids & Glycemic Labs ---
    tchol   = as.numeric(LBXTC),
    hdl     = as.numeric(LBDHDD),
    non_hdl = ifelse(!is.na(tchol) & !is.na(hdl), tchol - hdl, NA_real_),
    hba1c   = as.numeric(LBXGH),
    
    # --- E. Renal Function: 2021 Race-Free CKD-EPI eGFR from Serum Creatinine ---
    scr   = ifelse(!is.na(LBXSCR) & LBXSCR > 0, as.numeric(LBXSCR), NA_real_),
    kappa = ifelse(sex == "Female", 0.7, 0.9),
    alpha = ifelse(sex == "Female", -0.241, -0.302),
    scr_k = scr / kappa,
    egfr  = 142 * (pmin(scr_k, 1)^alpha) * (pmax(scr_k, 1)^(-1.200)) *
      (0.9938^age) * ifelse(sex == "Female", 1.012, 1.0),
    
    # --- F. Albuminuria: Urine Albumin-to-Creatinine Ratio (mg/g) ---
    uacr_calc = ifelse(!is.na(URXUMA) & !is.na(URXUCR) & URXUCR > 0,
                       (URXUMA / URXUCR) * 100, NA_real_),
    uacr     = coalesce(as.numeric(URDACT), uacr_calc),
    log_uacr = ifelse(!is.na(uacr) & uacr > 0, log(uacr), NA_real_),
    
    # --- G. Clean Questionnaire Codes (Recode 7/9 Refused/Don't Know -> NA) ---
    diq010_clean = ifelse(DIQ010 %in% c(1, 2, 3), DIQ010, NA_real_),
    dm_status = case_when(
      diq010_clean == 1 | (!is.na(hba1c) & hba1c >= 6.5) ~ "Diabetes",
      diq010_clean == 3 | (!is.na(hba1c) & hba1c >= 5.7 & hba1c < 6.5) ~ "Prediabetes",
      !is.na(diq010_clean) | !is.na(hba1c) ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    dm_status = factor(dm_status, levels = c("No Diabetes", "Prediabetes", "Diabetes")),
    
    # Hypertension history & Antihypertensive medication use
    htn_hx = case_when(
      BPQ020 == 1 ~ 1L,
      BPQ020 == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    tx_bp = case_when(
      BPQ050A == 1 ~ 1L,
      BPQ020 == 2 | BPQ050A == 2 ~ 0L,
      TRUE ~ 0L
    ),
    
    # Cholesterol / Statin medication use (BPQ100D in 2021-2023 / BPQ090D in 2017-2020)
    tx_statin = case_when(
      ("BPQ100D" %in% names(.) & BPQ100D == 1) ~ 1L,
      ("BPQ090D" %in% names(.) & BPQ090D == 1) ~ 1L,
      TRUE ~ 0L
    ),
    
    # Smoking status (SMQ020 = 100+ cigarettes in life; SMQ040 = current status)
    smoking_status = case_when(
      SMQ020 == 2 ~ "Never",
      SMQ020 == 1 & SMQ040 == 3 ~ "Former",
      SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
      TRUE ~ "Never"
    ),
    current_smoker = ifelse(smoking_status == "Current", 1L, 0L)
  ) %>%
  # Select only the clean variables needed for analysis & Shiny
  select(
    SEQN, wave, WT_POOL, STRA_POOL, SDMVPSU,
    age, age_grp, sex, bmi, bmi_tier,
    sbp, dbp, tchol, hdl, non_hdl, hba1c, scr, egfr, uacr, log_uacr,
    dm_status, htn_hx, tx_bp, tx_statin, smoking_status, current_smoker
  )

# Save the full individual-level cleaned analytic cohort (~11,200 rows)
saveRDS(nhanes_analytic_cohort, "nhanes_analytic_cohort.rds")
message(sprintf("Saved 'nhanes_analytic_cohort.rds' (%d adult participants, %d variables).",
                nrow(nhanes_analytic_cohort), ncol(nhanes_analytic_cohort)))

# -------------------------------------------------------------------------
# 4. PRE-COMPUTE SURVEY-WEIGHTED QUANTILE LOOKUP TABLE (DIRECT WEIGHTED ECDF)
# -------------------------------------------------------------------------
message("\n--- Computing Survey-Weighted 10th, 25th, 50th, 75th, 90th Percentiles ---")

# Pure base-R + dplyr survey-weighted quantile helper (zero svyby class issues)
calc_weighted_quantiles <- function(x, w, probs = c(0.10, 0.25, 0.50, 0.75, 0.90)) {
  valid <- !is.na(x) & !is.na(w) & w > 0
  x <- x[valid]
  w <- w[valid]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  # Cumulative survey weight proportion (0 to 1)
  cum_w <- cumsum(w) / sum(w)
  approx(cum_w, x, xout = probs, method = "linear", rule = 2, ties = "mean")$y
}

biomarkers_to_summarize <- c("bmi", "sbp", "non_hdl", "hba1c", "scr", "egfr", "uacr")

nhanes_quantiles <- nhanes_analytic_cohort %>%
  filter(!is.na(age_grp), !is.na(sex), !is.na(dm_status)) %>%
  pivot_longer(
    cols      = all_of(biomarkers_to_summarize),
    names_to  = "biomarker",
    values_to = "val"
  ) %>%
  group_by(biomarker, age_grp, sex, dm_status) %>%
  summarise(
    n_unweighted = sum(!is.na(val)),
    p10 = calc_weighted_quantiles(val, WT_POOL, 0.10),
    p25 = calc_weighted_quantiles(val, WT_POOL, 0.25),
    p50 = calc_weighted_quantiles(val, WT_POOL, 0.50),
    p75 = calc_weighted_quantiles(val, WT_POOL, 0.75),
    p90 = calc_weighted_quantiles(val, WT_POOL, 0.90),
    .groups = "drop"
  ) %>%
  mutate(across(c(p10, p25, p50, p75, p90), ~ round(.x, 2)))

saveRDS(nhanes_quantiles, "nhanes_quantiles.rds")
message(sprintf("Saved 'nhanes_quantiles.rds' (%d stratified quantile rows).", nrow(nhanes_quantiles)))

# Quick preview in console
print(head(nhanes_quantiles, 10))