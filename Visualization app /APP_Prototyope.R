# =========================================================================
# app.R - Cardiorenal Risk & Trajectory Explorer (Full Integrated Prototype)
# Requires: nhanes_quantiles.rds & nhanes_analytic_cohort.rds in working dir
# =========================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(ggiraph)
library(dplyr)
library(tidyr)
library(purrr)

# -------------------------------------------------------------------------
# 1. LOAD PRE-COMPUTED NHANES RDS FILES
# -------------------------------------------------------------------------
if (!file.exists("nhanes_quantiles.rds") || !file.exists("nhanes_analytic_cohort.rds")) {
  stop("Make sure both 'nhanes_quantiles.rds' and 'nhanes_analytic_cohort.rds' are in your working directory!")
}

nhanes_ref <- readRDS("nhanes_quantiles.rds") %>%
  mutate(across(c(age_grp, sex, dm_status), as.character))

nhanes_cohort <- readRDS("nhanes_analytic_cohort.rds") %>%
  mutate(across(c(age_grp, sex, dm_status), as.character))

has_preventr <- requireNamespace("preventr", quietly = TRUE)

# -------------------------------------------------------------------------
# 2. CLINICAL THRESHOLDS & NHANES DENSITY AUC HELPER
# -------------------------------------------------------------------------
biomarker_meta <- tibble(
  biomarker = c("bmi", "sbp", "non_hdl", "hba1c", "egfr", "uacr"),
  label     = factor(
    c("BMI (kg/m²)", "Systolic BP (mmHg)", "Non-HDL-C (mg/dL)",
      "HbA1c (%)", "eGFR (mL/min/1.73m²)", "uACR (mg/g)"),
    levels = c("BMI (kg/m²)", "Systolic BP (mmHg)", "Non-HDL-C (mg/dL)",
               "HbA1c (%)", "eGFR (mL/min/1.73m²)", "uACR (mg/g)")
  ),
  xmin      = c(15,  85,  50,  4.2, 15,  1),
  xmax      = c(50,  200, 260, 12.0, 125, 450),
  cut1      = c(25,  120, 130, 5.7, 60,  30),
  cut2      = c(30,  140, 160, 6.5, 90,  300),
  reverse   = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) # eGFR is healthier when higher
)

build_auc_density <- function(df_sub, bio_row) {
  var_col <- bio_row$biomarker
  df_valid <- df_sub %>%
    filter(!is.na(.data[[var_col]]), !is.na(WT_POOL), WT_POOL > 0)
  
  x_vals <- pmin(bio_row$xmax, pmax(bio_row$xmin, df_valid[[var_col]]))
  w_vals <- df_valid$WT_POOL / sum(df_valid$WT_POOL)
  if (length(x_vals) < 10) return(tibble())
  
  dens <- density(
    x_vals,
    weights = w_vals,
    from    = bio_row$xmin,
    to      = bio_row$xmax,
    n       = 256,
    adjust  = 1.15
  )
  
  tibble(
    biomarker = var_col,
    label     = bio_row$label,
    x         = dens$x,
    y         = dens$y / max(dens$y, na.rm = TRUE)
  ) %>%
    mutate(
      zone = case_when(
        !bio_row$reverse & x <  bio_row$cut1 ~ "Optimal (Green)",
        !bio_row$reverse & x <  bio_row$cut2 ~ "Borderline / Elevated (Yellow)",
        !bio_row$reverse                     ~ "High Clinical Risk (Red)",
        bio_row$reverse  & x <  bio_row$cut1 ~ "High Clinical Risk (Red)",
        bio_row$reverse  & x <  bio_row$cut2 ~ "Borderline / Elevated (Yellow)",
        TRUE                                 ~ "Optimal (Green)"
      ),
      zone = factor(
        zone,
        levels = c("Optimal (Green)", "Borderline / Elevated (Yellow)", "High Clinical Risk (Red)")
      )
    )
}

# -------------------------------------------------------------------------
# 3. CLINICAL PREDICTION MODELS (2021 CKD-EPI, AHA PREVENT, CKD-PC)
# -------------------------------------------------------------------------

# 2021 CKD-EPI Race-Free Creatinine eGFR Equation
calc_ckdepi_2021 <- function(scr, age, sex) {
  kappa <- ifelse(sex == "Female", 0.7, 0.9)
  alpha <- ifelse(sex == "Female", -0.241, -0.302)
  scr_k <- scr / kappa
  egfr  <- 142 * (pmin(scr_k, 1)^alpha) * (pmax(scr_k, 1)^(-1.200)) *
    (0.9938^age) * ifelse(sex == "Female", 1.012, 1.0)
  round(pmin(120, pmax(15, egfr)), 1)
}

# 10-Year AHA PREVENT ASCVD Risk Calculator (CRAN preventr + calibrated fallback)
calc_ascvd_10yr <- function(age, sex, bmi, sbp, non_hdl, hba1c, egfr, uacr,
                            dm, smoker, statin, bp_tx, sglt2, glp1) {
  prob_base <- NA_real_
  if (has_preventr) {
    try({
      res <- preventr::estimate_risk(
        age     = age,
        sex     = tolower(sex),
        sbp     = pmin(200, pmax(90, sbp)),
        bp_tx   = bp_tx,
        total_c = pmin(320, pmax(130, round(non_hdl + 48))),
        hdl_c   = 48,
        statin  = statin,
        dm      = (dm == "Diabetes"),
        smoking = smoker,
        egfr    = pmin(140, pmax(15, egfr)),
        bmi     = pmin(39.9, pmax(18.5, bmi)),
        hba1c   = pmin(15, pmax(4.5, hba1c)),
        uacr    = pmin(25000, pmax(0.1, uacr)),
        model   = "prevent",
        time    = "10yr",
        quiet   = TRUE
      )
      if (is.data.frame(res) && "ascvd" %in% names(res)) {
        prob_base <- as.numeric(res$ascvd[1])
      }
    }, silent = TRUE)
  }
  
  if (is.na(prob_base)) {
    age_c   <- (age - 55) / 10
    sbp_hi  <- pmax(0, sbp - 110) / 20
    egfr_lo <- pmax(0, 60 - egfr) / 15
    nhdl_c  <- (non_hdl - 130) / 40
    a1c_c   <- pmax(0, hba1c - 5.3)
    uacr_ln <- log(pmax(1, uacr) / 10)
    bmi_hi  <- pmax(0, bmi - 25) / 5
    b0      <- ifelse(sex == "Male", -2.95, -3.30)
    
    lp <- b0 +
      0.78 * age_c +
      0.34 * sbp_hi +
      0.22 * nhdl_c +
      0.19 * a1c_c +
      0.26 * egfr_lo +
      0.16 * uacr_ln +
      0.08 * bmi_hi +
      0.42 * (dm == "Diabetes") +
      0.48 * smoker -
      0.18 * statin -
      0.12 * bp_tx
    
    prob_base <- 1 / (1 + exp(-lp))
  }
  
  hr_mult <- ifelse(sglt2, 0.86, 1.0) * ifelse(glp1, 0.86, 1.0)
  pmin(0.95, pmax(0.01, prob_base * hr_mult))
}

# 20-Year CKD-PC Chronic eGFR Slope Calculator
calc_egfr_slopes <- function(sbp, hba1c, uacr, acei, sglt2, glp1) {
  base_aging <- 0.85
  pen_bp     <- pmax(0, (sbp - 120) * 0.022) * ifelse(acei, 0.72, 1.0)
  pen_a1c    <- pmax(0, (hba1c - 5.6) * 0.34) * ifelse(glp1, 0.78, 1.0)
  pen_uacr   <- pmax(0, log10(pmax(5, uacr) / 5) * 1.15) *
    ifelse(acei, 0.70, 1.0) * ifelse(sglt2, 0.58, 1.0)
  
  list(
    base       = base_aging,
    bp         = pen_bp,
    a1c        = pen_a1c,
    uacr       = pen_uacr,
    total_loss = base_aging + pen_bp + pen_a1c + pen_uacr
  )
}

# Helper to render FontAwesome Human Figure Breakout Pods with dynamic total figures
render_person_pod <- function(n_red, n_yellow, n_total = 10, cols = 6) {
  n_total  <- max(1, round(n_total))
  n_red    <- max(0, min(n_total, round(n_red)))
  n_yellow <- max(0, min(n_total - n_red, round(n_yellow)))
  n_green  <- max(0, n_total - n_red - n_yellow)
  states   <- c(rep("risk-mod", n_red), rep("risk-unknown", n_yellow), rep("risk-safe", n_green))
  
  div(
    style = sprintf("display:grid; grid-template-columns:repeat(%d, 1fr); gap:4px;", cols),
    lapply(seq_len(n_total), function(i) {
      ic <- if (states[i] == "risk-unknown") "person-circle-question" else "person"
      tags$span(class = paste("person-cell-sm", states[i]), shiny::icon(ic))
    })
  )
}

# -------------------------------------------------------------------------
# 4. USER INTERFACE (bslib)
# -------------------------------------------------------------------------
ui <- page_sidebar(
  title = "Cardiorenal Risk & Trajectory Explorer (AHA PREVENT™ + CKD-PC + NHANES)",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#1e293b"),
  
  tags$head(tags$style(HTML("
    .icon-grid-100 {
      display: grid; grid-template-columns: repeat(10, 1fr); gap: 5px;
      max-width: 340px; margin: 0 auto; padding: 10px;
      background: #f8fafc; border-radius: 8px; border: 1px solid #e2e8f0;
    }
    .person-cell { font-size: 1.35rem; text-align: center; transition: color 0.25s ease; }
    .person-cell-sm { font-size: 1.05rem; text-align: center; transition: color 0.25s ease; }
    .risk-locked  { color: #7f1d1d; }
    .risk-mod     { color: #ef4444; }
    .risk-unknown { color: #f59e0b; }
    .risk-safe    { color: #10b981; opacity: 0.78; }
    .pod-card {
      background: #f8fafc; border: 1px solid #e2e8f0;
      border-radius: 8px; padding: 6px 10px; margin-bottom: 6px;
    }
    .legend-dot { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 4px; }
  "))),
  
  sidebar = sidebar(
    width = 380,
    div(
      class = "d-flex justify-content-between align-items-center mb-2",
      tags$strong("Shared Clinical Controller"),
      actionButton("lock_baseline", "Lock Current as Baseline", class = "btn-sm btn-outline-primary")
    ),
    
    accordion(
      open = c("Tier 1: Demographics & Body Measures", "Tier 2 & 3: Biomarkers & Missing Flags"),
      
      accordion_panel(
        "Tier 1: Demographics & Body Measures",
        fluidRow(
          column(6, sliderInput("age", "Age (yrs)", 30, 79, 58, step = 1)),
          column(6, selectInput("sex", "Sex", c("Male", "Female"), selected = "Male"))
        ),
        fluidRow(
          column(7, selectInput("dm", "Glycemic Status", c("No Diabetes", "Prediabetes", "Diabetes"), selected = "Diabetes")),
          column(5, checkboxInput("smoker", "Smoker", value = FALSE))
        ),
        radioButtons("bmi_mode", "BMI Entry Mode:",
                     choices = c("Direct BMI Slider" = "direct", "Height & Weight" = "hw"), inline = TRUE),
        conditionalPanel(
          condition = "input.bmi_mode == 'direct'",
          sliderInput("bmi", "BMI (kg/m²)", 18, 48, 31.5, step = 0.5)
        ),
        conditionalPanel(
          condition = "input.bmi_mode == 'hw'",
          fluidRow(
            column(4, numericInput("ht_ft", "Feet", 5, min = 4, max = 7)),
            column(4, numericInput("ht_in", "Inches", 10, min = 0, max = 11)),
            column(4, numericInput("wt_lb", "Weight (lb)", 220, min = 80, max = 450))
          ),
          uiOutput("calc_bmi_text")
        )
      ),
      
      accordion_panel(
        "Tier 2 & 3: Biomarkers & Missing Flags",
        div(class = "d-flex justify-content-between",
            tags$label("Systolic BP (mmHg)"),
            checkboxInput("unk_sbp", "[? Unknown]", FALSE)),
        sliderInput("sbp", NULL, 90, 195, 148, step = 2),
        
        div(class = "d-flex justify-content-between",
            tags$label("Non-HDL Cholesterol (mg/dL)"),
            checkboxInput("unk_nhdl", "[? Unknown]", FALSE)),
        sliderInput("non_hdl", NULL, 60, 250, 165, step = 5),
        
        div(class = "d-flex justify-content-between",
            tags$label("HbA1c (%)"),
            checkboxInput("unk_a1c", "[? Unknown]", FALSE)),
        sliderInput("hba1c", NULL, 4.8, 11.5, 8.1, step = 0.1),
        
        div(class = "d-flex justify-content-between",
            tags$label("Renal Filtration (eGFR or Creatinine)"),
            checkboxInput("unk_egfr", "[? Unknown]", FALSE)),
        radioButtons("renal_mode", NULL,
                     choices = c("Direct eGFR" = "egfr", "Serum Creatinine (CKD-EPI)" = "scr"), inline = TRUE),
        conditionalPanel(
          condition = "input.renal_mode == 'egfr'",
          sliderInput("egfr", "Baseline eGFR (mL/min/1.73m²)", 15, 120, 58, step = 1)
        ),
        conditionalPanel(
          condition = "input.renal_mode == 'scr'",
          sliderInput("scr", "Serum Creatinine (mg/dL)", 0.4, 4.0, 1.35, step = 0.05),
          uiOutput("calc_egfr_text")
        ),
        
        div(class = "d-flex justify-content-between",
            tags$label("Urine Albumin-Creatinine Ratio (uACR, mg/g)"),
            checkboxInput("unk_uacr", "[? Unknown]", FALSE)),
        sliderInput("uacr", NULL, 5, 440, 180, step = 5)
      ),
      
      accordion_panel(
        "Target Therapies (What-If)",
        checkboxInput("statin", "Statin (Lipid-Lowering)", FALSE),
        checkboxInput("acei", "ACEi / ARB (BP & Glomerular Pressure)", FALSE),
        checkboxInput("sglt2", "SGLT2 Inhibitor (Renal & Heart Shield)", FALSE),
        checkboxInput("glp1", "GLP-1 RA (Cardiometabolic Shield)", FALSE)
      )
    )
  ),
  
  # TOP ROW: LEFT = ASCVD (5 cols), RIGHT = LINKED TRAJECTORY + CROSS-SECTION BAR (7 cols)
  layout_columns(
    col_widths = c(5, 7),
    
    # PANEL 1: 10x10 ICON ARRAY & DYNAMIC-DENOMINATOR BREAKOUT PODS
    card(
      full_screen = TRUE,
      card_header(
        class = "d-flex justify-content-between align-items-center",
        span("1. 10-Yr ASCVD Risk (AHA PREVENT™)"),
        uiOutput("ascvd_badge")
      ),
      card_body(
        layout_columns(
          col_widths = c(7, 5),
          div(
            uiOutput("ascvd_grid_ui"),
            div(
              class = "mt-2 small text-center",
              span(class = "legend-dot", style = "background:#7f1d1d;"), "Locked  ",
              span(class = "legend-dot", style = "background:#ef4444;"), "Modifiable  ",
              span(class = "legend-dot", style = "background:#f59e0b;"), "Unk (?)  ",
              span(class = "legend-dot", style = "background:#10b981;"), "Safe"
            )
          ),
          div(
            tags$strong(class = "small text-uppercase text-muted", "Dynamic Domain Pods (1 Fig = 1%)"),
            p(class = "small text-muted mb-1", style = "font-size:0.72rem;",
              "Total figures = domain ceiling; Red = active excess risk:"),
            uiOutput("breakout_pods_ui")
          )
        )
      )
    ),
    
    # PANEL 2: 20-YR TRAJECTORY + CROSS-SECTIONAL DECOMPOSITION BAR CHART
    card(
      full_screen = TRUE,
      card_header(
        class = "d-flex justify-content-between align-items-center",
        span("2. Renal Trajectory & Cross-Sectional Difference Composition (CKD-PC)"),
        div(
          class = "d-flex align-items-center gap-2",
          div(
            style = "width: 210px; margin-bottom: -14px;",
            sliderInput("slice_yr", NULL, min = 1, max = 20, value = 10, step = 1,
                        pre = "Cross-Section Year: ")
          ),
          uiOutput("ckd_badge")
        )
      ),
      card_body(
        layout_columns(
          col_widths = c(7, 5),
          girafeOutput("ckd_ribbon_plot", height = "365px"),
          girafeOutput("ckd_cross_section_bar", height = "365px")
        )
      )
    )
  ),
  
  # BOTTOM PANEL: GREEN / YELLOW / RED AUC NHANES DISTRIBUTIONS
  card(
    full_screen = TRUE,
    card_header(
      class = "d-flex justify-content-between align-items-center",
      span("3. Matched U.S. Peer Distributions (NHANES 2017–2023) with Clinical Risk AUC Zones"),
      uiOutput("nhanes_stratum_badge")
    ),
    card_body(
      girafeOutput("nhanes_auc_plot", height = "330px")
    )
  )
)

# -------------------------------------------------------------------------
# 5. SERVER
# -------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Click a year point on the trajectory plot to update the cross-section bar year
  observeEvent(input$ckd_ribbon_plot_selected, {
    sel_yr <- suppressWarnings(as.numeric(input$ckd_ribbon_plot_selected))
    if (!is.na(sel_yr) && sel_yr >= 1 && sel_yr <= 20) {
      updateSliderInput(session, "slice_yr", value = sel_yr)
    }
  })
  
  active_bmi <- reactive({
    if (input$bmi_mode == "hw") {
      ht_m  <- ((input$ht_ft * 12) + input$ht_in) * 0.0254
      wt_kg <- input$wt_lb * 0.453592
      round(wt_kg / (ht_m^2), 1)
    } else {
      input$bmi
    }
  })
  
  active_egfr <- reactive({
    if (input$renal_mode == "scr") {
      calc_ckdepi_2021(input$scr, input$age, input$sex)
    } else {
      input$egfr
    }
  })
  
  output$calc_bmi_text <- renderUI({
    div(class = "small text-primary fw-bold mb-2", sprintf("Calculated BMI: %0.1f kg/m²", active_bmi()))
  })
  
  output$calc_egfr_text <- renderUI({
    div(class = "small text-primary fw-bold mb-2", sprintf("2021 CKD-EPI eGFR: %0.1f mL/min/1.73m²", active_egfr()))
  })
  
  baseline_state <- reactiveVal(list(
    age = 58, sex = "Male", dm = "Diabetes", smoker = FALSE, bmi = 31.5,
    sbp = 148, non_hdl = 165, hba1c = 8.1, egfr = 58, uacr = 180,
    statin = FALSE, acei = FALSE, sglt2 = FALSE, glp1 = FALSE
  ))
  
  observeEvent(input$lock_baseline, {
    baseline_state(list(
      age = input$age, sex = input$sex, dm = input$dm, smoker = input$smoker,
      bmi = active_bmi(), sbp = input$sbp, non_hdl = input$non_hdl,
      hba1c = input$hba1c, egfr = active_egfr(), uacr = input$uacr,
      statin = input$statin, acei = input$acei, sglt2 = input$sglt2, glp1 = input$glp1
    ))
  })
  
  age_decade <- reactive({
    as.character(cut(input$age, breaks = c(29, 39, 49, 59, 69, 79),
                     labels = c("30-39", "40-49", "50-59", "60-69", "70-79")))
  })
  
  matched_nhanes <- reactive({
    nhanes_ref %>%
      filter(age_grp == age_decade(), sex == input$sex, dm_status == input$dm)
  })
  
  matched_cohort_microdata <- reactive({
    nhanes_cohort %>%
      filter(age_grp == age_decade(), sex == input$sex, dm_status == input$dm)
  })
  
  output$nhanes_stratum_badge <- renderUI({
    ref <- matched_nhanes()
    n_peers <- max(ref$n_unweighted, na.rm = TRUE)
    span(class = "badge bg-secondary",
         sprintf("Matched Stratum: %s | Age %s | %s (n = %d NHANES peers)",
                 input$sex, age_decade(), input$dm, n_peers))
  })
  
  resolved_inputs <- reactive({
    ref <- matched_nhanes()
    get_q <- function(var, col) {
      val <- ref %>% filter(biomarker == var) %>% pull(!!sym(col))
      if (length(val) == 0 || is.na(val[1])) 0 else val[1]
    }
    egfr_val <- active_egfr()
    
    list(
      bmi       = active_bmi(),
      sbp_med   = if (input$unk_sbp)  get_q("sbp", "p50")     else input$sbp,
      sbp_low   = if (input$unk_sbp)  get_q("sbp", "p10")     else input$sbp,
      sbp_high  = if (input$unk_sbp)  get_q("sbp", "p90")     else input$sbp,
      nhdl_med  = if (input$unk_nhdl) get_q("non_hdl", "p50") else input$non_hdl,
      nhdl_low  = if (input$unk_nhdl) get_q("non_hdl", "p10") else input$non_hdl,
      nhdl_high = if (input$unk_nhdl) get_q("non_hdl", "p90") else input$non_hdl,
      a1c_med   = if (input$unk_a1c)  get_q("hba1c", "p50")   else input$hba1c,
      a1c_low   = if (input$unk_a1c)  get_q("hba1c", "p10")   else input$hba1c,
      a1c_high  = if (input$unk_a1c)  get_q("hba1c", "p90")   else input$hba1c,
      egfr_med  = if (input$unk_egfr) get_q("egfr", "p50")    else egfr_val,
      egfr_low  = if (input$unk_egfr) get_q("egfr", "p90")    else egfr_val,
      egfr_high = if (input$unk_egfr) get_q("egfr", "p10")    else egfr_val,
      uacr_med  = if (input$unk_uacr) get_q("uacr", "p50")    else input$uacr,
      uacr_low  = if (input$unk_uacr) get_q("uacr", "p10")    else input$uacr,
      uacr_high = if (input$unk_uacr) get_q("uacr", "p90")    else input$uacr
    )
  })
  
  ascvd_metrics <- reactive({
    inp  <- resolved_inputs()
    base <- baseline_state()
    
    r_locked <- calc_ascvd_10yr(input$age, input$sex, 23, 110, 95, 5.1, 95, 5,
                                "No Diabetes", FALSE, FALSE, FALSE, FALSE, FALSE)
    r_low    <- calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_low, inp$nhdl_low,
                                inp$a1c_low, inp$egfr_low, inp$uacr_low, input$dm,
                                input$smoker, input$statin, input$acei, input$sglt2, input$glp1)
    r_med    <- calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med,
                                inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
                                input$smoker, input$statin, input$acei, input$sglt2, input$glp1)
    r_high   <- calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_high, inp$nhdl_high,
                                inp$a1c_high, inp$egfr_high, inp$uacr_high, input$dm,
                                input$smoker, input$statin, input$acei, input$sglt2, input$glp1)
    r_sq     <- calc_ascvd_10yr(base$age, base$sex, base$bmi, base$sbp, base$non_hdl,
                                base$hba1c, base$egfr, base$uacr, base$dm,
                                base$smoker, base$statin, base$acei, base$sglt2, base$glp1)
    
    n_locked  <- round(r_locked * 100)
    n_low     <- max(n_locked, round(r_low * 100))
    n_high    <- max(n_low, round(r_high * 100))
    n_mod     <- max(0, n_low - n_locked)
    n_unknown <- max(0, n_high - n_low)
    
    list(
      prob_med = round(r_med * 100, 1),
      prob_sq  = round(r_sq * 100, 1),
      n_locked = n_locked,
      n_mod    = n_mod,
      n_unk    = n_unknown,
      n_safe   = max(0, 100 - n_locked - n_mod - n_unknown)
    )
  })
  
  output$ascvd_badge <- renderUI({
    m <- ascvd_metrics()
    delta <- round(m$prob_med - m$prob_sq, 1)
    delta_txt <- if (delta < 0) sprintf(" (%0.1f%%)", delta) else if (delta > 0) sprintf(" (+%0.1f%%)", delta) else ""
    span(class = "badge bg-danger", sprintf("10-Yr: %0.1f%%%s", m$prob_med, delta_txt))
  })
  
  output$ascvd_grid_ui <- renderUI({
    m <- ascvd_metrics()
    states <- c(
      rep("risk-locked",  m$n_locked),
      rep("risk-mod",     m$n_mod),
      rep("risk-unknown", m$n_unk),
      rep("risk-safe",    m$n_safe)
    )
    div(
      class = "icon-grid-100",
      lapply(seq_len(100), function(i) {
        ic <- if (states[i] == "risk-unknown") "person-circle-question" else "person"
        tags$span(class = paste("person-cell", states[i]), shiny::icon(ic))
      })
    )
  })
  
  # Dynamic-Denominator Marginal Attribution Pods
  # Total figures (denominator) = Max physiological risk ceiling for that domain in this patient
  # Red figures (numerator)     = Current excess risk driven by the patient's active value
  output$breakout_pods_ui <- renderUI({
    inp <- resolved_inputs()
    
    # Baseline patient risk (using median/known values)
    r_full <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med,
      inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
      input$smoker, input$statin, input$acei, input$sglt2, input$glp1
    )
    
    # 1. OPTIMAL FLOOR: Risk if each domain were restored to ideal physiology
    r_opt_bp <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, 110, inp$nhdl_med,
      inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
      input$smoker, input$statin, FALSE, input$sglt2, input$glp1
    )
    r_opt_lip <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, inp$sbp_med, 95,
      inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
      input$smoker, FALSE, input$acei, input$sglt2, input$glp1
    )
    r_opt_gly <- calc_ascvd_10yr(
      input$age, input$sex, 23, inp$sbp_med, inp$nhdl_med,
      5.1, inp$egfr_med, inp$uacr_med, "No Diabetes",
      input$smoker, input$statin, input$acei, input$sglt2, FALSE
    )
    r_opt_ren <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med,
      inp$a1c_med, 95, 5, input$dm,
      FALSE, input$statin, input$acei, FALSE, input$glp1
    )
    
    # 2. CLINICAL CEILING (DYNAMIC DENOMINATOR): Untreated worst-case slider bound for each domain
    r_max_bp <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, 195, inp$nhdl_med,
      inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
      input$smoker, input$statin, FALSE, input$sglt2, input$glp1
    )
    r_max_lip <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, inp$sbp_med, 250,
      inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm,
      input$smoker, FALSE, input$acei, input$sglt2, input$glp1
    )
    r_max_gly <- calc_ascvd_10yr(
      input$age, input$sex, 45, inp$sbp_med, inp$nhdl_med,
      11.5, inp$egfr_med, inp$uacr_med, "Diabetes",
      input$smoker, input$statin, input$acei, input$sglt2, FALSE
    )
    r_max_ren <- calc_ascvd_10yr(
      input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med,
      inp$a1c_med, 20, 440, input$dm,
      TRUE, input$statin, input$acei, FALSE, input$glp1
    )
    
    # Current Marginal Risk Numerators (Red Icons)
    d_bp  <- max(0, (r_full - r_opt_bp)  * 100)
    d_lip <- max(0, (r_full - r_opt_lip) * 100)
    d_gly <- max(0, (r_full - r_opt_gly) * 100)
    d_ren <- max(0, (r_full - r_opt_ren) * 100)
    
    # Dynamic Domain Ceiling Denominators (Total Figures in Pod, minimum 1 figure)
    denom_bp  <- max(1, round(d_bp),  round((r_max_bp  - r_opt_bp)  * 100))
    denom_lip <- max(1, round(d_lip), round((r_max_lip - r_opt_lip) * 100))
    denom_gly <- max(1, round(d_gly), round((r_max_gly - r_opt_gly) * 100))
    denom_ren <- max(1, round(d_ren), round((r_max_ren - r_opt_ren) * 100))
    
    # 3. Exact 10th-to-90th NHANES uncertainty spread per domain when [? Unknown] is checked
    u_bp <- if (input$unk_sbp) {
      (calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_high, inp$nhdl_med, inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1) -
         calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_low,  inp$nhdl_med, inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1)) * 100
    } else 0
    
    u_lip <- if (input$unk_nhdl) {
      (calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_high, inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1) -
         calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_low,  inp$a1c_med, inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1)) * 100
    } else 0
    
    u_gly <- if (input$unk_a1c) {
      (calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med, inp$a1c_high, inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1) -
         calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med, inp$a1c_low,  inp$egfr_med, inp$uacr_med, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1)) * 100
    } else 0
    
    u_ren <- if (input$unk_egfr || input$unk_uacr) {
      (calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med, inp$a1c_med, inp$egfr_high, inp$uacr_high, input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1) -
         calc_ascvd_10yr(input$age, input$sex, inp$bmi, inp$sbp_med, inp$nhdl_med, inp$a1c_med, inp$egfr_low,  inp$uacr_low,  input$dm, input$smoker, input$statin, input$acei, input$sglt2, input$glp1)) * 100
    } else 0
    
    tagList(
      div(class = "pod-card",
          div(class = "d-flex justify-content-between small fw-bold mb-1",
              span("Blood Pressure"),
              span(class = "text-danger", sprintf("%d / %d figs (+%0.1f%%)", round(d_bp), denom_bp, d_bp))),
          render_person_pod(d_bp, u_bp, n_total = denom_bp, cols = 6)),
      
      div(class = "pod-card",
          div(class = "d-flex justify-content-between small fw-bold mb-1",
              span("Atherogenic Lipids"),
              span(class = "text-danger", sprintf("%d / %d figs (+%0.1f%%)", round(d_lip), denom_lip, d_lip))),
          render_person_pod(d_lip, u_lip, n_total = denom_lip, cols = 6)),
      
      div(class = "pod-card",
          div(class = "d-flex justify-content-between small fw-bold mb-1",
              span("Glycemic & BMI"),
              span(class = "text-danger", sprintf("%d / %d figs (+%0.1f%%)", round(d_gly), denom_gly, d_gly))),
          render_person_pod(d_gly, u_gly, n_total = denom_gly, cols = 6)),
      
      div(class = "pod-card",
          div(class = "d-flex justify-content-between small fw-bold mb-1",
              span("Renal & Smoking"),
              span(class = "text-danger", sprintf("%d / %d figs (+%0.1f%%)", round(d_ren), denom_ren, d_ren))),
          render_person_pod(d_ren, u_ren, n_total = denom_ren, cols = 6))
    )
  })
  
  output$ckd_badge <- renderUI({
    inp <- resolved_inputs()
    slp <- calc_egfr_slopes(inp$sbp_med, inp$a1c_med, inp$uacr_med,
                            input$acei, input$sglt2, input$glp1)
    span(class = "badge bg-primary", sprintf("Slope: -%0.2f/yr", slp$total_loss))
  })
  
  # -----------------------------------------------------------------------
  # PANEL 2A: 20-YEAR KIDNEY TRAJECTORY WITH VERTICAL CROSS-SECTION LINE
  # -----------------------------------------------------------------------
  output$ckd_ribbon_plot <- renderGirafe({
    inp      <- resolved_inputs()
    base     <- baseline_state()
    yr_slice <- input$slice_yr
    
    slp_curr <- calc_egfr_slopes(inp$sbp_med, inp$a1c_med, inp$uacr_med,
                                 input$acei, input$sglt2, input$glp1)
    slp_base <- calc_egfr_slopes(base$sbp, base$hba1c, base$uacr,
                                 base$acei, base$sglt2, base$glp1)
    slp_low  <- calc_egfr_slopes(inp$sbp_low, inp$a1c_low, inp$uacr_low,
                                 input$acei, input$sglt2, input$glp1)
    slp_high <- calc_egfr_slopes(inp$sbp_high, inp$a1c_high, inp$uacr_high,
                                 input$acei, input$sglt2, input$glp1)
    
    df_time <- tibble(year = 0:20) %>%
      mutate(
        e0              = inp$egfr_med,
        healthy_line    = pmax(5, e0 - slp_curr$base * year),
        wedge_bp_bot    = pmax(5, healthy_line - slp_curr$bp * year),
        wedge_a1c_bot   = pmax(5, wedge_bp_bot - slp_curr$a1c * year),
        target_line     = pmax(5, wedge_a1c_bot - slp_curr$uacr * year),
        status_quo_line = pmax(5, base$egfr - slp_base$total_loss * year),
        unc_upper       = pmax(5, inp$egfr_low  - slp_low$total_loss * year),
        unc_lower       = pmax(5, inp$egfr_high - slp_high$total_loss * year)
      )
    
    wedges_df <- bind_rows(
      df_time %>% transmute(year, ymax = healthy_line, ymin = wedge_bp_bot,
                            Wedge = "1. Excess BP Penalty",
                            tip = sprintf("BP Penalty: -%0.2f mL/min/yr", slp_curr$bp)),
      df_time %>% transmute(year, ymax = wedge_bp_bot, ymin = wedge_a1c_bot,
                            Wedge = "2. Glycemic (HbA1c) Penalty",
                            tip = sprintf("HbA1c Penalty: -%0.2f mL/min/yr", slp_curr$a1c)),
      df_time %>% transmute(year, ymax = wedge_a1c_bot, ymin = target_line,
                            Wedge = "3. Albuminuria (uACR) Stress",
                            tip = sprintf("uACR Glomerular Penalty: -%0.2f mL/min/yr", slp_curr$uacr))
    )
    
    any_unk <- input$unk_sbp || input$unk_a1c || input$unk_egfr || input$unk_uacr
    
    p <- ggplot() +
      annotate("rect", xmin = 0, xmax = 20, ymin = 0, ymax = 15, fill = "#fee2e2", alpha = 0.6) +
      annotate("text", x = 0.8, y = 7.5, label = "Kidney Failure / Dialysis (<15)",
               hjust = 0, size = 3.0, fontface = "bold", color = "#991b1b") +
      geom_ribbon_interactive(
        data = wedges_df,
        aes(x = year, ymin = ymin, ymax = ymax, fill = Wedge, tooltip = tip),
        alpha = 0.55
      )
    
    if (any_unk) {
      p <- p + geom_ribbon_interactive(
        data = df_time,
        aes(x = year, ymin = unc_lower, ymax = unc_upper),
        fill = "#f59e0b", alpha = 0.28,
        tooltip = "Epistemic Uncertainty Cone (10th-90th NHANES Percentiles)"
      )
    }
    
    p <- p +
      geom_vline(xintercept = yr_slice, color = "#1e293b", linewidth = 1.0, linetype = "longdash") +
      annotate("label", x = yr_slice, y = 113, label = sprintf("Year %d Slice ->", yr_slice),
               size = 3.1, fontface = "bold", fill = "#f1f5f9", color = "#0f172a") +
      geom_line(data = df_time, aes(x = year, y = healthy_line),
                color = "#10b981", linewidth = 1.0, linetype = "dashed") +
      geom_line(data = df_time, aes(x = year, y = status_quo_line),
                color = "#64748b", linewidth = 1.0, linetype = "dotted") +
      geom_line(data = df_time, aes(x = year, y = target_line),
                color = "#0f172a", linewidth = 1.3) +
      geom_point_interactive(
        data = filter(df_time, year >= 1),
        aes(x = year, y = target_line, data_id = as.character(year),
            tooltip = sprintf("Click to slice Cross-Section at Year %d (eGFR: %0.1f)", year, target_line)),
        size = 2.4, color = "#0f172a"
      ) +
      scale_fill_manual(values = c(
        "1. Excess BP Penalty"         = "#fb7185",
        "2. Glycemic (HbA1c) Penalty"  = "#c084fc",
        "3. Albuminuria (uACR) Stress" = "#fbbf24"
      )) +
      scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 20)) +
      labs(
        title = "A. 20-Yr Trajectory (Click Any Year Point)",
        x = "Years from Baseline (t = 0 to 20)",
        y = "eGFR (mL/min/1.73m²)",
        fill = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 11.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        panel.grid.minor = element_blank()
      )
    
    girafe(
      ggobj = p, width_svg = 4.6, height_svg = 3.9,
      options = list(
        opts_hover(css = "opacity:0.9;stroke:black;stroke-width:1.5px;"),
        opts_selection(type = "single", only_shiny = TRUE)
      )
    )
  })
  
  # -----------------------------------------------------------------------
  # PANEL 2B: CROSS-SECTIONAL STACKED BAR GRAPH AT SELECTED YEAR (t = yr)
  # -----------------------------------------------------------------------
  output$ckd_cross_section_bar <- renderGirafe({
    inp      <- resolved_inputs()
    base     <- baseline_state()
    yr_slice <- input$slice_yr
    e0       <- inp$egfr_med
    
    slp_curr <- calc_egfr_slopes(inp$sbp_med, inp$a1c_med, inp$uacr_med,
                                 input$acei, input$sglt2, input$glp1)
    slp_base <- calc_egfr_slopes(base$sbp, base$hba1c, base$uacr,
                                 base$acei, base$sglt2, base$glp1)
    
    decompose_at_year <- function(e_start, slp, scenario_label) {
      raw_aging <- slp$base * yr_slice
      raw_bp    <- slp$bp   * yr_slice
      raw_a1c   <- slp$a1c  * yr_slice
      raw_uacr  <- slp$uacr * yr_slice
      raw_tot   <- raw_aging + raw_bp + raw_a1c + raw_uacr
      
      rem_egfr    <- max(5, e_start - raw_tot)
      actual_drop <- e_start - rem_egfr
      scale_f     <- if (raw_tot > 0) actual_drop / raw_tot else 0
      
      tibble(
        Scenario = scenario_label,
        Component = factor(
          c("1. Excess BP Loss", "2. Glycemic (HbA1c) Loss", "3. Albuminuria (uACR) Loss",
            "4. Normal Aging (-0.85/yr)", "5. Preserved eGFR at Year t"),
          levels = c("1. Excess BP Loss", "2. Glycemic (HbA1c) Loss", "3. Albuminuria (uACR) Loss",
                     "4. Normal Aging (-0.85/yr)", "5. Preserved eGFR at Year t")
        ),
        ml_min = c(
          raw_bp    * scale_f,
          raw_a1c   * scale_f,
          raw_uacr  * scale_f,
          raw_aging * scale_f,
          rem_egfr
        )
      )
    }
    
    slp_healthy <- list(base = 0.85, bp = 0, a1c = 0, uacr = 0)
    
    bar_df <- bind_rows(
      decompose_at_year(e0, slp_healthy, "Healthy\nAging"),
      decompose_at_year(e0, slp_base,    "Status\nQuo"),
      decompose_at_year(e0, slp_curr,    "What-If\nTarget")
    ) %>%
      mutate(
        Scenario = factor(Scenario, levels = c("Healthy\nAging", "Status\nQuo", "What-If\nTarget")),
        tip = sprintf("%s (%s at Year %d): %0.1f mL/min/1.73m²",
                      gsub("\n", " ", Scenario), Component, yr_slice, ml_min)
      )
    
    rem_sq   <- bar_df %>% filter(Scenario == "Status\nQuo",     Component == "5. Preserved eGFR at Year t") %>% pull(ml_min)
    rem_tgt  <- bar_df %>% filter(Scenario == "What-If\nTarget", Component == "5. Preserved eGFR at Year t") %>% pull(ml_min)
    net_gain <- round(rem_tgt - rem_sq, 1)
    gain_sub <- if (net_gain > 0) {
      sprintf("Net Kidney Function Saved at Yr %d: +%0.1f mL/min", yr_slice, net_gain)
    } else if (net_gain < 0) {
      sprintf("Net Difference vs. Status Quo: %0.1f mL/min", net_gain)
    } else {
      sprintf("Modify sliders/meds to see Year %d eGFR gain", yr_slice)
    }
    
    preserved_labels <- bar_df %>%
      filter(Component == "5. Preserved eGFR at Year t")
    
    p_bar <- ggplot(bar_df, aes(x = Scenario, y = ml_min, fill = Component)) +
      annotate("rect", xmin = 0.4, xmax = 3.6, ymin = 0, ymax = 15, fill = "#fee2e2", alpha = 0.5) +
      geom_col_interactive(aes(tooltip = tip, data_id = Component), width = 0.62, color = "#ffffff", linewidth = 0.3) +
      geom_hline(yintercept = e0, linetype = "dotted", color = "#475569", linewidth = 0.7) +
      geom_text(
        data = preserved_labels,
        aes(x = Scenario, y = pmax(8, ml_min / 2), label = sprintf("%0.1f", ml_min)),
        inherit.aes = FALSE, color = "#ffffff", fontface = "bold", size = 3.4
      ) +
      scale_fill_manual(
        values = c(
          "1. Excess BP Loss"           = "#fb7185",
          "2. Glycemic (HbA1c) Loss"    = "#c084fc",
          "3. Albuminuria (uACR) Loss"  = "#fbbf24",
          "4. Normal Aging (-0.85/yr)"  = "#cbd5e1",
          "5. Preserved eGFR at Year t" = "#0f766e"
        ),
        guide = guide_legend(nrow = 3)
      ) +
      scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 20)) +
      labs(
        title    = sprintf("B. Year %d Cross-Section Breakdown", yr_slice),
        subtitle = gain_sub,
        x        = NULL,
        y        = NULL,
        fill     = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 9.5, color = "#0f766e", face = "bold"),
        legend.position = "bottom",
        legend.text   = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    girafe(
      ggobj = p_bar, width_svg = 3.6, height_svg = 3.9,
      options = list(opts_hover(css = "opacity:0.9;stroke:#0f172a;stroke-width:1.2px;"))
    )
  })
  
  # -----------------------------------------------------------------------
  # PANEL 3: SURVEY-WEIGHTED NHANES DISTRIBUTIONS WITH GREEN/YELLOW/RED AUC
  # -----------------------------------------------------------------------
  output$nhanes_auc_plot <- renderGirafe({
    sub_df <- matched_cohort_microdata()
    ref    <- matched_nhanes()
    inp    <- resolved_inputs()
    
    dens_all <- map_dfr(seq_len(nrow(biomarker_meta)), function(i) {
      build_auc_density(sub_df, biomarker_meta[i, ])
    })
    
    user_markers <- biomarker_meta %>%
      mutate(
        val         = c(inp$bmi, inp$sbp_med, inp$nhdl_med, inp$a1c_med, inp$egfr_med, inp$uacr_med),
        val_clamped = pmin(xmax, pmax(xmin, val)),
        is_unk      = c(FALSE, input$unk_sbp, input$unk_nhdl, input$unk_a1c, input$unk_egfr, input$unk_uacr)
      ) %>%
      left_join(ref, by = "biomarker") %>%
      mutate(
        tip = sprintf(
          "%s\nPatient Value: %0.1f %s\nNHANES Peer Median: %0.1f (IQR: %0.1f–%0.1f | 10th–90th: %0.1f–%0.1f)",
          as.character(label), val, ifelse(is_unk, "[Imputed Peer Median]", ""),
          p50, p25, p75, p10, p90
        )
      )
    
    p <- ggplot() +
      geom_ribbon_interactive(
        data = dens_all,
        aes(x = x, ymin = 0, ymax = y, fill = zone, group = interaction(label, zone),
            tooltip = paste0(label, " — ", zone)),
        alpha = 0.68
      ) +
      geom_line(data = dens_all, aes(x = x, y = y, group = label),
                color = "#1e293b", linewidth = 0.55) +
      geom_vline(data = user_markers, aes(xintercept = p50),
                 color = "#334155", linetype = "dotted", linewidth = 0.75) +
      geom_vline_interactive(
        data = user_markers,
        aes(xintercept = val_clamped, tooltip = tip),
        color = "#0f172a", linewidth = 1.25
      ) +
      geom_point_interactive(
        data = user_markers,
        aes(x = val_clamped, y = 1.03, shape = is_unk, tooltip = tip),
        size = 3.6, color = "#0f172a", fill = "#ffffff", stroke = 1.2
      ) +
      facet_wrap(~ label, scales = "free_x", ncol = 3) +
      scale_fill_manual(
        name = "Clinical Risk AUC Zone:",
        values = c(
          "Optimal (Green)"                = "#10b981",
          "Borderline / Elevated (Yellow)" = "#f59e0b",
          "High Clinical Risk (Red)"       = "#ef4444"
        )
      ) +
      scale_shape_manual(
        name = "Patient Marker:",
        values = c("FALSE" = 19, "TRUE" = 21),
        labels = c("FALSE" = "Known Patient Input (Solid Line)", "TRUE" = "Imputed NHANES Median ([? Unknown])")
      ) +
      scale_y_continuous(limits = c(0, 1.15), breaks = NULL) +
      labs(
        x = NULL, y = "Relative Peer Frequency",
        caption = "Shaded AUC = Survey-weighted NHANES peer distribution colored by clinical guidelines | Dotted line = Peer Median (50th) | Solid line = Patient Value"
      ) +
      theme_minimal(base_size = 11.5) +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      )
    
    girafe(ggobj = p, width_svg = 11.5, height_svg = 3.5,
           options = list(opts_hover(css = "opacity:0.9;stroke:#0f172a;stroke-width:1.5px;")))
  })
}

shinyApp(ui, server)