# ==============================================================
# Study 1 — Equity-Adjusted Thresholds for Safe Maternal Delivery
# Lancet GLH–styled outputs (figures + tables to ./output/)
# ==============================================================

# ---- 0) Packages ----
req <- c(
  "tidyverse","data.table","fixest","modelsummary","splines",
  "broom","ggrepel","sandwich","clubSandwich","dagitty"
)
new <- req[!req %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

source("lancet_style_helpers.R")

# ---- 1) Paths & setup ----
panel_path <- "equity_thresholds_panel.csv"   # set to your merged CSV
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# Set global theme
theme_set(theme_lancet_glh())

# ---- 2) Utilities ----
complete_for <- function(df, vars){
  df %>% dplyr::filter(dplyr::if_all(all_of(vars), ~ !is.na(.x)))
}
winsor <- function(x, p=0.01){
  q <- stats::quantile(x, probs = c(p, 1-p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
# Threshold search (Hansen-style grid) with NA-safety
threshold_search <- function(df, grid_g = seq(10,95,by=1), outcome = "lmmr",
                             expvar = "institutional_birth_pct",
                             covars = c("anemia_prev_pct","adolescent_birth_rate"),
                             fe = c("country_code","year")){
  need <- c(outcome, expvar, covars, fe)
  d <- df %>% complete_for(need)
  if(nrow(d) < 50) return(list(gamma = NA_real_, prof = tibble(gamma=grid_g, SSR=NA_real_)))
  ssr <- sapply(grid_g, function(g){
    d2 <- d %>% mutate(x_low = pmin(.data[[expvar]], g),
                       x_high = pmax(.data[[expvar]] - g, 0))
    fit <- tryCatch(
      fixest::feols(as.formula(paste(outcome, "~ x_low + x_high +",
                                     paste(covars, collapse = " + "), "|",
                                     paste(fe, collapse = " + "))),
                    data = d2, vcov = ~ country_code),
      error = function(e) NULL
    )
    if(is.null(fit)) return(NA_real_)
    sum(residuals(fit)^2, na.rm=TRUE)
  })
  if(all(is.na(ssr))) return(list(gamma=NA_real_, prof=tibble(gamma=grid_g, SSR=ssr)))
  ghat <- grid_g[which.min(ssr)]
  list(gamma = ghat, prof = tibble(gamma = grid_g, SSR = ssr))
}

# Finite-difference marginal effect for region-interaction model
fd_me_region <- function(fit, region, x, med_anemia, med_adbr, data_ref, h=0.1){
  base_row <- tibble(
    institutional_birth_pct = c(x - h, x + h),
    anemia_prev_pct         = med_anemia,
    adolescent_birth_rate   = med_adbr,
    who_region              = region,
    country_code            = data_ref$country_code[1],
    year                    = data_ref$year[1]
  )
  p <- as.numeric(predict(fit, newdata = base_row))
  (p[2] - p[1]) / (2*h)  # derivative on log scale per 1 pt coverage
}

# ---- 3) Load & transform ----
dt <- data.table::fread(panel_path, na.strings = c("", "NA", "NaN")) |> as_tibble()
dt <- dt |> filter(year >= 2000, year <= 2022)

need_vars <- c("country_code","country_name","who_region","year",
               "maternal_mortality_ratio","institutional_birth_pct",
               "anemia_prev_pct","adolescent_birth_rate")
miss <- setdiff(need_vars, names(dt))
if(length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

dt <- dt |> mutate(
  mmr_w = winsor(maternal_mortality_ratio, 0.01),
  lmmr  = log(mmr_w + 1)
)

# ---- 4) Table 1: Overview ----
tab1 <- dt |>
  summarise(
    rows      = n(),
    countries = n_distinct(country_code),
    years_min = min(year, na.rm=TRUE),
    years_max = max(year, na.rm=TRUE),
    mmr_med   = median(maternal_mortality_ratio, na.rm=TRUE),
    inst_med  = median(institutional_birth_pct, na.rm=TRUE),
    anemia_med= median(anemia_prev_pct, na.rm=TRUE),
    adbr_med  = median(adolescent_birth_rate, na.rm=TRUE)
  )
save_table_lancet(tab1, "table01_overview", to_docx = TRUE)

# ---- 5) Table 2: Availability by region ----
tab2 <- dt |>
  group_by(who_region) |>
  summarise(
    maternal_mortality_ratio_nonmissing   = sum(!is.na(maternal_mortality_ratio)),
    institutional_birth_pct_nonmissing    = sum(!is.na(institutional_birth_pct)),
    anemia_prev_pct_nonmissing            = sum(!is.na(anemia_prev_pct)),
    adolescent_birth_rate_nonmissing      = sum(!is.na(adolescent_birth_rate)),
    n_rows = n(), .groups="drop"
  ) |>
  arrange(who_region)
save_table_lancet(tab2, "table02_availability", to_docx = TRUE)

# ---- 6) Fig 1: Binned association by region ----
fig1 <- dt |>
  filter(!is.na(institutional_birth_pct), !is.na(maternal_mortality_ratio)) |>
  mutate(bin = floor(institutional_birth_pct/5)*5) |>
  group_by(who_region, bin) |>
  summarise(mmr = median(maternal_mortality_ratio, na.rm=TRUE), .groups="drop") |>
  ggplot(aes(bin, mmr, group = who_region)) +
  geom_line(linewidth = 0.6, color = "#0B2C40") +
  facet_wrap(~ who_region, scales = "free_y") +
  labs(x = "Institutional birth coverage (5% bins, %)",
       y = "Median MMR (per 100,000)",
       title = "Institutional birth coverage vs MMR (binned by WHO region)") +
  theme_lancet_glh()
save_figure_lancet(fig1, "fig01_binned_by_region", width_mm = 180, height_mm = 120)

# ---- 7) Table 3: TWFE + spline (primary) ----
need_main <- c("lmmr","institutional_birth_pct","anemia_prev_pct","adolescent_birth_rate","country_code","year")
d_main <- complete_for(dt, need_main)
if(nrow(d_main) == 0) stop("No complete cases for primary model.")

m_spline <- feols(
  lmmr ~ splines::ns(institutional_birth_pct, df = 3) +
    anemia_prev_pct + adolescent_birth_rate |
    country_code + year,
  data = d_main,
  vcov = ~ country_code
)
save_table_lancet(list("TWFE + spline (log MMR)" = m_spline), "table03_twfe_spline", to_docx = TRUE)

# ---- 8) Fig 2: Exposure–response (spline) ----
grid <- tibble(
  institutional_birth_pct = seq(0, 100, by = 1),
  anemia_prev_pct         = median(d_main$anemia_prev_pct, na.rm=TRUE),
  adolescent_birth_rate   = median(d_main$adolescent_birth_rate, na.rm=TRUE),
  country_code            = d_main$country_code[1],
  year                    = d_main$year[1]
)
grid$pred_lmmr <- predict(m_spline, newdata = grid, se.fit = FALSE)
grid$pred_mmr  <- exp(grid$pred_lmmr) - 1

fig2 <- ggplot(grid, aes(institutional_birth_pct, pred_mmr)) +
  geom_line(linewidth = 0.8, color = "#007F5F") +
  labs(x = "Institutional birth coverage (%)",
       y = "Predicted MMR (per 100,000)",
       title = "Predicted MMR across institutional birth coverage (TWFE + spline)") +
  theme_lancet_glh()
save_figure_lancet(fig2, "fig02_exposure_response", width_mm = 180, height_mm = 120)

# ---- 9) Table 4: Region interactions + finite-diff marginal effects ----
need_lin <- c("lmmr","institutional_birth_pct","anemia_prev_pct","adolescent_birth_rate","country_code","year","who_region")
d_lin <- complete_for(dt, need_lin)
if(nrow(d_lin) > 0){
  m_region_lin <- feols(
    lmmr ~ institutional_birth_pct * i(who_region) +
      anemia_prev_pct + adolescent_birth_rate |
      country_code + year,
    data = d_lin,
    vcov = ~ country_code
  )
  save_table_lancet(list("TWFE + region interactions (linear)" = m_region_lin),
                    "table04_region_interactions", to_docx = TRUE)
  
  rep_levels <- c(40, 60, 80)
  med_anemia <- median(d_lin$anemia_prev_pct, na.rm = TRUE)
  med_adbr   <- median(d_lin$adolescent_birth_rate, na.rm = TRUE)
  regions    <- sort(unique(d_lin$who_region))
  
  ME_list <- list()
  for (r in regions){
    for (xc in rep_levels){
      me <- tryCatch(fd_me_region(m_region_lin, r, xc, med_anemia, med_adbr, d_lin), error = function(e) NA_real_)
      ME_list[[length(ME_list)+1]] <- tibble(who_region = r, at_coverage = xc, dydx_logMMR_per1pct = me)
    }
  }
  ME_out <- bind_rows(ME_list) |>
    arrange(who_region, at_coverage)
  save_table_lancet(ME_out, "table04b_marginal_effects_finite_diff", to_docx = TRUE)
} else {
  message("Skipping region interactions: insufficient complete cases.")
}

# ---- 10) Tables 5–6 + Fig 3–4: Thresholds ----
need_thr <- c("lmmr","institutional_birth_pct","anemia_prev_pct","adolescent_birth_rate","country_code","year")
d_thr <- complete_for(dt, need_thr)

thr_res  <- threshold_search(d_thr)
thr_hat  <- thr_res$gamma
ssr_prof <- thr_res$prof

fig3 <- ssr_prof |>
  ggplot(aes(gamma, SSR)) +
  geom_line(linewidth = 0.7, color = "#6C1D5F", na.rm=TRUE) +
  {if(!is.na(thr_hat)) geom_vline(xintercept = thr_hat, linetype = 2, linewidth = 0.6)} +
  labs(x = "Candidate threshold (institutional birth, %)",
       y = "Sum of squared residuals",
       title = "Hansen-style threshold search (global)") +
  theme_lancet_glh()
save_figure_lancet(fig3, "fig03_threshold_ssr", width_mm = 140, height_mm = 100)

if(!is.na(thr_hat)){
  dt_thr <- d_thr |>
    mutate(x_low = pmin(institutional_birth_pct, thr_hat),
           x_high = pmax(institutional_birth_pct - thr_hat, 0))
  m_thr <- feols(
    lmmr ~ x_low + x_high + anemia_prev_pct + adolescent_birth_rate |
      country_code + year,
    data = dt_thr,
    vcov = ~ country_code
  )
  save_table_lancet(list("TWFE + threshold (log MMR)" = m_thr), "table05_threshold_model", to_docx = TRUE)
} else {
  message("No valid global threshold estimated — skipping Table 5.")
}

region_thr <- dt |>
  group_by(who_region) |>
  group_modify(~{
    d <- complete_for(.x, need_thr)
    if(nrow(d) < 50) return(tibble(threshold_pct = NA_real_, SSR_min = NA_real_))
    tr <- threshold_search(d)
    tibble(threshold_pct = tr$gamma, SSR_min = suppressWarnings(min(tr$prof$SSR, na.rm=TRUE)))
  }) |>
  ungroup()
save_table_lancet(region_thr, "table06_region_specific_thresholds", to_docx = TRUE)

if(nrow(region_thr) > 0 && any(!is.na(region_thr$threshold_pct))){
  fig4 <- region_thr |>
    ggplot(aes(reorder(who_region, threshold_pct), threshold_pct)) +
    geom_point(size = 2, color = "#C14C36", na.rm = TRUE) +
    coord_flip() +
    labs(x = "WHO region", y = "Estimated threshold (%)",
         title = "Region-specific thresholds for institutional birth coverage") +
    theme_lancet_glh()
  save_figure_lancet(fig4, "fig04_region_thresholds", width_mm = 140, height_mm = 100)
}

# ---- 11) Table 7: Sensitivity on raw MMR ----
need_raw <- c("maternal_mortality_ratio","institutional_birth_pct","anemia_prev_pct","adolescent_birth_rate","country_code","year")
d_raw <- complete_for(dt, need_raw)
if(nrow(d_raw) > 0){
  m_raw <- feols(
    maternal_mortality_ratio ~ splines::ns(institutional_birth_pct, df = 3) +
      anemia_prev_pct + adolescent_birth_rate |
      country_code + year,
    data = d_raw, vcov = ~ country_code
  )
  save_table_lancet(list("TWFE + spline (raw MMR)" = m_raw), "table07_sensitivity_raw_mmr", to_docx = TRUE)
}

# ---- 12) Table 8: Sensitivity + SBA ----
if("skilled_birth_attendance_pct" %in% names(dt)){
  need_sba <- c("lmmr","institutional_birth_pct","skilled_birth_attendance_pct",
                "anemia_prev_pct","adolescent_birth_rate","country_code","year")
  d_sba <- complete_for(dt, need_sba)
  if(nrow(d_sba) > 0){
    m_sba <- feols(
      lmmr ~ splines::ns(institutional_birth_pct, df = 3) +
        skilled_birth_attendance_pct + anemia_prev_pct + adolescent_birth_rate |
        country_code + year,
      data = d_sba, vcov = ~ country_code
    )
    save_table_lancet(list("TWFE + spline + SBA (log MMR)" = m_sba), "table08_sensitivity_add_sba", to_docx = TRUE)
  }
}

# ---- 13) Fig 5: DAG ----
dag_txt <- "
dag {
  MMR [outcome];
  InstBirth [exposure];
  Anemia;
  AgeStruct;
  Region;
  Income;

  InstBirth -> MMR
  Anemia -> MMR
  AgeStruct -> MMR

  Region -> InstBirth
  Region -> MMR

  Income -> InstBirth
  Income -> MMR
}
"
g <- dagitty(dag_txt)
png(file.path("output","figures","fig05_dag.png"), width = 1800, height = 1200, res = 300)
plot(g)
dev.off()

cat("\nLancet-styled outputs written to output/tables/ and output/figures/\n")
