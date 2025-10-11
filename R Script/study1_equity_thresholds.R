# ===========================
# Study 1: Equity-Adjusted Thresholds for Safe Maternal Delivery
# RStudio end-to-end script
# ===========================

# ---- 0) Packages ----
req <- c(
  "tidyverse","data.table","fixest","modelsummary","splines",
  "broom","broom.helpers","ggplot2","ggrepel","sandwich",
  "clubSandwich","marginaleffects","dagitty"
)
new <- req[!req %in% installed.packages()[,"Package"]]
if(length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ---- 1) Load data ----
# Set your path to the CSV (downloaded from ChatGPT)
# Option A (interactive): panel_path <- file.choose()
# Option B: manual path
panel_path <- "equity_thresholds_panel.csv"

dt <- fread(panel_path, na.strings = c("", "NA", "NaN")) |>
  as_tibble()

# Expect columns:
# country_code, country_name, region_code, who_region, year,
# maternal_mortality_ratio, institutional_birth_pct,
# skilled_birth_attendance_pct, anemia_prev_pct, adolescent_birth_rate

# Keep analysis window
dt <- dt |> filter(year >= 2000, year <= 2022)

# Basic sanity checks
skim_vars <- c("institutional_birth_pct","skilled_birth_attendance_pct",
               "anemia_prev_pct","adolescent_birth_rate",
               "maternal_mortality_ratio")
summary(dt[skim_vars])

# Optional: winsorize extreme MMR values (robustness)
winsor <- function(x, p=0.01){
  q <- quantile(x, probs = c(p, 1-p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
dt <- dt |> mutate(mmr_w = winsor(maternal_mortality_ratio, 0.01))

# ---- 2) Core TWFE model with nonlinear exposure (splines) ----
# We’ll log-transform MMR for scale; add 3-knot natural spline for inst. births
dt <- dt |>
  mutate(lmmr = log(mmr_w + 1)) # +1 guard in case of zeros (rare for MMR)

# Base covariates (can expand): anemia, adolescent fertility (proxy for age structure)
# Cluster-robust SE at country level
m_spline <- feols(
  lmmr ~ ns(institutional_birth_pct, df = 3) +
    anemia_prev_pct + adolescent_birth_rate |
    country_code + year,
  data = dt,
  vcov = ~ country_code
)

msummary(list("TWFE + spline (log MMR)" = m_spline), gof_omit = "IC|Log|Adj|AIC|BIC")

# Marginal effect curve over coverage
grid <- tibble(
  institutional_birth_pct = seq(0, 100, by = 1),
  anemia_prev_pct = median(dt$anemia_prev_pct, na.rm=TRUE),
  adolescent_birth_rate = median(dt$adolescent_birth_rate, na.rm=TRUE),
  # pick a baseline FE via reference country/year (dropped in prediction)
  country_code = dt$country_code[1],
  year = dt$year[1]
)

# Predict on the feols model (fixest handles FEs if supplied with reference levels)
grid$pred_lmmr <- predict(m_spline, newdata = grid, se.fit = FALSE)
grid$pred_mmr <- exp(grid$pred_lmmr) - 1

ggplot(grid, aes(institutional_birth_pct, pred_mmr)) +
  geom_line() +
  labs(x = "Institutional birth coverage (%)",
       y = "Predicted Maternal Mortality Ratio (per 100,000)",
       title = "Nonlinear relationship between institutional births and MMR") +
  theme_minimal()

# ---- 3) Region-specific marginal effects (via interaction) ----
# Allow slope to vary by region using interactions in a flexible way (binspline)
# Simpler: linear-by-region first for readable marginal effects
m_region_lin <- feols(
  lmmr ~ institutional_birth_pct * i(who_region) +
    anemia_prev_pct + adolescent_birth_rate |
    country_code + year,
  data = dt,
  vcov = ~ country_code
)
msummary(list("TWFE + region interactions (linear)" = m_region_lin),
         gof_omit = "IC|Log|Adj|AIC|BIC")

# Compute region-specific marginal effects at representative coverage levels
rep_levels <- c(40, 60, 80)
ME <- map_df(rep_levels, function(xc){
  marginaleffects::marginaleffects(
    m_region_lin,
    variables = "institutional_birth_pct",
    newdata = datagrid(
      institutional_birth_pct = xc,
      anemia_prev_pct = median(dt$anemia_prev_pct, na.rm=TRUE),
      adolescent_birth_rate = median(dt$adolescent_birth_rate, na.rm=TRUE),
      who_region = unique(dt$who_region)
    ),
    vcov = sandwich::vcovCL(m_region_lin, cluster = ~ country_code)
  ) |>
    as_tibble() |>
    mutate(at_coverage = xc)
})
ME |> arrange(who_region, at_coverage) |> select(who_region, at_coverage, dydx, conf.low, conf.high)

# ---- 4) Threshold search (Hansen-style, single threshold) ----
# We perform a simple grid search over candidate gamma in [10..95] percent coverage.
# For each gamma, fit a TWFE with piecewise linear terms:
#   x_low  = pmin(inst_cov, gamma)
#   x_high = pmax(inst_cov - gamma, 0)
# and record the sum of squared residuals (SSR). Pick gamma minimizing SSR.
grid_g <- seq(10, 95, by = 1)
ssr <- numeric(length(grid_g))

build_terms <- function(x, g){
  tibble(x_low = pmin(x, g), x_high = pmax(x - g, 0))
}

for(i in seq_along(grid_g)){
  g <- grid_g[i]
  tmp <- dt |> mutate(
    x_low = pmin(institutional_birth_pct, g),
    x_high = pmax(institutional_birth_pct - g, 0)
  )
  fit <- feols(
    lmmr ~ x_low + x_high + anemia_prev_pct + adolescent_birth_rate |
      country_code + year,
    data = tmp,
    vcov = ~ country_code
  )
  ssr[i] <- sum(residuals(fit)^2, na.rm = TRUE)
}

thresh_hat <- grid_g[which.min(ssr)]
cat("Estimated threshold (gamma*) at ~", thresh_hat, "% institutional birth coverage\n")

# Refit at threshold for interpretation
dt_thr <- dt |>
  mutate(x_low = pmin(institutional_birth_pct, thresh_hat),
         x_high = pmax(institutional_birth_pct - thresh_hat, 0))

m_thr <- feols(
  lmmr ~ x_low + x_high + anemia_prev_pct + adolescent_birth_rate |
    country_code + year,
  data = dt_thr,
  vcov = ~ country_code
)
msummary(list("TWFE + threshold (log MMR)" = m_thr), gof_omit = "IC|Log|Adj|AIC|BIC")

# Visualize SSR profile
data.frame(gamma = grid_g, SSR = ssr) |>
  ggplot(aes(gamma, SSR)) +
  geom_line() +
  geom_vline(xintercept = thresh_hat, linetype = 2) +
  labs(x = "Candidate threshold (institutional birth %)", y = "Sum of squared residuals",
       title = "Threshold search profile (Hansen-style grid search)") +
  theme_minimal()

# Optional (more rigorous): nonparametric bootstrap to form CI for gamma
# (Note: can take time; reduce B for quick checks)
# set.seed(123)
# B <- 200
# gamma_boot <- numeric(B)
# countries <- unique(dt$country_code)
# for(b in 1:B){
#   # block bootstrap by country
#   samp <- sample(countries, replace = TRUE)
#   bt <- map_dfr(samp, ~ dt |> filter(country_code == .x))
#   ssr_b <- sapply(grid_g, function(g){
#     dtmp <- bt |> mutate(
#       x_low = pmin(institutional_birth_pct, g),
#       x_high = pmax(institutional_birth_pct - g, 0)
#     )
#     fit <- feols(lmmr ~ x_low + x_high + anemia_prev_pct + adolescent_birth_rate |
#                    country_code + year, data = dtmp, vcov = ~ country_code)
#     sum(residuals(fit)^2, na.rm = TRUE)
#   })
#   gamma_boot[b] <- grid_g[which.min(ssr_b)]
# }
# quantile(gamma_boot, c(.025,.5,.975))

# ---- 5) Region-specific thresholds (optional extension)
# Repeat grid search within each WHO region
region_thr <- dt |>
  group_by(who_region) |>
  group_modify(~{
    d <- .x
    if(n_distinct(d$country_code) < 5) return(tibble(gamma = NA_real_, SSR = NA_real_))
    ssr_r <- sapply(grid_g, function(g){
      d2 <- d |> mutate(
        x_low = pmin(institutional_birth_pct, g),
        x_high = pmax(institutional_birth_pct - g, 0)
      )
      fit <- feols(lmmr ~ x_low + x_high + anemia_prev_pct + adolescent_birth_rate |
                     country_code + year, data = d2, vcov = ~ country_code)
      sum(residuals(fit)^2, na.rm = TRUE)
    })
    tibble(gamma = grid_g, SSR = ssr_r)
  }) |>
  group_by(who_region) |>
  slice_min(SSR, with_ties = FALSE) |>
  ungroup() |>
  rename(threshold_pct = gamma)

region_thr

# ---- 6) Simple DAG (for manuscript transparency)
dag_txt <- "
dag {
  MMR [outcome];
  InstBirth [exposure];
  Anemia;
  AgeStruct;
  Region;
  Income;  # If/when you add WB income mapping

  InstBirth -> MMR
  Anemia -> MMR
  AgeStruct -> MMR

  Region -> InstBirth
  Region -> MMR

  Income -> InstBirth
  Income -> MMR

  # Unobserved shocks (time FE handle global secular changes)
}
"
g <- dagitty(dag_txt)
plot(g)
adjustmentSets(g, exposure = "InstBirth", outcome = "MMR")
