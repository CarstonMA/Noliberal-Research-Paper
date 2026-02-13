# Neoliberal Policy Event Study Analysis
# DiD event study: effects of sharp policy score increases on economic outcomes
# - Standard errors clustered by country (ISO3)
# - Parallel trends tested via joint F-test on pre-treatment coefficients
# Uses merged_full_data.csv (run merge_full_data.R first)

library(tidyverse)
library(fixest)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
if (!file.exists("merged_full_data.csv")) {
  source("merge_full_data.R")
}
df <- read_csv("merged_full_data.csv", show_col_types = FALSE)
df <- df %>% filter(Year >= 1995, Year <= 2022, !is.na(ISO3), !is.na(Year))

# Panel: keep countries with at least 15 observations (avoids losing Gini/Poverty)
n_years <- length(unique(df$Year))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)
cat("Panel: ", length(panel_iso3), " countries (min 15 obs), ", n_years, " years\n", sep = "")

# ============================================================================
# 2. COUNTRY GROUPS (same as Thesis Data/neoliberal_event_study.R)
# ============================================================================
post_soviet <- c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", "LTU",
                 "MDA", "RUS", "TJK", "TKM", "UKR", "UZB")
latin_america <- c("ARG", "BOL", "BRA", "BRB", "BLZ", "CHL", "COL", "CRI", "CUB",
                  "DOM", "ECU", "SLV", "GTM", "GUY", "HND", "JAM", "MEX", "NIC",
                  "PAN", "PRY", "PER", "TTO", "URY", "VEN")
sub_saharan_africa <- c("AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CAF", "TCD",
                        "COM", "COG", "COD", "CIV", "GNQ", "ERI", "SWZ", "ETH",
                        "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", "LBR",
                        "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM", "NER",
                        "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM", "ZAF",
                        "SSD", "SDN", "TZA", "TGO", "UGA", "ZMB", "ZWE")
east_asia_pacific <- c("AUS", "BGD", "BTN", "KHM", "CHN", "FJI", "HKG", "IDN",
                      "JPN", "KIR", "LAO", "MYS", "MMR", "MNG", "NPL", "NZL",
                      "PHL", "PNG", "KOR", "WSM", "SGP", "SLB", "LKA", "THA",
                      "TON", "VUT", "VNM")
middle_east_na <- c("DZA", "BHR", "EGY", "IRN", "IRQ", "ISR", "JOR", "KWT",
                    "LBN", "LBY", "MAR", "OMN", "QAT", "SAU", "SYR", "TUN",
                    "ARE", "YEM", "PSE")
eastern_europe <- c("ALB", "BIH", "BGR", "HRV", "CZE", "HUN", "MKD", "POL",
                   "ROU", "SRB", "SVK", "SVN")
western_europe_na <- c("AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "DEU", "GRC",
                      "ISL", "IRL", "ITA", "LUX", "MLT", "NLD", "NOR", "PRT",
                      "ESP", "SWE", "CHE", "GBR", "USA")

df <- df %>%
  mutate(
    country_group = case_when(
      ISO3 %in% post_soviet ~ "Post-Soviet",
      ISO3 %in% latin_america ~ "Latin America",
      ISO3 %in% sub_saharan_africa ~ "Sub-Saharan Africa",
      ISO3 %in% east_asia_pacific ~ "East Asia & Pacific",
      ISO3 %in% middle_east_na ~ "Middle East & N. Africa",
      ISO3 %in% eastern_europe ~ "Eastern Europe",
      ISO3 %in% western_europe_na ~ "Western Europe & N. America",
      TRUE ~ "Other"
    ),
    region_year = paste0(country_group, "_", Year)
  )

# ============================================================================
# 3. TREATMENT DEFINITION
# ============================================================================
# Use fraser_Summary (or heritage_Overall as fallback)
policy_var <- if (sum(!is.na(df$fraser_Summary)) >= 1000) "fraser_Summary" else "heritage_Overall"
cat("Using policy variable:", policy_var, "\n")

# Identify treatment: top-decile jump (top 10% of year-over-year increases)
identify_treatment <- function(data, pvar, pct = 0.90) {
  d <- data %>% filter(!is.na(.data[[pvar]])) %>% arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(lag_p = dplyr::lag(.data[[pvar]]), jump = .data[[pvar]] - lag_p) %>%
    ungroup() %>% filter(!is.na(jump))
  if (nrow(d) == 0) return(tibble(ISO3 = character(), treated = logical(), event_year = numeric()))
  jump_threshold <- quantile(d$jump, probs = pct, na.rm = TRUE)
  d %>%
    filter(jump >= jump_threshold) %>%
    group_by(ISO3) %>% slice(1) %>% ungroup() %>%
    transmute(ISO3, treated = TRUE, event_year = Year)
}

treatment_info <- identify_treatment(df, policy_var)
treated_iso3 <- treatment_info$ISO3
policy_label <- if (policy_var == "fraser_Summary") "Fraser Economic Freedom" else "Heritage Economic Freedom"
treatment_def <- paste0("Treatment: top-decile jump in ", policy_label, " (top 10% year-over-year increase)")
n_treated <- nrow(treatment_info)
cat("Number of treated countries:", n_treated, "\n")
if (n_treated < 10) stop("Need at least 10 treated countries. Found ", n_treated, ".")
if (n_treated > 0) {
  print(treatment_info)
  tibble(policy_var = policy_var, n_treated = n_treated, treated_ISO3 = paste(treatment_info$ISO3, collapse = ", ")) %>%
    write_csv("event_study_treatment_count.csv")
  cat("Saved event_study_treatment_count.csv\n")
  treatment_info %>%
    left_join(distinct(df, ISO3, Country), by = "ISO3") %>%
    select(ISO3, Country, event_year) %>%
    arrange(Country) %>%
    write_csv("fraser_treated_country_year.csv")
  cat("Saved fraser_treated_country_year.csv\n")
}

# Merge treatment into main data
df <- df %>%
  left_join(treatment_info, by = "ISO3") %>%
  mutate(
    treated = replace_na(treated, FALSE),
    event_year = if_else(treated, event_year, NA_real_),
    post = if_else(treated, as.numeric(Year > event_year), 0)
  )

# ============================================================================
# 4. EVENT TIME
# ============================================================================
event_times <- c(-5:-2, 0:5)   # -1 omitted as reference, 5 years post-treatment
df <- df %>%
  mutate(event_time = if_else(treated, Year - event_year, NA_real_)) %>%
  mutate(event_time = if_else(event_time >= -5 & event_time <= 5, event_time, NA_real_))

for (k in event_times) {
  vname <- paste0("event_time_", k)
  df[[vname]] <- as.numeric(df$event_time == k & df$treated)
  df[[vname]] <- replace_na(df[[vname]], 0)
}

# ============================================================================
# 5. DEFINE OUTCOMES
# ============================================================================
outcome_cols <- c()
if ("log_GDP_per_capita" %in% names(df)) outcome_cols <- c(outcome_cols, "log_GDP_per_capita")
if ("GDP_per_capita_growth" %in% names(df)) outcome_cols <- c(outcome_cols, "GDP_per_capita_growth")
gini <- names(df)[grepl("Gini|GINI", names(df), ignore.case = TRUE)]
if (length(gini) > 0) outcome_cols <- c(outcome_cols, gini[1])
pov <- names(df)[grepl("Poverty|SI.POV", names(df), ignore.case = TRUE)]
if (length(pov) > 0) outcome_cols <- c(outcome_cols, pov[1])
inc <- names(df)[grepl("Income.share|SI.DST.FRST", names(df), ignore.case = TRUE)]
if (length(inc) > 0) outcome_cols <- c(outcome_cols, inc[1])
vuln <- names(df)[grepl("Vulnerable.employment|SL.EMP.VULN", names(df), ignore.case = TRUE)]
if (length(vuln) > 0) outcome_cols <- c(outcome_cols, vuln[1])
if (file.exists("component_impacts_significant.csv")) {
  sig_outcomes <- read_csv("component_impacts_significant.csv", show_col_types = FALSE) %>%
    pull(outcome) %>% unique()
  outcome_cols <- intersect(outcome_cols, sig_outcomes)
  if (length(outcome_cols) > 0) cat("Restricted to significant outcomes from component analysis:\n")
}
cat("Outcomes:", paste(outcome_cols, collapse = ", "), "\n")

# Winsorize DVs at 1% and 99% to reduce variance from extreme outliers
winsorize <- function(x, p_low = 0.01, p_high = 0.99) {
  q <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
df <- df %>% mutate(across(all_of(outcome_cols), winsorize, .names = "{.col}"))

# No lagged DV (avoids Nickell bias in FE)
control_var <- "log_GDP_per_capita"  # structural control (omit when outcome is log GDP)
if (!control_var %in% names(df)) control_var <- NULL

# ============================================================================
# 6. RUN EVENT STUDY (control: log GDP pc; no lagged DV to avoid Nickell bias)
# ============================================================================
run_event_study <- function(data, outcome, event_times) {
  use_control <- !is.null(control_var) && control_var %in% names(data) && control_var != outcome
  data_sub <- data %>% filter(!is.na(.data[[outcome]]))
  if (use_control) data_sub <- data_sub %>% filter(!is.na(.data[[control_var]]))
  if (nrow(data_sub) < 100 || sum(data_sub$treated) < 10) return(NULL)

  event_terms <- paste0("`event_time_", event_times, "`", collapse = " + ")
  out_escaped <- if (grepl("[^a-zA-Z0-9._]", outcome)) paste0("`", outcome, "`") else outcome
  rhs <- event_terms
  if (use_control) rhs <- paste(rhs, paste0("`", control_var, "`"), sep = " + ")
  # Baseline: ISO3 + Year. For robustness, replace Year with region_year.
  f <- as.formula(paste0(out_escaped, " ~ ", rhs, " | ISO3 + Year"))
  tryCatch(feols(f, data = data_sub, cluster = ~ ISO3), error = function(e) NULL)
}

results <- list()
for (out in outcome_cols) {
  m <- run_event_study(df, out, event_times)
  if (!is.null(m)) results[[out]] <- m
}

# ============================================================================
# 6b. PARALLEL TRENDS TEST (joint F-test: pre-treatment coefs = 0)
# ============================================================================
test_pre_trends <- function(model, event_times) {
  pre_periods <- event_times[event_times < 0]
  if (length(pre_periods) == 0) return(NULL)
  pre_vars <- paste0("event_time_", pre_periods)
  pre_vars_bt <- paste0("`event_time_", pre_periods, "`")
  model_coefs <- names(coef(model))
  available_vars <- pre_vars_bt[pre_vars_bt %in% model_coefs]
  if (length(available_vars) == 0) available_vars <- pre_vars[pre_vars %in% model_coefs]
  if (length(available_vars) == 0) return(NULL)
  tryCatch({
    wt <- wald(model, keep = available_vars,
               print = FALSE)
    list(test_stat = as.numeric(wt["stat"]), p_value = as.numeric(wt["p"]), pre_periods = list(pre_periods))
  }, error = function(e) {
    coefs <- coef(model)[available_vars]
    vcv <- vcov(model)[available_vars, available_vars, drop = FALSE]
    if (length(coefs) == 1) {
      test_stat <- (coefs / sqrt(vcv[1,1]))^2
      p_val <- 1 - pchisq(test_stat, 1)
    } else {
      test_stat <- as.numeric(t(coefs) %*% solve(vcv) %*% coefs)
      p_val <- 1 - pchisq(test_stat, length(coefs))
    }
    list(test_stat = test_stat, p_value = p_val, pre_periods = list(pre_periods))
  })
}

pre_trend_tests <- list()
for (out in names(results)) {
  pt <- test_pre_trends(results[[out]], event_times)
  if (!is.null(pt)) pre_trend_tests[[out]] <- pt
}

cat("\n=== PARALLEL TRENDS TEST (H0: pre-treatment coefs = 0) ===\n")
cat("p > 0.10 supports parallel trends assumption\n\n")
for (out in names(pre_trend_tests)) {
  pt <- pre_trend_tests[[out]]
  cat(sprintf("  %s: F = %.2f, p = %.3f %s\n", out, pt$test_stat, pt$p_value,
              if (pt$p_value > 0.10) "(OK)" else "(pre-trends detected)"))
}

# ============================================================================
# 6c. SIMPLE DiD (Treated * Post) with SE clustered by country
# ============================================================================
run_simple_did <- function(data, outcome) {
  use_control <- control_var %in% names(data) && control_var != outcome
  d <- data %>% filter(!is.na(.data[[outcome]]))
  if (use_control) d <- d %>% filter(!is.na(.data[[control_var]]))
  if (nrow(d) < 100 || sum(d$treated) < 10) return(NULL)
  out_esc <- if (grepl("[^a-zA-Z0-9._]", outcome)) paste0("`", outcome, "`") else outcome
  rhs <- "treated * post"
  if (use_control) rhs <- paste(rhs, paste0("`", control_var, "`"), sep = " + ")
  f <- as.formula(paste0(out_esc, " ~ ", rhs, " | ISO3 + Year"))
  tryCatch(feols(f, data = d, cluster = ~ ISO3), error = function(e) NULL)
}

did_models <- list()
for (out in outcome_cols) {
  m <- run_simple_did(df, out)
  if (!is.null(m)) did_models[[out]] <- m
}

cat("\n=== SIMPLE DiD (Treated * Post, SE clustered by ISO3) ===\n")
for (out in names(did_models)) {
  b <- coef(did_models[[out]])["treated:post"]
  se_b <- se(did_models[[out]])["treated:post"]
  if (!is.na(b)) cat(sprintf("  %s: DiD coef = %.4f (SE = %.4f)\n", out, b, se_b))
}

# ============================================================================
# 7. ATE CALCULATION
# ============================================================================
calculate_ate <- function(model, event_times) {
  post <- event_times[event_times >= 0]
  vnames <- paste0("event_time_", post)
  coefs <- coef(model)
  ses <- setNames(se(model), names(coefs))
  av <- vnames[vnames %in% names(coefs)]
  if (length(av) == 0) return(NULL)
  ate <- mean(coefs[av], na.rm = TRUE)
  ate_se <- sqrt(mean(ses[av]^2, na.rm = TRUE) / length(av))
  list(ate = ate, se = ate_se, ci_lower = ate - 1.96 * ate_se, ci_upper = ate + 1.96 * ate_se)
}

cat("\n=== AVERAGE TREATMENT EFFECTS (post-period) ===\n")
for (nm in names(results)) {
  ate <- calculate_ate(results[[nm]], event_times)
  if (!is.null(ate)) {
    cat(sprintf("  %s: ATE = %.4f, CI [%.4f, %.4f]\n", nm, ate$ate, ate$ci_lower, ate$ci_upper))
  }
}

# ============================================================================
# 8. EVENT-STUDY PLOT
# ============================================================================
plot_event_study <- function(model, outcome_name, event_times) {
  if (is.null(model)) return(NULL)
  coefs <- coef(model)
  ses <- setNames(se(model), names(coefs))
  pre_t <- event_times[event_times < 0]
  post_t <- event_times[event_times >= 0]
  plot_times <- c(pre_t, -1, post_t)
  est <- numeric(length(plot_times))
  se_est <- numeric(length(plot_times))
  est[] <- NA_real_
  se_est[] <- NA_real_
  all_names <- names(coefs)
  for (i in seq_along(plot_times)) {
    k <- plot_times[i]
    if (k == -1) { est[i] <- 0; se_est[i] <- 0; next }  # reference period
    vn <- paste0("event_time_", k)
    vn_bt <- paste0("`event_time_", k, "`")  # fixest uses backticks for neg event times
    name_use <- if (vn_bt %in% all_names) vn_bt else if (vn %in% all_names) vn else NULL
    if (!is.null(name_use)) {
      est[i] <- coefs[name_use]
      se_est[i] <- ses[name_use]
    }
  }
  tibble(
    event_time = plot_times,
    estimate = est,
    se = se_est,
    ci_lower = est - 1.96 * se_est,
    ci_upper = est + 1.96 * se_est
  ) %>%
    filter(!is.na(estimate) | event_time == -1)
}

if (length(results) > 0) {
  dir.create("event_study_plots", showWarnings = FALSE)
  for (out in names(results)) {
    pd <- plot_event_study(results[[out]], out, event_times)
    if (nrow(pd) > 0) {
      pd <- pd %>% mutate(period = if_else(event_time < 0, "Pre", "Post"))
      pt_label <- if (out %in% names(pre_trend_tests))
        sprintf("Pre-trend p = %.3f", pre_trend_tests[[out]]$p_value) else ""
      y_min <- min(pd$ci_lower, na.rm = TRUE)
      y_max <- max(pd$ci_upper, na.rm = TRUE)
      y_span <- y_max - y_min
      if (y_span < 0.1) y_span <- 0.1

      p <- ggplot() +
        annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "gray90", alpha = 0.5) +
        annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.2) +
        annotate("text", x = -2.5, y = y_max - 0.02 * y_span, label = "Pre", size = 4, fontface = "bold") +
        annotate("text", x = 5, y = y_max - 0.02 * y_span, label = "Post", size = 4, fontface = "bold") +
        geom_vline(xintercept = 0, linetype = "dashed", size = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_point(data = pd, aes(x = event_time, y = estimate, color = period), size = 3) +
        geom_errorbar(data = pd, aes(x = event_time, y = estimate, ymin = ci_lower, ymax = ci_upper, color = period), width = 0.2, size = 0.8) +
        scale_color_manual(values = c("Pre" = "#2166ac", "Post" = "#b2182b"),
                          breaks = c("Pre", "Post"),
                          na.translate = FALSE) +
        labs(x = "Years relative to treatment", y = "Coefficient / Normalized outcome",
             title = paste("Event study:", out),
             subtitle = paste0(treatment_def, ". Blue = pre, red = post. ", pt_label),
             color = "Period") +
        expand_limits(y = 0) +
        theme_minimal() +
        theme(legend.position = "bottom")
      fname <- paste0("event_study_plots/event_study_", gsub("[^a-zA-Z0-9]", "_", out), ".png")
      ggsave(fname, p, width = 7, height = 5, dpi = 150)
      cat("Saved", fname, "\n")
    }
  }
}

# ============================================================================
# 9. SAVE RESULTS
# ============================================================================
ate_tbl <- map_dfr(names(results), function(out) {
  ate <- calculate_ate(results[[out]], event_times)
  if (is.null(ate)) return(tibble())
  tibble(outcome = out, ate = ate$ate, se = ate$se, ci_lower = ate$ci_lower, ci_upper = ate$ci_upper)
})
if (nrow(ate_tbl) > 0) {
  write_csv(ate_tbl, "event_study_ate_results.csv")
  cat("\nSaved event_study_ate_results.csv\n")
}

# Pre-trend test results
if (length(pre_trend_tests) > 0) {
  pre_trend_df <- tibble(
    outcome = names(pre_trend_tests),
    f_statistic = map_dbl(pre_trend_tests, "test_stat"),
    p_value = map_dbl(pre_trend_tests, "p_value")
  )
  write_csv(pre_trend_df, "pre_trend_tests.csv")
  cat("Saved pre_trend_tests.csv\n")
}

cat("\nDone.\n")
