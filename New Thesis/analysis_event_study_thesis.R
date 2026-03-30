# Thesis event studies: one shock (Year 0), two outcome blocks
# 1) Shock: top-decile YoY jump in a chosen Fraser component (default: fraser_FreeTrade)
# 2) Analysis A — Macro: GDP growth + vulnerable employment
# 3) Analysis B — Micro: extreme poverty + income share bottom 20%
# Event time: reference t = -1 (omitted). Pre-trend test: coefs at t = -3, -2.
# Post: t = +1, ..., +5 (t = 0 included in model as separate dummy)
#
# Requires: tidyverse, fixest. Run from project root after data/merged_full_data.csv exists.

if (!requireNamespace("fixest", quietly = TRUE)) {
  install.packages("fixest", repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(fixest)

# --- Settings ---
policy_shock_var <- "fraser_FreeTrade" # Fraser component defining "year zero"
jump_quantile <- 0.90                # top 10% of all historical policy changes
min_year <- 1995L
max_year <- 2022L
min_obs_country <- 15L
min_treated_countries <- 10L

# Event dummies: -1 omitted (reference). Include -3, -2, 0, +1..+5
event_times <- c(-3L, -2L, 0L, 1L, 2L, 3L, 4L, 5L)

# --- Load data ---
if (!file.exists("data/merged_full_data.csv")) {
  source("cleaning/merge_full_data.R")
}
df <- read_csv("data/merged_full_data.csv", show_col_types = FALSE)
df <- df %>%
  filter(Year >= min_year, Year <= max_year, !is.na(ISO3), !is.na(Year))

panel_iso3 <- df %>% count(ISO3) %>% filter(n >= min_obs_country) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)
cat("Panel:", length(unique(df$ISO3)), "countries,", length(unique(df$Year)), "years\n")

# --- Resolve outcome column names (make.names-safe) ---
resolve_col <- function(d, exact, patterns) {
  if (exact %in% names(d)) return(exact)
  for (p in patterns) {
    hit <- names(d)[grepl(p, names(d), ignore.case = TRUE)][1]
    if (!is.na(hit)) return(hit)
  }
  NA_character_
}

col_gdp_growth <- resolve_col(df, "NY.GDP.MKTP.KD.ZG", c("^NY\\.GDP\\.MKTP\\.KD\\.ZG$", "MKTP", "GDP.*MKTP"))
col_gdp_pc_growth <- resolve_col(df, "GDP_per_capita_growth", c("^GDP_per_capita_growth$"))
col_vuln <- resolve_col(df, "SL.EMP.VULN.ZS", c("SL\\.EMP\\.VULN", "Vulnerable\\.employment"))
col_pov_extreme <- resolve_col(df, "SI.POV.DDAY", c("^SI\\.POV\\.DDAY$", "DDAY"))
col_share_bottom <- resolve_col(df, "SI.DST.FRST.20", c("SI\\.DST\\.FRST", "Income\\.share.*lowest"))

# Analysis A: national GDP growth (WDI) + vulnerable employment; fallback to GDP_per_capita_growth if MKTP missing
outcomes_macro <- c(col_gdp_growth, col_vuln)
if (is.na(col_gdp_growth) && !is.na(col_gdp_pc_growth)) {
  outcomes_macro <- c(col_gdp_pc_growth, col_vuln)
}
outcomes_macro <- outcomes_macro[!is.na(outcomes_macro)]

outcomes_micro <- c(col_pov_extreme, col_share_bottom)
outcomes_micro <- outcomes_micro[!is.na(outcomes_micro)]

if (length(outcomes_macro) == 0) stop("No macro outcome columns found (need NY.GDP.MKTP.KD.ZG and/or GDP_per_capita_growth, SL.EMP.VULN.ZS). Re-run merge_full_data.R after updating wdi_clean.")
if (length(outcomes_micro) < 1) {
  warning("Micro outcomes partly missing (SI.POV.DDAY and/or SI.DST.FRST.20). Re-run merge after wdi_clean includes those series.")
}

cat("Macro outcomes:", paste(outcomes_macro, collapse = ", "), "\n")
cat("Micro outcomes:", paste(outcomes_micro, collapse = ", "), "\n")

# --- Treatment (Year 0), same policy variable for all analyses ---
if (!policy_shock_var %in% names(df)) {
  stop("policy_shock_var not in data: ", policy_shock_var)
}

identify_treatment <- function(data, pvar, pct = jump_quantile) {
  d <- data %>%
    filter(!is.na(.data[[pvar]])) %>%
    arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(lag_p = dplyr::lag(.data[[pvar]]), jump = .data[[pvar]] - lag_p) %>%
    ungroup() %>%
    filter(!is.na(jump))
  if (nrow(d) == 0) {
    return(tibble(ISO3 = character(), treated = logical(), event_year = numeric()))
  }
  jump_threshold <- quantile(d$jump, probs = pct, na.rm = TRUE)
  d %>%
    filter(jump >= jump_threshold) %>%
    group_by(ISO3) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(ISO3, treated = TRUE, event_year = Year)
}

treatment_info <- identify_treatment(df, policy_shock_var)
n_treated <- nrow(treatment_info)
cat("Treated countries (", policy_shock_var, ", top ", (1 - jump_quantile) * 100, "% jumps): ", n_treated, "\n", sep = "")
if (n_treated < min_treated_countries) {
  stop("Need at least ", min_treated_countries, " treated countries. Found ", n_treated, ".")
}
print(treatment_info)

df <- df %>%
  left_join(treatment_info, by = "ISO3") %>%
  mutate(
    treated = replace_na(treated, FALSE),
    event_year = if_else(treated, event_year, NA_real_),
    post = if_else(treated, as.numeric(Year > event_year), 0)
  )

# Event time (years relative to shock); only for treated
df <- df %>%
  mutate(event_time = if_else(treated, as.integer(Year - event_year), NA_integer_))

# Dummies: omit t = -1 (reference)
for (k in event_times) {
  vn <- paste0("event_time_", k)
  df[[vn]] <- as.numeric(df$event_time == k & df$treated)
  df[[vn]] <- replace_na(df[[vn]], 0)
}

control_var <- if ("log_GDP_per_capita" %in% names(df)) "log_GDP_per_capita" else NULL

winsorize <- function(x, p_low = 0.01, p_high = 0.99) {
  q <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

all_outcomes <- unique(c(outcomes_macro, outcomes_micro))
all_outcomes <- all_outcomes[all_outcomes %in% names(df)]
df <- df %>% mutate(across(all_of(all_outcomes), winsorize))

# --- Model: ISO3 + Year FE, SE clustered by country ---
build_rhs <- function(outcome, event_times_vec, control_var) {
  terms <- paste0("`event_time_", event_times_vec, "`")
  rhs <- paste(terms, collapse = " + ")
  if (!is.null(control_var) && control_var %in% names(df) && control_var != outcome) {
    rhs <- paste(rhs, "+", paste0("`", control_var, "`"))
  }
  rhs
}

run_es <- function(outcome) {
  if (!outcome %in% names(df)) return(NULL)
  d <- df %>% filter(!is.na(.data[[outcome]]))
  if (!is.null(control_var) && control_var %in% names(d) && control_var != outcome) {
    d <- d %>% filter(!is.na(.data[[control_var]]))
  }
  if (nrow(d) < 100L || sum(d$treated) < min_treated_countries) return(NULL)
  out_esc <- if (grepl("[^a-zA-Z0-9._]", outcome)) paste0("`", outcome, "`") else outcome
  rhs <- build_rhs(outcome, event_times, control_var)
  f <- as.formula(paste0(out_esc, " ~ ", rhs, " | ISO3 + Year"))
  tryCatch(feols(f, data = d, cluster = ~ ISO3), error = function(e) NULL)
}

# Pre-trend: H0: coef(event_time_-3) = coef(event_time_-2) = 0
joint_pre_trend <- function(model) {
  if (is.null(model)) return(NULL)
  cn <- names(coef(model))
  available <- cn[grepl("event_time_-3|event_time_-2", cn)]
  if (length(available) == 0) return(NULL)
  tryCatch(wald(model, keep = available), error = function(e) NULL)
}

collect_results <- function(outcome, model, block) {
  if (is.null(model)) return(tibble())
  pt <- joint_pre_trend(model)
  pre_p <- NA_real_
  if (!is.null(pt) && is.matrix(pt)) {
    if ("p" %in% colnames(pt)) pre_p <- as.numeric(pt[, "p"])
    else if ("Pr(>Chisq)" %in% colnames(pt)) pre_p <- as.numeric(pt[, "Pr(>Chisq)"])
    else if (ncol(pt) >= 2) pre_p <- as.numeric(pt[1, ncol(pt)])
  }
  tibble(
    block = block,
    outcome = outcome,
    pre_trend_wald_p = pre_p,
    n_obs = as.integer(nobs(model))
  )
}

run_block <- function(block_name, outcomes) {
  res <- list()
  models <- list()
  for (y in outcomes) {
    m <- run_es(y)
    models[[y]] <- m
    res[[length(res) + 1]] <- collect_results(y, m, block_name)
    if (!is.null(m)) {
      cat("\n=== ", block_name, ": ", y, " ===\n", sep = "")
      print(summary(m))
      wpt <- joint_pre_trend(m)
      if (!is.null(wpt)) {
        cat("Pre-trend joint Wald (t=-3, -2; ref t=-1):\n")
        print(wpt)
      }
    }
  }
  list(models = models, tbl = bind_rows(res))
}

r_macro <- run_block("Analysis_A_Macro", outcomes_macro)
r_micro <- if (length(outcomes_micro) > 0) {
  run_block("Analysis_B_Micro", outcomes_micro)
} else {
  list(models = list(), tbl = tibble())
}

out_tbl <- bind_rows(r_macro$tbl, r_micro$tbl)
if (nrow(out_tbl) > 0) {
  write_csv(out_tbl, "event_study_thesis_summary.csv")
  cat("\nSaved event_study_thesis_summary.csv\n")
}

# Save treatment list for documentation
treatment_info %>%
  left_join(distinct(df, ISO3, Country), by = "ISO3") %>%
  select(ISO3, Country, event_year) %>%
  arrange(Country) %>%
  write_csv("thesis_shock_year_zero.csv")
cat("Saved thesis_shock_year_zero.csv\n")

cat("\nDone. Reference period for event time is t = -1 (omitted). Pre-trends: t = -3, -2. Post effects: t = 0, +1, ..., +5.\n")
