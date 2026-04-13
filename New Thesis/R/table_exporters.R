# Pooled OLS / TWFE extracts and CS-DiD ATT + BH–FDR consolidations for thesis tables.
# Expects `results/intermediate/*.csv` from the main pipeline. Writes to `results/final_tables/`.
# Normally sourced from R/analysis_methodology.R; can be run standalone from New Thesis after methodology.

for (.p in c(
  file.path(getwd(), "R", "proj_paths.R"),
  file.path(dirname(getwd()), "R", "proj_paths.R"),
  file.path(dirname(dirname(getwd())), "R", "proj_paths.R")
)) {
  if (file.exists(.p)) {
    source(.p, local = FALSE)
    break
  }
}
if (!exists("find_proj_root", mode = "function")) {
  stop("Cannot find R/proj_paths.R")
}

library(tidyverse)
library(fixest)

proj_root <- find_proj_root()
path_int <- file.path(proj_root, "results", "intermediate")
path_ft <- file.path(proj_root, "results", "final_tables")
if (!dir.exists(path_ft)) {
  dir.create(path_ft, recursive = TRUE)
}

# =============================================================================
# Pooled OLS: Y ~ policy score for all thesis outcomes (Level 1 style)
# =============================================================================

min_year <- 1995L
max_year <- 2022L
min_obs_country <- 10L
policy_components <- c(
  "fraser_SizeGov", "fraser_LegalSystem", "fraser_SoundMoney", "fraser_FreeTrade", "fraser_Regulation",
  "heritage_Property.Rights", "heritage_Government.Integrity", "heritage_Judicial.Effectiveness",
  "heritage_Tax.Burden", "heritage_Government.Spending", "heritage_Fiscal.Health",
  "heritage_Business.Freedom", "heritage_Labor.Freedom", "heritage_Monetary.Freedom",
  "heritage_Trade.Freedom", "heritage_Investment.Freedom", "heritage_Financial.Freedom"
)

path_merged <- file.path(proj_root, "data", "merged_full_data.csv")
df <- read_csv(path_merged, show_col_types = FALSE) %>%
  filter(Year >= min_year, Year <= max_year, !is.na(ISO3), !is.na(Year))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= min_obs_country) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)

resolve_outcome <- function(d) {
  tibble(
    label = c(
      "GDP_per_capita_growth",
      "Gini_SWIID",
      "Income_share_bottom50_WID",
      "log_under5_mort"
    ),
    col = c(
      names(d)[names(d) == "GDP_per_capita_growth"][1],
      names(d)[names(d) == "Gini_SWIID"][1],
      names(d)[names(d) == "Income_share_bottom50_WID"][1],
      names(d)[names(d) == "log_under5_mort"][1]
    )
  ) %>% filter(!is.na(col))
}

outcomes_tbl <- resolve_outcome(df)
policy_components <- intersect(policy_components, names(df))

p_stars <- function(p) {
  if (!is.finite(p)) {
    return("")
  }
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

coef_latex <- function(coef, stars) {
  x <- formatC(round(coef, 3), format = "f", digits = 3)
  if (nzchar(stars)) {
    sprintf("$%s^{%s}$", x, stars)
  } else {
    sprintf("$%s$", x)
  }
}

rows_ols <- list()
for (oc in outcomes_tbl$col) {
  ol <- outcomes_tbl$label[outcomes_tbl$col == oc][1]
  for (comp in policy_components) {
    d <- df %>% filter(!is.na(.data[[oc]]), !is.na(.data[[comp]]))
    m <- lm(as.formula(paste0("`", oc, "` ~ `", comp, "`")), data = d)
    sm <- summary(m)
    ct <- sm$coefficients
    est <- as.numeric(ct[2, 1])
    pv <- as.numeric(ct[2, 4])
    st <- p_stars(pv)
    rows_ols[[paste(ol, comp, sep = "||")]] <- tibble(
      outcome_label = ol,
      component = comp,
      coef = est,
      se = as.numeric(ct[2, 2]),
      p_value = pv,
      stars = st,
      coef_latex = coef_latex(est, st),
      n = nobs(m),
      r2 = sm$r.squared
    )
  }
}

out_ols <- bind_rows(rows_ols)
write_csv(out_ols, file.path(path_ft, "ols_pooled_all_outcomes.csv"))
cat("Wrote results/final_tables/ols_pooled_all_outcomes.csv\n")

# =============================================================================
# TWFE: Y ~ Score | ISO3 + Year, cluster SE ~ ISO3 (Level 2 style)
# =============================================================================

run_twfe <- function(data, outcome, comp_col) {
  d <- data %>% filter(!is.na(.data[[outcome]]), !is.na(.data[[comp_col]]))
  if (nrow(d) < 50L) {
    return(NULL)
  }
  out_esc <- paste0("`", outcome, "`")
  comp_esc <- paste0("`", comp_col, "`")
  f <- as.formula(paste0(out_esc, " ~ ", comp_esc, " | ISO3 + Year"))
  tryCatch(
    feols(f, data = d, cluster = ~ ISO3),
    error = function(e) NULL
  )
}

within_r2 <- function(m) {
  if (is.null(m)) {
    return(NA_real_)
  }
  tryCatch(
    as.numeric(fitstat(m, "wr2", simplify = TRUE)),
    error = function(e) {
      tryCatch(
        as.numeric(r2(m, type = "wr2")),
        error = function(e2) NA_real_
      )
    }
  )
}

sd_mult_main <- 1.5
compute_first_treat <- function(data, comp_col, sd_mult) {
  d <- data %>%
    arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(
      score = .data[[comp_col]],
      delta = score - dplyr::lag(score)
    ) %>%
    ungroup()
  sig <- stats::sd(d$delta, na.rm = TRUE) * sd_mult
  if (!is.finite(sig) || sig <= 0) {
    return(tibble(ISO3 = unique(data$ISO3), first_treat = 0L))
  }
  shocks <- d %>%
    filter(!is.na(delta), delta > sig) %>%
    group_by(ISO3) %>%
    slice_min(Year, n = 1L, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(ISO3, first_treat = as.integer(Year))
  all_iso <- tibble(ISO3 = unique(data$ISO3))
  all_iso %>%
    left_join(shocks, by = "ISO3") %>%
    mutate(first_treat = ifelse(is.na(first_treat), 0L, first_treat))
}

rows_twfe <- list()
for (comp in policy_components) {
  ft <- compute_first_treat(df, comp, sd_mult_main)
  d <- df %>% left_join(ft, by = "ISO3")
  for (oc in outcomes_tbl$col) {
    ol <- outcomes_tbl$label[outcomes_tbl$col == oc][1]
    m <- run_twfe(d, oc, comp)
    if (is.null(m)) {
      rows_twfe[[paste(ol, comp, sep = "||")]] <- tibble(
        outcome_label = ol,
        component = comp,
        coef = NA_real_,
        se = NA_real_,
        p_value = NA_real_,
        stars = NA_character_,
        coef_latex = NA_character_,
        n = NA_integer_,
        r2_within = NA_real_
      )
      next
    }
    est <- unname(coef(m)[1])
    pv <- unname(fixest::pvalue(m)[1])
    sev <- unname(fixest::se(m)[1])
    st <- p_stars(pv)
    rows_twfe[[paste(ol, comp, sep = "||")]] <- tibble(
      outcome_label = ol,
      component = comp,
      coef = est,
      se = sev,
      p_value = pv,
      stars = st,
      coef_latex = coef_latex(est, st),
      n = nobs(m),
      r2_within = within_r2(m)
    )
  }
}

out_twfe <- bind_rows(rows_twfe)
write_csv(out_twfe, file.path(path_ft, "twfe_pooled_all_outcomes.csv"))
cat("Wrote results/final_tables/twfe_pooled_all_outcomes.csv\n")

# =============================================================================
# CS-DiD: simple ATT dumps + BH–FDR (reads intermediate Level 3/4 CSVs)
# =============================================================================

main <- read_csv(file.path(path_int, "results_level3_cs_att.csv"), show_col_types = FALSE) %>%
  mutate(
    sd_mult_shock = 1.5,
    p_value_raw = p_value
  ) %>%
  select(sd_mult_shock, component, outcome, outcome_col, att, se, p_value_raw)

rob <- read_csv(file.path(path_int, "results_level3_cs_att_robust10sd.csv"), show_col_types = FALSE) %>%
  mutate(
    sd_mult_shock = 1.0,
    p_value_raw = p_value
  ) %>%
  select(sd_mult_shock, component, outcome, outcome_col, att, se, p_value_raw)

combined <- bind_rows(main, rob) %>%
  arrange(sd_mult_shock, outcome, component)

write_csv(main, file.path(path_ft, "cs_did_simple_att_main_sd15.csv"))
write_csv(rob, file.path(path_ft, "cs_did_simple_att_robust_sd10.csv"))
write_csv(combined, file.path(path_ft, "cs_did_simple_att_all_specs.csv"))

main_fdr <- read_csv(file.path(path_int, "results_level4_fdr.csv"), show_col_types = FALSE) %>%
  mutate(
    sd_mult_shock = 1.5,
    p_value_raw = p_value
  ) %>%
  select(sd_mult_shock, component, outcome, outcome_col, att, se, p_value_raw, fdr_bh, reject_bh_05)

rob_fdr <- read_csv(file.path(path_int, "results_level3_cs_att_robust10sd_fdr.csv"), show_col_types = FALSE) %>%
  mutate(
    sd_mult_shock = 1.0,
    p_value_raw = p_value
  ) %>%
  select(sd_mult_shock, component, outcome, outcome_col, att, se, p_value_raw, fdr_bh, reject_bh_05)

combined_fdr <- bind_rows(main_fdr, rob_fdr) %>%
  arrange(sd_mult_shock, outcome, component)

write_csv(main_fdr, file.path(path_ft, "did_fdr_main_15sd.csv"))
write_csv(rob_fdr, file.path(path_ft, "did_fdr_robust_10sd.csv"))
write_csv(combined_fdr, file.path(path_ft, "cs_did_bh_fdr_all_specs.csv"))

cat(
  "Wrote results/final_tables/cs_did_simple_att_*.csv; did_fdr_main_15sd.csv; did_fdr_robust_10sd.csv; cs_did_bh_fdr_all_specs.csv\n",
  sep = ""
)
