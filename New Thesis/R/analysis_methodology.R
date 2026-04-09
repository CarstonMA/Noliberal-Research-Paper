# Four-stage methodology (thesis Methods chapter):
#   Level 1 — OLS: Y ~ Score
#   Level 2 — TWFE: Y ~ Score | ISO3 + Year, cluster ISO3
#   Level 3 — Callaway & Sant'Anna: att_gt + aggte (simple ATT), cluster bootstrap
#   Level 4 — Benjamini–Hochberg FDR on Level-3 p-values
#
# Shock: first year in which YoY change in component > sd_mult * SD(Δ) over the full panel.
#        first_treat = 0 if never exceeds threshold (never-treated).
# Components: 17 sub-indices (Fraser without Summary; Heritage without Overall).
# Outcomes (4): GDP per capita growth; SWIID Gini; WID bottom-50% income share; WDI under-5 mortality.
#
# Requires: tidyverse, fixest, did.
# Run from the thesis folder (New Thesis):  Rscript R/analysis_methodology.R
#
# Also writes thesis-extraction outputs:
#   - Event study (dynamic aggte + ggdid), pre-trend coefficients
#   - Distributional outcome failure counts + survivor table
#   - Cohort size vs robustness threshold + Table 3 (GDP @ 1.0 SD)

req <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
req("tidyverse")
req("fixest")
req("did")

library(tidyverse)
library(fixest)
library(did)
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
proj_root <- find_proj_root()
path_results <- file.path(proj_root, "results")
if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}
path_figs <- file.path(path_results, "figures")
if (!dir.exists(path_figs)) {
  dir.create(path_figs, recursive = TRUE)
}

# --- Settings ---
min_year <- 1995L
max_year <- 2022L
min_obs_country <- 10L
sd_mult_main <- 1.5
sd_mult_robust <- 1.0
bootstrap_iters <- 199L

# 17 sub-components (exclude aggregate headline scores per thesis)
policy_components <- c(
  "fraser_SizeGov", "fraser_LegalSystem", "fraser_SoundMoney", "fraser_FreeTrade", "fraser_Regulation",
  "heritage_Property.Rights", "heritage_Government.Integrity", "heritage_Judicial.Effectiveness",
  "heritage_Tax.Burden", "heritage_Government.Spending", "heritage_Fiscal.Health",
  "heritage_Business.Freedom", "heritage_Labor.Freedom", "heritage_Monetary.Freedom",
  "heritage_Trade.Freedom", "heritage_Investment.Freedom", "heritage_Financial.Freedom"
)

path_merged <- file.path(proj_root, "data", "merged_full_data.csv")

if (!file.exists(path_merged)) {
  source(file.path(proj_root, "R", "cleaning", "merge_full_data.R"))
}

df <- read_csv(path_merged, show_col_types = FALSE) %>%
  filter(Year >= min_year, Year <= max_year, !is.na(ISO3), !is.na(Year))

panel_iso3 <- df %>% count(ISO3) %>% filter(n >= min_obs_country) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)

# Four thesis outcomes (resolve make.names)
resolve_outcome <- function(d) {
  tibble(
    label = c(
      "GDP_per_capita_growth",
      "Gini_SWIID",
      "Income_share_bottom50_WID",
      "Under5_mortality_rate"
    ),
    col = c(
      names(d)[names(d) == "GDP_per_capita_growth"][1],
      names(d)[names(d) == "Gini_SWIID"][1],
      names(d)[names(d) == "Income_share_bottom50_WID"][1],
      names(d)[names(d) == "Under5_mortality_rate"][1]
    )
  ) %>% filter(!is.na(col))
}

outcomes_tbl <- resolve_outcome(df)
outcome_cols <- outcomes_tbl$col
names(outcome_cols) <- outcomes_tbl$label
cat("Outcomes:", paste(names(outcome_cols), "=", outcome_cols, collapse = "; "), "\n")
expected <- c("GDP_per_capita_growth", "Gini_SWIID", "Income_share_bottom50_WID", "Under5_mortality_rate")
missing_labs <- setdiff(expected, names(outcome_cols))
if (length(missing_labs) > 0L) {
  warning(
    "Outcome column(s) missing from merged_full_data.csv — re-run R/cleaning/clean_wdi.R, ",
    "then merge_full_data.R (and add data/raw/swiid_summary.csv + WID via import_wid.R). Missing: ",
    paste(missing_labs, collapse = ", "),
    call. = FALSE
  )
}

policy_components <- intersect(policy_components, names(df))
if (length(policy_components) != 17L) {
  warning("Expected 17 policy columns; found ", length(policy_components), ": check merged_full_data.")
}

# --- Shock: first_treat (0 = never treated) ---
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

# --- Level 2: TWFE ---
run_twfe <- function(data, outcome, comp_col) {
  d <- data %>% filter(!is.na(.data[[outcome]]), !is.na(.data[[comp_col]]))
  if (nrow(d) < 50L) return(NULL)
  out_esc <- paste0("`", outcome, "`")
  comp_esc <- paste0("`", comp_col, "`")
  f <- as.formula(paste0(out_esc, " ~ ", comp_esc, " | ISO3 + Year"))
  tryCatch(
    feols(f, data = d, cluster = ~ ISO3),
    error = function(e) NULL
  )
}

# did::att_gt requires numeric idname (character ISO3 errors in pre_process_did).
att_gt_panel <- function(d, outcome, control_group = "nevertreated") {
  d <- as.data.frame(d)
  d$unit_id <- as.integer(factor(d$ISO3, levels = sort(unique(as.character(d$ISO3)))))
  att_gt(
    yname = outcome,
    tname = "Year",
    idname = "unit_id",
    gname = "first_treat",
    xformla = ~1,
    data = d,
    control_group = control_group,
    est_method = "dr",
    allow_unbalanced_panel = TRUE,
    clustervars = "unit_id",
    bstrap = TRUE,
    biters = bootstrap_iters,
    print_details = FALSE
  )
}

# --- Level 3: CS (did) ---
run_cs <- function(data, outcome) {
  d <- data %>%
    filter(!is.na(.data[[outcome]])) %>%
    mutate(Year = as.integer(Year))
  if (nrow(d) < 100L) return(NULL)
  if (length(unique(d$first_treat[d$first_treat > 0])) < 1L) return(NULL)
  tryCatch(
    {
      att <- att_gt_panel(d, outcome, control_group = "nevertreated")
      agg <- aggte(att, type = "simple", na.rm = TRUE)
      list(att = att, agg = agg)
    },
    error = function(e) {
      tryCatch(
        {
          att <- att_gt_panel(d, outcome, control_group = "notyettreated")
          agg <- aggte(att, type = "simple", na.rm = TRUE)
          list(att = att, agg = agg)
        },
        error = function(e2) NULL
      )
    }
  )
}

# att_gt + dynamic aggte (event study; GDP-only loop in extractions)
run_att_gt_dynamic <- function(data, outcome) {
  d <- data %>%
    filter(!is.na(.data[[outcome]])) %>%
    mutate(Year = as.integer(Year))
  if (nrow(d) < 100L) return(NULL)
  if (length(unique(d$first_treat[d$first_treat > 0])) < 1L) return(NULL)
  tryCatch(
    {
      att <- att_gt_panel(d, outcome, control_group = "nevertreated")
      agg_dyn <- aggte(att, type = "dynamic", na.rm = TRUE)
      list(att = att, agg_dynamic = agg_dyn)
    },
    error = function(e) {
      tryCatch(
        {
          att <- att_gt_panel(d, outcome, control_group = "notyettreated")
          agg_dyn <- aggte(att, type = "dynamic", na.rm = TRUE)
          list(att = att, agg_dynamic = agg_dyn)
        },
        error = function(e2) NULL
      )
    }
  )
}

dynamic_coefs_tibble <- function(agg_dyn, component) {
  if (is.null(agg_dyn) || is.null(agg_dyn$egt)) {
    return(NULL)
  }
  tibble(
    component = component,
    event_time = as.numeric(agg_dyn$egt),
    att = as.numeric(agg_dyn$att.egt),
    se = as.numeric(agg_dyn$se.egt)
  ) %>%
    mutate(
      ci_lower = att - 1.96 * se,
      ci_upper = att + 1.96 * se,
      ci_crosses_zero = ci_lower <= 0 & ci_upper >= 0
    )
}

pval_twfe <- function(m) {
  if (is.null(m)) return(NA_real_)
  tryCatch(as.numeric(fixest::pvalue(m)[1]), error = function(e) NA_real_)
}

pval_cs_simple <- function(agg) {
  if (is.null(agg)) return(NA_real_)
  att <- agg$overall.att
  se <- agg$overall.se
  if (!is.finite(att) || !is.finite(se) || se <= 0) return(NA_real_)
  2 * stats::pnorm(-abs(att / se))
}

# --- Main loop ---
results_ols <- list()
results_twfe <- list()
results_cs <- list()

for (comp in policy_components) {
  ft <- compute_first_treat(df, comp, sd_mult_main)
  d <- df %>% left_join(ft, by = "ISO3")

  for (i in seq_along(outcome_cols)) {
    oc <- outcome_cols[i]
    ol <- names(outcome_cols)[i]
    if (!oc %in% names(d)) next
    key <- paste0(comp, "||", ol)

    # Level 1
    m1 <- tryCatch(
      stats::lm(as.formula(paste0("`", oc, "` ~ `", comp, "`")), data = d %>% filter(!is.na(.data[[oc]]), !is.na(.data[[comp]]))),
      error = function(e) NULL
    )
    p1 <- if (!is.null(m1)) {
      tryCatch(summary(m1)$coefficients[2, 4], error = function(e) NA_real_)
    } else {
      NA_real_
    }
    results_ols[[key]] <- tibble(component = comp, outcome = ol, outcome_col = oc, p_value = p1)

    # Level 2
    m2 <- run_twfe(d, oc, comp)
    results_twfe[[key]] <- tibble(component = comp, outcome = ol, outcome_col = oc, p_value = pval_twfe(m2))

    # Level 3
    cs <- run_cs(d, oc)
    p3 <- if (!is.null(cs)) pval_cs_simple(cs$agg) else NA_real_
    att_hat <- if (!is.null(cs)) cs$agg$overall.att else NA_real_
    att_se <- if (!is.null(cs)) cs$agg$overall.se else NA_real_
    results_cs[[key]] <- tibble(
      component = comp,
      outcome = ol,
      outcome_col = oc,
      att = att_hat,
      se = att_se,
      p_value = p3
    )
  }
}

bind_rows(results_ols) %>% write_csv(file.path(path_results, "results_level1_ols.csv"))
bind_rows(results_twfe) %>% write_csv(file.path(path_results, "results_level2_twfe.csv"))
tbl_cs <- bind_rows(results_cs)
write_csv(tbl_cs, file.path(path_results, "results_level3_cs_att.csv"))

# Level 4: Benjamini–Hochberg (FDR) on Level-3 p-values
ok <- is.finite(tbl_cs$p_value) & !is.na(tbl_cs$p_value)
tbl_cs$fdr_bh <- NA_real_
tbl_cs$reject_bh_05 <- NA
if (sum(ok) > 0) {
  tbl_cs$fdr_bh[ok] <- stats::p.adjust(tbl_cs$p_value[ok], method = "fdr")
  tbl_cs$reject_bh_05 <- tbl_cs$fdr_bh < 0.05 & !is.na(tbl_cs$fdr_bh)
}
write_csv(tbl_cs, file.path(path_results, "results_level4_fdr.csv"))

cat("\nSaved under ", path_results, ": results_level1_ols.csv, results_level2_twfe.csv, results_level3_cs_att.csv, results_level4_fdr.csv\n", sep = "")
cat("Tests (CS):", nrow(tbl_cs), "\n")

# --- Robustness: sd_mult = 1.0 ---
results_cs_robust <- list()
for (comp in policy_components) {
  ft <- compute_first_treat(df, comp, sd_mult_robust)
  d <- df %>% left_join(ft, by = "ISO3")
  for (i in seq_along(outcome_cols)) {
    oc <- outcome_cols[i]
    ol <- names(outcome_cols)[i]
    if (!oc %in% names(d)) next
    key <- paste0(comp, "||", ol)
    cs <- run_cs(d, oc)
    p3 <- if (!is.null(cs)) pval_cs_simple(cs$agg) else NA_real_
    att_hat <- if (!is.null(cs)) cs$agg$overall.att else NA_real_
    att_se <- if (!is.null(cs)) cs$agg$overall.se else NA_real_
    results_cs_robust[[key]] <- tibble(
      component = comp,
      outcome = ol,
      outcome_col = oc,
      att = att_hat,
      se = att_se,
      p_value = p3
    )
  }
}
tbl_robust <- bind_rows(results_cs_robust)
write_csv(tbl_robust, file.path(path_results, "results_level3_cs_att_robust10sd.csv"))
cat("Saved: ", file.path(path_results, "results_level3_cs_att_robust10sd.csv"), " (1.0 * SD)\n", sep = "")

ok_r <- is.finite(tbl_robust$p_value) & !is.na(tbl_robust$p_value)
tbl_robust$fdr_bh <- NA_real_
tbl_robust$reject_bh_05 <- NA
if (sum(ok_r) > 0) {
  tbl_robust$fdr_bh[ok_r] <- stats::p.adjust(tbl_robust$p_value[ok_r], method = "fdr")
  tbl_robust$reject_bh_05 <- tbl_robust$fdr_bh < 0.05 & !is.na(tbl_robust$fdr_bh)
}
write_csv(tbl_robust, file.path(path_results, "results_level3_cs_att_robust10sd_fdr.csv"))

# =============================================================================
# Thesis extractions (event study, distributional, robustness tables)
# =============================================================================

oc_gdp <- outcome_cols[names(outcome_cols) == "GDP_per_capita_growth"][1]

# --- 1) Event study: dynamic coefficients + ggdid plots (GDP outcome) ---
dyn_rows <- list()
agg_dyn_for_plot <- list()

for (comp in policy_components) {
  ft <- compute_first_treat(df, comp, sd_mult_main)
  d <- df %>% left_join(ft, by = "ISO3")
  if (!oc_gdp %in% names(d)) next
  dyn <- run_att_gt_dynamic(d, oc_gdp)
  if (is.null(dyn)) {
    next
  }
  row <- dynamic_coefs_tibble(dyn$agg_dynamic, comp)
  if (!is.null(row)) {
    dyn_rows[[comp]] <- row
  }
  agg_dyn_for_plot[[comp]] <- dyn$agg_dynamic
}

tbl_dyn <- bind_rows(dyn_rows)
if (nrow(tbl_dyn) > 0) {
  write_csv(tbl_dyn, file.path(path_results, "component_dynamic_coefs_GDP_per_capita_growth.csv"))
  pretrend <- tbl_dyn %>% filter(event_time < 0)
  write_csv(pretrend, file.path(path_results, "component_dynamic_coefs_GDP_per_capita_growth_pretrend_only.csv"))
  n_pt <- nrow(pretrend)
  n_cross <- if ("ci_crosses_zero" %in% names(pretrend)) sum(pretrend$ci_crosses_zero, na.rm = TRUE) else NA_integer_
  cat(
    "\n[Event study / pre-trend] GDP growth: ", n_pt,
    " pre-treatment coefficient rows (event_time < 0). Of these, ",
    n_cross, " have 95% CIs that include zero (parallel trends diagnostic).\n",
    sep = ""
  )
}

tbl_gdp <- tbl_cs %>% filter(outcome == "GDP_per_capita_growth", is.finite(p_value))
best_name <- if (nrow(tbl_gdp) > 0) {
  tbl_gdp %>% dplyr::arrange(p_value) %>% dplyr::slice(1) %>% pull(component)
} else {
  NA_character_
}
second_name <- if (nrow(tbl_gdp) > 1) {
  tbl_gdp %>% dplyr::arrange(p_value) %>% dplyr::slice(2) %>% pull(component)
} else {
  NA_character_
}

plot_event <- function(comp_name, suffix) {
  if (is.na(comp_name) || !comp_name %in% names(agg_dyn_for_plot)) {
    return(invisible(NULL))
  }
  ad <- agg_dyn_for_plot[[comp_name]]
  if (is.null(ad)) {
    return(invisible(NULL))
  }
  g <- tryCatch(
    did::ggdid(ad),
    error = function(e) NULL
  )
  if (is.null(g)) {
    return(invisible(NULL))
  }
  fn <- file.path(path_figs, paste0("event_study_gdp_", suffix, ".png"))
  ggplot2::ggsave(fn, g, width = 8, height = 5, dpi = 150)
  cat("Saved event-study figure:", fn, "\n")
}

if (!is.na(best_name)) {
  plot_event(best_name, "main_lowest_p")
}
if (!is.na(second_name) && !is.na(best_name) && second_name != best_name) {
  plot_event(second_name, "appendix_second_lowest_p")
}

# --- 2) Distributional outcomes: failures + survivor table ---
dist_outcomes <- c("Gini_SWIID", "Under5_mortality_rate")
tbl_dist <- tbl_cs %>% filter(outcome %in% dist_outcomes)

for (ol in dist_outcomes) {
  sub <- tbl_dist %>% filter(outcome == ol)
  n_tot <- nrow(sub)
  n_fail <- sum(is.na(sub$att) | is.na(sub$p_value))
  n_ok <- n_tot - n_fail
  cat(
    "\n[Distributional] ", ol, ": ", n_fail, " of ", n_tot,
    " component-level ATT models failed (missing / non-converged).\n",
    sep = ""
  )
  if (n_ok == 0L) {
    cat("  All components failed for this outcome — results are data-limited.\n")
  }
}

tbl_survivors <- tbl_dist %>%
  filter(is.finite(att), is.finite(p_value)) %>%
  select(component, outcome, att, se, p_value, fdr_bh, reject_bh_05)
write_csv(tbl_survivors, file.path(path_results, "table2_distributional_survivors.csv"))
cat("Saved: ", file.path(path_results, "table2_distributional_survivors.csv"), "\n", sep = "")

# --- 3) Cohort size (1.5 vs 1.0 SD) + Table 3 robustness (GDP @ 1.0 SD) ---
cohort_tbl <- purrr::map_dfr(policy_components, function(comp) {
  ft15 <- compute_first_treat(df, comp, sd_mult_main)
  ft10 <- compute_first_treat(df, comp, sd_mult_robust)
  tibble(
    component = comp,
    n_treated_1_5_sd = sum(ft15$first_treat > 0L, na.rm = TRUE),
    n_treated_1_0_sd = sum(ft10$first_treat > 0L, na.rm = TRUE)
  )
})
write_csv(cohort_tbl, file.path(path_results, "cohort_size_by_component_15_vs_10sd.csv"))

mean15 <- mean(cohort_tbl$n_treated_1_5_sd)
mean10 <- mean(cohort_tbl$n_treated_1_0_sd)
pct_increase <- if (is.finite(mean15) && mean15 > 0) {
  100 * (mean10 - mean15) / mean15
} else {
  NA_real_
}
cat(
  "\n[Cohort sizes] Mean treated countries per component @ ", sd_mult_main, " SD: ",
  round(mean15, 2), "; @ ", sd_mult_robust, " SD: ", round(mean10, 2),
  " (avg. increase ", round(pct_increase, 1), "%).\n",
  sep = ""
)

tbl_t3 <- tbl_robust %>%
  filter(outcome == "GDP_per_capita_growth", is.finite(p_value)) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:4) %>%
  select(component, att, se, p_value, fdr_bh, reject_bh_05)
write_csv(tbl_t3, file.path(path_results, "table3_robustness_gdp_10sd_top4.csv"))
cat("Saved: ", file.path(path_results, "table3_robustness_gdp_10sd_top4.csv"), "\n", sep = "")

cat("\nDone.\n")
