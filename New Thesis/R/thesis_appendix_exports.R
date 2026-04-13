# Thesis appendix exports (event studies, dynamic coef tables, robustness CSV, cohort lists).
# Requires merged data and prior successful methodology setup (same as analysis_methodology.R).
#
# Run from New Thesis:
#   Rscript R/thesis_appendix_exports.R
#
# Optional: only flagship plot + tables (skip 68 event-study re-estimation):
#   set THESIS_APPENDIX_LIGHT=1  (Windows: $env:THESIS_APPENDIX_LIGHT="1")
#
# Outputs:
#   results/figures/appendix/A_event_studies/individual/*.png (one per outcome×component); flagship PNG (600 dpi)
#   results/appendix/ — B, D CSVs + manifest (pooled ATT + FDR tables: use results/final_tables/did_fdr_*.csv)

req <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
req("tidyverse")
req("fixest")
req("did")

library(tidyverse)
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
source(file.path(proj_root, "R", "event_study_plot.R"))
path_results <- file.path(proj_root, "results")
path_appendix <- file.path(path_results, "appendix")
path_figs_appendix <- file.path(path_results, "figures", "appendix")
path_figs_A <- file.path(path_figs_appendix, "A_event_studies")
path_figs_A_indiv <- file.path(path_figs_A, "individual")
for (p in c(path_appendix, path_figs_appendix, path_figs_A, path_figs_A_indiv)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

min_year <- 1995L
max_year <- 2022L
min_obs_country <- 10L
sd_mult_main <- 1.5
bootstrap_iters <- 199L

policy_components <- c(
  "fraser_SizeGov", "fraser_LegalSystem", "fraser_SoundMoney", "fraser_FreeTrade", "fraser_Regulation",
  "heritage_Property.Rights", "heritage_Government.Integrity", "heritage_Judicial.Effectiveness",
  "heritage_Tax.Burden", "heritage_Government.Spending", "heritage_Fiscal.Health",
  "heritage_Business.Freedom", "heritage_Labor.Freedom", "heritage_Monetary.Freedom",
  "heritage_Trade.Freedom", "heritage_Investment.Freedom", "heritage_Financial.Freedom"
)

path_merged <- file.path(proj_root, "data", "merged_full_data.csv")
if (!file.exists(path_merged)) {
  stop("Missing data/merged_full_data.csv — run cleaning/merge first.")
}

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
outcome_cols <- outcomes_tbl$col
names(outcome_cols) <- outcomes_tbl$label
policy_components <- intersect(policy_components, names(df))

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

run_att_gt_dynamic <- function(data, outcome) {
  d <- data %>%
    filter(!is.na(.data[[outcome]])) %>%
    mutate(Year = as.integer(Year))
  if (nrow(d) < 100L) {
    return(NULL)
  }
  if (length(unique(d$first_treat[d$first_treat > 0])) < 1L) {
    return(NULL)
  }
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

dynamic_coefs_tibble <- function(agg_dyn, component, outcome_label) {
  if (is.null(agg_dyn) || is.null(agg_dyn$egt)) {
    return(NULL)
  }
  tibble(
    outcome = outcome_label,
    component = component,
    event_time = as.numeric(agg_dyn$egt),
    att = as.numeric(agg_dyn$att.egt),
    se = as.numeric(agg_dyn$se.egt)
  ) %>%
    mutate(
      ci_lower = att - 1.96 * se,
      ci_upper = att + 1.96 * se,
      ci_crosses_zero = ci_lower <= 0 & ci_upper >= 0,
      period = if_else(event_time < 0, "pre_treatment", "post_treatment")
    )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# ggdid + thesis layers: trimmed window (cell counts / SE), plain-English title, no Pre/Post legend
save_event_study_plots <- function(agg_dyn, outcome_label, component, path_base_no_ext, d,
                                    title_override = NULL, dpi = 300, w = 10, h = 6) {
  g <- tryCatch(
    thesis_ggdid(agg_dyn, outcome_label, component, d = d, title_override = title_override),
    error = function(e) NULL
  )
  if (is.null(g)) {
    return(FALSE)
  }
  ggplot2::ggsave(paste0(path_base_no_ext, ".png"), g, width = w, height = h, dpi = dpi, bg = "white")
  TRUE
}

light_mode <- Sys.getenv("THESIS_APPENDIX_LIGHT", "") == "1"

# --- Appendix D: cohort lists @ 1.5 SD ---
cohort_long <- purrr::map_dfr(policy_components, function(comp) {
  ft <- compute_first_treat(df, comp, sd_mult_main)
  ft %>%
    filter(first_treat > 0L) %>%
    transmute(
      sd_mult_shock = sd_mult_main,
      policy_component = comp,
      ISO3,
      first_treat_year = first_treat
    )
})
write_csv(cohort_long, file.path(path_appendix, "D_cohort_treatment_year_sd15.csv"))
cat("Wrote ", path_appendix, "/D_cohort_treatment_year_sd15.csv\n", sep = "")

# --- Appendix C (pooled ATT + FDR): no duplicate CSVs here — use results/final_tables/did_fdr_main_15sd.csv
#     and did_fdr_robust_10sd.csv (from R/analysis_methodology.R → table_exporters.R).

# --- Event studies + Appendix B (dynamic coefs) ---
manifest <- list()
dyn_all <- list()

flagship_comp <- "heritage_Tax.Burden"
flagship_oc <- outcome_cols[["GDP_per_capita_growth"]]
flagship_label <- "GDP_per_capita_growth"

run_one_cell <- function(outcome_label, oc, comp, save_plots = TRUE) {
  ft <- compute_first_treat(df, comp, sd_mult_main)
  d <- df %>% left_join(ft, by = "ISO3")
  if (!oc %in% names(d)) {
    return(list(ok = FALSE, reason = "missing outcome column"))
  }
  dyn <- run_att_gt_dynamic(d, oc)
  if (is.null(dyn)) {
    return(list(ok = FALSE, reason = "att_gt/aggte dynamic failed"))
  }
  row <- dynamic_coefs_tibble(dyn$agg_dynamic, comp, outcome_label)
  if (is.null(row)) {
    return(list(ok = FALSE, reason = "no egt"))
  }
  base_name <- safe_filename(paste(outcome_label, comp, sep = "__"))
  if (save_plots) {
    indiv <- file.path(path_figs_A_indiv, base_name)
    save_event_study_plots(dyn$agg_dynamic, outcome_label, comp, indiv, d = d, dpi = 300, w = 10, h = 6)
  }
  list(ok = TRUE, agg_dynamic = dyn$agg_dynamic, coefs = row, base_name = base_name, d = d)
}

if (!light_mode) {
  for (ol in names(outcome_cols)) {
    oc <- outcome_cols[[ol]]
    cat("\nOutcome:", ol, "\n")
    for (comp in policy_components) {
      res <- tryCatch(run_one_cell(ol, oc, comp, save_plots = TRUE), error = function(e) {
        list(ok = FALSE, reason = conditionMessage(e))
      })
      manifest[[paste(ol, comp, sep = "||")]] <- tibble(
        outcome = ol,
        component = comp,
        success = isTRUE(res$ok),
        reason = if (isTRUE(res$ok)) {
          NA_character_
        } else {
          as.character(if (is.null(res$reason)) "failed" else res$reason)
        }
      )
      if (isTRUE(res$ok)) {
        dyn_all[[paste(ol, comp, sep = "||")]] <- res$coefs
      }
    }
    cat("Finished outcome:", ol, "(PNGs in A_event_studies/individual/)\n")
  }
}

# Flagship: Heritage Tax Burden x GDP (always run if data exist)
if (flagship_oc %in% names(df) && flagship_comp %in% names(df)) {
  ft <- compute_first_treat(df, flagship_comp, sd_mult_main)
  d <- df %>% left_join(ft, by = "ISO3")
  dyn_f <- run_att_gt_dynamic(d, flagship_oc)
  if (!is.null(dyn_f)) {
    fp_base <- file.path(path_figs_appendix, "flagship_event_study_heritage_tax_burden_gdp")
    save_event_study_plots(
      dyn_f$agg_dynamic,
      flagship_label,
      flagship_comp,
      fp_base,
      d = d,
      title_override = "Heritage Tax Burden and GDP per Capita Growth (Dynamic ATT)",
      dpi = 600,
      w = 10,
      h = 6.5
    )
    cat("Flagship (600 dpi PNG): ", fp_base, ".png\n", sep = "")
  } else {
    warning("Flagship Heritage Tax Burden × GDP: dynamic estimation failed.")
  }
}

# Appendix B: bind all dynamic coef rows
if (!light_mode && length(dyn_all) > 0) {
  tbl_b <- bind_rows(dyn_all) %>%
    filter(!is.na(event_time)) %>%
    arrange(outcome, component, event_time)
  write_csv(tbl_b, file.path(path_appendix, "B_dynamic_att_event_time_all_pairs.csv"))
  pret <- tbl_b %>% filter(event_time < 0)
  post <- tbl_b %>% filter(event_time >= 0)
  write_csv(pret, file.path(path_appendix, "B_dynamic_att_pre_treatment_only.csv"))
  write_csv(post, file.path(path_appendix, "B_dynamic_att_post_treatment_only.csv"))
  cat("Wrote appendix B CSVs (pre/post splits + full).\n")
}

if (!light_mode && length(manifest) > 0) {
  bind_rows(manifest) %>%
    write_csv(file.path(path_appendix, "A_event_study_manifest.csv"))
  cat("Wrote A_event_study_manifest.csv\n")
}

cat("\nDone thesis_appendix_exports.R\n")
