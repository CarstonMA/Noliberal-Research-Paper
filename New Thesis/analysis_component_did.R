# Component-level DiD event studies
# Treatment = sharp jump in each specific component (not aggregate)
# Produces event study plots for each component-outcome pair with significant FE effect
# Run analysis_components.R first to create component_impacts_significant.csv

library(tidyverse)
library(fixest)

if (!file.exists("merged_full_data.csv")) source("merge_full_data.R")
if (!file.exists("component_impacts_significant.csv")) {
  stop("Run analysis_components.R first to create component_impacts_significant.csv")
}

# Load significant component-outcome pairs
sig_pairs <- read_csv("component_impacts_significant.csv", show_col_types = FALSE)
df_full <- read_csv("merged_full_data.csv", show_col_types = FALSE)
df_full <- df_full %>% filter(Year >= 1995, Year <= 2022, !is.na(ISO3), !is.na(Year))

# Panel: min 15 observations (avoids losing Gini/Poverty)
n_years <- length(unique(df_full$Year))
panel_iso3 <- df_full %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)
df_full <- df_full %>% filter(ISO3 %in% panel_iso3)
cat("Panel:", length(panel_iso3), "countries (min 15 obs),", n_years, "years\n")

# Region definitions (match analysis_event_study.R)
post_soviet <- c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", "LTU", "MDA", "RUS", "TJK", "TKM", "UKR", "UZB")
latin_america <- c("ARG", "BOL", "BRA", "BRB", "BLZ", "CHL", "COL", "CRI", "CUB", "DOM", "ECU", "SLV", "GTM", "GUY", "HND", "JAM", "MEX", "NIC", "PAN", "PRY", "PER", "TTO", "URY", "VEN")
sub_saharan_africa <- c("AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CAF", "TCD", "COM", "COG", "COD", "CIV", "GNQ", "ERI", "SWZ", "ETH", "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TZA", "TGO", "UGA", "ZMB", "ZWE")
east_asia_pacific <- c("AUS", "BGD", "BTN", "KHM", "CHN", "FJI", "HKG", "IDN", "JPN", "KIR", "LAO", "MYS", "MMR", "MNG", "NPL", "NZL", "PHL", "PNG", "KOR", "WSM", "SGP", "SLB", "LKA", "THA", "TON", "VUT", "VNM")
middle_east_na <- c("DZA", "BHR", "EGY", "IRN", "IRQ", "ISR", "JOR", "KWT", "LBN", "LBY", "MAR", "OMN", "QAT", "SAU", "SYR", "TUN", "ARE", "YEM", "PSE")
eastern_europe <- c("ALB", "BIH", "BGR", "HRV", "CZE", "HUN", "MKD", "POL", "ROU", "SRB", "SVK", "SVN")
western_europe_na <- c("AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "DEU", "GRC", "ISL", "IRL", "ITA", "LUX", "MLT", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "GBR", "USA")
df_full <- df_full %>% mutate(
  country_group = case_when(
    ISO3 %in% post_soviet ~ "Post-Soviet", ISO3 %in% latin_america ~ "Latin America",
    ISO3 %in% sub_saharan_africa ~ "Sub-Saharan Africa", ISO3 %in% east_asia_pacific ~ "East Asia & Pacific",
    ISO3 %in% middle_east_na ~ "Middle East & N. Africa", ISO3 %in% eastern_europe ~ "Eastern Europe",
    ISO3 %in% western_europe_na ~ "Western Europe & N. America", TRUE ~ "Other"
  ),
  region_year = paste0(country_group, "_", Year)
)

# Get unique components and outcomes (exclude aggregate indices)
components <- unique(sig_pairs$component)
components <- components[components %in% names(df_full)]
components <- components[!components %in% c("fraser_Summary", "heritage_Overall")]
outcomes <- unique(sig_pairs$outcome)
outcomes <- outcomes[outcomes %in% names(df_full)]

# Winsorize outcomes at 1% and 99%
winsorize <- function(x, p_low = 0.01, p_high = 0.99) {
  q <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
df_full <- df_full %>% mutate(across(all_of(outcomes), winsorize, .names = "{.col}"))

# No lagged DV (avoids Nickell bias)
event_times <- c(-5:-2, 0:5)   # 5 years post-treatment
identify_treatment <- function(data, pvar, pct = 0.90) {
  d <- data %>% filter(!is.na(.data[[pvar]])) %>% arrange(ISO3, Year) %>%
    group_by(ISO3) %>% mutate(lag_p = dplyr::lag(.data[[pvar]]), jump = .data[[pvar]] - lag_p) %>%
    ungroup() %>% filter(!is.na(jump))
  if (nrow(d) == 0) return(tibble(ISO3 = character(), treated = logical(), event_year = numeric()))
  jump_threshold <- quantile(d$jump, probs = pct, na.rm = TRUE)
  d %>% filter(jump >= jump_threshold) %>%
    group_by(ISO3) %>% slice(1) %>% ungroup() %>%
    transmute(ISO3, treated = TRUE, event_year = Year)
}

run_event_study <- function(data, outcome, event_times, control_var = NULL) {
  data_sub <- data %>% filter(!is.na(.data[[outcome]]))
  if (!is.null(control_var) && control_var %in% names(data) && control_var != outcome)
    data_sub <- data_sub %>% filter(!is.na(.data[[control_var]]))
  if (nrow(data_sub) < 100 || sum(data_sub$treated) < 10) return(NULL)
  event_terms <- paste0("`event_time_", event_times, "`", collapse = " + ")
  out_esc <- if (grepl("[^a-zA-Z0-9._]", outcome)) paste0("`", outcome, "`") else outcome
  rhs <- event_terms
  if (!is.null(control_var) && control_var != outcome && control_var %in% names(data_sub))
    rhs <- paste(rhs, paste0("`", control_var, "`"), sep = " + ")
  f <- as.formula(paste0(out_esc, " ~ ", rhs, " | ISO3 + Year"))
  tryCatch(feols(f, data = data_sub, cluster = ~ ISO3), error = function(e) NULL)
}

plot_event_study <- function(model, event_times) {
  if (is.null(model)) return(NULL)
  coefs <- coef(model)
  ses <- setNames(se(model), names(coefs))
  pre_t <- event_times[event_times < 0]
  post_t <- event_times[event_times >= 0]
  plot_times <- c(pre_t, -1, post_t)
  all_names <- names(coefs)
  est <- numeric(length(plot_times))
  se_est <- numeric(length(plot_times))
  est[] <- NA_real_
  se_est[] <- NA_real_
  for (i in seq_along(plot_times)) {
    k <- plot_times[i]
    if (k == -1) { est[i] <- 0; se_est[i] <- 0; next }
    vn <- paste0("event_time_", k)
    vn_bt <- paste0("`event_time_", k, "`")
    name_use <- if (vn_bt %in% all_names) vn_bt else if (vn %in% all_names) vn else NULL
    if (!is.null(name_use)) { est[i] <- coefs[name_use]; se_est[i] <- ses[name_use] }
  }
  tibble(event_time = plot_times, estimate = est, se = se_est,
         ci_lower = est - 1.96 * se_est, ci_upper = est + 1.96 * se_est) %>%
    filter(!is.na(estimate) | event_time == -1)
}

control_var <- if ("log_GDP_per_capita" %in% names(df_full)) "log_GDP_per_capita" else NULL
med_ev_default <- 2000

dir.create("component_did_plots", showWarnings = FALSE)
comp_labels <- c(
  fraser_Summary = "Fraser Summary", fraser_SizeGov = "Size of Government",
  fraser_LegalSystem = "Legal System", fraser_SoundMoney = "Sound Money",
  fraser_FreeTrade = "Free Trade", fraser_Regulation = "Regulation",
  heritage_Overall = "Heritage Overall", heritage_Property.Rights = "Property Rights",
  heritage_Government.Integrity = "Government Integrity",
  heritage_Judicial.Effectiveness = "Judicial Effectiveness",
  heritage_Tax.Burden = "Tax Burden", heritage_Government.Spending = "Government Spending",
  heritage_Fiscal.Health = "Fiscal Health", heritage_Business.Freedom = "Business Freedom",
  heritage_Labor.Freedom = "Labor Freedom", heritage_Monetary.Freedom = "Monetary Freedom",
  heritage_Trade.Freedom = "Trade Freedom", heritage_Investment.Freedom = "Investment Freedom",
  heritage_Financial.Freedom = "Financial Freedom"
)

treatment_counts <- tibble()
for (comp in components) {
  cat("\n=== Component:", comp, "===\n")
  treatment_info <- identify_treatment(df_full, comp)
  n_t <- nrow(treatment_info)
  treatment_counts <- bind_rows(treatment_counts, tibble(component = comp, n_treated = n_t))
  if (n_t < 10) {
    cat("  Skipping: only", n_t, "treated countries\n")
    next
  }
  cat("  Treated countries:", n_t, "\n")
  df <- df_full %>%
    left_join(treatment_info, by = "ISO3") %>%
    mutate(
      treated = replace_na(treated, FALSE),
      event_year = if_else(treated, event_year, NA_real_),
      event_time = if_else(treated, Year - event_year, NA_real_),
      event_time = if_else(event_time >= -5 & event_time <= 5, event_time, NA_real_)
    )
  for (k in event_times) {
    vn <- paste0("event_time_", k)
    df[[vn]] <- as.numeric(df$event_time == k & df$treated)
    df[[vn]] <- replace_na(df[[vn]], 0)
  }
  comp_label <- comp_labels[comp]
  if (is.na(comp_label)) comp_label <- gsub("^fraser_|^heritage_", "", comp)
  treatment_def <- paste0("Treatment: top-decile jump in ", comp_label)

  outcomes_for_comp <- sig_pairs %>% filter(component == comp) %>% pull(outcome) %>% unique()
  outcomes_for_comp <- intersect(outcomes_for_comp, outcomes)

  med_ev <- median(df$event_year[df$treated], na.rm = TRUE)
  if (is.na(med_ev)) med_ev <- med_ev_default

  for (out in outcomes_for_comp) {
    m <- run_event_study(df, out, event_times, control_var)
    if (is.null(m)) next
    pd <- plot_event_study(m, event_times)
    if (is.null(pd) || nrow(pd) == 0) next
    pd <- pd %>% mutate(period = if_else(event_time < 0, "Pre", "Post"))

    y_min <- min(pd$ci_lower, na.rm = TRUE)
    y_max <- max(pd$ci_upper, na.rm = TRUE)
    y_span <- max(y_max - y_min, 0.1)
    out_short <- gsub("\\..*", "", out)
    if (nchar(out_short) < 5) out_short <- out

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
                        breaks = c("Pre", "Post"), na.translate = FALSE) +
      labs(x = "Years relative to treatment", y = "Coefficient / Normalized outcome",
           title = paste0("DiD: ", comp_label, " → ", out_short),
           subtitle = paste0(treatment_def, ". Blue = pre, red = post."),
           color = "Period") +
      expand_limits(y = 0) +
      theme_minimal() +
      theme(legend.position = "bottom")
    fname <- paste0("component_did_plots/did_", gsub("[^a-zA-Z0-9]", "_", comp), "_", gsub("[^a-zA-Z0-9]", "_", out), ".png")
    ggsave(fname, p, width = 7, height = 5, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
}

if (nrow(treatment_counts) > 0) {
  write_csv(treatment_counts, "component_did_treatment_counts.csv")
  cat("\n=== Treatment counts by component ===\n")
  print(treatment_counts)
  cat("\nSaved component_did_treatment_counts.csv\n")
}
cat("\nDone. Component DiD plots saved to component_did_plots/\n")
