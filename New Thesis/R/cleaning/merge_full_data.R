# Merge WDI economic outcomes with Fraser/Heritage policy indices
# Input:  data/wdi_clean.csv, data/merged_fraser_heritage.csv
# Output: data/merged_full_data.csv
# Creates growth variables and harmonized column names for analysis

library(tidyverse)

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
path_wdi <- file.path(proj_root, "data", "wdi_clean.csv")
path_policy <- file.path(proj_root, "data", "merged_fraser_heritage.csv")
path_out <- file.path(proj_root, "data", "merged_full_data.csv")

# 1. Load data
cat("Loading data...\n")
wdi <- read_csv(path_wdi, show_col_types = FALSE)
policy <- read_csv(path_policy, show_col_types = FALSE)

# 2. Harmonize WDI: ensure Year is integer
wdi <- wdi %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(ISO3) & !is.na(Year))

# 3. Policy: ensure Year is integer
policy <- policy %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(ISO3) & !is.na(Year))

# 4. Restrict to 1995-2022 (policy coverage + pre-COVID/partial-year)
wdi <- wdi %>% filter(Year >= 1995, Year <= 2022)
policy <- policy %>% filter(Year >= 1995, Year <= 2022)

# 5. Inner join on ISO3 + Year (keeps only overlapping country-years)
merged <- inner_join(
  policy,
  wdi,
  by = c("ISO3", "Year"),
  relationship = "one-to-one"
)

# 6. Harmonize country name (prefer policy's Country)
merged <- merged %>%
  mutate(Country = coalesce(Country, `Country Name`)) %>%
  select(-`Country Name`)

# 7. Safe column names for R (dots for spaces, etc.)
names(merged) <- make.names(names(merged), unique = TRUE)

# 8. Find GDP per capita column (PPP constant 2021)
gdp_pc_col <- names(merged)[grepl("GDP.*capita.*PPP|NY.GDP.PCAP.PP", names(merged), ignore.case = TRUE)][1]
if (is.na(gdp_pc_col)) {
  gdp_pc_col <- names(merged)[grepl("GDP.*capita|GDPPCAP", names(merged), ignore.case = TRUE)][1]
}

# 9. Create growth variables
if (!is.null(gdp_pc_col) && gdp_pc_col %in% names(merged)) {
  merged <- merged %>%
    arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(
      log_GDP_per_capita = log(.data[[gdp_pc_col]] + 1),
      GDP_per_capita_growth = 100 * (log(.data[[gdp_pc_col]] + 1) - log(dplyr::lag(.data[[gdp_pc_col]]) + 1))
    ) %>%
    ungroup()
  cat("Created log_GDP_per_capita, GDP_per_capita_growth from", gdp_pc_col, "\n")
} else {
  cat("Warning: GDP per capita column not found. Skipping growth vars.\n")
  merged <- merged %>%
    mutate(log_GDP_per_capita = NA_real_, GDP_per_capita_growth = NA_real_)
}

# 10. Time-invariant binary splits (baseline cross section; use AT MOST ONE for sub-samples)
#     Literature alignment: "high voice" ≈ stronger democratic institutions (WGI);
#     "high GDP" ≈ richer at baseline (development level). Median splits at first
#     observed year in 1995–2000 per country (avoids missing 1995 only).
#     Re-run clean_wdi.R first so VA.EST is present in wdi_clean.csv.
va_col <- names(merged)[names(merged) %in% c("VA.EST", "VA_EST")][1]
gdp_for_split <- names(merged)[grepl("^NY\\.GDP\\.PCAP\\.PP", names(merged))][1]

if (!is.na(gdp_for_split)) {
  ref_gdp <- merged %>%
    filter(Year >= 1995L, Year <= 2000L, !is.na(.data[[gdp_for_split]])) %>%
    group_by(ISO3) %>%
    slice_min(Year, n = 1L, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(ISO3, gdp_baseline = .data[[gdp_for_split]])
  med_g <- stats::median(ref_gdp$gdp_baseline, na.rm = TRUE)
  ref_gdp <- ref_gdp %>%
    mutate(high_gdp_pc_ppp_baseline = as.integer(gdp_baseline > med_g))
  merged <- merged %>% left_join(ref_gdp %>% select(ISO3, high_gdp_pc_ppp_baseline), by = "ISO3")
} else {
  merged$high_gdp_pc_ppp_baseline <- NA_integer_
}

if (!is.na(va_col) && va_col %in% names(merged)) {
  ref_va <- merged %>%
    filter(Year >= 1995L, Year <= 2000L, !is.na(.data[[va_col]])) %>%
    group_by(ISO3) %>%
    slice_min(Year, n = 1L, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(ISO3, va_baseline = .data[[va_col]])
  med_v <- stats::median(ref_va$va_baseline, na.rm = TRUE)
  ref_va <- ref_va %>%
    mutate(high_voice_accountability_baseline = as.integer(va_baseline > med_v))
  merged <- merged %>% left_join(ref_va %>% select(ISO3, high_voice_accountability_baseline), by = "ISO3")
} else {
  merged$high_voice_accountability_baseline <- NA_integer_
  cat("Note: VA.EST not in panel — run R/cleaning/clean_wdi.R to add WGI Voice & Accountability.\n")
}

# 11. Relocate key columns
merged <- merged %>%
  relocate(
    ISO3, Country, Year,
    any_of(c("high_gdp_pc_ppp_baseline", "high_voice_accountability_baseline")),
    fraser_Summary, heritage_Overall
  )

# 12. Save
write_csv(merged, path_out, na = "")

# 13. Summary
cat("\n=== MERGE SUMMARY ===\n")
cat("Rows:", nrow(merged), "\n")
cat("Years:", min(merged$Year, na.rm = TRUE), "to", max(merged$Year, na.rm = TRUE), "\n")
cat("Countries:", n_distinct(merged$ISO3), "\n")
cat("Non-missing fraser_Summary:", sum(!is.na(merged$fraser_Summary)), "\n")
cat("Non-missing heritage_Overall:", sum(!is.na(merged$heritage_Overall)), "\n")
if (exists("gdp_pc_col") && gdp_pc_col %in% names(merged)) {
  cat("Non-missing", gdp_pc_col, ":", sum(!is.na(merged[[gdp_pc_col]])), "\n")
}
gini_col <- names(merged)[grepl("Gini|GINI", names(merged))][1]
if (!is.na(gini_col)) {
  cat("Non-missing", gini_col, ":", sum(!is.na(merged[[gini_col]])), "\n")
}
if ("high_gdp_pc_ppp_baseline" %in% names(merged)) {
  cat(
    "Binary high_gdp_pc_ppp_baseline (1=above median baseline GDP pc PPP):",
    " n1=", sum(merged$high_gdp_pc_ppp_baseline == 1L, na.rm = TRUE),
    " n0=", sum(merged$high_gdp_pc_ppp_baseline == 0L, na.rm = TRUE),
    " NA=", sum(is.na(merged$high_gdp_pc_ppp_baseline)), "\n"
  )
}
if ("high_voice_accountability_baseline" %in% names(merged)) {
  cat(
    "Binary high_voice_accountability_baseline (1=above median WGI Voice):",
    " n1=", sum(merged$high_voice_accountability_baseline == 1L, na.rm = TRUE),
    " n0=", sum(merged$high_voice_accountability_baseline == 0L, na.rm = TRUE),
    " NA=", sum(is.na(merged$high_voice_accountability_baseline)), "\n"
  )
}
cat("Correlation fraser_Summary vs heritage_Overall:",
    round(cor(merged$fraser_Summary, merged$heritage_Overall, use = "pairwise.complete.obs"), 3), "\n")
cat("\nSaved merged_full_data.csv -> ", path_out, "\n", sep = "")
