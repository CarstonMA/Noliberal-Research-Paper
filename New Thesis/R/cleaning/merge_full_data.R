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

# 10. Relocate key columns
merged <- merged %>%
  relocate(ISO3, Country, Year, fraser_Summary, heritage_Overall)

# 11. Save
write_csv(merged, path_out, na = "")

# 12. Summary
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
cat("Correlation fraser_Summary vs heritage_Overall:",
    round(cor(merged$fraser_Summary, merged$heritage_Overall, use = "pairwise.complete.obs"), 3), "\n")
cat("\nSaved merged_full_data.csv -> ", path_out, "\n", sep = "")
