# Build data/swiid_clean.csv from SWIID summary (disposable-income Gini, 0–100).
# Raw file: data/raw/swiid9_92/swiid9_92/swiid9_92_summary.csv (or data/raw/swiid_summary.csv).
# Keeps: ISO3, Year, SWIID_gini_disp (from gini_disp).

if (!requireNamespace("countrycode", quietly = TRUE)) {
  install.packages("countrycode", repos = "https://cloud.r-project.org")
}
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
path_nested <- file.path(proj_root, "data", "raw", "swiid9_92", "swiid9_92", "swiid9_92_summary.csv")
path_flat <- file.path(proj_root, "data", "raw", "swiid_summary.csv")
path_out <- file.path(proj_root, "data", "swiid_clean.csv")

src <- if (file.exists(path_nested)) path_nested else if (file.exists(path_flat)) path_flat else NA_character_
if (is.na(src)) {
  stop("SWIID summary not found. Expected ", path_nested, " or ", path_flat)
}

raw <- read_csv(src, show_col_types = FALSE, locale = locale(decimal_mark = "."))

# Names SWIID uses that countrycode treats as ambiguous or multi-match
swiid_country_iso3 <- c(
  "Czechoslovakia" = "CSK",
  "Kosovo" = "XKX",
  "Micronesia" = "FSM",
  "Yugoslavia" = "YUG"
)

swiid_clean <- raw %>%
  transmute(
    country = as.character(.data$country),
    Year = as.integer(.data$year),
    SWIID_gini_disp = as.numeric(.data$gini_disp)
  ) %>%
  mutate(
    ISO3 = dplyr::coalesce(
      unname(swiid_country_iso3[.data$country]),
      countrycode::countrycode(.data$country, "country.name", "iso3c", warn = FALSE)
    )
  ) %>%
  select(ISO3, Year, SWIID_gini_disp) %>%
  filter(!is.na(ISO3), nchar(ISO3) == 3L, !is.na(Year)) %>%
  distinct(ISO3, Year, .keep_all = TRUE)

write_csv(swiid_clean, path_out, na = "")
message("Saved ", path_out, " (", nrow(swiid_clean), " rows).")
