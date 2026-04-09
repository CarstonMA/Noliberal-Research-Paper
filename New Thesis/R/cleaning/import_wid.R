# Build data/wid_clean.csv: pre-tax national income share of bottom 50% (World Inequality Database).
# Preferred: wid::download_wid() from CRAN (requires network).
# Fallback: place a CSV at data/raw/wid_bottom50_ptinc.csv with columns ISO3, Year, value
#           where value is the WID share as a fraction 0–1 (same as API).

req <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
req("tidyverse")
req("countrycode")

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
path_out <- file.path(proj_root, "data", "wid_clean.csv")
path_fallback <- file.path(proj_root, "data", "raw", "wid_bottom50_ptinc.csv")

year_min <- 1995L
year_max <- 2022L

read_fallback_csv <- function() {
  if (!file.exists(path_fallback)) {
    return(NULL)
  }
  x <- read_csv(path_fallback, show_col_types = FALSE)
  iso_c <- intersect(names(x), c("ISO3", "iso3"))[1]
  yr_c <- intersect(names(x), c("Year", "year"))[1]
  val_c <- intersect(names(x), c("value", "share", "WID_income_share_bottom50", "bottom50"))[1]
  if (any(is.na(c(iso_c, yr_c, val_c)))) {
    warning("wid_bottom50_ptinc.csv must have ISO3, Year, and value (or share).")
    return(NULL)
  }
  x %>%
    transmute(
      ISO3 = as.character(.data[[iso_c]]),
      Year = as.integer(.data[[yr_c]]),
      WID_income_share_bottom50 = as.numeric(.data[[val_c]])
    ) %>%
    filter(!is.na(ISO3), nchar(ISO3) == 3L, !is.na(Year))
}

download_wid_bottom50 <- function() {
  if (!requireNamespace("wid", quietly = TRUE)) {
    install.packages("wid", repos = "https://cloud.r-project.org")
  }
  path_wdi <- file.path(proj_root, "data", "wdi_clean.csv")
  # Use areas = "all" for one bulk API pull. Passing a long ISO2 vector makes wid
  # hit the API per area and can take a very long time; we filter to WDI countries after.
  iso3_keep <- NULL
  if (file.exists(path_wdi)) {
    iso3_keep <- unique(read_csv(path_wdi, show_col_types = FALSE)$ISO3)
    iso3_keep <- as.character(iso3_keep[!is.na(iso3_keep)])
  }
  d <- tryCatch(
    wid::download_wid(
      indicators = "sptinc",
      areas = "all",
      years = year_min:year_max,
      perc = "p0p50",
      ages = 999,
      pop = "i",
      include_extrapolations = TRUE,
      verbose = TRUE
    ),
    error = function(e) {
      warning("WID API download failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(d) || nrow(d) == 0L) {
    return(NULL)
  }
  # API returns 2-letter area codes; convert to ISO3
  d <- as_tibble(d)
  if ("variable" %in% names(d)) {
    d <- d %>% filter(grepl("^sptinc", as.character(.data$variable)))
  }
  d <- d %>%
    mutate(percentile = trimws(as.character(.data$percentile))) %>%
    filter(.data$percentile == "p0p50") %>%
    mutate(
      ISO3 = countrycode::countrycode(.data$country, origin = "iso2c", destination = "iso3c"),
      Year = as.integer(.data$year),
      # WID shares are fractions of 1; store as percent 0–100 for comparability with SWIID Gini scale
      WID_income_share_bottom50 = as.numeric(.data$value) * 100
    ) %>%
    filter(!is.na(ISO3), !is.na(Year)) %>%
    distinct(ISO3, Year, .keep_all = TRUE) %>%
    select(ISO3, Year, WID_income_share_bottom50)
  if (!is.null(iso3_keep) && length(iso3_keep) >= 10L) {
    d <- d %>% filter(.data$ISO3 %in% iso3_keep)
  }
  d
}

if (file.exists(path_out)) {
  message("wid_clean.csv already exists; delete it to rebuild.")
} else {
  out <- download_wid_bottom50()
  if (is.null(out)) {
    message("Trying fallback file: ", path_fallback)
    fb <- read_fallback_csv()
    if (!is.null(fb)) {
      mx <- suppressWarnings(max(fb$WID_income_share_bottom50, na.rm = TRUE))
      out <- if (is.finite(mx) && mx <= 1) {
        fb %>% mutate(WID_income_share_bottom50 = WID_income_share_bottom50 * 100)
      } else {
        fb
      }
    }
  }
  if (is.null(out) || nrow(out) == 0L) {
    stop(
      "Could not build WID panel. Install package 'wid' and run online, or add ",
      "data/raw/wid_bottom50_ptinc.csv with columns ISO3, Year, value (share 0–1 or 0–100)."
    )
  }
  write_csv(out, path_out, na = "")
  message("Saved ", path_out, " (", nrow(out), " rows).")
}
