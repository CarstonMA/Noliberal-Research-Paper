# Build data/wid_clean.csv: pre-tax national income share of bottom 50% (World Inequality Database).
# Preferred: wid::download_wid() from CRAN (requires network).
# Fallback: place a CSV at data/raw/wid_bottom50_ptinc.csv with columns ISO3, Year, value
#           where value is the WID share as a fraction 0–1 (same as API).
# Panel CS needs many country-years: analysis_methodology.R requires >=100 non-missing outcome
# rows; a single-country stub (e.g. USA only) will merge but all WID CS cells stay NA.

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

# WID area codes: usually ISO2 ("US"); sometimes ISO3-like 3 letters — map to ISO3
iso3_from_wid_area <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  idx3 <- !is.na(x) & nchar(x) == 3L & grepl("^[A-Za-z]{3}$", x)
  out[idx3] <- toupper(x[idx3])
  idx2 <- is.na(out) & !is.na(x) & nchar(x) == 2L
  if (any(idx2)) {
    out[idx2] <- countrycode::countrycode(x[idx2], origin = "iso2c", destination = "iso3c")
  }
  out
}

download_wid_bottom50 <- function() {
  if (!requireNamespace("wid", quietly = TRUE)) {
    install.packages("wid", repos = "https://cloud.r-project.org")
  }
  path_wdi <- file.path(proj_root, "data", "wdi_clean.csv")
  iso3_keep <- NULL
  iso2_list <- NULL
  if (file.exists(path_wdi)) {
    w <- read_csv(path_wdi, show_col_types = FALSE)
    iso3_keep <- unique(as.character(w$ISO3[!is.na(w$ISO3)]))
    iso2_list <- unique(countrycode::countrycode(iso3_keep, origin = "iso3c", destination = "iso2c"))
    iso2_list <- sort(iso2_list[!is.na(iso2_list) & nchar(iso2_list) == 2L])
  }
  # pop = "j" (equal-split adults): WID publishes sptinc bottom-50 for many countries.
  # pop = "i" (individuals) largely returns US-only for this series — do not use.
  wid_batch_size <- 25L
  parts <- list()
  if (!is.null(iso2_list) && length(iso2_list) >= 5L) {
    chunks <- split(iso2_list, ceiling(seq_along(iso2_list) / wid_batch_size))
    message(
      "Downloading WID sptinc p0p50 in ", length(chunks), " batch(es) of ~",
      wid_batch_size, " countries (pop = equal-split adults)..."
    )
    for (k in seq_along(chunks)) {
      chunk <- chunks[[k]]
      d_batch <- tryCatch(
        wid::download_wid(
          indicators = "sptinc",
          areas = chunk,
          years = year_min:year_max,
          perc = "p0p50",
          ages = 999L,
          pop = "j",
          include_extrapolations = TRUE,
          verbose = TRUE
        ),
        error = function(e) {
          warning("WID API batch ", k, "/", length(chunks), " failed: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(d_batch) && nrow(d_batch) > 0L) {
        parts[[length(parts) + 1L]] <- d_batch
      }
      Sys.sleep(0.35)
    }
  }
  if (length(parts) == 0L) {
    message("Falling back to areas = 'all' (single request; may be incomplete)...")
    d <- tryCatch(
      wid::download_wid(
        indicators = "sptinc",
        areas = "all",
        years = year_min:year_max,
        perc = "p0p50",
        ages = 999L,
        pop = "j",
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
  } else {
    d <- dplyr::bind_rows(parts)
  }
  d <- as_tibble(d)
  if ("variable" %in% names(d)) {
    d <- d %>% filter(grepl("^sptinc", as.character(.data$variable)))
  }
  country_col <- if ("country" %in% names(d)) {
    "country"
  } else if ("area" %in% names(d)) {
    "area"
  } else {
    stop("WID download: no 'country' or 'area' column in API result.")
  }
  d <- d %>%
    mutate(percentile = trimws(as.character(.data$percentile))) %>%
    filter(.data$percentile == "p0p50") %>%
    mutate(
      ISO3 = iso3_from_wid_area(.data[[country_col]]),
      Year = as.integer(.data$year),
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

# If wid_clean.csv is a stub (e.g. USA-only), re-download when API works — otherwise CS step stays NA.
force_rebuild <- Sys.getenv("WID_FORCE_REBUILD", "") == "1"
stub <- FALSE
if (file.exists(path_out)) {
  prev <- tryCatch(
    read_csv(path_out, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(prev) && nrow(prev) > 0L && "ISO3" %in% names(prev)) {
    n_iso <- length(unique(prev$ISO3[!is.na(prev$ISO3)]))
    stub <- nrow(prev) < 200L || n_iso < 10L
  }
}
run_download <- !file.exists(path_out) || force_rebuild || stub
if (!run_download) {
  message(
    "wid_clean.csv already exists (full panel). Delete data/wid_clean.csv or set ",
    "WID_FORCE_REBUILD=1 to rebuild."
  )
} else {
  if (file.exists(path_out)) {
    message(
      if (stub) {
        "Attempting to replace stub wid_clean.csv (< 200 rows or < 10 countries).\n"
      } else if (force_rebuild) {
        "WID_FORCE_REBUILD=1: rebuilding wid_clean.csv.\n"
      } else {
        "Building wid_clean.csv.\n"
      }
    )
  }
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
    warning(
      paste0(
        "Could not build WID panel (API unreachable or no fallback file).\n",
        "  • Retry later with network, or install CRAN package 'wid' if missing.\n",
        "  • Offline: add\n    ", path_fallback, "\n",
        "    with columns ISO3, Year, value (bottom-50% pre-tax income share; 0–1 or 0–100).\n",
        "  • merge_full_data.R will still run; Income_share_bottom50_WID will be missing.",
        if (stub) "\n  • Stub wid_clean.csv was left unchanged." else ""
      ),
      call. = FALSE
    )
  } else {
    write_csv(out, path_out, na = "")
    message(
      "Saved ", path_out, " (", nrow(out), " rows; ",
      length(unique(out$ISO3)), " countries)."
    )
  }
}
