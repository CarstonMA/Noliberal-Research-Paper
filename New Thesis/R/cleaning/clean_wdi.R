# Download WDI from the World Bank API and build the country-year panel
# Output: data/wdi_clean.csv (same schema as before: ISO3, Country Name, Year, indicators)

if (!requireNamespace("WDI", quietly = TRUE)) {
  install.packages("WDI", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo", repos = "https://cloud.r-project.org")
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
path_out <- file.path(proj_root, "data", "wdi_clean.csv")

# Years: align with policy merge (merge_full_data.R uses 1995–2022). Pulling 1995–latest
# avoids a long 1960s–80s stretch where most series are blank in the API.
wdi_year_start <- 1995L
wdi_year_end <- as.integer(format(Sys.Date(), "%Y"))

# Same regional/income aggregates as the old CSV cleaner (exclude from panel)
agg_codes <- c(
  "WLD", "EAS", "ECS", "LCN", "LAC", "SAS", "SSA", "ARB", "CEB",
  "EAR", "EAP", "EMU", "EUU", "FCS", "HPC", "IBD", "IDB", "IDX",
  "IDA", "LDC", "LMY", "LIC", "LMC", "MEA", "MNA", "MIC", "NAC",
  "INX", "OED", "OSS", "PSS", "PST", "PRE", "SST", "TSS", "UMC",
  "HIC", "TSA", "SSF"
)

# ISO3 codes that are actual countries/economies in WDI (excludes regional & income aggregates)
iso3_members <- WDI::WDI_data$country %>%
  as_tibble() %>%
  filter(
    !is.na(iso3c),
    nchar(iso3c) == 3L,
    !is.na(region),
    region != "Aggregates"
  ) %>%
  pull(iso3c) %>%
  unique()

# Indicators used downstream (GDP, poverty/inequality, vulnerable employment)
wdi_indicators <- c(
  "NY.GDP.PCAP.PP.KD",  # GDP per capita
  "NY.GDP.MKTP.KD.ZG",  # GDP growth %
  "SI.POV.DDAY",        # Extreme poverty $2.15/day
  "SI.POV.LMIC",        # Lower-mid poverty $3.65/day
  "SI.POV.UMIC",        # Upper-mid poverty $8.30/day
  "SI.POV.UMIC.GP",     # Poverty gap at $8.30
  "SI.POV.GINI",        # Gini index
  "SI.DST.FRST.20",     # Income share lowest 20%
  "SL.EMP.VULN.ZS"      # Vulnerable employment
)

# Survey-based poverty / inequality only: linear interpolation inside gaps of at most
# `interp_max_gap` consecutive missing years, between observed points only (no extrapolation).
# GDP and employment are left as-is (not imputed from other variables).
interp_survey_vars <- c(
  "SI.POV.DDAY",
  "SI.POV.LMIC",
  "SI.POV.UMIC",
  "SI.POV.UMIC.GP",
  "SI.POV.GINI",
  "SI.DST.FRST.20"
)
interp_max_gap <- 2L

interpolate_survey_column <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) {
    return(x)
  }
  # rule = 1: NA outside first/last observed value (no extrapolation)
  zoo::na.approx(x, na.rm = FALSE, maxgap = interp_max_gap, rule = 1)
}

message("Fetching WDI from World Bank API (this may take a minute)...")
# WDI::WDI uses https://api.worldbank.org (documented JSON API)
raw <- WDI::WDI(
  country = "all",
  indicator = wdi_indicators,
  start = wdi_year_start,
  end = wdi_year_end,
  extra = FALSE
)

stopifnot(all(c("iso3c", "country", "year") %in% names(raw)))

# Drop income groups (empty iso3c), known aggregate codes, and non-country economies
wdi_clean <- raw %>%
  mutate(iso3c = na_if(trimws(as.character(iso3c)), "")) %>%
  filter(
    !is.na(iso3c),
    nchar(iso3c) == 3L,
    !iso3c %in% agg_codes,
    iso3c %in% iso3_members
  ) %>%
  rename(ISO3 = iso3c, `Country Name` = country, Year = year) %>%
  select(-any_of(c("iso2c"))) %>%
  arrange(ISO3, Year) %>%
  relocate(ISO3, `Country Name`, Year)

interp_present <- intersect(interp_survey_vars, names(wdi_clean))
if (length(interp_present) > 0) {
  n_before <- sum(!is.na(as.matrix(wdi_clean[interp_present])))
  wdi_clean <- wdi_clean %>%
    group_by(ISO3) %>%
    mutate(across(all_of(interp_present), interpolate_survey_column)) %>%
    ungroup()
  n_after <- sum(!is.na(as.matrix(wdi_clean[interp_present])))
  message(
    "Constrained linear interpolation (max gap ", interp_max_gap,
    " yrs, no extrapolation) on: ",
    paste(interp_present, collapse = ", "),
    "\n  Non-missing cells in these columns: ", n_before, " -> ", n_after
  )
}

write_csv(wdi_clean, path_out, na = "")

ind_cols <- setdiff(names(wdi_clean), c("ISO3", "Country Name", "Year"))
nmiss <- map_int(wdi_clean[ind_cols], ~ sum(!is.na(.x)))
message(
  "Saved ", path_out, ": ",
  nrow(wdi_clean), " rows, ", ncol(wdi_clean), " cols (source: World Bank API via WDI package)"
)
message(
  "Non-missing counts by indicator (sparse series are normal in WDI; poverty/Gini are not annual):\n",
  paste(sprintf("  %s: %s", names(nmiss), nmiss), collapse = "\n")
)
