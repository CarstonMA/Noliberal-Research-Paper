# Merge all datasets: WDI, Heritage, and Fraser
# This script merges all neoliberal policy index data with WDI outcome data

library(tidyverse)
library(readxl)

# ============================================
# LOAD ALL DATASETS
# ============================================

cat("Loading datasets...\n")

# Load WDI data (outcomes)
wdi <- read_csv("wdi_tidy.csv")

# Load Heritage Foundation data (economic freedom index)
heritage <- read_csv("Heritage.csv")

# Load Fraser Institute data (economic freedom index)
fraser <- read_excel("Fraser.xlsx")

cat("WDI: ", nrow(wdi), " observations\n")
cat("Heritage: ", nrow(heritage), " observations\n")
cat("Fraser: ", nrow(fraser), " observations\n\n")

# ============================================
# EXPLORE DATA STRUCTURES
# ============================================

cat("=== WDI Structure ===\n")
cat("Columns:", paste(head(names(wdi), 5), collapse = ", "), "...\n")
cat("Key identifiers: Country Name, ISO3, Year\n\n")

cat("=== Heritage Structure ===\n")
cat("Columns:", paste(names(heritage), collapse = ", "), "\n")
cat("Key identifiers: name_web, Year\n\n")

cat("=== Fraser Structure ===\n")
cat("Columns:", paste(names(fraser), collapse = ", "), "\n")
cat("Key identifiers:", paste(names(fraser)[1:3], collapse = ", "), "\n\n")

# ============================================
# STANDARDIZE COUNTRY NAMES FOR MERGING
# ============================================

cat("Standardizing country names for merging...\n\n")

# Function to clean country names for matching
clean_country_name <- function(x) {
  x %>%
    tolower() %>%
    str_trim() %>%
    str_replace_all("[^a-z0-9]", "") %>%
    str_replace_all("the", "")
}

# Prepare WDI data
wdi_clean <- wdi %>%
  mutate(
    country_clean = clean_country_name(`Country Name`),
    country_original = `Country Name`
  )

# Prepare Heritage data
heritage_clean <- heritage %>%
  mutate(
    country_clean = clean_country_name(name_web),
    country_original = name_web
  ) %>%
  # Rename columns to avoid conflicts
  rename_with(~paste0("heritage_", .x), .cols = -c(name_web, Year, country_clean, country_original)) %>%
  rename(heritage_country = name_web)

# Prepare Fraser data
# First, identify the country and year columns
fraser_cols <- names(fraser)
cat("Fraser column names:\n")
print(fraser_cols)
cat("\nFraser data preview:\n")
print(head(fraser, 3))
cat("\n")

# Try to identify country and year columns
# Common patterns: "Country", "Year", "ISO", "Code", "Jurisdiction"
country_col <- fraser_cols[str_detect(tolower(fraser_cols), "country|name|jurisdiction|iso.name")]
year_col <- fraser_cols[str_detect(tolower(fraser_cols), "^year$|^yr$")]

# If not found, try common positions
if(length(country_col) == 0) {
  # Check first few columns
  if(any(str_detect(tolower(fraser_cols[1:3]), "country|name"))) {
    country_col <- fraser_cols[1:3][str_detect(tolower(fraser_cols[1:3]), "country|name")][1]
  } else {
    country_col <- fraser_cols[1]  # Default to first column
  }
}

if(length(year_col) == 0) {
  # Check for year in first few columns
  if(any(str_detect(tolower(fraser_cols[1:5]), "year"))) {
    year_col <- fraser_cols[1:5][str_detect(tolower(fraser_cols[1:5]), "year")][1]
  } else if("Year" %in% fraser_cols) {
    year_col <- "Year"
  } else {
    year_col <- fraser_cols[2]  # Default to second column
  }
}

cat("Using Fraser columns:\n")
cat("  Country:", country_col[1], "\n")
cat("  Year:", year_col[1], "\n\n")

# Get column names to exclude from renaming
exclude_cols <- c(country_col[1], year_col[1])

fraser_clean <- fraser %>%
  rename(
    fraser_country = !!sym(country_col[1]),
    Year = !!sym(year_col[1])
  ) %>%
  mutate(
    country_clean = clean_country_name(fraser_country),
    country_original = fraser_country,
    Year = as.numeric(Year)  # Ensure Year is numeric
  ) %>%
  # Rename all other columns with "fraser_" prefix (except country and year)
  rename_with(~paste0("fraser_", .x), .cols = -c(fraser_country, Year, country_clean, country_original))

# ============================================
# MERGE DATASETS
# ============================================

cat("Merging datasets...\n\n")

# Step 1: Merge WDI with Heritage
cat("Step 1: Merging WDI with Heritage...\n")
merged_data <- wdi_clean %>%
  left_join(
    heritage_clean %>%
      select(-country_original, -heritage_country),
    by = c("country_clean", "Year")
  )

cat("  After Heritage merge: ", nrow(merged_data), " observations\n")
cat("  Heritage matches: ", sum(!is.na(merged_data$heritage_Overall)), " observations\n\n")

# Step 2: Merge with Fraser
cat("Step 2: Merging with Fraser...\n")
merged_data <- merged_data %>%
  left_join(
    fraser_clean %>%
      select(-country_original, -fraser_country),
    by = c("country_clean", "Year")
  )

cat("  After Fraser merge: ", nrow(merged_data), " observations\n")
# Check for any Fraser column to count matches
fraser_match_col <- names(merged_data)[str_detect(names(merged_data), "^fraser_") & 
                                       !str_detect(names(merged_data), "country|Year")][1]
if(!is.na(fraser_match_col)) {
  cat("  Fraser matches: ", sum(!is.na(merged_data[[fraser_match_col]])), " observations\n")
} else {
  cat("  Fraser matches: Check manually - no Fraser columns found\n")
}
cat("\n")

# ============================================
# CHECK MERGE QUALITY
# ============================================

cat("=== Merge Quality Check ===\n\n")

# Check country matching
wdi_countries <- unique(wdi_clean$country_clean)
heritage_countries <- unique(heritage_clean$country_clean)
fraser_countries <- unique(fraser_clean$country_clean)

cat("Unique countries:\n")
cat("  WDI: ", length(wdi_countries), "\n")
cat("  Heritage: ", length(heritage_countries), "\n")
cat("  Fraser: ", length(fraser_countries), "\n")
cat("  Overlap (WDI & Heritage): ", length(intersect(wdi_countries, heritage_countries)), "\n")
cat("  Overlap (WDI & Fraser): ", length(intersect(wdi_countries, fraser_countries)), "\n")
cat("  Overlap (All three): ", length(intersect(intersect(wdi_countries, heritage_countries), fraser_countries)), "\n\n")

# Check year overlap
wdi_years <- unique(wdi_clean$Year)
heritage_years <- unique(heritage_clean$Year)
fraser_years <- unique(fraser_clean$Year)

cat("Year ranges:\n")
cat("  WDI: ", min(wdi_years), " to ", max(wdi_years), "\n")
cat("  Heritage: ", min(heritage_years), " to ", max(heritage_years), "\n")
cat("  Fraser: ", min(fraser_years), " to ", max(fraser_years), "\n")
cat("  Overlap (WDI & Heritage): ", length(intersect(wdi_years, heritage_years)), " years\n")
cat("  Overlap (WDI & Fraser): ", length(intersect(wdi_years, fraser_years)), " years\n\n")

# Check data availability by year
cat("Data availability by year (sample):\n")
# Find a Fraser column for matching
fraser_check_col <- names(merged_data)[str_detect(names(merged_data), "^fraser_") & 
                                      !str_detect(names(merged_data), "country|Year")][1]

if(!is.na(fraser_check_col)) {
  availability <- merged_data %>%
    group_by(Year) %>%
    summarise(
      wdi = sum(!is.na(`GDP per capita (constant 2015 US$)`)),
      heritage = sum(!is.na(heritage_Overall)),
      fraser = sum(!is.na(!!sym(fraser_check_col))),
      .groups = "drop"
    ) %>%
    filter(Year >= 1980 & Year <= 2020) %>%
    head(10)
} else {
  availability <- merged_data %>%
    group_by(Year) %>%
    summarise(
      wdi = sum(!is.na(`GDP per capita (constant 2015 US$)`)),
      heritage = sum(!is.na(heritage_Overall)),
      fraser = 0,
      .groups = "drop"
    ) %>%
    filter(Year >= 1980 & Year <= 2020) %>%
    head(10)
}

print(availability)
cat("\n")

# ============================================
# CLEAN UP AND SAVE
# ============================================

# Remove the temporary country_clean column (keep original country names)
final_data <- merged_data %>%
  select(-country_clean) %>%
  # Keep original country name from WDI
  rename(Country = country_original) %>%
  # Remove duplicate Country Name column (same as Country)
  select(-`Country Name`)

cat("=== Final Merged Dataset ===\n")
cat("Total observations: ", nrow(final_data), "\n")
cat("Total columns: ", ncol(final_data), "\n")
cat("Countries: ", n_distinct(final_data$Country), "\n")
cat("Years: ", min(final_data$Year, na.rm = TRUE), " to ", max(final_data$Year, na.rm = TRUE), "\n\n")

# Save merged dataset
write_csv(final_data, "wdi_merged_all_indices.csv")

cat("Merged dataset saved to: wdi_merged_all_indices.csv\n\n")

# ============================================
# SUMMARY OF AVAILABLE VARIABLES
# ============================================

cat("=== Available Variables ===\n\n")

cat("WDI Outcome Variables (sample):\n")
wdi_vars <- names(wdi)[4:min(10, ncol(wdi))]
cat(paste("  -", wdi_vars, collapse = "\n"))
cat("\n... and ", ncol(wdi) - 3, " more WDI variables\n\n")

cat("Heritage Index Variables:\n")
heritage_vars <- names(heritage_clean)[str_detect(names(heritage_clean), "^heritage_")]
cat(paste("  -", heritage_vars, collapse = "\n"))
cat("\n\n")

cat("Fraser Index Variables:\n")
fraser_vars <- names(fraser_clean)[str_detect(names(fraser_clean), "^fraser_")]
cat(paste("  -", head(fraser_vars, 10), collapse = "\n"))
if(length(fraser_vars) > 10) cat("\n  ... and ", length(fraser_vars) - 10, " more Fraser variables\n")
cat("\n")

cat("=== Merge Complete! ===\n")
cat("You can now use 'wdi_merged_all_indices.csv' for your analysis.\n")
cat("The dataset contains:\n")
cat("  - All WDI outcome variables\n")
cat("  - Heritage Foundation economic freedom indices\n")
cat("  - Fraser Institute economic freedom indices\n")
cat("  - All matched by country and year\n")
