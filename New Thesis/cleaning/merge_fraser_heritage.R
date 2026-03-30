# Merge Fraser Time Series and Heritage Foundation data
# Input:  data/raw/Fraser Time Series.csv, data/raw/Heritage.csv
# Output: data/merged_fraser_heritage.csv

library(tidyverse)
if (!requireNamespace("countrycode", quietly = TRUE)) install.packages("countrycode")

proj_root <- if (dir.exists("data/raw")) {
  getwd()
} else if (dir.exists("../data/raw")) {
  normalizePath("..")
} else {
  getwd()
}
path_fraser <- file.path(proj_root, "data", "raw", "Fraser Time Series.csv")
path_heritage <- file.path(proj_root, "data", "raw", "Heritage.csv")
path_out <- file.path(proj_root, "data", "merged_fraser_heritage.csv")

# 1. Load data
fraser <- read_csv(path_fraser, show_col_types = FALSE)
heritage <- read_csv(path_heritage, show_col_types = FALSE)

# 2. Harmonize Heritage: map name_web -> ISO3
# Manual mapping for hyphenated/special cases not in countrycode
name_to_iso <- c(
  "bahamas" = "BHS", "bosnia-herzegovina" = "BIH", "brunei-darussalam" = "BRN",
  "burkina-faso" = "BFA", "burma" = "MMR", "cabo-verde" = "CPV",
  "central-african-rep" = "CAF", "congo-dem-rep" = "COD", "congo-rep" = "COG",
  "costa-rica" = "CRI", "cote-divoire" = "CIV", "czech-rep" = "CZE",
  "dominican-rep" = "DOM", "el-salvador" = "SLV", "equatorial-guinea" = "GNQ",
  "eswatini" = "SWZ", "gambia" = "GMB", "guinea-bissau" = "GNB", "hong-kong" = "HKG",
  "kosovo" = "XKX",
  "korea-south" = "KOR", "korea-north" = "PRK", "kyrgyz-rep" = "KGZ",
  "laos" = "LAO", "macau" = "MAC", "micronesia" = "FSM",
  "new-zealand" = "NZL", "north-macedonia" = "MKD", "papua-new-guinea" = "PNG",
  "russia" = "RUS", "sao-tome-principe" = "STP", "saudi-arabia" = "SAU",
  "sierra-leone" = "SLE", "solomon-islands" = "SLB", "south-africa" = "ZAF",
  "sri-lanka" = "LKA", "st-lucia" = "LCA", "st-vincent-grenadines" = "VCT",
  "taiwan" = "TWN", "timor-leste" = "TLS", "trinidad-tobago" = "TTO",
  "turkey" = "TUR", "turkiye" = "TUR", "turkmenistan" = "TKM",
  "united-arab-emirates" = "ARE", "united-kingdom" = "GBR", "united-states" = "USA",
  "venezuela" = "VEN", "vietnam" = "VNM", "yemen" = "YEM"
)

# Map name_web to ISO3
heritage <- heritage %>%
  filter(name_web != "name_web") %>%  # drop header if duplicated
  mutate(
    name_proper = str_to_title(str_replace_all(name_web, "-", " ")),
    ISO3 = if_else(
      name_web %in% names(name_to_iso),
      unname(name_to_iso[name_web]),
      countrycode::countrycode(
        replace(name_proper, name_web %in% names(name_to_iso), NA_character_),
        "country.name", "iso3c"
      )
    )
  ) %>%
  select(-name_proper)

# 3. Prepare Fraser (rename columns with fraser_ prefix to avoid clashes)
fraser_clean <- fraser %>%
  select(
    ISO3 = ISO_Code,
    Year,
    Country = Countries,
    fraser_Summary = Summary,
    fraser_SizeGov = `Size of Government`,
    fraser_LegalSystem = `Legal System and Property Rights`,
    fraser_SoundMoney = `Sound Money`,
    fraser_FreeTrade = `Freedom to trade internationally`,
    fraser_Regulation = Regulation
  )

# 4. Prepare Heritage (prefix Heritage columns)
heritage_ind <- c("Overall", "Property Rights", "Government Integrity", "Judicial Effectiveness",
  "Tax Burden", "Government Spending", "Fiscal Health", "Business Freedom",
  "Labor Freedom", "Monetary Freedom", "Trade Freedom", "Investment Freedom", "Financial Freedom")

heritage_clean <- heritage %>%
  select(ISO3, Year, name_web, all_of(intersect(heritage_ind, names(heritage)))) %>%
  mutate(across(-c(ISO3, Year, name_web), ~ suppressWarnings(as.numeric(.x)))) %>%
  rename_with(~ paste0("heritage_", .x), .cols = -c(ISO3, Year, name_web))

# 5. Merge on ISO3 + Year (full join keeps all observations from both)
merged <- full_join(
  fraser_clean,
  heritage_clean,
  by = c("ISO3", "Year"),
  relationship = "many-to-one"
) %>%
  mutate(Country = coalesce(Country, str_to_title(str_replace_all(name_web, "-", " ")))) %>%
  select(-name_web)

# 6. Save
merged %>%
  relocate(ISO3, Country, Year) %>%
  arrange(ISO3, Year) %>%
  write_csv(path_out, na = "")

message("Saved ", path_out, ": ", nrow(merged), " rows, ", ncol(merged), " cols")
