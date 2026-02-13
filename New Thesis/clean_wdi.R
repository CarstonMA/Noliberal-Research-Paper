# Clean WDI data: pivot wide-year format to long, remove metadata, exclude aggregates
# Output: wdi_clean.csv

library(tidyverse)

# 1. Read raw WDI (n_max avoids parsing metadata section with different structure)
raw <- read_csv("WDI_Data.csv", show_col_types = FALSE, n_max = 2500)

# Drop metadata rows
wdi <- raw %>%
  filter(!is.na(`Series Name`), `Series Name` != "", `Series Name` != "Code",
         !str_detect(`Series Name`, "^Data from database|^Last Updated"))

# 2. Exclude regional/income aggregates (keep country-level only)
agg_codes <- c("WLD", "EAS", "ECS", "LCN", "LAC", "SAS", "SSA", "ARB", "CEB",
               "EAR", "EAP", "EMU", "EUU", "FCS", "HPC", "IBD", "IDB", "IDX",
               "IDA", "LDC", "LMY", "LIC", "LMC", "MEA", "MNA", "MIC", "NAC",
               "INX", "OED", "OSS", "PSS", "PST", "PRE", "SST", "TSS", "UMC",
               "HIC")
wdi <- wdi %>% filter(!`Country Code` %in% agg_codes)

# 3. Pivot longer: year columns -> Year
wdi_long <- wdi %>%
  pivot_longer(
    cols = matches("^\\d{4} \\[YR"),
    names_to = "Year",
    values_to = "Value"
  ) %>%
  mutate(
    Year = as.integer(str_extract(Year, "\\d{4}")),
    Value = na_if(Value, ".."),
    Value = as.numeric(Value)
  ) %>%
  select(`Series Name`, `Country Name`, `Country Code`, Year, Value)

# 4. Pivot wider: one column per indicator (use Series Name for readability)
wdi_clean <- wdi_long %>%
  pivot_wider(
    names_from = `Series Name`,
    values_from = Value
  ) %>%
  rename(ISO3 = `Country Code`) %>%
  relocate(ISO3, `Country Name`, Year)

# 5. Save
write_csv(wdi_clean, "wdi_clean.csv", na = "")

message("Saved wdi_clean.csv: ", nrow(wdi_clean), " rows, ", ncol(wdi_clean), " cols")
