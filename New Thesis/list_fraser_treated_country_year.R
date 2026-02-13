# Output Fraser treated country-year pairs to CSV
library(tidyverse)
df <- read_csv("merged_full_data.csv", show_col_types = FALSE) %>%
  filter(Year >= 1995, Year <= 2022, !is.na(ISO3))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)
d <- df %>% filter(!is.na(fraser_Summary)) %>% arrange(ISO3, Year) %>%
  group_by(ISO3) %>%
  mutate(lag_p = dplyr::lag(fraser_Summary), jump = fraser_Summary - lag_p) %>%
  ungroup() %>% filter(!is.na(jump))
thresh <- quantile(d$jump, 0.90, na.rm = TRUE)
tr <- d %>% filter(jump >= thresh) %>% group_by(ISO3) %>% slice(1) %>% ungroup() %>%
  transmute(ISO3, event_year = Year)
out <- tr %>% left_join(distinct(df, ISO3, Country), by = "ISO3") %>%
  select(ISO3, Country, event_year) %>% arrange(Country)
write_csv(out, "fraser_treated_country_year.csv")
cat("Saved fraser_treated_country_year.csv\n")
print(out, n = Inf)
