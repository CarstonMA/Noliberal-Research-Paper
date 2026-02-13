# Map: Countries colored by treatment year (event study)
# Replicates treatment logic from analysis_event_study.R
# Requires: rnaturalearth, sf, tidyverse

library(tidyverse)
library(sf)
library(rnaturalearth)

# 1. Load data and identify treatment (same logic as event study)
if (!file.exists("merged_full_data.csv")) stop("Run merge_full_data.R first.")
df <- read_csv("merged_full_data.csv", show_col_types = FALSE) %>%
  filter(Year >= 1995, Year <= 2022, !is.na(ISO3))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)

policy_var <- if (sum(!is.na(df$fraser_Summary)) >= 1000) "fraser_Summary" else "heritage_Overall"
identify_treatment <- function(data, pvar, pct = 0.90) {
  d <- data %>% filter(!is.na(.data[[pvar]])) %>% arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(lag_p = dplyr::lag(.data[[pvar]]), jump = .data[[pvar]] - lag_p) %>%
    ungroup() %>% filter(!is.na(jump))
  if (nrow(d) == 0) return(tibble(ISO3 = character(), event_year = numeric()))
  jump_threshold <- quantile(d$jump, probs = pct, na.rm = TRUE)
  d %>% filter(jump >= jump_threshold) %>%
    group_by(ISO3) %>% slice(1) %>% ungroup() %>%
    transmute(ISO3, event_year = as.integer(Year))
}

treatment_info <- identify_treatment(df, policy_var)

# 2. Build status: control vs treated by year (binned into 5-year periods)
status <- tibble(ISO3 = panel_iso3) %>%
  left_join(treatment_info, by = "ISO3") %>%
  mutate(
    period = case_when(
      is.na(event_year) ~ "Control",
      event_year < 2000 ~ "1996–1999",
      event_year < 2005 ~ "2000–2004",
      event_year < 2010 ~ "2005–2009",
      event_year < 2015 ~ "2010–2014",
      TRUE ~ "2015+"
    ),
    group = factor(period, levels = c("Control", "1996–1999", "2000–2004", "2005–2009", "2010–2014", "2015+"))
  )

# 4. World map (Natural Earth)
# Note: France (and sometimes Norway) have missing iso_a3 in Natural Earth; use iso_a3_eh as fallback
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world %>%
  mutate(join_iso3 = coalesce(
    na_if(na_if(iso_a3, "-99"), ""),
    iso_a3_eh
  )) %>%
  filter(!is.na(join_iso3)) %>%
  left_join(
    status %>% rename(join_iso3 = ISO3),
    by = "join_iso3"
  ) %>%
  mutate(
    group = factor(
      replace_na(as.character(group), "Not in sample"),
      levels = c("Control", "1996–1999", "2000–2004", "2005–2009", "2010–2014", "2015+", "Not in sample")
    )
  )

# 5. Plot
p <- ggplot(world) +
  geom_sf(aes(fill = group), color = "gray40", linewidth = 0.1) +
  scale_fill_manual(
    values = c(
      "Control" = "#deebf7",
      "1996–1999" = "#2166ac",
      "2000–2004" = "#4393c3",
      "2005–2009" = "#92c5de",
      "2010–2014" = "#67a9cf",
      "2015+" = "#0571b0",
      "Not in sample" = "gray90"
    ),
    name = "Treatment year",
    na.value = "gray90",
    drop = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  labs(
    title = "Event Study: Countries by Treatment Year",
    subtitle = paste0("Treatment: top-decile jump in ", policy_var, ". ",
                     nrow(treatment_info), " treated, ", sum(is.na(status$event_year)), " control.")
  )

dir.create("maps", showWarnings = FALSE)
ggsave("maps/treatment_map.png", p, width = 12, height = 7, dpi = 150)
cat("Saved maps/treatment_map.png\n")
