# Map: Treated countries by regional/typological group with wave subdivisions
# Groups and waves follow refined narrative
# Requires: rnaturalearth, sf, tidyverse

library(tidyverse)
library(sf)
library(rnaturalearth)

# 1. Regional groups with wave subdivisions (from refined typology)
group_lookup <- tribble(
  ~ISO3, ~region, ~wave,
  # Post-Soviet & Eastern Bloc
  "ARM", "Post-Soviet & Eastern Bloc", "2001",
  "BIH", "Post-Soviet & Eastern Bloc", "2001",
  "BGR", "Post-Soviet & Eastern Bloc", "2001",
  "GEO", "Post-Soviet & Eastern Bloc", "2001",
  "HUN", "Post-Soviet & Eastern Bloc", "2001",
  "ROU", "Post-Soviet & Eastern Bloc", "2001",
  "SVN", "Post-Soviet & Eastern Bloc", "2001",
  "TJK", "Post-Soviet & Eastern Bloc", "2001",
  "UKR", "Post-Soviet & Eastern Bloc", "2001",
  "LVA", "Post-Soviet & Eastern Bloc", "2002-03",
  "LTU", "Post-Soviet & Eastern Bloc", "2002-03",
  "MDA", "Post-Soviet & Eastern Bloc", "2002-03",
  "MNE", "Post-Soviet & Eastern Bloc", "2002-03",
  "MKD", "Post-Soviet & Eastern Bloc", "2002-03",
  "POL", "Post-Soviet & Eastern Bloc", "2002-03",
  "RUS", "Post-Soviet & Eastern Bloc", "2002-03",
  "SRB", "Post-Soviet & Eastern Bloc", "2002-03",
  "AZE", "Post-Soviet & Eastern Bloc", "2004",
  "BLR", "Post-Soviet & Eastern Bloc", "2004",
  "KGZ", "Post-Soviet & Eastern Bloc", "2004",
  "SVK", "Post-Soviet & Eastern Bloc", "2004",
  # Latin America & Caribbean
  "BLZ", "Latin America & Caribbean", "2001-02",
  "BRA", "Latin America & Caribbean", "2001-02",
  "ECU", "Latin America & Caribbean", "2001-02",
  "HND", "Latin America & Caribbean", "2001-02",
  "NIC", "Latin America & Caribbean", "2001-02",
  "PER", "Latin America & Caribbean", "2001-02",
  "VEN", "Latin America & Caribbean", "2001-02",
  "ARG", "Latin America & Caribbean", "2004-05",
  "CRI", "Latin America & Caribbean", "2004-05",
  "DOM", "Latin America & Caribbean", "2004-05",
  "GUY", "Latin America & Caribbean", "2004-05",
  "URY", "Latin America & Caribbean", "2004-05",
  "BHS", "Latin America & Caribbean", "Recent",
  "BRB", "Latin America & Caribbean", "Recent",
  "HTI", "Latin America & Caribbean", "Recent",
  "PAN", "Latin America & Caribbean", "Recent",
  "PRY", "Latin America & Caribbean", "Recent",
  "SUR", "Latin America & Caribbean", "Recent",
  "TTO", "Latin America & Caribbean", "Recent",
  # Sub-Saharan Africa
  "BFA", "Sub-Saharan Africa", "2001",
  "BDI", "Sub-Saharan Africa", "2001",
  "CPV", "Sub-Saharan Africa", "2001",
  "GNB", "Sub-Saharan Africa", "2001",
  "MDG", "Sub-Saharan Africa", "2001",
  "MWI", "Sub-Saharan Africa", "2001",
  "MOZ", "Sub-Saharan Africa", "2001",
  "TCD", "Sub-Saharan Africa", "2002-04",
  "COD", "Sub-Saharan Africa", "2002-04",
  "ETH", "Sub-Saharan Africa", "2002-04",
  "GMB", "Sub-Saharan Africa", "2002-04",
  "GHA", "Sub-Saharan Africa", "2002-04",
  "KEN", "Sub-Saharan Africa", "2002-04",
  "LBR", "Sub-Saharan Africa", "2002-04",
  "NGA", "Sub-Saharan Africa", "2002-04",
  "RWA", "Sub-Saharan Africa", "2002-04",
  "SYC", "Sub-Saharan Africa", "2002-04",
  "SLE", "Sub-Saharan Africa", "2002-04",
  "AGO", "Sub-Saharan Africa", "2005-10",
  "BWA", "Sub-Saharan Africa", "2005-10",
  "CMR", "Sub-Saharan Africa", "2005-10",
  "CIV", "Sub-Saharan Africa", "2005-10",
  "GIN", "Sub-Saharan Africa", "2005-10",
  "MUS", "Sub-Saharan Africa", "2005-10",
  "SEN", "Sub-Saharan Africa", "2005-10",
  "ZAF", "Sub-Saharan Africa", "2005-10",
  "TZA", "Sub-Saharan Africa", "2005-10",
  "UGA", "Sub-Saharan Africa", "2005-10",
  "ZWE", "Sub-Saharan Africa", "2005-10",
  "BEN", "Sub-Saharan Africa", "Late",
  "COG", "Sub-Saharan Africa", "Late",
  "TGO", "Sub-Saharan Africa", "Late",
  "ZMB", "Sub-Saharan Africa", "Late",
  # MENA
  "IRN", "Middle East & North Africa", "2001-02",
  "KWT", "Middle East & North Africa", "2001-02",
  "LBN", "Middle East & North Africa", "2001-02",
  "QAT", "Middle East & North Africa", "2001-02",
  "SYR", "Middle East & North Africa", "2001-02",
  "ARE", "Middle East & North Africa", "2001-02",
  "YEM", "Middle East & North Africa", "2001-02",
  "DZA", "Middle East & North Africa", "2003-06",
  "EGY", "Middle East & North Africa", "2003-06",
  "IRQ", "Middle East & North Africa", "2003-06",
  "LBY", "Middle East & North Africa", "2003-06",
  "SAU", "Middle East & North Africa", "2003-06",
  "JOR", "Middle East & North Africa", "2009-10",
  "SDN", "Middle East & North Africa", "2009-10",
  # Asia & Pacific
  "IDN", "Asia & Pacific", "2001-02",
  "MNG", "Asia & Pacific", "2001-02",
  "NPL", "Asia & Pacific", "2001-02",
  "PAK", "Asia & Pacific", "2001-02",
  "TLS", "Asia & Pacific", "2001-02",
  "BGD", "Asia & Pacific", "2003-05",
  "BTN", "Asia & Pacific", "2003-05",
  "KHM", "Asia & Pacific", "2003-05",
  "FJI", "Asia & Pacific", "2003-05",
  "JPN", "Asia & Pacific", "2003-05",
  "LAO", "Asia & Pacific", "2003-05",
  "MMR", "Asia & Pacific", "2003-05",
  "PNG", "Asia & Pacific", "2003-05",
  "VNM", "Asia & Pacific", "2003-05",
  "IND", "Asia & Pacific", "Late",
  "MYS", "Asia & Pacific", "Late",
  "PHL", "Asia & Pacific", "Late",
  "LKA", "Asia & Pacific", "Late",
  # High-Income & Western
  "CYP", "High-Income & Western", "Early 2000s",
  "TUR", "High-Income & Western", "Early 2000s",
  "GRC", "High-Income & Western", "Crisis/Policy",
  "ISL", "High-Income & Western", "Crisis/Policy",
  "NOR", "High-Income & Western", "Crisis/Policy",
  "GBR", "High-Income & Western", "Crisis/Policy"
)

# 3. Load panel and treatment
if (!file.exists("merged_full_data.csv")) stop("Run merge_full_data.R first.")
df <- read_csv("merged_full_data.csv", show_col_types = FALSE) %>%
  filter(Year >= 1995, Year <= 2022, !is.na(ISO3))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)

identify_treatment <- function(data, pvar, pct = 0.90) {
  d <- data %>% filter(!is.na(.data[[pvar]])) %>% arrange(ISO3, Year) %>%
    group_by(ISO3) %>%
    mutate(lag_p = dplyr::lag(.data[[pvar]]), jump = .data[[pvar]] - lag_p) %>%
    ungroup() %>% filter(!is.na(jump))
  if (nrow(d) == 0) return(tibble(ISO3 = character(), event_year = numeric()))
  thresh <- quantile(d$jump, probs = pct, na.rm = TRUE)
  d %>% filter(jump >= thresh) %>% group_by(ISO3) %>% slice(1) %>% ungroup() %>%
    transmute(ISO3, event_year = as.integer(Year))
}
treatment_info <- identify_treatment(df %>% filter(ISO3 %in% panel_iso3), "fraser_Summary")

# 4. Status: control vs treated by region and wave
group_levels <- c(
  "Post-Soviet & Eastern Bloc (2001)", "Post-Soviet & Eastern Bloc (2002-03)", "Post-Soviet & Eastern Bloc (2004)",
  "Latin America & Caribbean (2001-02)", "Latin America & Caribbean (2004-05)", "Latin America & Caribbean (Recent)",
  "Sub-Saharan Africa (2001)", "Sub-Saharan Africa (2002-04)", "Sub-Saharan Africa (2005-10)", "Sub-Saharan Africa (Late)",
  "Middle East & North Africa (2001-02)", "Middle East & North Africa (2003-06)", "Middle East & North Africa (2009-10)",
  "Asia & Pacific (2001-02)", "Asia & Pacific (2003-05)", "Asia & Pacific (Late)",
  "High-Income & Western (Early 2000s)", "High-Income & Western (Crisis/Policy)",
  "Other treated", "Control", "Not in sample"
)

status <- tibble(ISO3 = panel_iso3) %>%
  left_join(treatment_info, by = "ISO3") %>%
  left_join(group_lookup, by = "ISO3") %>%
  mutate(
    group = case_when(
      is.na(event_year) ~ "Control",
      !is.na(region) ~ paste0(region, " (", wave, ")"),
      TRUE ~ "Other treated"
    ),
    group = factor(group, levels = group_levels)
  )

# 5. World map
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world %>%
  mutate(join_iso3 = coalesce(na_if(na_if(iso_a3, "-99"), ""), iso_a3_eh)) %>%
  filter(!is.na(join_iso3)) %>%
  left_join(status %>% select(ISO3, group) %>% rename(join_iso3 = ISO3), by = "join_iso3") %>%
  mutate(
    group = factor(
      replace_na(as.character(group), "Not in sample"),
      levels = group_levels
    )
  )

# 6. Plot - colors by region (shades for waves)
group_colors <- c(
  "Post-Soviet & Eastern Bloc (2001)" = "#08519c",
  "Post-Soviet & Eastern Bloc (2002-03)" = "#2171b5",
  "Post-Soviet & Eastern Bloc (2004)" = "#6baed6",
  "Latin America & Caribbean (2001-02)" = "#a50f15",
  "Latin America & Caribbean (2004-05)" = "#de2d26",
  "Latin America & Caribbean (Recent)" = "#fb6a4a",
  "Sub-Saharan Africa (2001)" = "#cc4c02",
  "Sub-Saharan Africa (2002-04)" = "#ec7014",
  "Sub-Saharan Africa (2005-10)" = "#fc8d59",
  "Sub-Saharan Africa (Late)" = "#fdbb84",
  "Middle East & North Africa (2001-02)" = "#016c59",
  "Middle East & North Africa (2003-06)" = "#3690c0",
  "Middle East & North Africa (2009-10)" = "#67a9cf",
  "Asia & Pacific (2001-02)" = "#006d2c",
  "Asia & Pacific (2003-05)" = "#31a354",
  "Asia & Pacific (Late)" = "#74c476",
  "High-Income & Western (Early 2000s)" = "#54278f",
  "High-Income & Western (Crisis/Policy)" = "#756bb1",
  "Other treated" = "#e6ab02",
  "Control" = "#deebf7",
  "Not in sample" = "gray90"
)

p <- ggplot(world) +
  geom_sf(aes(fill = group), color = "gray40", linewidth = 0.1) +
  scale_fill_manual(values = group_colors, name = "Region (Wave)", na.value = "gray90", drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    panel.grid = element_blank(),
    axis.text = element_blank()
  ) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +
  labs(
    title = "Treated Countries by Regional Group and Wave",
    subtitle = "Treatment: top-decile jump in Fraser Economic Freedom. Wave = timing of shock."
  )

dir.create("maps", showWarnings = FALSE)
ggsave("maps/treatment_map_by_region.png", p, width = 16, height = 9, dpi = 150)
cat("Saved maps/treatment_map_by_region.png\n")
