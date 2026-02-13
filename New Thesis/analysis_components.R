# Component-level impact analysis: Fraser & Heritage indices
# Which components drive outcomes? Graphs significant effects.
# Uses fixed-effects panel: outcome ~ component | ISO3 + Year, cluster(ISO3)

library(tidyverse)
library(fixest)
library(ggplot2)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
if (!file.exists("merged_full_data.csv")) source("merge_full_data.R")
df <- read_csv("merged_full_data.csv", show_col_types = FALSE)
df <- df %>% filter(Year >= 1995, Year <= 2022, !is.na(ISO3), !is.na(Year))

# Panel: min 15 observations
n_years <- length(unique(df$Year))
panel_iso3 <- df %>% count(ISO3) %>% filter(n >= 15) %>% pull(ISO3)
df <- df %>% filter(ISO3 %in% panel_iso3)
cat("Panel:", length(panel_iso3), "countries (min 15 obs)\n")

# Region definitions (for optional robustness: ISO3 + region_year)
post_soviet <- c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", "LTU", "MDA", "RUS", "TJK", "TKM", "UKR", "UZB")
latin_america <- c("ARG", "BOL", "BRA", "BRB", "BLZ", "CHL", "COL", "CRI", "CUB", "DOM", "ECU", "SLV", "GTM", "GUY", "HND", "JAM", "MEX", "NIC", "PAN", "PRY", "PER", "TTO", "URY", "VEN")
sub_saharan_africa <- c("AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CAF", "TCD", "COM", "COG", "COD", "CIV", "GNQ", "ERI", "SWZ", "ETH", "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TZA", "TGO", "UGA", "ZMB", "ZWE")
east_asia_pacific <- c("AUS", "BGD", "BTN", "KHM", "CHN", "FJI", "HKG", "IDN", "JPN", "KIR", "LAO", "MYS", "MMR", "MNG", "NPL", "NZL", "PHL", "PNG", "KOR", "WSM", "SGP", "SLB", "LKA", "THA", "TON", "VUT", "VNM")
middle_east_na <- c("DZA", "BHR", "EGY", "IRN", "IRQ", "ISR", "JOR", "KWT", "LBN", "LBY", "MAR", "OMN", "QAT", "SAU", "SYR", "TUN", "ARE", "YEM", "PSE")
eastern_europe <- c("ALB", "BIH", "BGR", "HRV", "CZE", "HUN", "MKD", "POL", "ROU", "SRB", "SVK", "SVN")
western_europe_na <- c("AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "DEU", "GRC", "ISL", "IRL", "ITA", "LUX", "MLT", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "GBR", "USA")
df <- df %>% mutate(
  country_group = case_when(
    ISO3 %in% post_soviet ~ "Post-Soviet", ISO3 %in% latin_america ~ "Latin America",
    ISO3 %in% sub_saharan_africa ~ "Sub-Saharan Africa", ISO3 %in% east_asia_pacific ~ "East Asia & Pacific",
    ISO3 %in% middle_east_na ~ "Middle East & N. Africa", ISO3 %in% eastern_europe ~ "Eastern Europe",
    ISO3 %in% western_europe_na ~ "Western Europe & N. America", TRUE ~ "Other"
  ),
  region_year = paste0(country_group, "_", Year)
)

# ============================================================================
# 2. DEFINE COMPONENTS & OUTCOMES
# ============================================================================
fraser_components <- c("fraser_SizeGov", "fraser_LegalSystem", "fraser_SoundMoney",
                      "fraser_FreeTrade", "fraser_Regulation", "fraser_Summary")

heritage_components <- c("heritage_Property.Rights", "heritage_Government.Integrity",
                        "heritage_Judicial.Effectiveness", "heritage_Tax.Burden",
                        "heritage_Government.Spending", "heritage_Fiscal.Health",
                        "heritage_Business.Freedom", "heritage_Labor.Freedom",
                        "heritage_Monetary.Freedom", "heritage_Trade.Freedom",
                        "heritage_Investment.Freedom", "heritage_Financial.Freedom",
                        "heritage_Overall")

outcomes <- c()
if ("log_GDP_per_capita" %in% names(df)) outcomes <- c(outcomes, "log_GDP_per_capita")
if ("GDP_per_capita_growth" %in% names(df)) outcomes <- c(outcomes, "GDP_per_capita_growth")
gini <- names(df)[grepl("Gini|GINI", names(df), ignore.case = TRUE)]
if (length(gini) > 0) outcomes <- c(outcomes, gini[1])
pov <- names(df)[grepl("Poverty|SI.POV.UMIC", names(df), ignore.case = TRUE)]
if (length(pov) > 0) outcomes <- c(outcomes, pov[1])
inc <- names(df)[grepl("Income.share|SI.DST.FRST", names(df), ignore.case = TRUE)]
if (length(inc) > 0) outcomes <- c(outcomes, inc[1])
vuln <- names(df)[grepl("Vulnerable.employment|SL.EMP.VULN", names(df), ignore.case = TRUE)]
if (length(vuln) > 0) outcomes <- c(outcomes, vuln[1])

# Winsorize outcomes at 1% and 99%
winsorize <- function(x, p_low = 0.01, p_high = 0.99) {
  q <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
df <- df %>% mutate(across(all_of(outcomes), winsorize, .names = "{.col}"))

# ============================================================================
# 3. RUN FE REGRESSIONS (ISO3 + Year; use region_year for robustness check)
# ============================================================================
run_fe <- function(data, outcome, component) {
  d <- data %>% filter(!is.na(.data[[outcome]]) & !is.na(.data[[component]]))
  if (nrow(d) < 200) return(NULL)
  out_esc <- if (grepl("[^a-zA-Z0-9._]", outcome)) paste0("`", outcome, "`") else outcome
  comp_esc <- if (grepl("[^a-zA-Z0-9._]", component)) paste0("`", component, "`") else component
  f <- as.formula(paste0(out_esc, " ~ ", comp_esc, " | ISO3 + Year"))
  tryCatch(
    feols(f, data = d, cluster = ~ ISO3),
    error = function(e) NULL
  )
}

all_components <- c(fraser_components, heritage_components)
results <- list()

for (comp in all_components) {
  if (comp %in% names(df)) {
    for (out in outcomes) {
      m <- run_fe(df, out, comp)
      if (!is.null(m)) {
        b <- coef(m)[1]
        se_b <- se(m)[1]
        ct <- coeftable(m)
        pcol <- which(colnames(ct) %in% c("Pr(>|t|)", "pvalue", "Pr(>|z|)"))[1]
        pval <- if (nrow(ct) >= 1 && length(pcol) > 0) as.numeric(ct[1, pcol]) else NA_real_
        results[[length(results) + 1]] <- tibble(
          index = if_else(grepl("^fraser_", comp), "Fraser", "Heritage"),
          component = comp,
          outcome = out,
          coef = b,
          se = se_b,
          pvalue = pval,
          significant = pval < 0.05
        )
      }
    }
  }
}

res_df <- bind_rows(results)

# ============================================================================
# 4. SUMMARY: IMPACTFUL COMPONENTS
# ============================================================================
cat("\n=== COMPONENT IMPACT ANALYSIS (FE panel, SE clustered by country) ===\n")
cat("Coefficient = change in outcome per 1-unit increase in component score\n\n")

# Significant results only
sig <- res_df %>% filter(significant) %>% arrange(index, outcome, pvalue)
if (nrow(sig) > 0) {
  cat("SIGNIFICANT EFFECTS (p < 0.05):\n")
  print(sig %>% select(index, component, outcome, coef, se, pvalue))
  write_csv(sig, "component_impacts_significant.csv")
  cat("\nSaved: component_impacts_significant.csv\n\n")
}

# Full results
write_csv(res_df, "component_impacts_full.csv")
cat("Saved: component_impacts_full.csv (all component-outcome pairs)\n\n")

# Compact summary: components with most impactful effects
impact_count <- res_df %>% filter(significant) %>%
  count(index, component, name = "n_sig") %>%
  arrange(desc(n_sig))
if (nrow(impact_count) > 0) {
  cat("COMPONENTS BY NUMBER OF SIGNIFICANT OUTCOME EFFECTS:\n")
  print(impact_count)
}

cat("\nDone.\n")
