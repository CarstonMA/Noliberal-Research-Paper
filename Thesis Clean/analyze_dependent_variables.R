# Analyze Dependent Variables for Neoliberal Policy Impact Study
# Lists all outcome variables and proposes criteria for selection

library(tidyverse)

cat("=== Loading Merged Data ===\n")
df <- read_csv("wdi_merged_all_indices.csv", show_col_types = FALSE)

# Identify dependent variables (WDI outcomes, excluding identifiers and policy indices)
# Exclude: Country, ISO3, Year, heritage_*, fraser_*
dependent_vars <- names(df)[!str_detect(names(df), "^heritage_|^fraser_|^Country$|^ISO3$|^Year$|^Country Name$")]

cat("\n=== ALL DEPENDENT VARIABLES (", length(dependent_vars), " total) ===\n\n")

# Categorize variables by domain
categorize_variable <- function(var_name) {
  var_lower <- tolower(var_name)
  
  if(str_detect(var_lower, "gdp|gni|income|consumption|prosperity")) {
    return("Economic Growth/Income")
  } else if(str_detect(var_lower, "poverty|poor")) {
    return("Poverty")
  } else if(str_detect(var_lower, "gini|inequality|income share|lowest|highest")) {
    return("Inequality")
  } else if(str_detect(var_lower, "unemployment|employment|labor")) {
    return("Labor Market")
  } else if(str_detect(var_lower, "education|school|enrollment|attainment")) {
    return("Education")
  } else if(str_detect(var_lower, "health|expenditure|out-of-pocket")) {
    return("Health")
  } else if(str_detect(var_lower, "social|insurance|safety|benefit|coverage")) {
    return("Social Protection")
  } else if(str_detect(var_lower, "carbon|co2|emission|environment")) {
    return("Environment")
  } else if(str_detect(var_lower, "account|financial|bank")) {
    return("Financial Inclusion")
  } else if(str_detect(var_lower, "price|inflation|cpi")) {
    return("Inflation")
  } else if(str_detect(var_lower, "export|trade")) {
    return("Trade")
  } else {
    return("Other")
  }
}

# Create categorized list
var_categories <- tibble(
  Variable = dependent_vars,
  Category = map_chr(dependent_vars, categorize_variable)
) %>%
  arrange(Category, Variable)

# Print by category
cat("=== VARIABLES BY CATEGORY ===\n\n")
for(cat in unique(var_categories$Category)) {
  cat_vars <- var_categories %>% filter(Category == cat)
  cat(sprintf("%s (%d variables):\n", cat, nrow(cat_vars)))
  for(i in 1:nrow(cat_vars)) {
    cat(sprintf("  %d. %s\n", i, cat_vars$Variable[i]))
  }
  cat("\n")
}

# Analyze data availability for each variable
cat("=== DATA AVAILABILITY ANALYSIS ===\n\n")

var_availability <- var_categories %>%
  mutate(
    Non_Missing = map_int(Variable, ~sum(!is.na(df[[.x]]))),
    Total = nrow(df),
    Pct_Available = round(100 * Non_Missing / Total, 1),
    Countries = map_int(Variable, ~n_distinct(df$Country[!is.na(df[[.x]])])),
    Years = map_int(Variable, ~n_distinct(df$Year[!is.na(df[[.x]])]))
  ) %>%
  arrange(desc(Non_Missing))

cat("Top 20 variables by data availability:\n")
print(head(var_availability, 20) %>% 
      select(Variable, Category, Non_Missing, Pct_Available, Countries, Years))
cat("\n")

# Key outcome variables for neoliberal policy analysis
cat("=== RECOMMENDED KEY OUTCOME VARIABLES ===\n\n")

key_outcomes <- c(
  # Economic growth
  "GDP per capita (constant 2015 US$)",
  "GDP growth (annual %)",
  "GDP per capita growth (annual %)",
  
  # Inequality
  "Gini index",
  "Income share held by highest 10%",
  "Income share held by lowest 10%",
  
  # Poverty
  "Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population)",
  "Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population)",
  "Poverty headcount ratio at national poverty lines (% of population)",
  
  # Labor
  "Unemployment, total (% of total labor force) (modeled ILO estimate)",
  
  # Social protection
  "Coverage of social insurance programs (% of population)",
  "Coverage of social safety net programs (% of population)",
  
  # Health
  "Proportion of population spending more than 10% of household consumption or income on out-of-pocket health care expenditure (%)",
  
  # Education
  "School enrollment, primary (% net)",
  "Educational attainment, at least completed lower secondary, population 25+, total (%) (cumulative)"
)

cat("Recommended key outcomes for neoliberal policy analysis:\n")
for(i in seq_along(key_outcomes)) {
  if(key_outcomes[i] %in% dependent_vars) {
    avail <- var_availability %>% filter(Variable == key_outcomes[i])
    cat(sprintf("%d. %s\n", i, key_outcomes[i]))
    cat(sprintf("   Availability: %d obs (%.1f%%), %d countries, %d years\n\n", 
                avail$Non_Missing, avail$Pct_Available, avail$Countries, avail$Years))
  }
}

# Save full list
write_csv(var_availability, "dependent_variables_list.csv")
cat("Full variable list saved to: dependent_variables_list.csv\n\n")

# ============================================
# CRITERIA FOR NARROWING DOWN VARIABLES
# ============================================

cat("=== CRITERIA FOR SELECTING DEPENDENT VARIABLES ===\n\n")

cat("1. THEORETICAL RELEVANCE\n")
cat("   - Which outcomes does neoliberal theory predict will be affected?\n")
cat("   - Focus on variables directly related to your research question\n")
cat("   - Consider: growth, inequality, poverty, social protection\n\n")

cat("2. DATA QUALITY\n")
cat("   - Minimum threshold: At least 30-40% non-missing values\n")
cat("   - Coverage: Data for multiple countries and years\n")
cat("   - Temporal coverage: Variables with data spanning policy shift periods\n\n")

cat("3. MEASUREMENT QUALITY\n")
cat("   - Prefer standardized measures (e.g., PPP-adjusted poverty lines)\n")
cat("   - Avoid highly correlated redundant variables\n")
cat("   - Consider measurement error and reliability\n\n")

cat("4. POLICY RELEVANCE\n")
cat("   - Variables that policymakers care about\n")
cat("   - Outcomes that can be influenced by neoliberal policies\n")
cat("   - Variables with clear policy implications\n\n")

cat("5. STATISTICAL CONSIDERATIONS\n")
cat("   - Sufficient variation for regression analysis\n")
cat("   - Variables that can be log-transformed if needed\n")
cat("   - Consider distribution (normal vs. skewed)\n\n")

# Suggest a focused set
cat("=== SUGGESTED FOCUSED SET (10-15 variables) ===\n\n")

suggested_vars <- var_availability %>%
  filter(
    Pct_Available >= 20,  # At least 20% data availability
    Countries >= 50,      # At least 50 countries
    Years >= 10           # At least 10 years of data
  ) %>%
  filter(Category %in% c("Economic Growth/Income", "Inequality", "Poverty", 
                         "Labor Market", "Social Protection", "Health", "Education")) %>%
  arrange(Category, desc(Non_Missing)) %>%
  group_by(Category) %>%
  slice_head(n = 2) %>%  # Top 2 per category
  ungroup()

cat("Suggested focused set based on data quality and relevance:\n")
print(suggested_vars %>% select(Variable, Category, Pct_Available, Countries, Years))
cat("\n")

cat("=== NEXT STEPS ===\n\n")
cat("1. Review the categorized variable list above\n")
cat("2. Select 5-10 variables based on:\n")
cat("   - Your research question and theoretical framework\n")
cat("   - Data availability (use dependent_variables_list.csv)\n")
cat("   - Policy relevance\n")
cat("3. Check for multicollinearity among selected variables\n")
cat("4. Consider creating composite indices for related outcomes\n")
cat("5. Test robustness across different outcome measures\n\n")

cat("=== ANALYSIS COMPLETE ===\n")
