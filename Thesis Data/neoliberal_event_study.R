# Neoliberal Policy Event Study Analysis
# Difference-in-Differences Event Study Specification
# Examining effects of neoliberal policy changes on aggregate and distributional outcomes

# ============================================================================
# 1. DATA PREPARATION
# ============================================================================

# Load required packages
library(tidyverse)
library(fixest)
library(ggplot2)
library(data.table)
library(car)  # For linearHypothesis function

# Try to load modelsummary, but make it optional
modelsummary_available <- require(modelsummary, quietly = TRUE, warn.conflicts = FALSE)
if (!modelsummary_available) {
  cat("Warning: modelsummary package not available. Regression tables will be simplified.\n")
  cat("To install: update.packages('data.table') then install.packages('modelsummary')\n")
}

# Load data
data <- read.csv("merged_data.csv", stringsAsFactors = FALSE)

# Convert to data.table for efficiency
data <- as.data.table(data)

# Clean column names (remove special characters that might cause issues)
# First, let's see what columns we actually have
cat("Checking column names...\n")

# Define mapping of old to new names
col_mapping <- list(
  "Account ownership at a financial institution or with a mobile-money-service provider, poorest 40% (% of population ages 15+)" = "account_ownership_poorest40",
  "Annualized average growth rate in per capita real survey mean consumption or income, bottom 40% of population (%)" = "bottom40_growth",
  "GDP (constant 2015 US$)" = "GDP",
  "GDP growth (annual %)" = "GDP_growth",
  "GDP per capita (constant 2015 US$)" = "GDP_per_capita",
  "GDP per capita growth (annual %)" = "GDP_per_capita_growth",
  "GNI (constant 2015 US$)" = "GNI",
  "GNI growth (annual %)" = "GNI_growth",
  "GNI per capita (constant 2015 US$)" = "GNI_per_capita",
  "GNI per capita growth (annual %)" = "GNI_per_capita_growth",
  "Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population)" = "poverty_365",
  "Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population)" = "poverty_685"
)

# Try to rename columns that exist
renamed_count <- 0
for (old_name in names(col_mapping)) {
  if (old_name %in% names(data)) {
    setnames(data, old = old_name, new = col_mapping[[old_name]], skip_absent = TRUE)
    renamed_count <- renamed_count + 1
  }
}

if (renamed_count > 0) {
  cat("Renamed", renamed_count, "columns\n")
} else {
  cat("Note: Column renaming skipped - columns may have different names. Using original names.\n")
  # Create aliases for easier access later
  # We'll reference columns by their original names in the analysis
}

# Ensure Year is numeric
data[, Year := as.numeric(Year)]

# Sort by country and year
setorder(data, ISO3, Year)

# Create World Bank income classifications based on GDP per capita
# Using 2020 thresholds (will be time-varying)
# Low: < $1,045
# Lower-middle: $1,045 - $4,095
# Upper-middle: $4,095 - $12,695
# High: > $12,695

# Find GDP per capita column 
# Also search for partial matches
gdp_pc_col <- NULL
if ("GDP_per_capita" %in% names(data)) {
  gdp_pc_col <- "GDP_per_capita"
} else if ("GDP per capita (constant 2015 US$)" %in% names(data)) {
  gdp_pc_col <- "GDP per capita (constant 2015 US$)"
} else {
  # Try to find any column with "GDP" and "capita" or "per capita"
  gdp_cols <- names(data)[grepl("GDP.*capita|GDP.*per.*capita", names(data), ignore.case = TRUE)]
  if (length(gdp_cols) > 0) {
    gdp_pc_col <- gdp_cols[1]
    cat("Found GDP per capita column:", gdp_pc_col, "\n")
  }
}

if (!is.null(gdp_pc_col)) {
  data[, income_group := cut(get(gdp_pc_col),
                             breaks = c(-Inf, 1045, 4095, 12695, Inf),
                             labels = c("Low", "Lower-middle", "Upper-middle", "High"),
                             include.lowest = TRUE)]
} else {
  cat("Warning: GDP per capita column not found. Income groups cannot be created.\n")
  data[, income_group := NA_character_]
}

# Create binary indicators for income groups
data[, low_income := as.numeric(income_group == "Low")]
data[, middle_income := as.numeric(income_group %in% c("Lower-middle", "Upper-middle"))]
data[, high_income := as.numeric(income_group == "High")]

# ============================================================================
# 1.5. COUNTRY GROUPING BY REGION/HISTORICAL CHARACTERISTICS
# ============================================================================

# Define country groups based on geographic and historical characteristics
cat("\n=== CREATING COUNTRY GROUPS ===\n")

# Post-Soviet countries (former USSR)
post_soviet <- c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", "LTU", 
                 "MDA", "RUS", "TJK", "TKM", "UKR", "UZB")

# Latin America and Caribbean
latin_america <- c("ARG", "BOL", "BRA", "BRB", "BLZ", "CHL", "COL", "CRI", "CUB", 
                   "DOM", "ECU", "SLV", "GTM", "GUY", "HND", "JAM", "MEX", "NIC", 
                   "PAN", "PRY", "PER", "TTO", "URY", "VEN")

# Sub-Saharan Africa
sub_saharan_africa <- c("AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CAF", "TCD", 
                        "COM", "COG", "COD", "CIV", "GNQ", "ERI", "SWZ", "ETH", 
                        "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", "LBR", 
                        "MDG", "MWI", "MLI", "MRT", "MUS", "MOZ", "NAM", "NER", 
                        "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", 
                        "SSD", "SDN", "TZA", "TGO", "UGA", "ZMB", "ZWE")

# East Asia and Pacific
east_asia_pacific <- c("AUS", "BGD", "BTN", "KHM", "CHN", "FJI", "HKG", "IDN", 
                       "JPN", "KIR", "LAO", "MYS", "MMR", "MNG", "NPL", "NZL", 
                       "PHL", "PNG", "KOR", "WSM", "SGP", "SLB", "LKA", "THA", 
                       "TON", "VUT", "VNM")

# Middle East and North Africa
middle_east_na <- c("DZA", "BHR", "EGY", "IRN", "IRQ", "ISR", "JOR", "KWT", 
                    "LBN", "LBY", "MAR", "OMN", "QAT", "SAU", "SYR", "TUN", 
                    "ARE", "YEM", "PSE")

# Eastern Europe (non-Soviet)
eastern_europe <- c("ALB", "BIH", "BGR", "HRV", "CZE", "HUN", "MKD", "POL", 
                    "ROU", "SRB", "SVK", "SVN")

# Western Europe and North America
western_europe_na <- c("AUT", "BEL", "CAN", "DNK", "FIN", "FRA", "DEU", "GRC", 
                       "ISL", "IRL", "ITA", "LUX", "MLT", "NLD", "NOR", "PRT", 
                       "ESP", "SWE", "CHE", "GBR", "USA")

# South Asia
south_asia <- c("AFG", "BGD", "BTN", "IND", "MDV", "NPL", "PAK", "LKA")

# Create country group variable
data[, country_group := "Other"]
data[ISO3 %in% post_soviet, country_group := "Post-Soviet"]
data[ISO3 %in% latin_america, country_group := "Latin America"]
data[ISO3 %in% sub_saharan_africa, country_group := "Sub-Saharan Africa"]
data[ISO3 %in% east_asia_pacific, country_group := "East Asia & Pacific"]
data[ISO3 %in% middle_east_na, country_group := "Middle East & N. Africa"]
data[ISO3 %in% eastern_europe, country_group := "Eastern Europe"]
data[ISO3 %in% western_europe_na, country_group := "Western Europe & N. America"]
data[ISO3 %in% south_asia, country_group := "South Asia"]

# Print group counts
cat("Country groups created:\n")
group_counts <- data[, .N, by = country_group][order(-N)]
print(group_counts)
cat("\n")

# Log transform GDP variables for better interpretation
# Find GDP columns (reuse gdp_pc_col from above if found)
# For total GDP
gdp_col <- NULL
if ("GDP" %in% names(data)) {
  gdp_col <- "GDP"
} else if ("GDP (constant 2015 US$)" %in% names(data)) {
  gdp_col <- "GDP (constant 2015 US$)"
} else {
  # Try to find any column with "GDP" but not "capita" or "per"
  gdp_total_cols <- names(data)[grepl("^GDP[^p]|GDP.*constant.*US", names(data), ignore.case = TRUE) & 
                                 !grepl("capita|per", names(data), ignore.case = TRUE)]
  if (length(gdp_total_cols) > 0) {
    gdp_col <- gdp_total_cols[1]
    cat("Found GDP column:", gdp_col, "\n")
  }
}

if (!is.null(gdp_pc_col)) {
  data[, log_GDP_per_capita := log(get(gdp_pc_col) + 1)]
} else {
  cat("Warning: GDP per capita column not found for log transformation.\n")
  data[, log_GDP_per_capita := NA_real_]
}

if (!is.null(gdp_col)) {
  data[, log_GDP := log(get(gdp_col) + 1)]
} else {
  cat("Warning: GDP column not found for log transformation.\n")
  data[, log_GDP := NA_real_]
}

# ============================================================================
# 2. TREATMENT DEFINITION
# ============================================================================

# Define neoliberal policy indicators to analyze
# FOCUS ON SPECIFIC POLICY COMPONENTS instead of aggregate scores
# These individual policies may have clearer, more consistent effects
policy_vars <- c("Property Rights", "Trade Freedom", "Investment Freedom",
                 "Business Freedom", "Financial Freedom", "Monetary Freedom")

# Also keep aggregate scores for comparison, but prioritize components
aggregate_scores <- c("Heritage_Score", "fraser_Score")

# Function to identify treatment events (significant policy score increases)
# RELAXED TREATMENT DEFINITION to get more countries:
# 1. Moderate jump: Increase of at least min_jump_sd standard deviations (default 0.5)
# 2. Optional stability requirement: Can be turned off
# 3. Optional sustained increase: Can be turned off
# 4. Threshold: Must cross from bottom 50% to top 50% (using median)
identify_treatment <- function(data, policy_var, country_var = "ISO3", year_var = "Year",
                               min_jump_sd = 0.5, require_stability = FALSE, require_sustained = FALSE) {
  
  # Check if policy variable exists
  if (!policy_var %in% names(data)) {
    cat("Warning: Policy variable", policy_var, "not found in data\n")
    return(data.table(ISO3 = character(), treated = logical(), event_year = numeric()))
  }
  
  # Get non-missing values for the policy variable
  policy_data <- data[!is.na(get(policy_var)) & get(policy_var) != "", 
                      .(ISO3 = get(country_var), 
                        Year = as.numeric(get(year_var)), 
                        policy = as.numeric(get(policy_var)))]
  
  # Remove any remaining NAs
  policy_data <- policy_data[!is.na(policy) & !is.na(Year)]
  
  if (nrow(policy_data) == 0) {
    cat("Warning: No valid data for policy variable", policy_var, "\n")
    return(data.table(ISO3 = character(), treated = logical(), event_year = numeric()))
  }
  
  # Calculate thresholds: use median (50th percentile) for relaxed criteria
  # This allows crossing from bottom half to top half (more countries)
  q50_low <- quantile(policy_data$policy, 0.50, na.rm = TRUE)
  q50_high <- quantile(policy_data$policy, 0.50, na.rm = TRUE)  # Same as median
  # For crossing threshold, use: below median to above median
  median_policy <- median(policy_data$policy, na.rm = TRUE)
  sd_policy <- sd(policy_data$policy, na.rm = TRUE)
  
  if (is.na(sd_policy) || sd_policy == 0) {
    cat("Warning: Insufficient variation in", policy_var, "\n")
    return(data.table(ISO3 = character(), treated = logical(), event_year = numeric()))
  }
  
  # For each country, identify first year when policy crosses threshold
  # Sort the entire dataset first
  setorder(policy_data, ISO3, Year)
  
  # Process each country separately
  treatment_years <- policy_data[, {
    treated <- FALSE
    event_year <- NA_real_
    
    # Need at least 2 years of data (3 if stability check required)
    min_years <- if (require_stability) 3 else 2
    
    if (.N >= min_years) {
      policy_vals <- .SD$policy
      year_vals <- .SD$Year
      
      # Start checking from year 2 (or 3 if stability check required)
      start_year <- if (require_stability) 3 else 2
      for (i in start_year:.N) {
        curr_policy <- policy_vals[i]
        curr_year <- year_vals[i]
        prev_policy <- policy_vals[i-1]
        prev2_policy <- if (i >= 3) policy_vals[i-2] else NA
        
        # Check if we have next year for sustained check
        has_next <- i < .N
        next_policy <- if (has_next) policy_vals[i+1] else NA
        
        # Check data availability based on requirements
        has_prev_data <- !is.na(prev_policy) && !is.na(curr_policy)
        has_stability_data <- if (require_stability) (!is.na(prev2_policy)) else TRUE
        
        if (has_prev_data && has_stability_data) {
          
          # Criterion 1: Large jump - increase of at least min_jump_sd standard deviations
          jump_size <- curr_policy - prev_policy
          large_jump <- jump_size >= min_jump_sd * sd_policy
          
          # Criterion 2: Cross from bottom 50% to top 50% (using median - more relaxed)
          # Changed from bottom 40% to top 40% to median crossing for more countries
          crosses_threshold <- prev_policy <= median_policy && curr_policy > median_policy
          
          # Criterion 3: Stability requirement - policy should be relatively stable in 2 years before
          # (not already on a strong upward trend)
          if (require_stability && !is.na(prev2_policy)) {
            prev_change <- prev_policy - prev2_policy
            # Allow small increases but not large ones (pre-trend check)
            stable_pre <- abs(prev_change) < 0.5 * sd_policy
          } else {
            stable_pre <- TRUE
          }
          
          # Criterion 4: Sustained increase - policy should remain high in next year
          if (require_sustained && has_next && !is.na(next_policy)) {
            sustained <- next_policy > median_policy  # Remains above median
          } else if (require_sustained) {
            sustained <- FALSE  # Can't check, so don't treat
          } else {
            sustained <- TRUE
          }
          
          # All criteria must be met
          if (large_jump && crosses_threshold && stable_pre && sustained) {
            treated <- TRUE
            event_year <- curr_year
            break
          }
        }
      }
    }
    
    list(treated = treated, event_year = event_year)
  }, by = ISO3]
  
  return(treatment_years)
}

# Identify treatments for each policy component
# Using RELAXED criteria to get more countries while still focusing on components
cat("\n=== IDENTIFYING TREATMENT EVENTS FOR SPECIFIC POLICY COMPONENTS ===\n")
cat("Focusing on individual policy components for clearer effects:\n")
cat("  - Property Rights (privatization)\n")
cat("  - Trade Freedom (trade liberalization)\n")
cat("  - Investment Freedom (capital account liberalization)\n")
cat("  - Business Freedom (deregulation)\n")
cat("  - Financial Freedom (financial liberalization)\n")
cat("  - Monetary Freedom (monetary policy)\n\n")
cat("Using RELAXED treatment definition to get more countries:\n")
cat("  - Minimum jump: 0.5 standard deviation (relaxed from 1.0)\n")
cat("  - Must cross from bottom 50% to top 50% (relaxed from 40%/60%)\n")
cat("  - No stability requirement (allows some pre-trends)\n")
cat("  - No sustained increase requirement\n\n")

treatments_list <- list()
for (policy in policy_vars) {
  if (policy %in% names(data)) {
    cat("Identifying treatment for:", policy, "\n")
    # Start with relaxed criteria to get more countries
    treatments_list[[policy]] <- identify_treatment(data, policy, 
                                                      min_jump_sd = 0.5,
                                                      require_stability = FALSE,
                                                      require_sustained = FALSE)
    
    # Count treated countries
    n_treated <- sum(treatments_list[[policy]]$treated, na.rm = TRUE)
    cat("  -> Found", n_treated, "treated countries\n")
    
    # If still too few, try even more relaxed
    if (n_treated < 10) {
      cat("  -> Trying more relaxed criteria (0.4 SD, bottom 50% to top 50%)...\n")
      # Temporarily modify thresholds in the function call
      # We'll need to adjust the function to allow different percentile thresholds
      treatments_list[[policy]] <- identify_treatment(data, policy,
                                                        min_jump_sd = 0.4,
                                                        require_stability = FALSE,
                                                        require_sustained = FALSE)
      n_treated <- sum(treatments_list[[policy]]$treated, na.rm = TRUE)
      cat("  -> Found", n_treated, "treated countries with more relaxed criteria\n")
    }
  } else {
    cat("  -> Policy", policy, "not found in data\n")
  }
}

# Also check aggregate scores for comparison
cat("\n=== CHECKING AGGREGATE SCORES FOR COMPARISON ===\n")
for (policy in aggregate_scores) {
  if (policy %in% names(data)) {
    cat("Identifying treatment for:", policy, "\n")
    treatments_list[[policy]] <- identify_treatment(data, policy,
                                                      min_jump_sd = 0.5,
                                                      require_stability = FALSE,
                                                      require_sustained = FALSE)
    n_treated <- sum(treatments_list[[policy]]$treated, na.rm = TRUE)
    cat("  -> Found", n_treated, "treated countries\n")
  }
}

# Merge treatment indicators back to main dataset
for (policy in names(treatments_list)) {
  treatment_col <- paste0("treated_", gsub(" ", "_", policy))
  event_col <- paste0("event_year_", gsub(" ", "_", policy))
  
  # Merge treatment indicators
  data <- merge(data, 
                treatments_list[[policy]][, .(ISO3, 
                                               treated = as.numeric(treated),
                                               event_year)],
                by = "ISO3", 
                all.x = TRUE, 
                suffixes = c("", paste0("_", gsub(" ", "_", policy))))
  
  # Rename columns
  setnames(data, 
           old = c("treated", "event_year"),
           new = c(treatment_col, event_col))
}

# Select primary policy based on which has most treated countries
# Prioritize specific components over aggregate scores
cat("\n=== SELECTING PRIMARY POLICY FOR ANALYSIS ===\n")

# Count treated countries for each policy
policy_counts <- data.frame(
  policy = character(),
  n_treated = numeric(),
  stringsAsFactors = FALSE
)

for (policy in names(treatments_list)) {
  n_treated <- sum(treatments_list[[policy]]$treated, na.rm = TRUE)
  if (n_treated > 0) {
    policy_counts <- rbind(policy_counts,
                          data.frame(policy = policy, n_treated = n_treated))
  }
}

# Sort by number of treated countries
if (nrow(policy_counts) > 0) {
  policy_counts <- policy_counts[order(-policy_counts$n_treated), ]
  cat("\nPolicies ranked by number of treated countries:\n")
  print(policy_counts)
  
  # Select the policy with most treated countries (prioritize components)
  # First try to find a component with sufficient countries (lowered threshold to 10)
  component_policies <- policy_counts[policy_counts$policy %in% policy_vars, ]
  if (nrow(component_policies) > 0 && max(component_policies$n_treated) >= 10) {
    primary_policy <- component_policies$policy[1]
    cat("\nSelected component policy:", primary_policy, 
        "with", component_policies$n_treated[1], "treated countries\n")
  } else if (nrow(component_policies) > 0) {
    # Use component even if < 10, as long as it's the best option
    primary_policy <- component_policies$policy[1]
    cat("\nSelected component policy:", primary_policy, 
        "with", component_policies$n_treated[1], "treated countries\n")
    cat("(Note: Lower than ideal, but best component option)\n")
  } else {
    # Fall back to aggregate score if components don't have enough
    if (nrow(policy_counts) > 0) {
      primary_policy <- policy_counts$policy[1]
      cat("\nSelected policy:", primary_policy, 
          "with", policy_counts$n_treated[1], "treated countries\n")
      cat("(Note: Using aggregate score because components have insufficient treatment)\n")
    } else {
      primary_policy <- "Trade Freedom"  # Default fallback
      cat("\nWarning: No policies found with treated countries. Using default:", primary_policy, "\n")
    }
  }
  
  treatment_var <- paste0("treated_", gsub(" ", "_", primary_policy))
  event_year_var <- paste0("event_year_", gsub(" ", "_", primary_policy))
} else {
  cat("\nWarning: No policies found with treated countries!\n")
  primary_policy <- "Trade Freedom"  # Default fallback
  treatment_var <- paste0("treated_", gsub(" ", "_", primary_policy))
  event_year_var <- paste0("event_year_", gsub(" ", "_", primary_policy))
}

# Create main treatment indicator
data[, treated := get(treatment_var)]
data[, event_year := get(event_year_var)]
data[is.na(treated), treated := 0]

# Extract and display treatment countries and years
cat("\n=== TREATMENT COUNTRIES AND YEARS ===\n")
if (primary_policy %in% names(treatments_list)) {
  treatment_info <- treatments_list[[primary_policy]]
  treated_countries <- treatment_info[treated == TRUE & !is.na(ISO3) & ISO3 != "", .(ISO3, event_year)]
  
  if (nrow(treated_countries) > 0) {
    # Sort by country code
    setorder(treated_countries, ISO3)
    
    cat("\nCountries with treatment events (", primary_policy, "):\n", sep = "")
    print(treated_countries)
    cat("\nTotal number of treated countries:", nrow(treated_countries), "\n")
    
    # Also save to CSV
    write.csv(treated_countries, "treatment_countries_years.csv", row.names = FALSE)
    cat("\nTreatment list saved to: treatment_countries_years.csv\n")
  } else {
    cat("\nNo countries identified as treated.\n")
  }
} else {
  cat("\nWarning:", primary_policy, "treatment not found in treatments_list.\n")
}

# ============================================================================
# 2.5. RUN SEPARATE ANALYSES FOR EACH POLICY COMPONENT
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RUNNING SEPARATE EVENT STUDIES FOR EACH POLICY COMPONENT\n")
cat(rep("=", 70), "\n", sep = "")

# Store results for each policy component
component_results_summary <- data.frame(
  policy_component = character(),
  outcome = character(),
  ate = numeric(),
  se = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  n_treated = numeric(),
  stringsAsFactors = FALSE
)

# Run analysis for each policy component that has sufficient treatment
for (policy_comp in policy_vars) {
  if (policy_comp %in% names(treatments_list)) {
    n_treated_comp <- sum(treatments_list[[policy_comp]]$treated, na.rm = TRUE)
    
    if (n_treated_comp >= 3) {  # Lowered threshold to get more components analyzed
      cat("\n", rep("-", 70), "\n", sep = "")
      cat("Analyzing:", policy_comp, "(", n_treated_comp, "treated countries)\n")
      cat(rep("-", 70), "\n", sep = "")
      
      # Create treatment indicators for this component
      treatment_col_comp <- paste0("treated_", gsub(" ", "_", policy_comp))
      event_col_comp <- paste0("event_year_", gsub(" ", "_", policy_comp))
      
      # Merge treatment indicators
      data_comp <- copy(data)
      data_comp <- merge(data_comp,
                        treatments_list[[policy_comp]][, .(ISO3,
                                                           treated_comp = as.numeric(treated),
                                                           event_year_comp = event_year)],
                        by = "ISO3", all.x = TRUE)
      
      data_comp[is.na(treated_comp), treated_comp := 0]
      
      # Create event time for this component
      data_comp[, event_time_comp := Year - event_year_comp]
      data_comp[is.na(event_time_comp) | treated_comp == 0, event_time_comp := NA]
      data_comp[event_time_comp < -5 | event_time_comp > 10, event_time_comp := NA]
      
      # Create event time indicators
      for (k in event_times) {
        var_name <- paste0("event_time_", k, "_comp")
        data_comp[, (var_name) := as.numeric(event_time_comp == k & treated_comp == 1)]
        data_comp[is.na(get(var_name)), (var_name) := 0]
      }
      
      # Run event study for key outcomes
      key_outcomes <- c("log_GDP_per_capita")
      
      # Add GDP growth if available
      gdp_growth_cols <- names(data_comp)[grepl("GDP.*growth.*annual", names(data_comp), ignore.case = TRUE)]
      if (length(gdp_growth_cols) > 0) {
        key_outcomes <- c(key_outcomes, gdp_growth_cols[1])
      }
      
      # Add Gini if available
      if ("Gini.index" %in% names(data_comp)) {
        key_outcomes <- c(key_outcomes, "Gini.index")
      }
      
      # Add poverty measure if available
      poverty_cols <- names(data_comp)[grepl("Poverty.*headcount.*3\\.65", names(data_comp), ignore.case = TRUE)]
      if (length(poverty_cols) > 0) {
        key_outcomes <- c(key_outcomes, poverty_cols[1])
      }
      
      for (outcome in key_outcomes) {
        if (outcome %in% names(data_comp)) {
          # Create formula
          event_var_names_comp <- paste0("`event_time_", event_times, "_comp`")
          event_terms_comp <- paste(event_var_names_comp, collapse = " + ")
          
          outcome_formula <- if (grepl(" ", outcome) || grepl("[()]", outcome)) {
            paste0("`", outcome, "`")
          } else {
            outcome
          }
          
          formula_str <- paste0(outcome_formula, " ~ ", event_terms_comp, " | ISO3 + Year")
          
          tryCatch({
            model_comp <- feols(as.formula(formula_str),
                               data = data_comp[!is.na(get(outcome))],
                               cluster = ~ ISO3)
            
            # Calculate ATE
            ate_result <- calculate_ate(model_comp, event_times)
            if (!is.null(ate_result)) {
              component_results_summary <- rbind(component_results_summary,
                                                data.frame(
                                                  policy_component = policy_comp,
                                                  outcome = outcome,
                                                  ate = ate_result$ate,
                                                  se = ate_result$se,
                                                  ci_lower = ate_result$ci_lower,
                                                  ci_upper = ate_result$ci_upper,
                                                  n_treated = n_treated_comp
                                                ))
              
              # Check if significant
              is_sig <- (ate_result$ci_lower > 0 && ate_result$ci_upper > 0) || 
                       (ate_result$ci_lower < 0 && ate_result$ci_upper < 0)
              sig_marker <- if (is_sig) "*** SIGNIFICANT ***" else ""
              
              cat("  ", outcome, ": ATE =", round(ate_result$ate, 4), 
                  "CI: [", round(ate_result$ci_lower, 4), ",", 
                  round(ate_result$ci_upper, 4), "]", sig_marker, "\n")
            }
          }, error = function(e) {
            cat("  ", outcome, ": Could not estimate -", e$message, "\n")
          })
        }
      }
    } else {
      cat("\nSkipping", policy_comp, "- only", n_treated_comp, "treated countries (need >= 5)\n")
    }
  }
}

# Save component results summary
if (nrow(component_results_summary) > 0) {
  write.csv(component_results_summary, "policy_component_results_summary.csv", row.names = FALSE)
  cat("\n\nComponent results saved to: policy_component_results_summary.csv\n")
  
  # Print summary table
  cat("\n=== SUMMARY OF POLICY COMPONENT EFFECTS ===\n")
  print(component_results_summary)
}

# ============================================================================
# 3. EVENT STUDY SPECIFICATION
# ============================================================================

# Create event time variable (years relative to treatment)
data[, event_time := Year - event_year]
data[is.na(event_time) | treated == 0, event_time := NA]

# Restrict to reasonable window around treatment (-5 to +10 years)
data[event_time < -5 | event_time > 10, event_time := NA]

# Create lead/lag indicators (normalize to -1, so -1 is omitted)
# We'll create indicators for -5, -4, -3, -2, 0, 1, 2, ..., 10
event_times <- c(-5:-2, 0:10)  # Exclude -1 as reference

for (k in event_times) {
  var_name <- paste0("event_time_", k)
  data[, (var_name) := as.numeric(event_time == k & treated == 1)]
  data[is.na(get(var_name)), (var_name) := 0]
}

# ============================================================================
# 3.5. HELPER FUNCTIONS
# ============================================================================

# Function to calculate average treatment effect
calculate_ate <- function(model, event_times) {
  post_periods <- event_times[event_times >= 0]
  post_vars <- paste0("event_time_", post_periods)
  
  coefs <- coef(model)
  ses <- se(model)
  
  available_vars <- post_vars[post_vars %in% names(coefs)]
  
  if (length(available_vars) > 0) {
    post_coefs <- coefs[available_vars]
    ate <- mean(post_coefs, na.rm = TRUE)
    
    post_ses <- ses[available_vars]
    ate_se <- sqrt(mean(post_ses^2, na.rm = TRUE) / length(available_vars))
    
    return(list(ate = ate, se = ate_se, 
                ci_lower = ate - 1.96 * ate_se,
                ci_upper = ate + 1.96 * ate_se))
  }
  return(NULL)
}

# ============================================================================
# 4. REGRESSION ANALYSIS - BY COUNTRY GROUP
# ============================================================================

# Function to run event study for a specific country group
run_event_study_for_group <- function(data_subset, group_name, event_times, 
                                       aggregate_outcomes, distributional_outcomes) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Running event study for:", group_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Count treated countries in this group
  treated_in_group <- unique(data_subset[treated == 1, .(ISO3)])
  n_treated_group <- nrow(treated_in_group)
  cat("Treated countries in", group_name, ":", n_treated_group, "\n")
  
  if (n_treated_group < 3) {
    cat("Warning: Too few treated countries (", n_treated_group, 
        ") in ", group_name, ". Skipping this group.\n", sep = "")
    return(list(aggregate_results = list(), distributional_results = list()))
  }
  
  aggregate_results_group <- list()
  distributional_results_group <- list()
  
  # Run aggregate outcomes
  for (outcome in aggregate_outcomes) {
    if (outcome %in% names(data_subset)) {
      cat("\nRunning event study for:", outcome, "\n")
      
      event_var_names <- paste0("`event_time_", event_times, "`")
      event_terms <- paste(event_var_names, collapse = " + ")
      
      outcome_formula <- if (grepl(" ", outcome) || grepl("[()]", outcome)) {
        paste0("`", outcome, "`")
      } else {
        outcome
      }
      
      formula_str <- paste0(outcome_formula, " ~ ", event_terms, " | ISO3 + Year")
      
      tryCatch({
        model <- feols(as.formula(formula_str), 
                       data = data_subset[!is.na(get(outcome))],
                       cluster = ~ ISO3)
        
        aggregate_results_group[[outcome]] <- model
      }, error = function(e) {
        cat("Warning: Could not estimate model for", outcome, "-", e$message, "\n")
      })
    }
  }
  
  # Run distributional outcomes
  for (outcome in distributional_outcomes) {
    outcome_exists <- FALSE
    outcome_col <- NULL
    
    if (outcome %in% names(data_subset)) {
      outcome_exists <- TRUE
      outcome_col <- outcome
    } else {
      matching_cols <- names(data_subset)[grepl(gsub(" ", ".*", outcome), 
                                                 names(data_subset), ignore.case = TRUE)]
      if (length(matching_cols) > 0) {
        outcome_exists <- TRUE
        outcome_col <- matching_cols[1]
      }
    }
    
    if (outcome_exists && !is.null(outcome_col)) {
      cat("\nRunning event study for:", outcome, "\n")
      
      if (grepl(" ", outcome_col) || grepl("%", outcome_col)) {
        outcome_formula <- paste0("`", outcome_col, "`")
      } else {
        outcome_formula <- outcome_col
      }
      
      event_var_names <- paste0("`event_time_", event_times, "`")
      event_terms <- paste(event_var_names, collapse = " + ")
      
      formula_str <- paste0(outcome_formula, " ~ ", event_terms, " | ISO3 + Year")
      
      tryCatch({
        model <- feols(as.formula(formula_str),
                       data = data_subset[!is.na(get(outcome_col))],
                       cluster = ~ ISO3)
        
        distributional_results_group[[outcome]] <- model
      }, error = function(e) {
        cat("Warning: Could not estimate model for", outcome, "-", e$message, "\n")
      })
    }
  }
  
  return(list(aggregate_results = aggregate_results_group,
              distributional_results = distributional_results_group))
}

# Get list of groups with sufficient treated countries
groups_with_treatment <- data[treated == 1, .(n_treated = .N, 
                                                n_countries = length(unique(ISO3))), 
                               by = country_group]
groups_with_treatment <- groups_with_treatment[n_countries >= 3]  # Lowered from 2 to 3 for better quality

cat("\n=== COUNTRY GROUPS WITH SUFFICIENT TREATMENT ===\n")
print(groups_with_treatment)

# Store results by group
results_by_group <- list()

# Run event studies for each group
for (group in groups_with_treatment$country_group) {
  data_group <- data[country_group == group]
  
  # Get outcomes (defined later in original code, but we'll define them here)
  aggregate_outcomes <- c()
  if ("log_GDP_per_capita" %in% names(data_group)) {
    aggregate_outcomes <- c(aggregate_outcomes, "log_GDP_per_capita")
  }
  
  gdp_growth_cols <- names(data_group)[grepl("GDP.*growth|GDP.*Growth", names(data_group), ignore.case = TRUE) & 
                                       grepl("annual|%", names(data_group), ignore.case = TRUE)]
  if (length(gdp_growth_cols) > 0) {
    aggregate_outcomes <- c(aggregate_outcomes, gdp_growth_cols)
  }
  
  gdp_pc_growth_cols <- names(data_group)[grepl("GDP.*per.*capita.*growth|GDP.*capita.*Growth", names(data_group), ignore.case = TRUE)]
  if (length(gdp_pc_growth_cols) > 0) {
    aggregate_outcomes <- c(aggregate_outcomes, gdp_pc_growth_cols)
  }
  
  gni_cols <- names(data_group)[grepl("^GNI", names(data_group), ignore.case = TRUE)]
  if (length(gni_cols) > 0) {
    gni_pc_cols <- gni_cols[grepl("per.*capita|capita", gni_cols, ignore.case = TRUE) & 
                            !grepl("growth", gni_cols, ignore.case = TRUE)]
    gni_growth_cols <- gni_cols[grepl("growth", gni_cols, ignore.case = TRUE)]
    aggregate_outcomes <- c(aggregate_outcomes, gni_pc_cols, gni_growth_cols)
  }
  
  distributional_outcomes <- c()
  poverty_cols <- names(data_group)[grepl("poverty|Poverty", names(data_group), ignore.case = TRUE) & 
                                     grepl("headcount|3\\.65|6\\.85|\\$3|\\$6", names(data_group), ignore.case = TRUE)]
  if (length(poverty_cols) > 0) {
    distributional_outcomes <- c(distributional_outcomes, poverty_cols)
  }
  
  gini_cols <- names(data_group)[grepl("Gini|gini", names(data_group), ignore.case = TRUE)]
  if (length(gini_cols) > 0) {
    distributional_outcomes <- c(distributional_outcomes, gini_cols)
  }
  
  income_share_cols <- names(data_group)[grepl("Income.*share|income.*share", names(data_group), ignore.case = TRUE) & 
                                         grepl("10%|lowest|highest", names(data_group), ignore.case = TRUE)]
  if (length(income_share_cols) > 0) {
    distributional_outcomes <- c(distributional_outcomes, income_share_cols)
  }
  
  # Run event study for this group
  group_results <- run_event_study_for_group(data_group, group, event_times,
                                             aggregate_outcomes, distributional_outcomes)
  results_by_group[[group]] <- group_results
  
  # Create plots for this group
  group_name_clean <- gsub("[^A-Za-z0-9]", "_", group)
  
  # Aggregate plots
  for (outcome in names(group_results$aggregate_results)) {
    title <- paste("Event Study:", outcome, "-", group)
    p <- create_event_plot(group_results$aggregate_results[[outcome]], title, event_times, debug = FALSE)
    ggsave(paste0("event_study_", outcome, "_", group_name_clean, ".png"), 
           plot = p, width = 10, height = 6, dpi = 300)
  }
  
  # Distributional plots
  for (outcome in names(group_results$distributional_results)) {
    title <- paste("Event Study:", outcome, "-", group)
    outcome_clean <- gsub(" ", "_", outcome)
    p <- create_event_plot(group_results$distributional_results[[outcome]], title, event_times, debug = FALSE)
    ggsave(paste0("event_study_", outcome_clean, "_", group_name_clean, ".png"),
           plot = p, width = 10, height = 6, dpi = 300)
  }
  
  # Save results summary for this group
  results_summary_group <- data.frame(
    outcome = character(),
    ate = numeric(),
    se = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (outcome in names(group_results$aggregate_results)) {
    ate_result <- calculate_ate(group_results$aggregate_results[[outcome]], event_times)
    if (!is.null(ate_result)) {
      results_summary_group <- rbind(results_summary_group,
                                     data.frame(outcome = outcome,
                                               ate = ate_result$ate,
                                               se = ate_result$se,
                                               ci_lower = ate_result$ci_lower,
                                               ci_upper = ate_result$ci_upper))
    }
  }
  
  for (outcome in names(group_results$distributional_results)) {
    ate_result <- calculate_ate(group_results$distributional_results[[outcome]], event_times)
    if (!is.null(ate_result)) {
      results_summary_group <- rbind(results_summary_group,
                                     data.frame(outcome = outcome,
                                               ate = ate_result$ate,
                                               se = ate_result$se,
                                               ci_lower = ate_result$ci_lower,
                                               ci_upper = ate_result$ci_upper))
    }
  }
  
  if (nrow(results_summary_group) > 0) {
    results_summary_group$country_group <- group
    write.csv(results_summary_group, 
              paste0("treatment_effects_summary_", group_name_clean, ".csv"), 
              row.names = FALSE)
  }
}

# ============================================================================
# 4. REGRESSION ANALYSIS - AGGREGATE OUTCOMES (OVERALL)
# ============================================================================

# Aggregate outcomes - try renamed first, then original names
aggregate_outcomes <- c()

# Always include log_GDP_per_capita if it exists
if ("log_GDP_per_capita" %in% names(data)) {
  aggregate_outcomes <- c(aggregate_outcomes, "log_GDP_per_capita")
}

# Try to find GDP growth columns - search more flexibly
gdp_growth_cols <- names(data)[grepl("GDP.*growth|GDP.*Growth", names(data), ignore.case = TRUE) & 
                                grepl("annual|%", names(data), ignore.case = TRUE)]
if (length(gdp_growth_cols) > 0) {
  aggregate_outcomes <- c(aggregate_outcomes, gdp_growth_cols)
  cat("Found GDP growth columns:", paste(gdp_growth_cols, collapse = ", "), "\n")
}

# Try to find GDP per capita growth
gdp_pc_growth_cols <- names(data)[grepl("GDP.*per.*capita.*growth|GDP.*capita.*Growth", names(data), ignore.case = TRUE)]
if (length(gdp_pc_growth_cols) > 0) {
  aggregate_outcomes <- c(aggregate_outcomes, gdp_pc_growth_cols)
  cat("Found GDP per capita growth columns:", paste(gdp_pc_growth_cols, collapse = ", "), "\n")
}

# Try to find GNI columns
gni_cols <- names(data)[grepl("^GNI", names(data), ignore.case = TRUE)]
if (length(gni_cols) > 0) {
  # Filter to per capita and growth
  gni_pc_cols <- gni_cols[grepl("per.*capita|capita", gni_cols, ignore.case = TRUE) & 
                          !grepl("growth", gni_cols, ignore.case = TRUE)]
  gni_growth_cols <- gni_cols[grepl("growth", gni_cols, ignore.case = TRUE)]
  aggregate_outcomes <- c(aggregate_outcomes, gni_pc_cols, gni_growth_cols)
  if (length(gni_pc_cols) > 0 || length(gni_growth_cols) > 0) {
    cat("Found GNI columns:", paste(c(gni_pc_cols, gni_growth_cols), collapse = ", "), "\n")
  }
}

cat("Aggregate outcomes to analyze:", paste(aggregate_outcomes, collapse = ", "), "\n")

# Store regression results
aggregate_results <- list()

# Run event study for each aggregate outcome
for (outcome in aggregate_outcomes) {
  if (outcome %in% names(data)) {
    cat("\nRunning event study for:", outcome, "\n")
    
    # Create formula with all event time indicators
    # Handle negative numbers in variable names by using backticks
    event_var_names <- paste0("`event_time_", event_times, "`")
    event_terms <- paste(event_var_names, collapse = " + ")
    
    # Base specification: country and year fixed effects
    # Use backticks for outcome if it has special characters
    outcome_formula <- if (grepl(" ", outcome) || grepl("[()]", outcome)) {
      paste0("`", outcome, "`")
    } else {
      outcome
    }
    
    formula_str <- paste0(outcome_formula, " ~ ", event_terms, " | ISO3 + Year")
    
    # Run regression with clustered standard errors
    tryCatch({
      model <- feols(as.formula(formula_str), 
                     data = data[!is.na(get(outcome))],
                     cluster = ~ ISO3)
      
      aggregate_results[[outcome]] <- model
    }, error = function(e) {
      cat("Warning: Could not estimate model for", outcome, "-", e$message, "\n")
      cat("This may be due to insufficient variation (all observations are singletons).\n")
    })
    
    # Also run with income group interactions
    # Create interaction terms (only if income group data exists)
    if (sum(!is.na(data$low_income)) > 0) {
      # Create interaction terms
      for (k in event_times) {
        var_name <- paste0("event_time_", k)
        low_name <- paste0("event_time_", k, "_low")
        mid_name <- paste0("event_time_", k, "_mid")
        high_name <- paste0("event_time_", k, "_high")
        
        if (var_name %in% names(data)) {
          data[, (low_name) := get(var_name) * low_income]
          data[, (mid_name) := get(var_name) * middle_income]
          data[, (high_name) := get(var_name) * high_income]
        }
      }
      
      # Interaction specification
      event_terms_low <- paste0("`event_time_", event_times, "_low`", collapse = " + ")
      event_terms_mid <- paste0("`event_time_", event_times, "_mid`", collapse = " + ")
      event_terms_high <- paste0("`event_time_", event_times, "_high`", collapse = " + ")
      
      formula_interact <- paste0(outcome_formula, " ~ ", event_terms_low, " + ", 
                                 event_terms_mid, " + ", event_terms_high, 
                                 " | ISO3 + Year")
      
      # Check if formula is valid
      tryCatch({
        model_interact <- feols(as.formula(formula_interact),
                               data = data[!is.na(get(outcome))],
                               cluster = ~ ISO3)
        aggregate_results[[paste0(outcome, "_interact")]] <- model_interact
      }, error = function(e) {
        cat("Warning: Could not estimate interaction model for", outcome, "\n")
      })
    }
  }
}

# ============================================================================
# 5. REGRESSION ANALYSIS - DISTRIBUTIONAL OUTCOMES
# ============================================================================

# Distributional outcomes - try renamed first, then original names
distributional_outcomes <- c()

# Poverty measures - search flexibly
poverty_cols <- names(data)[grepl("poverty|Poverty", names(data), ignore.case = TRUE) & 
                            grepl("headcount|3\\.65|6\\.85|\\$3|\\$6", names(data), ignore.case = TRUE)]
if (length(poverty_cols) > 0) {
  distributional_outcomes <- c(distributional_outcomes, poverty_cols)
  cat("Found poverty columns:", paste(poverty_cols, collapse = ", "), "\n")
}

# Inequality measures - Gini index
gini_cols <- names(data)[grepl("Gini|gini", names(data), ignore.case = TRUE)]
if (length(gini_cols) > 0) {
  distributional_outcomes <- c(distributional_outcomes, gini_cols)
  cat("Found Gini columns:", paste(gini_cols, collapse = ", "), "\n")
}

# Income share measures
income_share_cols <- names(data)[grepl("Income.*share|income.*share", names(data), ignore.case = TRUE) & 
                                  grepl("10%|lowest|highest", names(data), ignore.case = TRUE)]
if (length(income_share_cols) > 0) {
  distributional_outcomes <- c(distributional_outcomes, income_share_cols)
  cat("Found income share columns:", paste(income_share_cols, collapse = ", "), "\n")
}

# Bottom 40% growth
bottom40_cols <- names(data)[grepl("bottom.*40|40.*bottom|poorest.*40|40.*poorest", names(data), ignore.case = TRUE) & 
                              grepl("growth|Growth", names(data), ignore.case = TRUE)]
if (length(bottom40_cols) > 0) {
  distributional_outcomes <- c(distributional_outcomes, bottom40_cols)
  cat("Found bottom 40% growth columns:", paste(bottom40_cols, collapse = ", "), "\n")
}

cat("Distributional outcomes to analyze:", paste(distributional_outcomes, collapse = ", "), "\n")

# Note: Some column names have spaces and special characters
# We'll handle them carefully in the regression loops

# Store regression results
distributional_results <- list()

# Run event study for each distributional outcome
for (outcome in distributional_outcomes) {
  # Check if outcome exists in data (handle spaces in column names)
  outcome_exists <- FALSE
  outcome_col <- NULL
  
  # Try exact match first
  if (outcome %in% names(data)) {
    outcome_exists <- TRUE
    outcome_col <- outcome
  } else {
    # Try to find similar column name
    matching_cols <- names(data)[grepl(gsub(" ", ".*", outcome), names(data), ignore.case = TRUE)]
    if (length(matching_cols) > 0) {
      outcome_exists <- TRUE
      outcome_col <- matching_cols[1]
      cat("Note: Using column '", outcome_col, "' for outcome '", outcome, "'\n", sep = "")
    }
  }
  
  if (outcome_exists && !is.null(outcome_col)) {
    cat("\nRunning event study for:", outcome, "\n")
    
    # Create formula with all event time indicators
    # Use backticks for column names with spaces
    if (grepl(" ", outcome_col) || grepl("%", outcome_col)) {
      outcome_formula <- paste0("`", outcome_col, "`")
    } else {
      outcome_formula <- outcome_col
    }
    
    # Handle negative numbers in variable names by using backticks
    event_var_names <- paste0("`event_time_", event_times, "`")
    event_terms <- paste(event_var_names, collapse = " + ")
    
    # Base specification
    formula_str <- paste0(outcome_formula, " ~ ", event_terms, " | ISO3 + Year")
    
    # Run regression with clustered standard errors
    tryCatch({
      model <- feols(as.formula(formula_str),
                     data = data[!is.na(get(outcome_col))],
                     cluster = ~ ISO3)
      
      distributional_results[[outcome]] <- model
    }, error = function(e) {
      cat("Warning: Could not estimate model for", outcome, "-", e$message, "\n")
      cat("This may be due to insufficient variation (all observations are singletons).\n")
    })
    
    # Interaction specification with income groups (if income data exists)
    if (sum(!is.na(data$low_income)) > 0) {
      event_terms_low <- paste0("`event_time_", event_times, "_low`", collapse = " + ")
      event_terms_mid <- paste0("`event_time_", event_times, "_mid`", collapse = " + ")
      event_terms_high <- paste0("`event_time_", event_times, "_high`", collapse = " + ")
      
      formula_interact <- paste0(outcome_formula, " ~ ", event_terms_low, " + ",
                                 event_terms_mid, " + ", event_terms_high,
                                 " | ISO3 + Year")
      
      tryCatch({
        model_interact <- feols(as.formula(formula_interact),
                               data = data[!is.na(get(outcome_col))],
                               cluster = ~ ISO3)
        distributional_results[[paste0(outcome, "_interact")]] <- model_interact
      }, error = function(e) {
        cat("Warning: Could not estimate interaction model for", outcome, "\n")
      })
    }
  } else {
    cat("Warning: Outcome '", outcome, "' not found in data\n", sep = "")
  }
}

# ============================================================================
# 6. PRE-TREND TESTING
# ============================================================================

# Extract coefficients and standard errors for pre-trend analysis
extract_event_coefs <- function(model, event_times, debug = FALSE) {
  coefs <- coef(model)
  ses <- se(model)
  coef_names <- names(coefs)
  
  # Debug: print available coefficients (only if debug=TRUE)
  if (debug) {
    cat("Available coefficients in model:\n")
    event_coef_names <- coef_names[grepl("event_time", coef_names)]
    if (length(event_coef_names) > 0) {
      cat(paste(event_coef_names, collapse = ", "), "\n")
    } else {
      cat("No event_time coefficients found!\n")
    }
  }
  
  # Create full sequence of event times from -5 to 10, including -1
  full_event_times <- c(-5:-1, 0:10)
  
  # Extract event time coefficients
  event_coefs <- data.frame(
    event_time = full_event_times,
    coefficient = NA,
    se = NA,
    ci_lower = NA,
    ci_upper = NA
  )
  
  for (i in seq_along(full_event_times)) {
    k <- full_event_times[i]
    
    # Reference period (-1) is always 0
    if (k == -1) {
      event_coefs$coefficient[i] <- 0
      event_coefs$se[i] <- 0
      event_coefs$ci_lower[i] <- 0
      event_coefs$ci_upper[i] <- 0
    } else {
      # Try different possible variable names
      var_name <- paste0("event_time_", k)
      var_name_alt1 <- paste0("`event_time_", k, "`")  # With backticks
      var_name_alt2 <- gsub("-", ".", var_name)  # Fixest might convert - to .
      
      # Check which name exists
      found_var <- NULL
      if (var_name %in% names(coefs)) {
        found_var <- var_name
      } else if (var_name_alt1 %in% names(coefs)) {
        found_var <- var_name_alt1
      } else if (var_name_alt2 %in% names(coefs)) {
        found_var <- var_name_alt2
      } else {
        # Try to find any coefficient that matches the pattern
        pattern <- paste0("event_time.*", abs(k))
        matches <- coef_names[grepl(pattern, coef_names, ignore.case = TRUE)]
        if (length(matches) > 0) {
          found_var <- matches[1]
          if (debug) {
            cat("Found coefficient with pattern matching:", found_var, "for event_time", k, "\n")
          }
        }
      }
      
      if (!is.null(found_var) && found_var %in% names(coefs)) {
        event_coefs$coefficient[i] <- coefs[found_var]
        if (found_var %in% names(ses)) {
          event_coefs$se[i] <- ses[found_var]
          event_coefs$ci_lower[i] <- coefs[found_var] - 1.96 * ses[found_var]
          event_coefs$ci_upper[i] <- coefs[found_var] + 1.96 * ses[found_var]
        }
      }
      # If coefficient doesn't exist, leave as NA (will be shown but not connected)
    }
  }
  
  # Sort by event_time to ensure proper ordering
  event_coefs <- event_coefs[order(event_coefs$event_time), ]
  
  # Normalize to average of pre-periods to reduce discontinuity at t=-1
  # This centers the pre-period average around 0, making the plot smoother
  pre_period_coefs <- event_coefs$coefficient[event_coefs$event_time >= -5 & 
                                               event_coefs$event_time <= -2 & 
                                               !is.na(event_coefs$coefficient)]
  
  if (length(pre_period_coefs) > 0) {
    avg_pre_coef <- mean(pre_period_coefs, na.rm = TRUE)
    
    # Subtract average from all coefficients (including t=-1)
    # This makes the average of pre-periods = 0, reducing visual discontinuity
    event_coefs$coefficient <- event_coefs$coefficient - avg_pre_coef
    event_coefs$ci_lower[!is.na(event_coefs$ci_lower)] <- event_coefs$ci_lower[!is.na(event_coefs$ci_lower)] - avg_pre_coef
    event_coefs$ci_upper[!is.na(event_coefs$ci_upper)] <- event_coefs$ci_upper[!is.na(event_coefs$ci_upper)] - avg_pre_coef
    
    # Now t=-1 will be at its relative position (likely close to 0 if pre-trends are flat)
    # and the pre-period average will be approximately 0
  }
  
  # Debug: print what we extracted (only if debug=TRUE)
  if (debug) {
    cat("Extracted coefficients:\n")
    print(event_coefs[, c("event_time", "coefficient")])
    if (length(pre_period_coefs) > 0) {
      cat("Average of pre-period coefficients (before normalization):", mean(pre_period_coefs, na.rm = TRUE), "\n")
    }
  }
  
  return(event_coefs)
}

# Test pre-trends: test if lead coefficients (pre-treatment) are jointly zero
test_pre_trends <- function(model, event_times) {
  # Pre-treatment periods are negative event times
  pre_periods <- event_times[event_times < 0]
  
  if (length(pre_periods) > 0) {
    # Get coefficient names
    pre_vars <- paste0("event_time_", pre_periods)
    
    # Check which are in the model
    model_coefs <- names(coef(model))
    available_vars <- pre_vars[pre_vars %in% model_coefs]
    
    if (length(available_vars) > 0) {
      # Joint F-test using fixest's wald function
      tryCatch({
        # Use fixest's built-in wald test
        wald_test <- wald(model, available_vars)
        return(list(
          test_stat = wald_test$stat,
          p_value = wald_test$p,
          pre_periods = pre_periods
        ))
      }, error = function(e) {
        # Fallback: manual calculation
        coefs <- coef(model)[available_vars]
        vcov_mat <- vcov(model)[available_vars, available_vars]
        if (length(coefs) == 1) {
          test_stat <- (coefs / sqrt(vcov_mat))^2
          p_val <- 1 - pchisq(test_stat, df = 1)
        } else {
          test_stat <- as.numeric(t(coefs) %*% solve(vcov_mat) %*% coefs)
          p_val <- 1 - pchisq(test_stat, df = length(coefs))
        }
        return(list(
          test_stat = test_stat,
          p_value = p_val,
          pre_periods = pre_periods
        ))
      })
    }
  }
  return(NULL)
}

# Apply pre-trend tests
pre_trend_tests <- list()
for (outcome in names(aggregate_results)) {
  test <- test_pre_trends(aggregate_results[[outcome]], event_times)
  if (!is.null(test)) {
    pre_trend_tests[[outcome]] <- test
  }
}

for (outcome in names(distributional_results)) {
  test <- test_pre_trends(distributional_results[[outcome]], event_times)
  if (!is.null(test)) {
    pre_trend_tests[[outcome]] <- test
  }
}

# ============================================================================
# 7. VISUALIZATION
# ============================================================================

# Function to create event study plot
create_event_plot <- function(model, title, event_times, debug = FALSE) {
  coefs_data <- extract_event_coefs(model, event_times, debug = debug)
  
  # Separate pre and post periods
  coefs_data$period <- ifelse(coefs_data$event_time < 0, "Pre", "Post")
  
  # Create data for plotting - separate complete and incomplete cases
  coefs_complete <- coefs_data[!is.na(coefs_data$coefficient), ]
  coefs_missing <- coefs_data[is.na(coefs_data$coefficient), ]
  
  # Ensure we have all event times from -5 to 10
  all_times <- seq(-5, 10, by = 1)
  missing_times <- setdiff(all_times, coefs_data$event_time)
  if (length(missing_times) > 0) {
    # Add missing times with NA values
    missing_df <- data.frame(
      event_time = missing_times,
      coefficient = NA,
      se = NA,
      ci_lower = NA,
      ci_upper = NA,
      period = ifelse(missing_times < 0, "Pre", "Post")
    )
    coefs_data <- rbind(coefs_data, missing_df)
    coefs_data <- coefs_data[order(coefs_data$event_time), ]
    # Recalculate complete/missing
    coefs_complete <- coefs_data[!is.na(coefs_data$coefficient), ]
    coefs_missing <- coefs_data[is.na(coefs_data$coefficient), ]
  }
  
  # Split data into pre and post periods for reference
  coefs_pre <- coefs_complete[coefs_complete$event_time < 0, ]
  coefs_post <- coefs_complete[coefs_complete$event_time >= 0, ]
  coefs_ref <- coefs_complete[coefs_complete$event_time == -1, ]
  
  # Create base plot with all data
  p <- ggplot(coefs_data, aes(x = event_time, y = coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.5)
  
  # Add confidence intervals - now we can connect through t=-1 since we normalized
  if (nrow(coefs_complete) > 0) {
    p <- p + geom_ribbon(data = coefs_complete, aes(ymin = ci_lower, ymax = ci_upper), 
                         alpha = 0.2, fill = "blue")
  }
  
  # Add line connecting all points (including through t=-1)
  # Normalization to average of pre-periods should make this smooth
  if (nrow(coefs_complete) > 0) {
    p <- p + geom_line(data = coefs_complete, color = "blue", size = 1)
  }
  
  # Add points for all coefficients
  if (nrow(coefs_complete) > 0) {
    p <- p + geom_point(data = coefs_complete, color = "blue", size = 2)
  }
  
  # Add error bars for all coefficients
  if (nrow(coefs_complete) > 0) {
    p <- p + geom_errorbar(data = coefs_complete, aes(ymin = ci_lower, ymax = ci_upper), 
                           width = 0.2, color = "blue")
  }
  
  # Highlight reference period (t=-1) with a different style
  # Note: After normalization, t=-1 may not be exactly 0, but shows its relative position
  if (nrow(coefs_ref) > 0) {
    p <- p + geom_point(data = coefs_ref, color = "darkred", size = 3, shape = 21, fill = "white", stroke = 2)
    # Add a subtle horizontal line at t=-1's value to show it's the reference
    if (!is.na(coefs_ref$coefficient[1])) {
      p <- p + geom_hline(yintercept = coefs_ref$coefficient[1], linetype = "dotted", 
                          color = "darkred", alpha = 0.3, linewidth = 0.5)
    }
  }
  
  # Show missing periods as empty points or gaps
  if (nrow(coefs_missing) > 0) {
    p <- p + geom_point(data = coefs_missing, color = "gray", size = 1, shape = 1, alpha = 0.5)
  }
  
  # Ensure all event times are shown on x-axis - force all breaks
  p <- p + scale_x_continuous(breaks = seq(-5, 10, by = 1), 
                               minor_breaks = NULL,
                               limits = c(-5.5, 10.5),
                               expand = c(0.05, 0.05)) +
    labs(title = title,
         x = "Years Relative to Treatment",
         y = "Coefficient Estimate",
         subtitle = "Normalized to average of pre-periods (t=-5 to t=-2). Reference: t=-1 omitted from regression.") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  return(p)
}

# Create plots for aggregate outcomes
aggregate_plots <- list()
for (i in seq_along(aggregate_outcomes)) {
  outcome <- aggregate_outcomes[i]
  if (outcome %in% names(aggregate_results)) {
    title <- paste("Event Study:", outcome)
    # Debug only for first plot
    p <- create_event_plot(aggregate_results[[outcome]], title, event_times, debug = (i == 1))
    aggregate_plots[[outcome]] <- p
    
    # Save plot
    ggsave(paste0("event_study_", outcome, ".png"), 
           plot = p, width = 10, height = 6, dpi = 300)
  }
}

# Create plots for distributional outcomes
distributional_plots <- list()
for (i in seq_along(distributional_outcomes)) {
  outcome <- distributional_outcomes[i]
  if (outcome %in% names(distributional_results)) {
    title <- paste("Event Study:", outcome)
    # Debug only for first plot
    p <- create_event_plot(distributional_results[[outcome]], title, event_times, debug = (i == 1))
    distributional_plots[[outcome]] <- p
    
    # Save plot
    outcome_clean <- gsub(" ", "_", outcome)
    ggsave(paste0("event_study_", outcome_clean, ".png"),
           plot = p, width = 10, height = 6, dpi = 300)
  }
}

# ============================================================================
# 8. RESULTS SUMMARY
# ============================================================================

# Create summary table of average treatment effects
# Average effect is the average of post-treatment coefficients (t=0 to t=10)

calculate_ate <- function(model, event_times) {
  post_periods <- event_times[event_times >= 0]
  post_vars <- paste0("event_time_", post_periods)
  
  coefs <- coef(model)
  ses <- se(model)
  
  available_vars <- post_vars[post_vars %in% names(coefs)]
  
  if (length(available_vars) > 0) {
    post_coefs <- coefs[available_vars]
    # Simple average (could weight by precision)
    ate <- mean(post_coefs, na.rm = TRUE)
    
    # Calculate SE of average (assuming independence for simplicity)
    post_ses <- ses[available_vars]
    ate_se <- sqrt(mean(post_ses^2, na.rm = TRUE) / length(available_vars))
    
    return(list(ate = ate, se = ate_se, 
                ci_lower = ate - 1.96 * ate_se,
                ci_upper = ate + 1.96 * ate_se))
  }
  return(NULL)
}

# Compile results summary
results_summary <- data.frame(
  outcome = character(),
  ate = numeric(),
  se = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  stringsAsFactors = FALSE
)

# Aggregate outcomes
for (outcome in aggregate_outcomes) {
  if (outcome %in% names(aggregate_results)) {
    ate_result <- calculate_ate(aggregate_results[[outcome]], event_times)
    if (!is.null(ate_result)) {
      results_summary <- rbind(results_summary,
                               data.frame(outcome = outcome,
                                         ate = ate_result$ate,
                                         se = ate_result$se,
                                         ci_lower = ate_result$ci_lower,
                                         ci_upper = ate_result$ci_upper))
    }
  }
}

# Distributional outcomes
for (outcome in distributional_outcomes) {
  if (outcome %in% names(distributional_results)) {
    ate_result <- calculate_ate(distributional_results[[outcome]], event_times)
    if (!is.null(ate_result)) {
      results_summary <- rbind(results_summary,
                               data.frame(outcome = outcome,
                                         ate = ate_result$ate,
                                         se = ate_result$se,
                                         ci_lower = ate_result$ci_lower,
                                         ci_upper = ate_result$ci_upper))
    }
  }
}

# Export results
write.csv(results_summary, "treatment_effects_summary.csv", row.names = FALSE)

# Export pre-trend test results
pre_trend_df <- data.frame(
  outcome = names(pre_trend_tests),
  f_statistic = sapply(pre_trend_tests, function(x) x$test_stat),
  p_value = sapply(pre_trend_tests, function(x) x$p_value),
  stringsAsFactors = FALSE
)
write.csv(pre_trend_df, "pre_trend_tests.csv", row.names = FALSE)

# Create regression tables using modelsummary (if available) or simple summary
if (length(aggregate_results) > 0) {
  if (modelsummary_available) {
    tryCatch({
      modelsummary(aggregate_results, 
                   output = "aggregate_results_table.html",
                   stars = TRUE,
                   coef_rename = function(x) gsub("event_time_", "t=", x))
      cat("Aggregate results table saved to aggregate_results_table.html\n")
    }, error = function(e) {
      cat("Warning: Could not create modelsummary table. Saving simple summaries instead.\n")
      # Save simple summaries
      for (i in 1:length(aggregate_results)) {
        sink(paste0("aggregate_results_", names(aggregate_results)[i], ".txt"))
        print(summary(aggregate_results[[i]]))
        sink()
      }
    })
  } else {
    # Save simple summaries
    for (i in 1:length(aggregate_results)) {
      sink(paste0("aggregate_results_", names(aggregate_results)[i], ".txt"))
      print(summary(aggregate_results[[i]]))
      sink()
    }
    cat("Aggregate results summaries saved as text files\n")
  }
}

if (length(distributional_results) > 0) {
  if (modelsummary_available) {
    tryCatch({
      modelsummary(distributional_results,
                   output = "distributional_results_table.html",
                   stars = TRUE,
                   coef_rename = function(x) gsub("event_time_", "t=", x))
      cat("Distributional results table saved to distributional_results_table.html\n")
    }, error = function(e) {
      cat("Warning: Could not create modelsummary table. Saving simple summaries instead.\n")
      # Save simple summaries
      for (i in 1:length(distributional_results)) {
        sink(paste0("distributional_results_", names(distributional_results)[i], ".txt"))
        print(summary(distributional_results[[i]]))
        sink()
      }
    })
  } else {
    # Save simple summaries
    for (i in 1:length(distributional_results)) {
      sink(paste0("distributional_results_", names(distributional_results)[i], ".txt"))
      print(summary(distributional_results[[i]]))
      sink()
    }
    cat("Distributional results summaries saved as text files\n")
  }
}

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:\n")
cat("- treatment_effects_summary.csv\n")
cat("- pre_trend_tests.csv\n")
cat("- aggregate_results_table.html\n")
cat("- distributional_results_table.html\n")
cat("- Event study plots saved as PNG files\n")

