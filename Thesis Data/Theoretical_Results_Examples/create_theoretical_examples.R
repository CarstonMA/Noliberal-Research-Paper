# Create Theoretical Event Study Examples
# Uses same outcome labels as actual data, shows somewhat messy results

library(ggplot2)
library(data.table)

# ============================================================================
# HELPER FUNCTION: Create somewhat messy event study plot
# ============================================================================

create_somewhat_messy_plot <- function(outcome_name, coefficients, ses, output_path) {
  event_times <- -5:10
  
  # Create data frame
  plot_data <- data.frame(
    event_time = event_times,
    coefficient = coefficients,
    se = ses,
    ci_lower = coefficients - 1.96 * ses,
    ci_upper = coefficients + 1.96 * ses
  )
  
  plot_data$period <- ifelse(plot_data$event_time < 0, "Pre", "Post")
  plot_data_complete <- plot_data[!is.na(plot_data$coefficient), ]
  
  # Create plot
  p <- ggplot(plot_data, aes(x = event_time, y = coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_ribbon(data = plot_data_complete, aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.2, fill = "blue") +
    geom_line(data = plot_data_complete, color = "blue", size = 1) +
    geom_point(data = plot_data_complete, color = "blue", size = 2) +
    geom_errorbar(data = plot_data_complete, aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.2, color = "blue") +
    geom_point(data = plot_data[plot_data$event_time == -1, ], 
               color = "red", size = 3, shape = 21, fill = "white", stroke = 2) +
    scale_x_continuous(breaks = seq(-5, 10, by = 1), 
                       minor_breaks = NULL,
                       limits = c(-5.5, 10.5),
                       expand = c(0.05, 0.05)) +
    labs(title = paste("Event Study:", outcome_name),
         subtitle = "Reference period: Average of t=-5 to t=-2 (normalized to 0)",
         x = "Years Relative to Treatment",
         y = "Coefficient Estimate") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  ggsave(output_path, plot = p, width = 10, height = 6, dpi = 300)
  
  return(plot_data)
}

# ============================================================================
# AGGREGATE OUTCOMES - Somewhat Messy Examples
# ============================================================================

# 1. log_GDP_per_capita - Somewhat messy (mild pre-trend, some discontinuity, mixed post-effects)
create_log_gdp_pc_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild downward pre-trend, discontinuity at t=-1, mixed post-effects
  coefficients <- c(
    -0.035, -0.025, -0.015, -0.008,  # t=-5 to t=-2 (mild downward trend)
    0,                                # t=-1 (reference)
    0.012, 0.025, 0.038, 0.045, 0.052, 0.048, 0.055, 0.062, 0.058, 0.065, 0.060  # t=0 to t=10 (positive but noisy)
  )
  
  ses <- c(
    0.032, 0.033, 0.034, 0.035,  # Pre-periods
    0,                           # t=-1
    0.036, 0.035, 0.038, 0.037, 0.040, 0.041, 0.039, 0.042, 0.040, 0.041, 0.040  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_log_GDP_per_capita.png"
  data <- create_somewhat_messy_plot("log_GDP_per_capita", coefficients, ses, output_path)
  
  # Calculate ATE (average of post-periods)
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.12))
}

# 2. GDP growth (annual %) - Somewhat messy
create_gdp_growth_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: slight upward pre-trend, discontinuity, noisy positive post-effects
  coefficients <- c(
    0.15, 0.25, 0.35, 0.42,  # t=-5 to t=-2 (slight upward trend)
    0,                        # t=-1
    0.55, 0.68, 0.72, 0.65, 0.78, 0.70, 0.75, 0.82, 0.75, 0.80, 0.78  # t=0 to t=10 (positive but noisy)
  )
  
  ses <- c(
    0.28, 0.30, 0.32, 0.33,  # Pre-periods
    0,                       # t=-1
    0.35, 0.34, 0.36, 0.37, 0.38, 0.39, 0.37, 0.40, 0.38, 0.39, 0.38  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_GDP.growth..annual.....png"
  data <- create_somewhat_messy_plot("GDP growth (annual %)", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.18))
}

# 3. GDP per capita growth (annual %) - Somewhat messy
create_gdp_pc_growth_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild pre-trend, discontinuity, mixed post-effects
  coefficients <- c(
    0.08, 0.12, 0.18, 0.22,  # t=-5 to t=-2 (mild upward)
    0,                        # t=-1
    0.28, 0.35, 0.32, 0.38, 0.40, 0.35, 0.42, 0.38, 0.40, 0.45, 0.38  # t=0 to t=10 (positive, somewhat noisy)
  )
  
  ses <- c(
    0.15, 0.16, 0.17, 0.18,  # Pre-periods
    0,                       # t=-1
    0.20, 0.19, 0.21, 0.20, 0.22, 0.21, 0.22, 0.23, 0.21, 0.22, 0.21  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_GDP.per.capita.growth..annual.....png"
  data <- create_somewhat_messy_plot("GDP per capita growth (annual %)", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.15))
}

# 4. GNI per capita (constant 2015 US$) - Somewhat messy
create_gni_pc_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild pre-trend, discontinuity, positive but noisy post-effects
  coefficients <- c(
    -120, -85, -50, -25,  # t=-5 to t=-2 (mild upward trend in levels)
    0,                     # t=-1
    35, 68, 92, 115, 138, 125, 142, 158, 145, 162, 150  # t=0 to t=10 (positive but noisy)
  )
  
  ses <- c(
    95, 98, 102, 105,  # Pre-periods
    0,                 # t=-1
    110, 108, 115, 112, 118, 120, 115, 122, 118, 120, 118  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_GNI.per.capita..constant.2015.US.....png"
  data <- create_somewhat_messy_plot("GNI per capita (constant 2015 US$)", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.14))
}

# 5. GNI per capita growth (annual %) - Somewhat messy
create_gni_pc_growth_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: slight pre-trend, discontinuity, mixed post-effects
  coefficients <- c(
    0.10, 0.18, 0.25, 0.32,  # t=-5 to t=-2 (slight upward)
    0,                        # t=-1
    0.42, 0.50, 0.45, 0.55, 0.58, 0.52, 0.60, 0.55, 0.58, 0.62, 0.55  # t=0 to t=10 (positive, noisy)
  )
  
  ses <- c(
    0.22, 0.24, 0.25, 0.26,  # Pre-periods
    0,                       # t=-1
    0.28, 0.27, 0.29, 0.28, 0.30, 0.29, 0.30, 0.31, 0.29, 0.30, 0.29  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_GNI.per.capita.growth..annual.....png"
  data <- create_somewhat_messy_plot("GNI per capita growth (annual %)", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.16))
}

# ============================================================================
# DISTRIBUTIONAL OUTCOMES - Somewhat Messy Examples
# ============================================================================

# 6. Gini.index - Somewhat messy
create_gini_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild downward pre-trend (inequality decreasing before), discontinuity, mixed post-effects
  coefficients <- c(
    -1.8, -1.2, -0.6, -0.3,  # t=-5 to t=-2 (mild downward - inequality decreasing)
    0,                        # t=-1
    0.4, 0.8, 0.5, 1.1, 1.3, 0.9, 1.2, 1.0, 1.1, 1.4, 0.95  # t=0 to t=10 (positive - inequality increasing, noisy)
  )
  
  ses <- c(
    1.5, 1.6, 1.4, 1.5,  # Pre-periods
    0,                   # t=-1
    1.1, 1.0, 1.2, 1.1, 1.3, 1.2, 1.3, 1.1, 1.2, 1.3, 1.1  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_Gini.index.png"
  data <- create_somewhat_messy_plot("Gini.index", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.11))
}

# 7. Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population) - Somewhat messy
create_poverty_365_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: slight upward pre-trend (poverty increasing), discontinuity, mixed post-effects
  coefficients <- c(
    0.8, 1.2, 1.6, 1.9,  # t=-5 to t=-2 (slight upward - poverty increasing)
    0,                    # t=-1
    -1.5, -2.2, -1.8, -2.8, -3.2, -2.5, -3.0, -2.7, -3.1, -3.5, -2.9  # t=0 to t=10 (negative - poverty decreasing, noisy)
  )
  
  ses <- c(
    0.9, 1.0, 1.1, 1.2,  # Pre-periods
    0,                   # t=-1
    1.3, 1.2, 1.4, 1.3, 1.5, 1.4, 1.5, 1.3, 1.4, 1.5, 1.3  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_Poverty.headcount.ratio.at..3.65.a.day..2017.PPP.....of.population..png"
  data <- create_somewhat_messy_plot("Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population)", 
                                     coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.13))
}

# 8. Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population) - Somewhat messy
create_poverty_685_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild pre-trend, discontinuity, mixed post-effects
  coefficients <- c(
    0.6, 0.9, 1.3, 1.5,  # t=-5 to t=-2 (mild upward)
    0,                    # t=-1
    -1.2, -1.8, -1.5, -2.2, -2.5, -2.0, -2.4, -2.1, -2.3, -2.6, -2.2  # t=0 to t=10 (negative, noisy)
  )
  
  ses <- c(
    0.8, 0.9, 1.0, 1.1,  # Pre-periods
    0,                   # t=-1
    1.2, 1.1, 1.3, 1.2, 1.4, 1.3, 1.4, 1.2, 1.3, 1.4, 1.2  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_Poverty.headcount.ratio.at..6.85.a.day..2017.PPP.....of.population..png"
  data <- create_somewhat_messy_plot("Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population)", 
                                     coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.17))
}

# 9. Income share held by highest 10% - Somewhat messy
create_income_high10_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild pre-trend, discontinuity, mixed post-effects
  coefficients <- c(
    -0.4, -0.25, -0.15, -0.08,  # t=-5 to t=-2 (mild upward - top 10% share increasing)
    0,                           # t=-1
    0.12, 0.25, 0.18, 0.32, 0.38, 0.28, 0.35, 0.30, 0.33, 0.40, 0.32  # t=0 to t=10 (positive, noisy)
  )
  
  ses <- c(
    0.35, 0.36, 0.37, 0.38,  # Pre-periods
    0,                        # t=-1
    0.40, 0.39, 0.41, 0.40, 0.42, 0.41, 0.42, 0.40, 0.41, 0.42, 0.40  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_Income.share.held.by.highest.10..png"
  data <- create_somewhat_messy_plot("Income share held by highest 10%", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.19))
}

# 10. Income share held by lowest 10% - Somewhat messy
create_income_low10_example <- function() {
  event_times <- -5:10
  
  # Somewhat messy: mild pre-trend, discontinuity, mixed post-effects
  coefficients <- c(
    0.08, 0.05, 0.03, 0.01,  # t=-5 to t=-2 (mild downward - bottom 10% share decreasing)
    0,                        # t=-1
    -0.02, -0.05, -0.03, -0.08, -0.10, -0.06, -0.09, -0.07, -0.08, -0.11, -0.08  # t=0 to t=10 (negative, noisy)
  )
  
  ses <- c(
    0.04, 0.045, 0.05, 0.055,  # Pre-periods
    0,                          # t=-1
    0.06, 0.055, 0.065, 0.06, 0.07, 0.065, 0.07, 0.06, 0.065, 0.07, 0.06  # Post-periods
  )
  
  output_path <- "Theoretical_Results_Examples/Messy_Results/event_study_Income.share.held.by.lowest.10..png"
  data <- create_somewhat_messy_plot("Income share held by lowest 10%", coefficients, ses, output_path)
  
  post_coefs <- coefficients[6:16]
  ate <- mean(post_coefs, na.rm = TRUE)
  ate_se <- sqrt(mean(ses[6:16]^2, na.rm = TRUE))
  
  return(list(data = data, ate = ate, ate_se = ate_se, pre_trend_p = 0.14))
}

# ============================================================================
# CREATE RESULTS TABLES
# ============================================================================

create_results_tables <- function(results_list) {
  # Aggregate outcomes table
  aggregate_table <- data.frame(
    Outcome = c("log_GDP_per_capita", 
                "GDP growth (annual %)", 
                "GDP per capita growth (annual %)",
                "GNI per capita (constant 2015 US$)",
                "GNI per capita growth (annual %)"),
    ATE = c(results_list$log_gdp_pc$ate,
            results_list$gdp_growth$ate,
            results_list$gdp_pc_growth$ate,
            results_list$gni_pc$ate,
            results_list$gni_pc_growth$ate),
    SE = c(results_list$log_gdp_pc$ate_se,
           results_list$gdp_growth$ate_se,
           results_list$gdp_pc_growth$ate_se,
           results_list$gni_pc$ate_se,
           results_list$gni_pc_growth$ate_se),
    CI_Lower = c(results_list$log_gdp_pc$ate - 1.96 * results_list$log_gdp_pc$ate_se,
                 results_list$gdp_growth$ate - 1.96 * results_list$gdp_growth$ate_se,
                 results_list$gdp_pc_growth$ate - 1.96 * results_list$gdp_pc_growth$ate_se,
                 results_list$gni_pc$ate - 1.96 * results_list$gni_pc$ate_se,
                 results_list$gni_pc_growth$ate - 1.96 * results_list$gni_pc_growth$ate_se),
    CI_Upper = c(results_list$log_gdp_pc$ate + 1.96 * results_list$log_gdp_pc$ate_se,
                 results_list$gdp_growth$ate + 1.96 * results_list$gdp_growth$ate_se,
                 results_list$gdp_pc_growth$ate + 1.96 * results_list$gdp_pc_growth$ate_se,
                 results_list$gni_pc$ate + 1.96 * results_list$gni_pc$ate_se,
                 results_list$gni_pc_growth$ate + 1.96 * results_list$gni_pc_growth$ate_se),
    Significant = c(ifelse(abs(results_list$log_gdp_pc$ate) > 1.96 * results_list$log_gdp_pc$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$gdp_growth$ate) > 1.96 * results_list$gdp_growth$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$gdp_pc_growth$ate) > 1.96 * results_list$gdp_pc_growth$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$gni_pc$ate) > 1.96 * results_list$gni_pc$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$gni_pc_growth$ate) > 1.96 * results_list$gni_pc_growth$ate_se, "Yes", "No")),
    Pre_Trend_P_Value = c(results_list$log_gdp_pc$pre_trend_p,
                          results_list$gdp_growth$pre_trend_p,
                          results_list$gdp_pc_growth$pre_trend_p,
                          results_list$gni_pc$pre_trend_p,
                          results_list$gni_pc_growth$pre_trend_p),
    stringsAsFactors = FALSE
  )
  
  # Distributional outcomes table
  distributional_table <- data.frame(
    Outcome = c("Gini.index",
                "Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population)",
                "Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population)",
                "Income share held by highest 10%",
                "Income share held by lowest 10%"),
    ATE = c(results_list$gini$ate,
            results_list$poverty_365$ate,
            results_list$poverty_685$ate,
            results_list$income_high10$ate,
            results_list$income_low10$ate),
    SE = c(results_list$gini$ate_se,
           results_list$poverty_365$ate_se,
           results_list$poverty_685$ate_se,
           results_list$income_high10$ate_se,
           results_list$income_low10$ate_se),
    CI_Lower = c(results_list$gini$ate - 1.96 * results_list$gini$ate_se,
                 results_list$poverty_365$ate - 1.96 * results_list$poverty_365$ate_se,
                 results_list$poverty_685$ate - 1.96 * results_list$poverty_685$ate_se,
                 results_list$income_high10$ate - 1.96 * results_list$income_high10$ate_se,
                 results_list$income_low10$ate - 1.96 * results_list$income_low10$ate_se),
    CI_Upper = c(results_list$gini$ate + 1.96 * results_list$gini$ate_se,
                 results_list$poverty_365$ate + 1.96 * results_list$poverty_365$ate_se,
                 results_list$poverty_685$ate + 1.96 * results_list$poverty_685$ate_se,
                 results_list$income_high10$ate + 1.96 * results_list$income_high10$ate_se,
                 results_list$income_low10$ate + 1.96 * results_list$income_low10$ate_se),
    Significant = c(ifelse(abs(results_list$gini$ate) > 1.96 * results_list$gini$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$poverty_365$ate) > 1.96 * results_list$poverty_365$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$poverty_685$ate) > 1.96 * results_list$poverty_685$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$income_high10$ate) > 1.96 * results_list$income_high10$ate_se, "Yes", "No"),
                    ifelse(abs(results_list$income_low10$ate) > 1.96 * results_list$income_low10$ate_se, "Yes", "No")),
    Pre_Trend_P_Value = c(results_list$gini$pre_trend_p,
                          results_list$poverty_365$pre_trend_p,
                          results_list$poverty_685$pre_trend_p,
                          results_list$income_high10$pre_trend_p,
                          results_list$income_low10$pre_trend_p),
    stringsAsFactors = FALSE
  )
  
  write.csv(aggregate_table, "Theoretical_Results_Examples/Tables/aggregate_results_theoretical.csv", row.names = FALSE)
  write.csv(distributional_table, "Theoretical_Results_Examples/Tables/distributional_results_theoretical.csv", row.names = FALSE)
  
  return(list(aggregate = aggregate_table, distributional = distributional_table))
}

# ============================================================================
# RUN ALL
# ============================================================================

cat("Creating somewhat messy theoretical examples with actual outcome labels...\n\n")

# Create aggregate outcome examples
cat("Creating aggregate outcome examples...\n")
results_list <- list()
results_list$log_gdp_pc <- create_log_gdp_pc_example()
cat("  Created log_GDP_per_capita\n")

results_list$gdp_growth <- create_gdp_growth_example()
cat("  Created GDP growth (annual %)\n")

results_list$gdp_pc_growth <- create_gdp_pc_growth_example()
cat("  Created GDP per capita growth (annual %)\n")

results_list$gni_pc <- create_gni_pc_example()
cat("  Created GNI per capita (constant 2015 US$)\n")

results_list$gni_pc_growth <- create_gni_pc_growth_example()
cat("  Created GNI per capita growth (annual %)\n\n")

# Create distributional outcome examples
cat("Creating distributional outcome examples...\n")
results_list$gini <- create_gini_example()
cat("  Created Gini.index\n")

results_list$poverty_365 <- create_poverty_365_example()
cat("  Created Poverty headcount ratio at $3.65 a day\n")

results_list$poverty_685 <- create_poverty_685_example()
cat("  Created Poverty headcount ratio at $6.85 a day\n")

results_list$income_high10 <- create_income_high10_example()
cat("  Created Income share held by highest 10%\n")

results_list$income_low10 <- create_income_low10_example()
cat("  Created Income share held by lowest 10%\n\n")

# Create results tables
cat("Creating results tables...\n")
tables <- create_results_tables(results_list)
cat("  Created aggregate_results_theoretical.csv\n")
cat("  Created distributional_results_theoretical.csv\n\n")

cat("=== All theoretical examples created ===\n")
cat("Files saved in: Theoretical_Results_Examples/Messy_Results/ and Tables/\n")
cat("\nNote: These are 'somewhat messy' examples showing:\n")
cat("  - Mild pre-trends (not perfectly flat)\n")
cat("  - Some discontinuity at t=-1\n")
cat("  - Noisy but somewhat patterned post-effects\n")
cat("  - Mixed statistical significance\n")
