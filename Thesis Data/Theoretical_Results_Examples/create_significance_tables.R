# Create Detailed Significance Tables
# Shows p-values, t-statistics, and significance levels for theoretical results

# Read the existing results
aggregate_results <- read.csv("Theoretical_Results_Examples/Tables/aggregate_results_theoretical.csv", stringsAsFactors = FALSE)
distributional_results <- read.csv("Theoretical_Results_Examples/Tables/distributional_results_theoretical.csv", stringsAsFactors = FALSE)

# Function to calculate t-statistic and p-value
calculate_significance <- function(ate, se) {
  if (is.na(ate) || is.na(se) || se == 0) {
    return(list(t_stat = NA, p_value = NA, sig_level = ""))
  }
  
  t_stat <- abs(ate / se)
  p_value <- 2 * (1 - pnorm(t_stat))  # Two-tailed test
  
  sig_level <- ""
  if (p_value < 0.001) {
    sig_level <- "***"
  } else if (p_value < 0.01) {
    sig_level <- "**"
  } else if (p_value < 0.05) {
    sig_level <- "*"
  } else if (p_value < 0.10) {
    sig_level <- "."
  }
  
  return(list(t_stat = t_stat, p_value = p_value, sig_level = sig_level))
}

# Add significance columns to aggregate results
aggregate_results$t_statistic <- NA
aggregate_results$p_value <- NA
aggregate_results$Significance_Level <- ""

for (i in 1:nrow(aggregate_results)) {
  sig <- calculate_significance(aggregate_results$ATE[i], aggregate_results$SE[i])
  aggregate_results$t_statistic[i] <- round(sig$t_stat, 3)
  aggregate_results$p_value[i] <- round(sig$p_value, 3)
  aggregate_results$Significance_Level[i] <- sig$sig_level
}

# Add interpretation
aggregate_results$Interpretation <- ifelse(
  aggregate_results$p_value < 0.05,
  "Significant at 5% level",
  ifelse(
    aggregate_results$p_value < 0.10,
    "Marginally significant (p<0.10)",
    "Not significant - CI includes 0"
  )
)

# Reorder columns
aggregate_results <- aggregate_results[, c("Outcome", "ATE", "SE", "t_statistic", "p_value", 
                                           "Significance_Level", "CI_Lower", "CI_Upper", 
                                           "Pre_Trend_P_Value", "Interpretation")]

# Round numeric columns
aggregate_results$ATE <- round(aggregate_results$ATE, 4)
aggregate_results$SE <- round(aggregate_results$SE, 4)
aggregate_results$CI_Lower <- round(aggregate_results$CI_Lower, 4)
aggregate_results$CI_Upper <- round(aggregate_results$CI_Upper, 4)

# Same for distributional results
distributional_results$t_statistic <- NA
distributional_results$p_value <- NA
distributional_results$Significance_Level <- ""

for (i in 1:nrow(distributional_results)) {
  sig <- calculate_significance(distributional_results$ATE[i], distributional_results$SE[i])
  distributional_results$t_statistic[i] <- round(sig$t_stat, 3)
  distributional_results$p_value[i] <- round(sig$p_value, 3)
  distributional_results$Significance_Level[i] <- sig$sig_level
}

distributional_results$Interpretation <- ifelse(
  distributional_results$p_value < 0.05,
  "Significant at 5% level",
  ifelse(
    distributional_results$p_value < 0.10,
    "Marginally significant (p<0.10)",
    "Not significant - CI includes 0"
  )
)

distributional_results <- distributional_results[, c("Outcome", "ATE", "SE", "t_statistic", "p_value", 
                                                    "Significance_Level", "CI_Lower", "CI_Upper", 
                                                    "Pre_Trend_P_Value", "Interpretation")]

distributional_results$ATE <- round(distributional_results$ATE, 4)
distributional_results$SE <- round(distributional_results$SE, 4)
distributional_results$CI_Lower <- round(distributional_results$CI_Lower, 4)
distributional_results$CI_Upper <- round(distributional_results$CI_Upper, 4)

# Write enhanced tables
write.csv(aggregate_results, "Theoretical_Results_Examples/Tables/aggregate_results_significance.csv", row.names = FALSE)
write.csv(distributional_results, "Theoretical_Results_Examples/Tables/distributional_results_significance.csv", row.names = FALSE)

# Create summary table
summary_table <- rbind(
  data.frame(
    Outcome = aggregate_results$Outcome,
    ATE = aggregate_results$ATE,
    Significant_5pct = ifelse(aggregate_results$p_value < 0.05, "Yes", "No"),
    Significant_10pct = ifelse(aggregate_results$p_value < 0.10, "Yes", "No"),
    Pre_Trend_Issue = paste0("Mild (p=", aggregate_results$Pre_Trend_P_Value, ")"),
    Overall_Assessment = aggregate_results$Interpretation,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Outcome = distributional_results$Outcome,
    ATE = distributional_results$ATE,
    Significant_5pct = ifelse(distributional_results$p_value < 0.05, "Yes", "No"),
    Significant_10pct = ifelse(distributional_results$p_value < 0.10, "Yes", "No"),
    Pre_Trend_Issue = paste0("Mild (p=", distributional_results$Pre_Trend_P_Value, ")"),
    Overall_Assessment = distributional_results$Interpretation,
    stringsAsFactors = FALSE
  )
)

write.csv(summary_table, "Theoretical_Results_Examples/Tables/significance_summary.csv", row.names = FALSE)

cat("Created significance tables:\n")
cat("  - aggregate_results_significance.csv\n")
cat("  - distributional_results_significance.csv\n")
cat("  - significance_summary.csv\n")
cat("\nSignificance levels: *** p<0.001, ** p<0.01, * p<0.05, . p<0.10\n")


