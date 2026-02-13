# Quick diagnostic to check if treatment countries meet strict criteria
library(data.table)

# Load data
data <- read.csv("merged_data.csv", stringsAsFactors = FALSE)
data <- as.data.table(data)

# Load treatment list
treatments <- read.csv("treatment_countries_years.csv", stringsAsFactors = FALSE)
treatments <- as.data.table(treatments)

# Check Heritage_Score for a sample of treated countries
cat("Checking treatment criteria for sample countries...\n\n")

# Get Heritage_Score data
if ("Heritage_Score" %in% names(data)) {
  # Calculate thresholds
  heritage_data <- data[!is.na(Heritage_Score) & Heritage_Score != "", 
                        .(ISO3, Year, Heritage_Score = as.numeric(Heritage_Score))]
  heritage_data <- heritage_data[!is.na(Heritage_Score)]
  
  q40 <- quantile(heritage_data$Heritage_Score, 0.40, na.rm = TRUE)
  q60 <- quantile(heritage_data$Heritage_Score, 0.60, na.rm = TRUE)
  sd_policy <- sd(heritage_data$Heritage_Score, na.rm = TRUE)
  
  cat("Thresholds:\n")
  cat("  Bottom 40% threshold (q40):", q40, "\n")
  cat("  Top 40% threshold (q60):", q60, "\n")
  cat("  Standard deviation:", sd_policy, "\n")
  cat("  Required jump (1.0 SD):", 1.0 * sd_policy, "\n\n")
  
  # Check a sample of countries
  sample_countries <- head(treatments, 10)
  
  for (i in 1:nrow(sample_countries)) {
    country <- sample_countries$ISO3[i]
    event_year <- sample_countries$event_year[i]
    
    country_data <- heritage_data[ISO3 == country & 
                                   Year >= (event_year - 3) & 
                                   Year <= (event_year + 1)]
    setorder(country_data, Year)
    
    if (nrow(country_data) >= 3) {
      cat("Country:", country, "Event Year:", event_year, "\n")
      
      # Get values around treatment
      prev2 <- country_data[Year == event_year - 2]$Heritage_Score
      prev1 <- country_data[Year == event_year - 1]$Heritage_Score
      curr <- country_data[Year == event_year]$Heritage_Score
      next1 <- country_data[Year == event_year + 1]$Heritage_Score
      
      if (length(prev2) > 0 && length(prev1) > 0 && length(curr) > 0) {
        jump <- curr - prev1
        prev_change <- prev1 - prev2
        
        cat("  t-2:", ifelse(length(prev2) > 0, prev2, "NA"), 
            "t-1:", prev1, 
            "t:", curr,
            "t+1:", ifelse(length(next1) > 0, next1, "NA"), "\n")
        cat("  Jump (t - t-1):", jump, 
            "Required:", 1.0 * sd_policy,
            "Meets:", jump >= 1.0 * sd_policy, "\n")
        cat("  Prev change (t-1 - t-2):", prev_change,
            "Stable (< 0.5 SD):", abs(prev_change) < 0.5 * sd_policy, "\n")
        cat("  Crosses threshold (prev <= q40 & curr >= q60):", 
            prev1 <= q40 && curr >= q60, "\n")
        if (length(next1) > 0) {
          cat("  Sustained (next >= q60):", next1 >= q60, "\n")
        }
        cat("\n")
      }
    }
  }
} else {
  cat("Heritage_Score not found in data\n")
}


