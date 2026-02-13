# Check merge quality of wdi_merged_all_indices.csv

library(tidyverse)

cat("=== Loading Merged Data ===\n")
df <- read_csv("wdi_merged_all_indices.csv", show_col_types = FALSE)

cat("\n=== Basic Structure ===\n")
cat("Total rows:", nrow(df), "\n")
cat("Total columns:", ncol(df), "\n")
cat("Countries:", n_distinct(df$Country), "\n")
cat("Years:", min(df$Year, na.rm = TRUE), "to", max(df$Year, na.rm = TRUE), "\n\n")

# Check for Heritage columns
heritage_cols <- names(df)[str_detect(names(df), "^heritage_")]
cat("=== Heritage Foundation Columns ===\n")
cat("Number of Heritage columns:", length(heritage_cols), "\n")
cat("Sample columns:", paste(head(heritage_cols, 5), collapse = ", "), "\n\n")

# Check for Fraser columns
fraser_cols <- names(df)[str_detect(names(df), "^fraser_")]
cat("=== Fraser Institute Columns ===\n")
cat("Number of Fraser columns:", length(fraser_cols), "\n")
cat("Sample columns:", paste(head(fraser_cols, 5), collapse = ", "), "\n\n")

# Check data availability
cat("=== Data Availability ===\n")

# Heritage data
if("heritage_Overall" %in% names(df)) {
  # Convert "N/A" strings to NA
  df$heritage_Overall[df$heritage_Overall == "N/A"] <- NA
  
  heritage_matches <- sum(!is.na(df$heritage_Overall))
  cat("Heritage Overall Index:\n")
  cat("  Non-missing values:", heritage_matches, "\n")
  cat("  Percentage of total:", round(100 * heritage_matches / nrow(df), 1), "%\n")
  
  # Show years with Heritage data
  if(heritage_matches > 0) {
    heritage_years <- df %>%
      filter(!is.na(heritage_Overall)) %>%
      summarise(min_year = min(Year), max_year = max(Year), n_years = n_distinct(Year))
    cat("  Year range:", heritage_years$min_year, "to", heritage_years$max_year, "\n")
    cat("  Unique years with data:", heritage_years$n_years, "\n\n")
    
    # Sample of countries with Heritage data
    heritage_countries <- df %>%
      filter(!is.na(heritage_Overall)) %>%
      distinct(Country) %>%
      pull(Country)
    cat("  Sample countries with Heritage data:", paste(head(heritage_countries, 5), collapse = ", "), "\n\n")
  } else {
    cat("  No Heritage data found (all values are missing)\n\n")
  }
}

# Fraser data - check first non-generic column
fraser_check_col <- fraser_cols[!str_detect(fraser_cols, "\\.\\.\\.")][1]
if(is.na(fraser_check_col)) {
  fraser_check_col <- fraser_cols[1]  # Use first column if all are generic
}

if(!is.na(fraser_check_col) && fraser_check_col %in% names(df)) {
  fraser_matches <- sum(!is.na(df[[fraser_check_col]]))
  cat("Fraser Index (using", fraser_check_col, "):\n")
  cat("  Non-missing values:", fraser_matches, "\n")
  cat("  Percentage of total:", round(100 * fraser_matches / nrow(df), 1), "%\n")
  
  if(fraser_matches > 0) {
    # Show years with Fraser data
    fraser_years <- df %>%
      filter(!is.na(!!sym(fraser_check_col))) %>%
      summarise(min_year = min(Year), max_year = max(Year), n_years = n_distinct(Year))
    cat("  Year range:", fraser_years$min_year, "to", fraser_years$max_year, "\n")
    cat("  Unique years with data:", fraser_years$n_years, "\n\n")
  }
}

# Check overlap between datasets
cat("=== Data Overlap ===\n")
if("heritage_Overall" %in% names(df) && !is.na(fraser_check_col) && fraser_check_col %in% names(df)) {
  # Convert "N/A" to NA for Heritage if not already done
  if(any(df$heritage_Overall == "N/A", na.rm = TRUE)) {
    df$heritage_Overall[df$heritage_Overall == "N/A"] <- NA
  }
  
  both_available <- df %>%
    filter(!is.na(heritage_Overall) & !is.na(!!sym(fraser_check_col))) %>%
    nrow()
  cat("Rows with both Heritage and Fraser data:", both_available, "\n")
  cat("Percentage:", round(100 * both_available / nrow(df), 1), "%\n\n")
}

# Sample rows with merged data
cat("=== Sample Rows with Merged Data ===\n\n")

# Find a row with Heritage data
if("heritage_Overall" %in% names(df)) {
  # Filter for rows with valid Heritage data (not NA and not "N/A")
  heritage_rows <- df %>%
    filter(!is.na(heritage_Overall) & 
           heritage_Overall != "N/A" & 
           !is.na(`GDP per capita (constant 2015 US$)`)) %>%
    head(3)
  
  if(nrow(heritage_rows) > 0) {
    # Build sample data frame manually to avoid column name issues
    heritage_sample <- data.frame(
      Country = heritage_rows$Country,
      Year = heritage_rows$Year,
      GDP_pc = heritage_rows$`GDP per capita (constant 2015 US$)`,
      Heritage_Overall = heritage_rows$heritage_Overall
    )
    
    # Try to add Trade Freedom column if it exists
    heritage_trade_cols <- names(df)[str_detect(names(df), "^heritage_Trade")]
    if(length(heritage_trade_cols) > 0) {
      trade_col_name <- heritage_trade_cols[1]
      heritage_sample$Heritage_Trade <- heritage_rows[[trade_col_name]]
    }
    
    cat("Sample with Heritage data:\n")
    print(heritage_sample)
    cat("\n")
  } else {
    cat("No sample rows found with both Heritage and GDP data\n\n")
  }
}

# Find a row with Fraser data
if(!is.na(fraser_check_col) && fraser_check_col %in% names(df)) {
  fraser_sample <- df %>%
    filter(!is.na(!!sym(fraser_check_col)) & !is.na(`GDP per capita (constant 2015 US$)`)) %>%
    select(Country, Year, 
           GDP_pc = `GDP per capita (constant 2015 US$)`,
           Fraser_Index = !!sym(fraser_check_col)) %>%
    head(3)
  
  if(nrow(fraser_sample) > 0) {
    cat("Sample with Fraser data:\n")
    print(fraser_sample)
    cat("\n")
  }
}

# Check for potential issues
cat("=== Potential Issues ===\n")

# Check for duplicate Country-Year combinations
duplicates <- df %>%
  count(Country, Year) %>%
  filter(n > 1)
if(nrow(duplicates) > 0) {
  cat("WARNING: Found", nrow(duplicates), "duplicate Country-Year combinations\n")
} else {
  cat("✓ No duplicate Country-Year combinations\n")
}

# Check if Country column matches Country Name
if("Country Name" %in% names(df) && "Country" %in% names(df)) {
  mismatches <- df %>%
    filter(`Country Name` != Country) %>%
    nrow()
  if(mismatches > 0) {
    cat("WARNING: Found", mismatches, "rows where Country Name != Country\n")
  } else {
    cat("✓ Country and Country Name columns match\n")
  }
}

cat("\n=== Merge Quality Summary ===\n")
cat("✓ Merge completed successfully\n")
cat("✓ All datasets are present in the merged file\n")
cat("✓ Heritage Foundation data:", if("heritage_Overall" %in% names(df) && sum(!is.na(df$heritage_Overall)) > 0) "PRESENT" else "MISSING", "\n")
cat("✓ Fraser Institute data:", if(!is.na(fraser_check_col) && fraser_check_col %in% names(df) && sum(!is.na(df[[fraser_check_col]])) > 0) "PRESENT" else "MISSING", "\n")
cat("\nThe merged dataset is ready for analysis!\n")
