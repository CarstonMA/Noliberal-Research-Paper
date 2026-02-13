# Data Merge Instructions

## Overview
This script merges all your neoliberal policy index data (Heritage Foundation and Fraser Institute) with WDI outcome data.

## Files Required
- `wdi_tidy.csv` - World Development Indicators (outcomes)
- `Heritage.csv` - Heritage Foundation Economic Freedom Index
- `Fraser.xlsx` - Fraser Institute Economic Freedom Index

## How to Run

1. Open R or RStudio
2. Set working directory to this folder
3. Run:
```r
source("merge_data.R")
```

## What It Does

1. **Loads all three datasets**
2. **Standardizes country names** for matching (handles different naming conventions)
3. **Merges by country and year**:
   - WDI data (base dataset with outcomes)
   - Heritage Foundation indices (prefixed with `heritage_`)
   - Fraser Institute indices (prefixed with `fraser_`)
4. **Checks merge quality** and reports:
   - Number of matched observations
   - Country overlap between datasets
   - Year overlap between datasets
5. **Saves merged dataset** as `wdi_merged_all_indices.csv`

## Output

The merged file `wdi_merged_all_indices.csv` contains:
- All original WDI variables (GDP, poverty, inequality, etc.)
- Heritage Foundation variables (prefixed with `heritage_`)
- Fraser Institute variables (prefixed with `fraser_`)
- All matched by `Country` and `Year`

## Notes

- Uses **left join**: keeps all WDI observations, adds index data where available
- Country name matching handles variations (e.g., "United States" vs "united states")
- Missing matches are expected (different countries/years in each dataset)
- The script will print diagnostic information about the merge quality

## Next Steps

After merging, you can:
- Use the merged dataset for regression analysis
- Compare Heritage vs Fraser indices
- Analyze relationships between policy indices and outcomes
- Create treatment variables based on policy shifts
