# Theoretical Results Examples - Using Your Actual Outcome Labels

This folder contains **somewhat messy** theoretical examples using the **exact same outcome variable names** as your actual analysis.

## What's Included

### Plots (in `Messy_Results/` folder)
Event study plots for all your outcome variables, showing **somewhat messy** patterns:
- Mild pre-trends (not perfectly flat)
- Some discontinuity at t=-1
- Noisy but somewhat patterned post-effects
- Mixed statistical significance

**Outcomes included:**
- `log_GDP_per_capita`
- `GDP growth (annual %)`
- `GDP per capita growth (annual %)`
- `GNI per capita (constant 2015 US$)`
- `GNI per capita growth (annual %)`
- `Gini.index`
- `Poverty headcount ratio at $3.65 a day (2017 PPP) (% of population)`
- `Poverty headcount ratio at $6.85 a day (2017 PPP) (% of population)`
- `Income share held by highest 10%`
- `Income share held by lowest 10%`

### Tables (in `Tables/` folder)
- `aggregate_results_theoretical.csv` - Results table for aggregate outcomes
- `distributional_results_theoretical.csv` - Results table for distributional outcomes

## How to Generate the Plots

Run the R script:
```r
source("Theoretical_Results_Examples/create_theoretical_examples.R")
```

Or open `create_theoretical_examples.R` in RStudio and run it.

## Characteristics of These "Somewhat Messy" Examples

### Pre-Periods (t=-5 to t=-2):
- **Mild trends**: Not perfectly flat, but not extreme either
- **Some coefficients significant**: A few pre-period coefficients may be marginally significant
- **Pre-trend p-values**: Between 0.11-0.19 (borderline, not extreme)

### Reference Period (t=-1):
- **Some discontinuity**: Visible jump, but not extreme
- **Normalization artifact**: Due to normalizing to average of pre-periods

### Post-Periods (t=0 to t=10):
- **Noisy but patterned**: Effects are visible but not perfectly smooth
- **Mixed significance**: Some periods significant, some not
- **Wider confidence intervals**: Less precise than ideal

### Statistical Properties:
- **Pre-trend test p-values**: 0.11-0.19 (borderline)
- **ATE significance**: Mixed (some significant, some not)
- **Confidence intervals**: Some include 0, some don't

## Comparison to Your Actual Results

These examples are designed to be **somewhat messy** - similar to what you might see in real data:
- Not perfectly clean (like ideal theoretical examples)
- Not completely chaotic (like worst-case scenarios)
- Realistic middle ground with some problems but interpretable patterns

## Use Cases

1. **Compare your actual plots** to these theoretical examples
2. **Understand what "somewhat messy" looks like** - not perfect, but not hopeless
3. **See how to interpret** borderline pre-trends and noisy post-effects
4. **Use tables** to see what realistic result tables look like


