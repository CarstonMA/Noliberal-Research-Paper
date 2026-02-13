# Theoretical Results Examples: Good vs. Messy Event Studies

This folder contains examples of what **good** vs. **messy** event study results look like, to help interpret your actual results.

## Folder Structure

- **Good_Results/**: Examples of clean, interpretable event study patterns
- **Messy_Results/**: Examples of problematic patterns (pre-trends, discontinuities)
- **Tables/**: Example result tables showing what to look for

## Key Characteristics of GOOD Results

### 1. Pre-Trends
- ✅ **Flat pre-trends**: Pre-period coefficients (t=-5 to t=-2) are close to 0
- ✅ **No significant pre-trends**: Pre-period coefficients are not statistically significant
- ✅ **Parallel trends satisfied**: Treated and control groups follow similar trends before treatment

### 2. Visual Pattern
- ✅ **Smooth transition**: No large discontinuity at t=-1
- ✅ **Clear post-treatment effect**: Post-period coefficients show consistent pattern
- ✅ **Confidence intervals**: Most post-period CIs don't include 0 (significant effects)

### 3. Statistical Properties
- ✅ **High within R²**: Event study variables explain meaningful variation
- ✅ **Significant effects**: Post-treatment coefficients are statistically significant
- ✅ **Pre-trend test p-value > 0.10**: Joint test of pre-period coefficients not significant

### 4. Example Countries (Hypothetical Good Cases)
- **Chile (1990)**: Clear trade liberalization, clean break
- **Poland (1990)**: Post-communist transition, clear policy change
- **Estonia (1992)**: Independence + market reforms, clear break
- **Czech Republic (1991)**: Clear financial reforms
- **Slovenia (1991)**: Clear capital account opening

## Key Characteristics of MESSY Results

### 1. Pre-Trends
- ❌ **Significant pre-trends**: Pre-period coefficients are statistically significant
- ❌ **Non-flat pre-trends**: Clear upward or downward trend before treatment
- ❌ **Parallel trends violated**: Treated group already on different trajectory

### 2. Visual Pattern
- ❌ **Large discontinuity at t=-1**: Artificial jump due to normalization
- ❌ **Noisy post-effects**: Inconsistent, alternating positive/negative effects
- ❌ **Wide confidence intervals**: Many CIs include 0 (non-significant effects)

### 3. Statistical Properties
- ❌ **Low within R²**: Event study variables explain little variation (< 0.01)
- ❌ **Non-significant effects**: Most post-treatment coefficients not significant
- ❌ **Pre-trend test p-value < 0.10**: Joint test of pre-period coefficients is significant

### 4. Example Countries (Hypothetical Messy Cases)
- **Argentina (1996)**: Gradual reforms, pre-existing trends
- **Brazil (1999)**: Multiple reforms over time, no clear break
- **Mexico (2001)**: NAFTA effects confounded with other factors
- **India (1996)**: Gradual liberalization, no clear break
- **China (various)**: Continuous reforms, no discrete treatment event

## Interpreting Your Results

### If Your Results Look Like "Good Results":
- ✅ Parallel trends assumption likely satisfied
- ✅ Treatment effects are credible
- ✅ Can interpret post-treatment coefficients as causal effects
- ✅ Results are publishable

### If Your Results Look Like "Messy Results":
- ⚠️ Parallel trends assumption likely violated
- ⚠️ Treatment effects may be confounded by pre-existing trends
- ⚠️ Need to:
  - Try different treatment definition
  - Use alternative identification strategy
  - Focus on regional heterogeneity (if that's clearer)
  - Consider that effects might be genuinely weak/null

## What to Report

### Good Results:
1. Show event study plot with flat pre-trends
2. Report pre-trend test p-value (should be > 0.10)
3. Report significant post-treatment effects
4. Interpret as causal effects

### Messy Results:
1. Acknowledge pre-trend problems
2. Report pre-trend test p-value (likely < 0.10)
3. Focus on regional heterogeneity if that's clearer
4. Consider alternative specifications
5. Be cautious in causal interpretation

## Files in This Folder

- `create_theoretical_examples.R`: Script to generate example plots
- `Tables/good_results_example.csv`: Example of good results table
- `Tables/messy_results_example.csv`: Example of messy results table
- `Tables/good_countries_example.csv`: Countries that might show good results
- `Tables/messy_countries_example.csv`: Countries that might show messy results


