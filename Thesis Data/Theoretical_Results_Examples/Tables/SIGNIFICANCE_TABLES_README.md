# Significance Tables Guide

## Files Created

1. **`aggregate_results_significance.csv`** - Detailed significance for aggregate outcomes
2. **`distributional_results_significance.csv`** - Detailed significance for distributional outcomes  
3. **`significance_summary.csv`** - Summary table with all outcomes

## Column Descriptions

### Main Columns:
- **Outcome**: Variable name
- **ATE**: Average Treatment Effect (post-period average)
- **SE**: Standard Error
- **t_statistic**: t-statistic (ATE/SE)
- **p_value**: Two-tailed p-value
- **Significance_Level**: Significance stars
  - `***` = p < 0.001
  - `**` = p < 0.01
  - `*` = p < 0.05
  - `.` = p < 0.10
  - (blank) = p ≥ 0.10
- **CI_Lower/CI_Upper**: 95% Confidence Interval
- **Pre_Trend_P_Value**: Pre-trend test p-value
- **Interpretation**: Plain language interpretation

## Significance Levels

### Standard Levels:
- **p < 0.001**: Highly significant (`***`)
- **p < 0.01**: Very significant (`**`)
- **p < 0.05**: Significant (`*`)
- **p < 0.10**: Marginally significant (`.`)
- **p ≥ 0.10**: Not significant (no star)

### In These Examples:
Most results are **not significant at 5%** but some are **marginally significant at 10%**:
- GDP growth: p=0.052 (marginally significant)
- GDP per capita growth: p=0.077 (marginally significant)
- GNI per capita growth: p=0.064 (marginally significant)
- Poverty $3.65: p=0.054 (marginally significant)

## How to Read

1. **Check p-value**: Lower is better (more significant)
2. **Check confidence interval**: If CI includes 0, not significant
3. **Check pre-trend p-value**: Should be > 0.10 (no pre-trends)
4. **Check interpretation**: Summary of what the result means

## Example Interpretation

**GDP growth (annual %)**: 
- ATE = 0.7255
- p-value = 0.052
- Significance: `*` (marginally significant at 10%)
- Interpretation: "Marginally significant positive effect (p<0.10)"
- Pre-trend p = 0.18 (acceptable, no significant pre-trends)

This means: The treatment increases GDP growth by about 0.73 percentage points on average, and this effect is marginally statistically significant (p<0.10), but not quite significant at the conventional 5% level.


