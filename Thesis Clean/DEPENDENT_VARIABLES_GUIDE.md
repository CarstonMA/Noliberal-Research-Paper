# Guide to Selecting Dependent Variables for Neoliberal Policy Impact Analysis

## Overview
This guide helps you narrow down which outcome variables to use in your analysis of neoliberal policy impacts.

## Run the Analysis Script
First, run the analysis script to see all available variables:
```r
source("analyze_dependent_variables.R")
```

This will:
- List all dependent variables by category
- Show data availability for each variable
- Suggest a focused set of variables
- Save a CSV with full variable details

## Selection Criteria

### 1. **Theoretical Relevance**
Focus on variables that neoliberal theory predicts will be affected:

**Core Outcomes:**
- **Economic Growth**: GDP per capita, GDP growth
- **Inequality**: Gini index, income shares (top 10%, bottom 10%)
- **Poverty**: Poverty headcount ratios at different thresholds
- **Social Protection**: Coverage of social insurance/safety nets
- **Labor Markets**: Unemployment rates

**Secondary Outcomes:**
- Health expenditures
- Education enrollment/attainment
- Financial inclusion
- Environmental indicators

### 2. **Data Quality Thresholds**
Minimum requirements:
- **Data availability**: ≥30-40% non-missing values
- **Country coverage**: ≥50 countries
- **Time coverage**: ≥10 years (ideally spanning policy shift periods)
- **Temporal consistency**: Data available before and after policy shifts

### 3. **Measurement Quality**
Prefer:
- Standardized measures (e.g., PPP-adjusted poverty lines)
- Internationally comparable indicators
- Variables with clear definitions
- Avoid highly correlated redundant variables

### 4. **Policy Relevance**
Select variables that:
- Policymakers care about
- Can be influenced by neoliberal policies
- Have clear policy implications
- Are commonly used in academic literature

### 5. **Statistical Considerations**
- Sufficient variation for regression analysis
- Variables that can be log-transformed if needed
- Consider distribution (normal vs. skewed)
- Check for multicollinearity

## Recommended Variable Categories

### **Primary Outcomes (Must-Have)**
1. **GDP per capita (constant 2015 US$)** - Economic growth
2. **Gini index** - Inequality
3. **Poverty headcount ratio at $3.65 a day** - Poverty
4. **Unemployment, total (% of total labor force)** - Labor market

### **Secondary Outcomes (Choose 2-3)**
5. **Income share held by highest 10%** - Inequality (alternative/complement)
6. **Coverage of social insurance programs** - Social protection
7. **Poverty headcount ratio at national poverty lines** - Poverty (alternative)
8. **GDP growth (annual %)** - Economic growth (alternative)

### **Tertiary Outcomes (Optional, for robustness)**
9. **School enrollment, primary (% net)** - Education
10. **Proportion spending >10% on health** - Health burden
11. **Coverage of social safety net programs** - Social protection (alternative)

## Narrowing Down Strategy

### Step 1: Start with Core Outcomes
Begin with 4-5 core variables covering:
- Economic growth (1 variable)
- Inequality (1 variable)
- Poverty (1 variable)
- Social protection (1 variable)
- Labor market (1 variable)

### Step 2: Check Data Availability
Run the analysis script and filter by:
- ≥30% data availability
- ≥50 countries
- ≥10 years of data

### Step 3: Test for Redundancy
Check correlations between similar variables:
- If two poverty measures are highly correlated (r > 0.9), pick one
- If GDP per capita and GNI per capita are highly correlated, pick one
- Keep variables that measure different aspects

### Step 4: Consider Your Research Question
- **If focusing on growth**: Emphasize GDP variables
- **If focusing on distribution**: Emphasize inequality and poverty
- **If focusing on social outcomes**: Emphasize social protection and health

### Step 5: Robustness Checks
Plan to test your main results with:
- Alternative measures of the same concept
- Different poverty thresholds
- Different inequality measures

## Example: Focused Set (10-12 variables)

Based on typical neoliberal policy research:

**Economic:**
1. GDP per capita (constant 2015 US$)
2. GDP growth (annual %)

**Inequality:**
3. Gini index
4. Income share held by highest 10%
5. Income share held by lowest 10%

**Poverty:**
6. Poverty headcount ratio at $3.65 a day
7. Poverty headcount ratio at national poverty lines

**Social Protection:**
8. Coverage of social insurance programs (% of population)
9. Coverage of social safety net programs (% of population)

**Labor:**
10. Unemployment, total (% of total labor force)

**Health (optional):**
11. Proportion spending >10% on out-of-pocket health expenditure

## Next Steps

1. **Run the analysis script** to see data availability
2. **Review the categorized list** in the output
3. **Select 5-10 variables** based on:
   - Your research question
   - Data availability
   - Theoretical relevance
4. **Check correlations** between selected variables
5. **Test robustness** with alternative measures

## Files Created

- `analyze_dependent_variables.R` - Analysis script
- `dependent_variables_list.csv` - Full list with availability stats
- This guide

## Tips

- **Start small**: Begin with 5-7 core variables, add more if needed
- **Balance breadth and depth**: Cover multiple domains but not too many variables
- **Consider composite indices**: Can combine related variables (e.g., multiple poverty measures)
- **Document your choices**: Explain why you selected each variable in your paper
