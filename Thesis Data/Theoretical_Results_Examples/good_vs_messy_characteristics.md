# Good vs. Messy Results: Detailed Characteristics

## GOOD RESULTS - What to Look For

### Event Study Plot Characteristics:
1. **Pre-Period (t=-5 to t=-2)**:
   - Coefficients close to 0
   - Confidence intervals include 0 (not significant)
   - Flat, horizontal pattern
   - No clear trend

2. **Reference Period (t=-1)**:
   - Normalized to 0 (or close to 0 after normalization)
   - Smooth transition from pre to post periods
   - No large discontinuity

3. **Post-Period (t=0 to t=10)**:
   - Clear, consistent pattern
   - Most coefficients statistically significant
   - Confidence intervals don't include 0
   - Effect builds or stabilizes over time

### Statistical Indicators:
- **Pre-trend test p-value**: > 0.10 (not significant)
- **Within R²**: > 0.01 (event study explains meaningful variation)
- **Post-treatment effects**: Most are significant (p < 0.05)
- **Confidence intervals**: Narrow, don't include 0

### Example Pattern:
```
Pre:  [0.01] [0.02] [-0.01] [0.01] → [0] → [0.05] [0.08] [0.10] [0.12] [0.15] ...
      (flat)                        (ref)  (clear positive effect)
```

---

## MESSY RESULTS - What to Avoid

### Event Study Plot Characteristics:
1. **Pre-Period (t=-5 to t=-2)**:
   - Coefficients significantly different from 0
   - Clear upward or downward trend
   - Confidence intervals don't include 0 (significant pre-trends)
   - Violates parallel trends assumption

2. **Reference Period (t=-1)**:
   - Large discontinuity/jump
   - Artificial break in the line
   - Visual "kink" at t=-1

3. **Post-Period (t=0 to t=10)**:
   - Noisy, inconsistent pattern
   - Alternating positive/negative effects
   - Many coefficients not significant
   - Wide confidence intervals that include 0

### Statistical Indicators:
- **Pre-trend test p-value**: < 0.10 (significant pre-trends)
- **Within R²**: < 0.01 (event study explains little variation)
- **Post-treatment effects**: Most not significant (p > 0.05)
- **Confidence intervals**: Wide, many include 0

### Example Pattern:
```
Pre:  [-0.08] [-0.06] [-0.04] [-0.02] → [0] → [0.02] [-0.01] [0.03] [0.01] [-0.02] ...
      (downward trend)                  (jump) (noisy, inconsistent)
```

---

## Your Current Results Assessment

Based on your actual results:

### Overall Assessment: **MESSY**
- ❌ Pre-trends present (t=-5, t=-4 marginally significant)
- ❌ Very low within R² (0.000455)
- ❌ Most aggregate outcomes not significant
- ⚠️ Regional heterogeneity is clearer than overall effects

### Post-Soviet: **MIXED** (Clear patterns but negative growth)
- ✅ Strong, significant poverty/inequality reduction
- ❌ Negative growth effects
- ⚠️ Trade-off between distribution and growth

### Latin America: **MIXED** (Clear patterns but increased poverty)
- ✅ Strong, significant growth effects
- ❌ Increased poverty/inequality
- ⚠️ Trade-off between growth and distribution

### Recommendation:
Focus on **regional heterogeneity** as the main finding rather than overall average effects.


