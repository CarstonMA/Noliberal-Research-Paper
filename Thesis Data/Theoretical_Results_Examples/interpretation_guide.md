# Event Study Interpretation Guide

## How to Read Event Study Plots

### X-Axis: Years Relative to Treatment
- **Negative values (t=-5 to t=-2)**: Pre-treatment periods
- **t=-1**: Reference period (omitted from regression, normalized)
- **Zero and positive (t=0 to t=10)**: Post-treatment periods

### Y-Axis: Coefficient Estimate
- **Value**: Effect relative to reference period (t=-1)
- **Zero line**: No effect (baseline)
- **Positive values**: Increase in outcome
- **Negative values**: Decrease in outcome

### Visual Elements:
- **Blue line**: Point estimates
- **Shaded area**: 95% confidence interval
- **Red circle at t=-1**: Reference period
- **Red vertical line**: Treatment boundary

---

## What Good Results Look Like

### Pattern 1: Positive Treatment Effect
```
Pre:  [≈0] [≈0] [≈0] [≈0] → [0] → [+0.05] [+0.08] [+0.10] [+0.12] ...
      (flat)                (ref)  (clear positive effect)
```
**Interpretation**: No pre-trends, clear positive effect after treatment

### Pattern 2: Negative Treatment Effect
```
Pre:  [≈0] [≈0] [≈0] [≈0] → [0] → [-0.05] [-0.08] [-0.10] [-0.12] ...
      (flat)                (ref)  (clear negative effect)
```
**Interpretation**: No pre-trends, clear negative effect after treatment

### Pattern 3: Delayed Effect
```
Pre:  [≈0] [≈0] [≈0] [≈0] → [0] → [≈0] [≈0] [+0.05] [+0.10] [+0.15] ...
      (flat)                (ref)  (no immediate effect, builds over time)
```
**Interpretation**: Effect takes time to materialize

---

## What Messy Results Look Like

### Pattern 1: Pre-Trends (Most Common Problem)
```
Pre:  [-0.08] [-0.06] [-0.04] [-0.02] → [0] → [0.02] [0.03] ...
      (downward trend)                  (jump) (effect unclear)
```
**Problem**: Treated group already declining before treatment
**Solution**: Different treatment definition or acknowledge limitation

### Pattern 2: Large Discontinuity at t=-1
```
Pre:  [4000] [4200] [4200] [4400] → [0] → [7000] [6500] ...
      (high values)                  (huge jump) (post effects)
```
**Problem**: Normalization artifact - pre-periods were high, forcing t=-1 to 0 creates jump
**Solution**: Normalize to average of pre-periods (already implemented)

### Pattern 3: Noisy Post-Effects
```
Pre:  [≈0] [≈0] [≈0] [≈0] → [0] → [+0.02] [-0.01] [+0.03] [-0.02] [+0.01] ...
      (flat)                (ref)  (alternating, inconsistent)
```
**Problem**: No clear pattern, effects inconsistent
**Solution**: May indicate weak/null effects or need better specification

### Pattern 4: Wide Confidence Intervals
```
Pre:  [≈0 ± 0.05] [≈0 ± 0.05] → [0] → [+0.02 ± 0.10] [+0.03 ± 0.12] ...
      (wide CIs)                  (ref)  (very wide CIs, include 0)
```
**Problem**: Imprecise estimates, can't distinguish from zero
**Solution**: Need more data or different specification

---

## Checklist for Evaluating Your Results

### Pre-Treatment Period (t=-5 to t=-2):
- [ ] Are coefficients close to 0?
- [ ] Are confidence intervals wide (include 0)?
- [ ] Is there a clear trend (upward or downward)?
- [ ] Are any pre-period coefficients significant?

**Good**: All close to 0, wide CIs, no trend, none significant
**Messy**: Non-zero, narrow CIs, clear trend, some significant

### Reference Period (t=-1):
- [ ] Is there a large discontinuity/jump?
- [ ] Does the line connect smoothly?
- [ ] Is t=-1 close to the pre-period average?

**Good**: Small discontinuity, smooth line, close to pre-period average
**Messy**: Large jump, broken line, far from pre-period average

### Post-Treatment Period (t=0 to t=10):
- [ ] Is there a clear pattern?
- [ ] Are most coefficients significant?
- [ ] Do confidence intervals exclude 0?
- [ ] Is the effect consistent over time?

**Good**: Clear pattern, most significant, CIs exclude 0, consistent
**Messy**: Noisy, few significant, CIs include 0, inconsistent

### Statistical Tests:
- [ ] Pre-trend test p-value > 0.10?
- [ ] Within R² > 0.01?
- [ ] Post-treatment effects significant?

**Good**: p > 0.10, R² > 0.01, effects significant
**Messy**: p < 0.10, R² < 0.01, effects not significant

---

## Your Results: What They Mean

Based on your analysis:

### Overall Results: **MESSY**
- Pre-trends present (violates parallel trends)
- Low within R² (weak identification)
- Most effects not significant
- **Conclusion**: Overall effects are weak/null, focus on regional differences

### Regional Results: **CLEARER**
- Post-Soviet: Clear pattern (poverty↓, inequality↓, growth↓)
- Latin America: Clear pattern (growth↑, poverty↑, inequality↑)
- **Conclusion**: Regional heterogeneity is the main finding

### Recommendation:
1. **Acknowledge** overall effects are weak
2. **Emphasize** regional heterogeneity as main finding
3. **Document** pre-trend problems
4. **Interpret** regional results cautiously given pre-trends


