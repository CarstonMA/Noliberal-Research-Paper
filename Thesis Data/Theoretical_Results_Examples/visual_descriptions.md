# Visual Descriptions of Good vs. Messy Event Study Plots

## To Generate the Actual Plots

Run the R script `create_theoretical_examples.R` in R or RStudio. It will create:
- `Good_Results/good_event_study_example.png`
- `Messy_Results/messy_event_study_example.png`
- `comparison_good_vs_messy.png`

---

## GOOD RESULT - Visual Description

### Plot Appearance:
- **X-axis**: Years from -5 to 10 (relative to treatment)
- **Y-axis**: Coefficient estimates
- **Pre-period (t=-5 to t=-2)**: 
  - Points clustered around 0 (e.g., -0.01, 0.02, -0.01, 0.01)
  - Confidence intervals wide and include 0
  - Blue line is flat, horizontal
- **Reference period (t=-1)**:
  - Red circle at (x=-1, y=0)
  - Line connects smoothly through this point
  - No visible jump or discontinuity
- **Post-period (t=0 to t=10)**:
  - Points show clear upward trend (e.g., 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.25, 0.25)
  - Confidence intervals narrow and mostly exclude 0
  - Blue line shows smooth, consistent increase
  - Shaded area (confidence band) is narrow

### Key Visual Features:
✅ Smooth, continuous line from pre to post
✅ No break or jump at t=-1
✅ Clear visual separation between flat pre-period and rising post-period
✅ Most post-period points clearly above zero line
✅ Narrow confidence intervals in post-period

---

## MESSY RESULT - Visual Description

### Plot Appearance:
- **X-axis**: Years from -5 to 10 (relative to treatment)
- **Y-axis**: Coefficient estimates
- **Pre-period (t=-5 to t=-2)**:
  - Points show clear downward trend (e.g., -0.08, -0.06, -0.04, -0.02)
  - Confidence intervals may not include 0 (significant pre-trends)
  - Orange line slopes downward
- **Reference period (t=-1)**:
  - Red circle at (x=-1, y=0)
  - **Large visible jump** from t=-2 (at -0.02) to t=-1 (at 0)
  - Line appears broken or has a sharp kink
- **Post-period (t=0 to t=10)**:
  - Points are noisy and inconsistent (e.g., 0.02, -0.01, 0.03, 0.01, -0.02, 0.04, 0.02, -0.01, 0.03, 0.01, 0.02)
  - Confidence intervals wide and many include 0
  - Orange line zigzags, no clear pattern
  - Shaded area (confidence band) is very wide

### Key Visual Features:
❌ Clear downward trend in pre-period (violates parallel trends)
❌ Large discontinuity/jump at t=-1 (normalization artifact)
❌ No clear pattern in post-period (noisy, inconsistent)
❌ Many post-period confidence intervals include 0
❌ Wide confidence intervals throughout

---

## Comparison Plot - Side by Side

When viewed side by side:
- **Left panel (Good)**: Clean, smooth pattern with flat pre-trends and clear post-effect
- **Right panel (Messy)**: Messy pattern with pre-trends, discontinuity, and noisy post-effects

The contrast makes it immediately clear which results are interpretable and which have problems.

---

## Numerical Values for Reference

### Good Result Coefficients:
```
t=-5:  0.01  (SE: 0.02, CI: [-0.03, 0.05])
t=-4:  0.02  (SE: 0.02, CI: [-0.02, 0.06])
t=-3: -0.01  (SE: 0.02, CI: [-0.05, 0.03])
t=-2:  0.01  (SE: 0.02, CI: [-0.03, 0.05])
t=-1:  0.00  (Reference, normalized)
t=0:   0.05  (SE: 0.03, CI: [-0.01, 0.11])
t=1:   0.08  (SE: 0.03, CI: [0.02, 0.14])
t=2:   0.10  (SE: 0.03, CI: [0.04, 0.16])
... (continues upward)
```

### Messy Result Coefficients:
```
t=-5: -0.08  (SE: 0.03, CI: [-0.14, -0.02])  *significant*
t=-4: -0.06  (SE: 0.03, CI: [-0.12, 0.00])   *significant*
t=-3: -0.04  (SE: 0.03, CI: [-0.10, 0.02])
t=-2: -0.02  (SE: 0.03, CI: [-0.08, 0.04])
t=-1:  0.00  (Reference, normalized) ← LARGE JUMP
t=0:   0.02  (SE: 0.04, CI: [-0.06, 0.10])
t=1:  -0.01  (SE: 0.04, CI: [-0.09, 0.07])
t=2:   0.03  (SE: 0.04, CI: [-0.05, 0.11])
... (noisy, inconsistent)
```

---

## How to Use These Examples

1. **Compare your actual plots** to these descriptions
2. **Identify which pattern** your results match (good or messy)
3. **Use the interpretation guides** to understand what your results mean
4. **Refer to the tables** to see what good vs. messy statistical results look like
5. **Check the country examples** to see which countries might show cleaner patterns


