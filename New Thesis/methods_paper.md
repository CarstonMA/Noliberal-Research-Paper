# Empirical Methods: Neoliberal Policy Shocks and Macro–Distributional Outcomes

This document describes the full estimation pipeline used to relate discrete **policy shocks** (large year-on-year moves in economic freedom sub-indices) to **macroeconomic and distributional outcomes** in a **staggered** panel (1995–2022). It is written to parallel the implementation in `R/analysis_methodology.R` and related data-cleaning scripts.

---

## 1. Conceptual framework

### 1.1 Causal question

The analysis asks whether **unusually large changes** in specific dimensions of the Fraser Institute’s *Economic Freedom of the World* (EFW) and the Heritage Foundation’s *Index of Economic Freedom* (IEF) coincide with subsequent changes in **GDP per capita growth**, **inequality** (SWIID Gini), **bottom-50% income share** (WID), and **child mortality** (WDI, log scale). The estimand is not the long-run association of *levels* of freedom with outcomes in a purely cross-sectional sense, but rather the **average treatment effect on the treated (ATT)** associated with **timing** of a sharp reform episode relative to a defined control group, under the identifying assumptions discussed below.

### 1.2 Why staggered DiD and Callaway–Sant’Anna

Classical two-way fixed effects (TWFE) with a binary treatment can **bias** estimates when treatment turns on at different calendar times (staggered adoption): already-treated units can serve as controls for later-treated units, which contaminates the estimand (see de Chaisemartin and D’Haultfoeuille, Goodman-Bacon, Sun and Abraham, among others). **Callaway and Sant’Anna (2021)** provide **group–time ATT** estimates `ATT(g,t)` for cohorts `g` (first treatment year) and periods `t`, built from **never-treated** or **not-yet-treated** comparison groups, then aggregate them (e.g. to a single **simple ATT** or **dynamic** event-study path). The implementation uses the `did` package in R (`att_gt` + `aggte`).

---

## 2. Sample and panel structure

### 2.1 Time window and country filter

- **Years:** `min_year = 1995`, `max_year = 2022`.
- **Unit:** country (ISO3), after merge of Fraser/Heritage scores with WDI, SWIID, WID, etc. (`merged_full_data.csv`).
- **Minimum series length:** countries with **fewer than 10** country-year observations in the merged panel are **dropped** (`min_obs_country = 10`) to avoid extremely sparse series driving estimation.

The resulting panel is **unbalanced**: not every country appears in every year; `allow_unbalanced_panel = TRUE` is passed to `att_gt`.

### 2.2 Effective sample size by outcome

Because outcomes come from different sources, **missingness differs by variable**. Any given `att_gt` run uses rows where the **outcome** is non-missing (and the package uses the constructed treatment and calendar year). Interpretation of “\(N\)” should be **outcome-specific**: distributional outcomes (Gini, WID, log mortality) may imply a slightly different set of country-years than GDP growth even when country counts are similar. Reporting **unique ISO3 counts** and non-missing country-years for each outcome is recommended in the main text.

---

## 3. Policy measures: Fraser and Heritage sub-indices

### 3.1 Components included

The analysis uses **17 sub-indices**, excluding aggregate headline scores:

- **Fraser (5):** Size of Government, Legal System and Property Rights, Sound Money, Freedom to Trade Internationally, Regulation (implementation uses internal column names such as `fraser_SizeGov`, …).
- **Heritage (12):** Property Rights, Government Integrity, Judicial Effectiveness, Tax Burden, Government Spending, Fiscal Health, Business Freedom, Labor Freedom, Monetary Freedom, Trade Freedom, Investment Freedom, Financial Freedom.

Each component is analyzed **separately**: treatment timing is **recomputed** for each column. This reflects the substantive view that “neoliberal” reforms may arrive in **different policy margins** (fiscal vs. trade vs. monetary) rather than as a single scalar.

### 3.2 Treatment of missing policy scores

Shock construction uses **within-country** year-on-year changes; missing scores produce missing `delta` and do not count as crossing the threshold until a valid comparison exists. This should be stated clearly so readers do not assume continuous coverage for every country-year.

---

## 4. Treatment assignment: “large jump” shocks and staggered adoption

### 4.1 First treatment year (`first_treat`)

For each component \(c\):

1. Sort by country and year; define **score** \(S_{it}\) and first difference \(\Delta S_{it} = S_{it} - S_{i,t-1}\).
2. Compute a **global** standard deviation of \(\Delta S_{it}\) across **all** country-years in the estimation sample (with `na.rm = TRUE`), denoted \(\mathrm{SD}(\Delta)\).
3. Define a **threshold** \( \tau = \texttt{sd\_mult} \times \mathrm{SD}(\Delta)\). The **main** specification uses **`sd_mult = 1.5`**; a **robustness** run uses **`sd_mult = 1.0`**, which lowers the bar for what counts as a “shock” and typically **increases** the number of treated countries per component.
4. The **first** calendar year \(t\) such that \(\Delta S_{it} > \tau\) is recorded as **`first_treat`**. If no year satisfies this, **`first_treat = 0`**, which encodes **never-treated** in the `did` package’s convention.

Thus treatment is **absorbing** in the sense of a single cohort label: the estimator uses the **first** large upward jump; subsequent jumps in later years do not redefine cohort membership in this implementation.

### 4.2 Interpretation

This rule operationalizes a **large positive reform episode** in the direction of higher reported economic freedom on that sub-index. It is **not** equivalent to a statutory reform date or an election; it is **data-driven** relative to the global distribution of annual changes. The thesis should discuss **anticipation** (effects visible at \(t=-1\)) and **measurement** in the indices.

### 4.3 Cohort overlap

Different components yield **different** `first_treat` vectors. Treated cohort sizes vary; `cohort_size_by_component_15_vs_10sd.csv` documents how many countries are treated under 1.5× vs 1.0× SD rules.

---

## 5. Outcomes

Four outcomes are used (column names in the merged data):

| Label | Description |
|--------|-------------|
| `GDP_per_capita_growth` | WDI-based growth of GDP per capita (definition as in cleaning scripts). |
| `Gini_SWIID` | Gini from SWIID (net or as merged). |
| `Income_share_bottom50_WID` | Pre-tax national income share of the bottom 50% (WID). |
| `log_under5_mort` | Natural log of under-five mortality rate (deaths per 1,000 live births, WDI); **log** stabilizes variance and avoids explosive standard errors in CS estimation. |

Each outcome is analyzed **separately** in CS: one `att_gt` per **(component × outcome)** cell after merging that outcome’s `first_treat` for the component.

---

## 6. Estimation stages (Levels 1–4)

### 6.1 Level 1 — Pooled OLS

For each **(component, outcome)** pair:

\[
Y_{it} = \beta_0 + \beta_1 \, \text{Score}_{it,c} + \varepsilon_{it},
\]

using all country-years with non-missing \(Y\) and the policy score. This is a **descriptive** association, not causal; no country or year fixed effects. Outputs: `results_level1_ols.csv` (p-values on \(\beta_1\)).

### 6.2 Level 2 — Two-way fixed effects (TWFE)

\[
Y_{it} = \alpha_i + \gamma_t + \beta \, \text{Score}_{it,c} + \varepsilon_{it},
\]

with **country** (`ISO3`) and **year** fixed effects, **cluster-robust** standard errors at the country level (`fixest::feols`, `cluster = ~ ISO3`). This partials out time-invariant country traits and common global shocks. TWFE remains vulnerable to **staggered-treatment bias** for causal interpretation of treatment effects; it is best viewed as a **bridge** specification between OLS and CS. Outputs: `results_level2_twfe.csv`.

### 6.3 Level 3 — Callaway–Sant’Anna (`att_gt` + `aggte`)

**Model (conceptual):** For each calendar year \(t\) and cohort \(g\) (first treatment year, or control),

- `yname`: outcome column.
- `tname`: `Year`.
- `idname`: **`unit_id`**, a **numeric** factor code for `ISO3` (required by `did` preprocessing; character IDs are not used directly).
- `gname`: `first_treat` (0 = never treated).
- `xformla = ~ 1`: no additional covariates (outcome and treatment are effectively adjusted via the DR framework’s structure and comparison sets).
- **`est_method = "dr"`**: **doubly robust** combination of outcome regression and propensity-type weighting (see Callaway and Sant’Anna, 2021).
- **`control_group`:** Primary **`"nevertreated"`** — units with `first_treat = 0` for the entire panel serve as the counterfactual. If `att_gt` **errors**, the script **retries** with **`"notyettreated"`**, which uses units not yet treated by \(t\) (including future-treated). The thesis must state which succeeded for the main results (or report the fallback rule honestly).
- **`allow_unbalanced_panel = TRUE`**: accommodates missing years.
- **Inference:** `bstrap = TRUE`, `biters = 199` (`bootstrap_iters`), **`clustervars = "unit_id"`** — cluster bootstrap at the **country** level.

**Aggregation:** `aggte(att, type = "simple", na.rm = TRUE)` produces a **single** pooled ATT (`overall.att`, `overall.se`) for that component–outcome run. **P-values** for the main tables are computed in script as two-sided normal tests: \(2\Phi(-|\widehat{\mathrm{ATT}}/\widehat{\mathrm{SE}}|)\).

Outputs: `results_level3_cs_att.csv` (ATT, SE, p-value per cell).

### 6.4 Level 4 — Multiple testing (Benjamini–Hochberg FDR)

Across all CS tests in a given table, **p-values are adjusted** with **`p.adjust(..., method = "fdr")`** (Benjamini–Hochberg). This controls the **false discovery rate** under independence or positive dependence among tests. **Adjusted** quantities are stored as `fdr_bh`; `reject_bh_05` flags `fdr_bh < 0.05`.

The **main** FDR-reported file for the 1.0 SD robustness run is `results_level3_cs_att_robust10sd_fdr.csv` (naming reflects the robust threshold; confirm which SD rule corresponds to your chapter’s “primary” vs “robustness” presentation).

### 6.5 “Methodological funnel” figure (Levels 1–4 side by side)

For the **flagship** illustration (one component × GDP), the script builds a **forest-style plot** that lines up **four estimates** on the same scale:

1. **Baseline OLS** — coefficient on the policy **score** \( \beta \) from Level 1 (homoskedastic standard errors).
2. **TWFE** — \( \beta \) from Level 2 with country and year fixed effects and **country-clustered** SEs.
3. **CS–DiD (unadjusted)** — Level 3 **simple ATT** with **nominal** 95% confidence interval: \( \widehat{\mathrm{ATT}} \pm z_{0.975}\,\widehat{\mathrm{SE}} \).
4. **CS–DiD (FDR / multiplicity-aware interval)** — **same point estimate and SE** as (3); the plotted interval uses a **Bonferroni** critical value \( z_{1-\alpha/(2M)} \) with \( M \) equal to the number of CS–DiD tests in the table (all component × outcome cells with valid p-values), holding \( \alpha = 0.05 \). This **does not** replace BH–FDR inference; it is a **conservative visual** for how uncertainty looks when one accounts for **many** hypotheses on the same scale. The figure subtitle reports the **BH–FDR \(q\)** for that cell from Level 4.

The **horizontal axis** lists those four models; the **vertical axis** is the point estimate (\( \beta \) or ATT). A **thick dashed line** at **zero** is drawn for reference. The component is **`best_name`**: the policy sub-index with the **smallest Level-3 CS p-value** for **GDP growth** in the main (1.5 SD) run (so the funnel matches the strongest CS signal on GDP, not necessarily the strongest OLS association).

**Output file:** `results/figures/methodological_funnel_<component>_gdp.png` (e.g. `methodological_funnel_fraser_Regulation_gdp.png`).

---

## 7. Event study, parallel trends, and anticipation

### 7.1 Dynamic aggregation

For **GDP growth only** (in the extraction block), the script also runs `aggte(att, type = "dynamic", na.rm = TRUE)`, producing **event-time** coefficients \(\widehat{\mathrm{ATT}}_e\) for leads and lags relative to treatment (`egt`). These are exported to `component_dynamic_coefs_GDP_per_capita_growth.csv`; pre-treatment rows (`event_time < 0`) are saved separately as `component_dynamic_coefs_GDP_per_capita_growth_pretrend_only.csv`.

### 7.2 Parallel trends (diagnostic)

Under the **parallel trends** assumption, **pre-treatment** coefficients (negative event times) should be **centered near zero**; formally, **null** pre-treatment effects are consistent with no systematic divergence before the shock. The script adds **95% intervals** as \(\widehat{\theta} \pm 1.96 \,\widehat{\mathrm{SE}}\) and a flag **`ci_crosses_zero`**. 

**Important nuance:** Very **distant** leads (e.g. 15–25 years before treatment) are often **imprecisely estimated**; many papers emphasize **short** pre-windows (e.g. \(t \in \{-5,\ldots,-1\}\) or \(\{-3,-2,-1\}\)). A failure at a single distant lead need not weigh equally with failure at \(t=-1\).

### 7.3 Anticipation

A large pre-treatment coefficient at **\(t=-1\)** may indicate **anticipation** (markets or policies moving before the coded shock) or misalignment between the **economic** shock and **legal** dates. The discussion section can acknowledge this without over-interpreting a single point estimate.

### 7.4 Figures (`ggdid`)

The script saves **`did::ggdid`** plots of the **dynamic** group–time ATT path for **GDP growth**, for the **lowest-p** (and optionally second-lowest-p) CS components among GDP tests:

- **`event_study_gdp_main_lowest_p.png`** — full event window returned by `aggte(..., type = "dynamic")`.
- **`event_study_gdp_appendix_second_lowest_p.png`** — same for the second-ranked component, when applicable.

For **parallel-trends** inspection over a **short** pre-window, a version **restricted in event time** to **\(t \in \{-5,\ldots,5\}\)** is saved as **`event_study_gdp_main_t5_to_t5.png`** (via `coord_cartesian(xlim = c(-5, 5))` on the main `ggdid` object). Readers can verify that **pre-treatment** periods (\(t < 0\)) have intervals that **overlap zero** in that window, supporting the **no pre-trend** diagnostic for the main GDP case.

---

## 8. Robustness: shock threshold

Lowering **`sd_mult`** from 1.5 to 1.0 **increases** the number of treated countries (see cohort summary output). Re-estimating CS at 1.0 SD checks whether results depend on defining “reform” very strictly. **GDP** top components by p-value under the 1.0 SD run are summarized in `table3_robustness_gdp_10sd_top4.csv`.

---

## 9. Distributional outcomes: estimation health

The script prints how many **(component, outcome)** CS models **fail** (missing ATT or p-value) for Gini and log mortality, and writes `table2_distributional_survivors.csv` for successfully estimated cells. This transparency matters when arguing that null results are **not** artifacts of universal non-convergence.

---

## 10. Identification: assumptions stated explicitly

1. **Parallel trends:** Between treated and control units, **counterfactual** outcome paths would have been parallel absent treatment, **conditional** on the CS setup (and, here, no extra covariates in `xformla`).
2. **No interference / SUTVA:** Outcomes in one country are not directly affected by treatment in another (a standard caveat in macro panels).
3. **Stable measurement:** Index revisions and WID/SWIID revisions are treated as comparable over time; sensitivity analysis is outside the core script.
4. **Treatment definition:** The shock is a **large positive \(\Delta\)** on one sub-index, not a comprehensive measure of “neoliberalism” in all dimensions simultaneously.

---

## 11. Reporting checklist (alignment with best practice)

- State **never-treated vs not-yet-treated** and any fallback.
- Report **sample sizes** (countries and country-years) **per outcome**.
- Show **pre-trend** diagnostics: CSV exports **and** **`ggdid`** figures; for GDP, include the **\(t_{-5}\)–\(t_{+5}\)** plot when discussing parallel trends near the shock.
- Present **nominal** CS p-values and **FDR-adjusted** (\(q\)) p-values; if using the **funnel** figure, explain that the **fourth** interval is a **Bonferroni-style** multiplicity display on the **same** ATT/SE, not a second estimator.
- Include the **methodological funnel** (OLS → TWFE → CS–DiD → multiplicity-aware interval) for at least one headline component–outcome pair when summarizing how conclusions **tighten** or **dissolve** across stages.
- Distinguish **Level 1–2** associations from **Level 3** ATT estimates.
- Discuss **null results after FDR** as a **valid** scientific outcome when the pipeline is pre-specified and diagnostics are reported.

---

## 12. Software

- **R**, with packages **`tidyverse`**, **`fixest`**, **`did`** (and **`ggplot2`** via `ggdid`).
- Reproducibility: run from the thesis project root, e.g. `Rscript R/analysis_methodology.R`, after `merged_full_data.csv` exists (via `R/cleaning/merge_full_data.R` and upstream imports). That run writes the CSV results, the funnel PNG, and the event-study PNGs under `results/figures/`.

---

*This methods paper is descriptive of the code as implemented; if any parameter (years, `sd_mult`, outcomes, or `biters`) changes, update the prose to match.*
