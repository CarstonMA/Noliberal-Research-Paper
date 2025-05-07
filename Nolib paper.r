########################################################################
##  Neolib Analysis – Master Script
##  -------------------------------------------------------------
##  1. Load & clean merged_data.csv
##  2. Run OLS + Two-way FE regressions            → regression_summary_clean.csv
##  3. Run Difference-in-Differences (±2 SD shock) → did_results_clean.csv
##  4. Produce illustrative scatter & DiD panels   → .png files
########################################################################

## ---------------- 0.  Packages -------------------------------------------------
required <- c(
  "tidyverse", "plm", "broom", "glue", "gridExtra", "grid",
  "sf", "visreg", "rnaturalearth", "rnaturalearthdata"
)

invisible(lapply(required, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

## ---------------- 1.  Load & clean data ---------------------------------------
data_raw <- read_csv("merged_data.csv", show_col_types = FALSE)
static_cols <- c("ISO3", "Year")

df <- data_raw %>%                       # remove rows w/ no info except ISO3,Year
  filter(!if_all(-all_of(static_cols), is.na)) %>%
  filter(!is.na(Heritage_Score) | !is.na(fraser_Score),
         Year >= 1995)

## ---------------- 2.  Variable labels -----------------------------------------
variable_labels <- tribble(
  ~raw, ~pretty,
  "Access.to.clean.fuels.and.technologies.for.cooking....of.population.",                       "Clean cooking access",
  "Account.ownership.at.a.financial.institution.or.with.a.mobile.money.service.provider..poorest.40.....of.population.ages.15..", "Account ownership (poorest 40 %)",
  "Annualized.average.growth.rate.in.per.capita.real.survey.mean.consumption.or.income..bottom.40..of.population....",           "Income growth (bottom 40 %)",
  "Carbon.dioxide..CO2..emissions..total..excluding.LULUCF..Mt.CO2e.",                       "Total CO₂ emissions",
  "Carbon.dioxide..CO2..emissions.excluding.LULUCF.per.capita..t.CO2e.capita.",              "CO₂ emissions per capita",
  "Children.in.employment..work.only....of.children.in.employment..ages.7.14.",              "Child labour (% 7–14)",
  "Consumer.price.index..2010...100.",                                                       "CPI (2010 = 100)",
  "Coverage.of.social.insurance.programs....of.population.",                                 "Social-insurance coverage",
  "Coverage.of.social.insurance.programs.in.poorest.quintile....of.population.",             "Social-insurance (poorest qnt.)",
  "Coverage.of.social.safety.net.programs....of.population.",                                "Safety-net coverage",
  "Coverage.of.social.safety.net.programs.in.poorest.quintile....of.population.",            "Safety-net (poorest qnt.)",
  "Educational.attainment..at.least.completed.lower.secondary..population.25...total......cumulative.", "Lower-secondary education (% 25+)",
  "Exports.of.goods.and.services..constant.2015.US..",                                       "Exports (constant USD)",
  "GDP..constant.2015.US..",                                                                 "GDP (constant USD)",
  "GDP.growth..annual...",                                                                   "GDP growth (%)",
  "GDP.per.capita..constant.2015.US..",                                                      "GDP per capita (constant USD)",
  "GDP.per.capita.growth..annual...",                                                        "GDP per-capita growth (%)",
  "Gini.index",                                                                              "Gini index (inequality)",
  "Income.share.held.by.highest.10.",                                                        "Income share (top 10 %)",
  "Income.share.held.by.lowest.10.",                                                         "Income share (bottom 10 %)",
  "Poverty.headcount.ratio.at..6.85.a.day..2017.PPP.....of.population.",                     "Poverty rate ($6.85/day)",
  "Carbon.dioxide..CO2..emissions.excluding.LULUCF.per.capita..t.CO2e.capita.",              "CO₂ emissions per capita",
  "Proportion.of.population.pushed.below.the..3.65....2017.PPP..poverty.line.by.out.of.pocket.health.care.expenditure....",       "Pushed below $3.65 by health costs"
  # … add or drop rows here as needed …
)

lookup_pretty  <- set_names(variable_labels$pretty, variable_labels$raw)
lookup_raw     <- set_names(variable_labels$raw,    variable_labels$pretty)

df <- df %>% rename(any_of(lookup_pretty))   # nice column names for modelling

## ---------------- 3.  Helper functions ----------------------------------------
run_ols <- function(outcome, score, data) {
  lm(glue("`{outcome}` ~ {score}"), data = data) |>
    tidy() |> slice(2) |> select(estimate, p.value, std.error)
}

run_fe <- function(outcome, score, data) {
  plm(
    glue("`{outcome}` ~ {score}"), data = data,
    index = c("ISO3", "Year"), model = "within", effect = "twoways"
  ) |>
    tidy() |> slice(1) |> select(estimate, p.value, std.error)
}

reg_summary <- function(outcome, pretty, score, data) {
  bind_cols(
    tibble(outcome = pretty, score = score),
    run_ols(outcome, score, data),
    run_fe(outcome, score, data) %>% rename_with(\(x) paste0("fe_", x))
  )
}

## ---------------- 4.  OLS & Fixed-Effects regressions --------------------------
scores <- c("Heritage_Score", "fraser_Score")

reg_tbl <- map_dfr(variable_labels$pretty, function(pretty) {
  map_dfr(scores, reg_summary, outcome = pretty,
          pretty = pretty, data = df)
})

write_csv(reg_tbl, "regression_summary_clean.csv")
print("Saved regression_summary_clean.csv")

## ---------------- 5.  Difference-in-Differences -------------------------------
df_did <- df %>%
  group_by(ISO3) %>%
  mutate(
    heri_diff = Heritage_Score - lag(Heritage_Score),
    fras_diff = fraser_Score   - lag(fraser_Score),
    big_h     = abs(heri_diff) >= 2 * sd(heri_diff, na.rm = TRUE),
    big_f     = abs(fras_diff) >= 2 * sd(fras_diff, na.rm = TRUE),
    treat_h   = as.numeric(cumsum(big_h) > 0),
    treat_f   = as.numeric(cumsum(big_f) > 0),
    post_h    = as.numeric(Year > min(Year[big_h],  na.rm = TRUE)),
    post_f    = as.numeric(Year > min(Year[big_f],  na.rm = TRUE))
  ) %>%
  ungroup()

run_did <- function(outcome, treat, post, data) {
  lm(glue("`{outcome}` ~ {treat} * {post} + factor(Year)"), data) %>%
    tidy() %>%
    filter(str_detect(term, ":")) %>%
    select(term, estimate, p.value, std.error)
}

did_tbl <- map_dfr(variable_labels$pretty, function(pretty) {
  bind_rows(
    run_did(pretty, "treat_h", "post_h", df_did) %>%
      mutate(outcome = pretty, score = "Heritage"),
    run_did(pretty, "treat_f", "post_f", df_did) %>%
      mutate(outcome = pretty, score = "Fraser")
  )
})

write_csv(did_tbl, "did_results_clean.csv")
print("Saved did_results_clean.csv")

## ---------------- 6.  Plot helpers --------------------------------------------
plot_scatter <- function(outcome, pretty, data) {
  ggplot(data, aes(x = Heritage_Score, y = .data[[outcome]])) +
    geom_point(alpha = .4) +
    geom_smooth(method = "lm", se = TRUE, colour = "steelblue") +
    labs(title = pretty, x = "Heritage Score", y = NULL) +
    theme_minimal(base_size = 13)
}

plot_did <- function(outcome, pretty, data) {
  data %>%
    group_by(treat_h, post_h) %>%
    summarise(y = mean(.data[[outcome]], na.rm = TRUE), .groups = "drop") %>%
    mutate(group  = ifelse(treat_h == 1, "Treatment", "Control"),
           period = ifelse(post_h  == 1, "Post", "Pre")) %>%
    ggplot(aes(period, y, colour = group, group = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(Control = "#377EB8", Treatment = "#E41A1C")) +
    labs(title = pretty, x = NULL, y = NULL, colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

## ---------------- 7.  Example panels ------------------------------------------
chosen <- c(
  "GDP per capita (constant USD)",
  "Gini index (inequality)",
  "Poverty rate ($6.85/day)",
  "Clean cooking access",
  "CO₂ emissions per capita",
  "Pushed below $3.65 by health costs"
)

## 7a  Scatter panel
scatter_plots <- map(chosen, ~ plot_scatter(.x, .x, df))
scatter_panel <- grid.arrange(
  grobs = scatter_plots, ncol = 3,
  top   = textGrob("OLS: Heritage Score vs. Outcomes",
                   gp = gpar(fontsize = 15, fontface = "bold"))
)
ggsave("panel_scatter_heritage.png", scatter_panel,
       width = 14, height = 9)

## 7b  DiD panel
did_plots <- map(chosen, ~ plot_did(.x, .x, df_did))
did_panel <- grid.arrange(
  grobs = did_plots, ncol = 3,
  top   = textGrob("Difference-in-Differences: Heritage Shocks",
                   gp = gpar(fontsize = 15, fontface = "bold"))
)
ggsave("panel_did_heritage.png", did_panel,
       width = 14, height = 9)

print("Panels saved: panel_scatter_heritage.png, panel_did_heritage.png")
