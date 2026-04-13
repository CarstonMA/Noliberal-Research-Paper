# Thesis-ready ggdid() wrappers: trimmed event window, plain-English labels, no Pre/Post legend.

#' Count treated unit-years at each event time (Year - first_treat), first_treat > 0 only.
event_time_cell_counts <- function(d) {
  if (!all(c("ISO3", "Year", "first_treat") %in% names(d))) {
    return(NULL)
  }
  d %>%
    dplyr::filter(first_treat > 0L, !is.na(Year)) %>%
    dplyr::mutate(egt = as.integer(Year) - as.integer(first_treat)) %>%
    dplyr::count(.data$egt, name = "n")
}

#' Choose x-axis limits: use cell counts when available; otherwise trim extreme bins by SE vs median.
#' Rule: keep event times with n >= max(3, min_frac * max(n)); cap |pre| and post span; require overlap with agg_dyn$egt.
event_study_xlim <- function(agg_dyn, d = NULL, min_frac_of_peak_n = 0.12, max_pre = 10L, max_post = 10L) {
  egt <- as.numeric(agg_dyn$egt)
  if (!length(egt) || !all(is.finite(egt))) {
    return(c(-5, 5))
  }
  in_agg <- function(lo, hi) {
    c(max(lo, min(egt)), min(hi, max(egt)))
  }
  ct <- event_time_cell_counts(d)
  if (!is.null(ct) && nrow(ct) > 0L) {
    peak <- max(ct$n, na.rm = TRUE)
    thr <- max(3L, floor(min_frac_of_peak_n * peak))
    ok <- ct$egt[ct$n >= thr & is.finite(ct$egt)]
    ok <- ok[ok %in% egt]
    if (length(ok)) {
      lo <- max(-max_pre, min(ok))
      hi <- min(max_post, max(ok))
      return(in_agg(lo, hi))
    }
  }
  se <- as.numeric(agg_dyn$se.egt)
  med <- stats::median(se[is.finite(se)], na.rm = TRUE)
  if (!is.finite(med) || med <= 0) {
    return(in_agg(max(-max_pre, min(egt)), min(max_post, max(egt))))
  }
  ok <- is.finite(egt) & is.finite(se) & (se <= 6 * med)
  if (!any(ok)) {
    return(in_agg(max(-max_pre, min(egt)), min(max_post, max(egt))))
  }
  rg <- range(egt[ok], na.rm = TRUE)
  in_agg(max(-max_pre, rg[1]), min(max_post, rg[2]))
}

pretty_outcome_label <- function(ol) {
  dict <- c(
    "GDP_per_capita_growth" = "GDP per Capita Growth",
    "Gini_SWIID" = "Gini Coefficient (SWIID)",
    "Income_share_bottom50_WID" = "Bottom-50% National Income Share (WID)",
    "log_under5_mort" = "Log Under-Five Mortality Rate"
  )
  if (ol %in% names(dict)) {
    dict[[ol]]
  } else {
    tools::toTitleCase(gsub("_", " ", ol))
  }
}

pretty_policy_component <- function(comp) {
  dict <- c(
    "fraser_SizeGov" = "Fraser: Size of Government",
    "fraser_LegalSystem" = "Fraser: Legal System and Property Rights",
    "fraser_SoundMoney" = "Fraser: Sound Money",
    "fraser_FreeTrade" = "Fraser: Freedom to Trade Internationally",
    "fraser_Regulation" = "Fraser: Regulation",
    "heritage_Property.Rights" = "Heritage: Property Rights",
    "heritage_Government.Integrity" = "Heritage: Government Integrity",
    "heritage_Judicial.Effectiveness" = "Heritage: Judicial Effectiveness",
    "heritage_Tax.Burden" = "Heritage: Tax Burden",
    "heritage_Government.Spending" = "Heritage: Government Spending",
    "heritage_Fiscal.Health" = "Heritage: Fiscal Health",
    "heritage_Business.Freedom" = "Heritage: Business Freedom",
    "heritage_Labor.Freedom" = "Heritage: Labor Freedom",
    "heritage_Monetary.Freedom" = "Heritage: Monetary Freedom",
    "heritage_Trade.Freedom" = "Heritage: Trade Freedom",
    "heritage_Investment.Freedom" = "Heritage: Investment Freedom",
    "heritage_Financial.Freedom" = "Heritage: Financial Freedom"
  )
  if (comp %in% names(dict)) {
    dict[[comp]]
  } else {
    tools::toTitleCase(gsub("[._]", " ", comp))
  }
}

#' Keep only event times in [lo, hi] so ggdid does not draw points outside the trimmed window
#' (coord_cartesian alone clips them and looks like "overflow").
subset_agg_dynamic <- function(agg, lo, hi) {
  e <- as.numeric(agg$egt)
  if (!length(e)) {
    return(agg)
  }
  k <- is.finite(e) & e >= lo & e <= hi
  if (!any(k)) {
    return(agg)
  }
  agg$egt <- agg$egt[k]
  agg$att.egt <- agg$att.egt[k]
  agg$se.egt <- agg$se.egt[k]
  if (!is.null(agg$crit.val.egt) && length(agg$crit.val.egt) == length(k)) {
    agg$crit.val.egt <- agg$crit.val.egt[k]
  }
  agg
}

#' @param title_override If non-NULL, used as the sole title (e.g. flagship one-liner).
thesis_ggdid <- function(agg_dyn, outcome_label, component, d = NULL, title_override = NULL,
                         ylab = "Average Treatment Effect",
                         xlab = "Event Time (Years Relative to First Shock)",
                         max_pre = 10L, max_post = 10L) {
  xl <- event_study_xlim(agg_dyn, d = d, max_pre = max_pre, max_post = max_post)
  agg_plot <- subset_agg_dynamic(agg_dyn, xl[1], xl[2])
  ttl <- if (!is.null(title_override)) {
    title_override
  } else {
    paste0(pretty_outcome_label(outcome_label), " — ", pretty_policy_component(component))
  }
  g0 <- did::ggdid(agg_plot, legend = FALSE, title = " ")
  eplot <- as.numeric(agg_plot$egt)
  pad <- 0.45
  xr <- if (length(eplot) && all(is.finite(range(eplot, na.rm = TRUE)))) {
    rg <- range(eplot, na.rm = TRUE)
    c(rg[1] - pad, rg[2] + pad)
  } else {
    NULL
  }
  g1 <- g0 + ggplot2::labs(title = ttl, subtitle = NULL, y = ylab, x = xlab, color = NULL)
  if (!is.null(xr)) {
    g1 <- g1 + ggplot2::coord_cartesian(xlim = xr, ylim = NULL, expand = TRUE, clip = "off")
  } else {
    g1 <- g1 + ggplot2::coord_cartesian(ylim = NULL, clip = "off")
  }
  g1 +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 12, 8, 12)
    )
}
