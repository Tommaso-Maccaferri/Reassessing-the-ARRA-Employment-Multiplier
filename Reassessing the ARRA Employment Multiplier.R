################################################################################
# TWFE_DiD_ChodorowReich.R  —  ARRA Employment Multiplier
#                               Chodorow-Reich (2019) Dataset
#
# STRUCTURE:
#   Step 1  — Data loading, first stage, treatment assignment
#   Step 2  — Staggering non-evidence
#   Step 3  — Pre-trends plots
#   Step 4  — Parallel trends tests
#   Step 5  — Chodorow-Reich baseline IV
#   Step 6  — TWFE median split
#   Step 7  — TWFE Q4 vs Q1
#   Step 8  — Callaway-Sant'Anna Q4 vs Q1
#   Step 9  — CGBS continuous dose-response
#   Step 10 — Summary table with Wald conversion to $/job
#
# OUTCOME:     ces_sa_pc  (CES employment / civilian population 16+)
# TREATMENT:   outlays_total_pc1208  (ARRA outlays per Dec 2008 capita)
# INSTRUMENTS: instrument_fmap_pc1208 + instrument_dot_pc1208
# CONTROLS:    RealAnnualPI_3yrMADifference, log_sa_fhfa, manshare
# SE:          clustered at state (FIPS) level throughout
################################################################################

library(haven)
library(dplyr)
library(ggplot2)
library(fixest)
library(did)
library(AER)
library(patchwork)

# ── Helpers ───────────────────────────────────────────────────────────────────
ym <- function(year, month) (year - 1960L) * 12L + (month - 1L)
stata_month_to_date <- function(m)
  as.Date(paste(1960L + m %/% 12L, m %% 12L + 1L, "01", sep = "-"))

# ── Key dates ─────────────────────────────────────────────────────────────────
start_date <- ym(2007, 12)
post_date  <- ym(2009,  2)
apr_date   <- ym(2009,  4)
jul_date   <- ym(2009,  7)
end_date   <- ym(2010, 12)
fs_date    <- ym(2010, 12)

arra_date  <- stata_month_to_date(post_date)
apr_ddate  <- stata_month_to_date(apr_date)
jul_ddate  <- stata_month_to_date(jul_date)

# ── Theme ─────────────────────────────────────────────────────────────────────
theme_cr <- function()
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(1.8, "cm"),
    legend.text      = element_text(size = 9),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.text        = element_text(color = "black"),
    axis.ticks       = element_line(color = "black"),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8,  color = "grey30"),
    strip.background = element_blank()
  )

adoption_vlines <- function()
  list(
    geom_vline(xintercept = as.numeric(arra_date),
               linetype = "dashed", color = "black", linewidth = 0.6),
    geom_vline(xintercept = as.numeric(apr_ddate),
               linetype = "dotted", color = "black", linewidth = 0.6),
    geom_vline(xintercept = as.numeric(jul_ddate),
               linetype = "longdash", color = "black", linewidth = 0.6),
    annotate("text", x = arra_date, y = Inf,
             label = "ARRA\nFeb 09", vjust = 1.5, hjust = -0.1,
             size = 2.8, color = "black", family = "serif"),
    annotate("text", x = apr_ddate, y = Inf,
             label = "Wave 1\nApr 09", vjust = 1.5, hjust = -0.1,
             size = 2.8, color = "black", family = "serif"),
    annotate("text", x = jul_ddate, y = Inf,
             label = "Wave 2\nJul 09", vjust = 1.5, hjust = -0.1,
             size = 2.8, color = "black", family = "serif")
  )


################################################################################
# STEP 1 — DATA + TREATMENT ASSIGNMENT
################################################################################

df_raw <- read_dta("table1-dataset.dta") %>%
  filter(FIPS != 11) %>%
  mutate(date = stata_month_to_date(monthly))

df_cross    <- df_raw %>% filter(monthly == fs_date)

first_stage <- lm(
  outlays_total_pc1208 ~ instrument_fmap_pc1208 + instrument_dot_pc1208 +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare,
  data = df_cross
)

df_treatment <- df_cross %>%
  select(FIPS) %>%
  mutate(
    predicted_outlays = predict(first_stage, newdata = df_cross),
    quartile          = ntile(predicted_outlays, 4),
    treated_median    = as.integer(predicted_outlays >=
                                     median(predicted_outlays, na.rm = TRUE)),
    treated_q         = case_when(quartile == 4L ~ 1L,
                                  quartile == 1L ~ 0L,
                                  TRUE           ~ NA_integer_)
  )

df_panel <- df_raw %>%
  filter(between(monthly, start_date, end_date)) %>%
  left_join(df_treatment, by = "FIPS") %>%
  mutate(
    post           = as.integer(monthly >= post_date),
    treated_post   = treated_median * post,
    treated_post_q = treated_q * post,
    time_trend     = as.integer(monthly) - min(as.integer(monthly)),
    state_id       = as.integer(factor(FIPS)),
    dose           = (predicted_outlays - min(predicted_outlays, na.rm = TRUE)) /
      (max(predicted_outlays, na.rm = TRUE) -
         min(predicted_outlays, na.rm = TRUE))
  )


################################################################################
# STEP 2 — STAGGERING NON-EVIDENCE
################################################################################

adopt_wave <- df_panel %>%
  group_by(FIPS) %>%
  mutate(final_pc = max(outlays_total_pc1208, na.rm = TRUE)) %>%
  filter(outlays_total_pc1208 >= 0.10 * final_pc) %>%
  summarise(adoption_date = min(date), .groups = "drop") %>%
  mutate(wave = case_when(
    adoption_date == apr_ddate ~ "Wave 1 \u2014 April 2009 (27 states)",
    adoption_date == jul_ddate ~ "Wave 2 \u2014 July 2009 (23 states)",
    TRUE                       ~ "Other"
  ))

p_stagger <- df_panel %>%
  left_join(adopt_wave %>% select(FIPS, wave), by = "FIPS") %>%
  mutate(wave = replace(wave, is.na(wave), "Other")) %>%
  group_by(date, wave) %>%
  summarise(mean_obl = mean(outlays_total_pc1208, na.rm = TRUE),
            .groups  = "drop") %>%
  ggplot(aes(date, mean_obl, linetype = wave)) +
  geom_line(linewidth = 0.7, color = "black") +
  adoption_vlines() +
  scale_linetype_manual(values = c(
    "Wave 1 \u2014 April 2009 (27 states)" = "solid",
    "Wave 2 \u2014 July 2009 (23 states)"  = "dashed",
    "Other"                                = "dotted")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Outlays per Dec 2008 Capita ($)", linetype = NULL) +
  theme_cr()
p_stagger


################################################################################
# STEP 3 — PRE-TRENDS PLOTS
################################################################################

# Panel A: median split
trends_median <- df_panel %>%
  mutate(group = if_else(treated_median == 1L,
                         "High Outlays (Treated)", "Low Outlays (Control)")) %>%
  group_by(date, group) %>%
  summarise(mean_emp = mean(ces_sa_pc, na.rm = TRUE), .groups = "drop")

p_trends_median <- ggplot(trends_median,
                          aes(date, mean_emp, linetype = group)) +
  geom_line(linewidth = 1.2, color = "black") +
  adoption_vlines() +
  scale_linetype_manual(values = c("High Outlays (Treated)" = "solid",
                                   "Low Outlays (Control)"  = "dashed")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Avg. CES Emp-to-Pop Ratio", linetype = NULL) +
  theme_cr() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))

# Panel B: quartile averages
p_trends_q <- df_panel %>%
  filter(!is.na(quartile)) %>%
  mutate(qlabel = c("1" = "Q1 (Lowest)", "2" = "Q2",
                    "3" = "Q3", "4" = "Q4 (Highest)")[as.character(quartile)]) %>%
  group_by(date, qlabel) %>%
  summarise(mean_emp = mean(ces_sa_pc, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(date, mean_emp, linetype = qlabel)) +
  geom_line(linewidth = 1.2, color = "black") +
  adoption_vlines() +
  scale_linetype_manual(values = c("Q1 (Lowest)"  = "dotted", "Q2" = "dashed",
                                   "Q3"           = "longdash",
                                   "Q4 (Highest)" = "solid")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "Avg. CES Emp-to-Pop Ratio", linetype = NULL) +
  theme_cr() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))

(p_trends_median | p_trends_q) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))


################################################################################
# STEP 4 — PARALLEL TRENDS TESTS
################################################################################

df_pre <- df_panel %>% filter(monthly < post_date)

# Test 1: Differential linear trend — median split
slope_median <- feols(
  ces_sa_pc ~ treated_median * time_trend | FIPS + monthly,
  data = df_pre, vcov = ~FIPS
)

# Test 2: Differential linear trend — Q4 vs Q1
slope_q <- feols(
  ces_sa_pc ~ treated_q * time_trend | FIPS + monthly,
  data = df_pre %>% filter(quartile %in% c(1L, 4L)),
  vcov = ~FIPS
)

# Test 3: Joint F-test on pre-period event-study leads
es_pretest <- feols(
  ces_sa_pc ~ i(monthly, treated_median, ref = post_date - 1L) | FIPS + monthly,
  data = df_panel, vcov = ~FIPS
)
pre_names <- names(coef(es_pretest))[
  grepl("^monthly::", names(coef(es_pretest))) &
    as.integer(gsub("monthly::|:treated_median", "",
                    names(coef(es_pretest)))) < post_date - 1L
]
wald_res <- if (length(pre_names) > 0)
  tryCatch(wald(es_pretest, keep = pre_names), error = function(e) NULL)

# Test 4: Cross-sectional placebo (Dec 2008 → Jan 2009)
df_t0 <- df_panel %>%
  filter(monthly == ym(2008, 12), !is.na(ces_sa_pc)) %>%
  select(FIPS, treated_median, emp_t0 = ces_sa_pc)
df_t1 <- df_panel %>%
  filter(monthly == ym(2009,  1), !is.na(ces_sa_pc)) %>%
  select(FIPS, emp_t1 = ces_sa_pc)
pl_ols <- inner_join(df_t0, df_t1, by = "FIPS") %>%
  mutate(emp_change = emp_t1 - emp_t0) %>%
  { lm(emp_change ~ treated_median, data = .) }

# Summary
etable(slope_median, slope_q,
       title    = "Parallel Trends \u2014 Differential Slope Tests",
       headers  = c("Median split", "Q4 vs Q1"),
       keep     = "%treated_median:time_trend|treated_q:time_trend",
       se.below = TRUE)



################################################################################
# STEP 5 — CHODOROW-REICH BASELINE IV
################################################################################

df_cr_cross <- df_raw %>%
  filter(monthly == fs_date, FIPS != 11) %>%
  left_join(df_treatment %>% select(FIPS, predicted_outlays), by = "FIPS")

cr_ols <- lm(
  Cces_sa_pc ~ outlays_total_pc1208 +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare,
  data = df_cr_cross
)

cr_iv <- ivreg(
  Cces_sa_pc ~ outlays_total_pc1208 +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare |
    instrument_fmap_pc1208 + instrument_dot_pc1208 +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare,
  data = df_cr_cross
)

summary(cr_iv, diagnostics = TRUE)


################################################################################
# STEP 6 — TWFE MEDIAN SPLIT
################################################################################

twfe_median <- feols(
  ces_sa_pc ~ treated_post | FIPS + monthly,
  data = df_panel, vcov = ~FIPS
)
twfe_median_trend <- feols(
  ces_sa_pc ~ treated_post + state_id:time_trend | FIPS + monthly,
  data = df_panel, vcov = ~FIPS
)
twfe_median_ctrl <- feols(
  ces_sa_pc ~ treated_post + state_id:time_trend +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare | FIPS + monthly,
  data = df_panel, vcov = ~FIPS
)

etable(twfe_median, twfe_median_trend, twfe_median_ctrl,
       title    = "TWFE \u2014 Median Split",
       headers  = c("(1) Baseline", "(2) + Trends", "(3) + Trends + Controls"),
       keep     = "treated_post",
       se.below = TRUE)


################################################################################
# STEP 7 — TWFE Q4 vs Q1
################################################################################

df_panel_q <- df_panel %>% filter(quartile %in% c(1L, 4L))

twfe_q <- feols(
  ces_sa_pc ~ treated_post_q | FIPS + monthly,
  data = df_panel_q, vcov = ~FIPS
)
twfe_q_trend <- feols(
  ces_sa_pc ~ treated_post_q + state_id:time_trend | FIPS + monthly,
  data = df_panel_q, vcov = ~FIPS
)
twfe_q_ctrl <- feols(
  ces_sa_pc ~ treated_post_q + state_id:time_trend +
    RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare | FIPS + monthly,
  data = df_panel_q, vcov = ~FIPS
)

etable(twfe_q, twfe_q_trend, twfe_q_ctrl,
       title    = "TWFE \u2014 Q4 vs Q1",
       headers  = c("(1) Baseline", "(2) + Trends", "(3) + Trends + Controls"),
       keep     = "treated_post_q",
       se.below = TRUE)


################################################################################
# STEP 8 — CALLAWAY-SANT'ANNA (2021): Q4 vs Q1
################################################################################

pre_slopes <- df_panel %>%
  filter(monthly < post_date, quartile %in% c(1L, 4L)) %>%
  group_by(state_id) %>%
  summarise(pre_slope = coef(lm(ces_sa_pc ~ time_trend))[["time_trend"]],
            .groups   = "drop")

df_cs_q <- df_panel_q %>%
  left_join(pre_slopes, by = "state_id") %>%
  mutate(gname = as.numeric(if_else(treated_q == 1L, post_date, 0L)))

cs_q <- tryCatch(
  att_gt(
    yname         = "ces_sa_pc",
    tname         = "monthly",
    idname        = "state_id",
    gname         = "gname",
    xformla       = ~pre_slope + RealAnnualPI_3yrMADifference +
      log_sa_fhfa + manshare,
    data          = df_cs_q,
    control_group = "nevertreated",
    anticipation  = 0,
    base_period   = "universal",
    est_method    = "reg",
    clustervars   = "state_id"
  ),
  error = function(e) { message("CS error: ", e$message); NULL }
)

if (!is.null(cs_q)) {
  
  es_q <- data.frame(
    event_time = cs_q$t - cs_q$group,
    ATT        = cs_q$att,
    SE         = cs_q$se
  ) %>% mutate(upper = ATT + 1.96*SE, lower = ATT - 1.96*SE,
               pre   = event_time < 0)
  
  cs_att_val <- mean(cs_q$att[cs_q$t >= cs_q$group], na.rm = TRUE)
  cs_att_se  <- mean(cs_q$se[cs_q$t  >= cs_q$group], na.rm = TRUE)
  
  # Reference: TWFE median with trends + controls (preferred specification)
  twfe_ref <- coef(twfe_median_ctrl)["treated_post"]
  
  p_cs_q <- ggplot(es_q %>% filter(between(event_time, -14, 22)),
                   aes(event_time, ATT)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.12, fill = "grey50", show.legend = FALSE) +
    geom_point(aes(shape = pre), size = 2, color = "black") +
    geom_line(aes(linetype = pre), linewidth = 0.5, color = "black") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               color = "black", linewidth = 0.5) +
    geom_hline(yintercept = twfe_ref,
               linetype = "dotted", color = "black", linewidth = 0.6) +
    annotate("text", x = 1, y = twfe_ref,
             label = paste0("TWFE (median, controlled) = ",
                            round(twfe_ref, 4)),
             vjust = -0.7, size = 3, family = "serif") +
    scale_shape_manual(values = c("TRUE" = 1,       "FALSE" = 16),
                       labels = c("TRUE" = "Pre-treatment",
                                  "FALSE" = "Post-treatment")) +
    scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                          labels = c("TRUE" = "Pre-treatment",
                                     "FALSE" = "Post-treatment")) +
    scale_x_continuous(breaks = seq(-14, 22, by = 2)) +
    labs(
      subtitle = paste0("Conditioned on pre-slope + controls | Clustered SE | ",
                        "ATT = ", round(cs_att_val, 4),
                        "  SE = ", round(cs_att_se, 4)),
      x = "Months Relative to Feb 2009",
      y = "ATT: Change in CES Emp-to-Pop Ratio",
      shape = NULL, linetype = NULL
    ) +
    theme_cr()
  p_cs_q
}

################################################################################
# STEP 9 — CGBS (2024): CONTINUOUS DOSE-RESPONSE
################################################################################

pre_slopes_all <- df_panel %>%
  filter(monthly < post_date) %>%
  group_by(state_id) %>%
  summarise(pre_slope = coef(lm(ces_sa_pc ~ time_trend))[["time_trend"]],
            .groups   = "drop")

df_cgbs     <- df_panel %>% left_join(pre_slopes_all, by = "state_id")
pre_baseline <- ym(2009, 1)
df_pre_base  <- df_cgbs %>%
  filter(monthly == pre_baseline) %>%
  select(state_id, emp_pre = ces_sa_pc)

# Dynamic event study
dynamic_cgbs <- lapply(sort(unique(df_cgbs$monthly)), function(t) {
  if (t == pre_baseline) return(NULL)
  df_t <- df_cgbs %>%
    filter(monthly == t) %>%
    left_join(df_pre_base, by = "state_id") %>%
    mutate(emp_change = ces_sa_pc - emp_pre) %>%
    filter(!is.na(emp_change), !is.na(dose))
  if (nrow(df_t) < 10) return(NULL)
  m <- tryCatch(
    lm(emp_change ~ dose + pre_slope +
         RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare, data = df_t),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  tibble(monthly    = t, event_time = t - post_date,
         ATT_slope  = coef(m)["dose"],
         SE         = sqrt(diag(vcov(m)))["dose"])
}) %>% bind_rows()

# Dose-response at 12 months post
dose_range <- range(df_panel$predicted_outlays, na.rm = TRUE)
dose_grid  <- seq(0.1, 0.9, by = 0.05)
bandwidth  <- 0.2

df_post12 <- df_cgbs %>%
  filter(monthly == post_date + 12L) %>%
  left_join(df_pre_base, by = "state_id") %>%
  mutate(emp_change = ces_sa_pc - emp_pre) %>%
  filter(!is.na(emp_change))

dose_response <- lapply(dose_grid, function(d0) {
  df_w <- df_post12 %>%
    mutate(w = dnorm((dose - d0) / bandwidth)) %>%
    filter(w > 0.001)
  if (nrow(df_w) < 5) return(NULL)
  m <- tryCatch(
    lm(emp_change ~ dose + pre_slope +
         RealAnnualPI_3yrMADifference + log_sa_fhfa + manshare,
       data = df_w, weights = df_w$w),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  tibble(
    dose_std    = d0,
    dose_actual = d0 * diff(dose_range) + dose_range[1],
    ATT         = predict(m, newdata = tibble(
      dose = d0,
      pre_slope                    = mean(df_w$pre_slope),
      RealAnnualPI_3yrMADifference = mean(df_w$RealAnnualPI_3yrMADifference,
                                          na.rm = TRUE),
      log_sa_fhfa = mean(df_w$log_sa_fhfa, na.rm = TRUE),
      manshare    = mean(df_w$manshare,    na.rm = TRUE))),
    SE = sqrt(diag(vcov(m)))["(Intercept)"])
}) %>% bind_rows()

# Dose-response plot
p_dose <- ggplot(dose_response, aes(dose_actual, ATT)) +
  geom_ribbon(aes(ymin = ATT - 1.96*SE, ymax = ATT + 1.96*SE),
              alpha = 0.12, fill = "grey50") +
  geom_line(color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3,
             linetype = "dotted") +
  labs(x = "Predicted ARRA Outlays Per Dec 2008 Capita ($)",
       y = "ATT: Change in CES Emp-to-Pop Ratio") +
  theme_cr()
p_dose

# Dynamic event study plot
p_cgbs_es <- ggplot(dynamic_cgbs %>% filter(between(event_time, -14, 22)),
                    aes(event_time, ATT_slope)) +
  geom_ribbon(aes(ymin = ATT_slope - 1.96*SE,
                  ymax = ATT_slope + 1.96*SE),
              alpha = 0.12, fill = "grey50", show.legend = FALSE) +
  geom_point(aes(shape = event_time < 0), size = 1.8, color = "black") +
  geom_line(aes(linetype = event_time < 0), linewidth = 0.5, color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = -0.5, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  scale_shape_manual(values = c("TRUE" = 1,       "FALSE" = 16),
                     labels = c("TRUE" = "Pre-treatment",
                                "FALSE" = "Post-treatment")) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        labels = c("TRUE" = "Pre-treatment",
                                   "FALSE" = "Post-treatment")) +
  scale_x_continuous(breaks = seq(-14, 22, by = 2)) +
  labs(x = "Months Relative to Feb 2009",
       y = "ATT Slope: Emp Effect per Unit Predicted Outlays",
       shape = NULL, linetype = NULL) +
  theme_cr()
p_cgbs_es


################################################################################
# STEP 10 — SUMMARY TABLE WITH WALD CONVERSION TO $/JOB
################################################################################

df_endline <- df_panel %>%
  filter(monthly == fs_date) %>%
  select(FIPS, quartile, treated_median, outlays_total_pc1208, predicted_outlays)

gap_median <- df_endline %>%
  group_by(treated_median) %>%
  summarise(m = mean(outlays_total_pc1208, na.rm = TRUE), .groups = "drop") %>%
  summarise(gap = diff(m)) %>% pull(gap)

gap_q <- df_endline %>%
  filter(quartile %in% c(1L, 4L)) %>%
  group_by(quartile) %>%
  summarise(m = mean(outlays_total_pc1208, na.rm = TRUE), .groups = "drop") %>%
  arrange(quartile) %>%
  summarise(gap = diff(m)) %>% pull(gap)

dose_range_dollars <- diff(range(df_panel$predicted_outlays, na.rm = TRUE))

cgbs_post_avg <- mean(dynamic_cgbs$ATT_slope[dynamic_cgbs$event_time > 0],
                      na.rm = TRUE)
cgbs_post_se  <- mean(dynamic_cgbs$SE[dynamic_cgbs$event_time > 0],
                      na.rm = TRUE)

cr_iv_est <- coef(cr_iv)["outlays_total_pc1208"]
cr_iv_se  <- sqrt(vcov(cr_iv)["outlays_total_pc1208", "outlays_total_pc1208"])

wald_cost <- function(beta, gap) round((gap / beta) / 1000, 1)

summary_table <- tibble(
  Estimator = c(
    "CR Baseline IV",
    "TWFE median \u2014 baseline",
    "TWFE median \u2014 trends + controls",
    "TWFE Q4vsQ1 \u2014 baseline",
    "TWFE Q4vsQ1 \u2014 trends + controls",
    "CS Q4vsQ1 (2021)",
    "CGBS continuous (2024)"
  ),
  Estimate = round(c(
    cr_iv_est,
    coef(twfe_median)["treated_post"],
    coef(twfe_median_ctrl)["treated_post"],
    coef(twfe_q)["treated_post_q"],
    coef(twfe_q_ctrl)["treated_post_q"],
    if (exists("cs_att_val")) cs_att_val else NA_real_,
    cgbs_post_avg
  ), 5),
  Clustered_SE = round(c(
    cr_iv_se,
    se(twfe_median)["treated_post"],
    se(twfe_median_ctrl)["treated_post"],
    se(twfe_q)["treated_post_q"],
    se(twfe_q_ctrl)["treated_post_q"],
    if (exists("cs_att_se")) cs_att_se else NA_real_,
    cgbs_post_se
  ), 5),
  p_value = round(c(
    summary(cr_iv)$coefficients["outlays_total_pc1208", 4],
    pvalue(twfe_median)["treated_post"],
    pvalue(twfe_median_ctrl)["treated_post"],
    pvalue(twfe_q)["treated_post_q"],
    pvalue(twfe_q_ctrl)["treated_post_q"],
    NA_real_,
    NA_real_
  ), 4),
  Cost_per_Job_k = c(
    round(100000 / cr_iv_est / 1000, 1),
    wald_cost(coef(twfe_median)["treated_post"],      gap_median),
    wald_cost(coef(twfe_median_ctrl)["treated_post"], gap_median),
    wald_cost(coef(twfe_q)["treated_post_q"],         gap_q),
    wald_cost(coef(twfe_q_ctrl)["treated_post_q"],    gap_q),
    if (exists("cs_att_val")) wald_cost(cs_att_val, gap_q) else NA_real_,
    round((dose_range_dollars / cgbs_post_avg) / 1000, 1)
  )
)

summary_table