# =============================================================================
# Q1: Does aggressive behaviour track Acropora cover in line with economic
#     theories of resource defence?
#
# STUDY DESIGN
#   Three reef sites (Nata, Unarizaki, Sonai), Iriomote Island, Japan.
#   Surveys: 2016 (pre-bleaching), 2017 (bleaching), 2018 (early recovery),
#            2022 (longer-term recovery). Gap of 4 years between 2018–2022.
#
# AGGRESSION DIMENSIONS
#   Frequency   : binary per-encounter outcome (1 = aggressive)
#                 → decision to engage: "is it worth defending?"
#   Chase distance : metres chased when aggressive
#                 → energetic investment: "how hard do I defend?"
#
# ANALYTICAL APPROACH
#   Coral cover is the primary predictor rather than year. Year is collinear
#   with coral by design (the temporal pattern is hypothesised to be mediated
#   by resource state). Including year alongside coral prevents detection of
#   that relationship — coral is the biological signal year represents.
#
#   Frequency models are fit separately for conspecific and heterospecific
#   encounters, which represent fundamentally different behavioural processes
#   (resource competition vs competitor recognition; Keith et al. 2023).
#
#   Chase distance models pool both encounter types with encounter_type as a
#   fixed effect, because heterospecific-only chase data are too sparse to
#   support a standalone model.
#
# MODEL STRUCTURE
#   aggressive ~ coral_std + site + (1 | unique_fish)   [freq, per encounter type]
#   chase_m    ~ coral_std + encounter_type + site + (1 | unique_fish)
#
#   site as FIXED effect (Bolker: < 5 levels unsuitable as random)
#   (1 | unique_fish) accounts for repeated measures within focal follows
#
# =============================================================================

rm(list = ls())


# 0. PACKAGES =================================================================

library(tidyverse)
library(lubridate)
library(brms)
library(tidybayes)
library(patchwork)
library(scales)

set.seed(42)



# 1. LOAD DATA ================================================================

load("SimplifiedAggFishCoral.RData")

cat("Objects loaded:", ls(), "\n")
cat("\n-- agg --\n");   glimpse(agg)
cat("\n-- coral --\n"); glimpse(coral)

cat("\nOutcome values in agg:\n")
print(table(agg$outcome, useNA = "always"))

agg <- agg[!is.na(agg$outcome), ]
agg <- agg[agg$outcome != "n", ]



# 2. CONSTANTS ================================================================
#
# All key decisions, MCMC settings, priors, and visual parameters are defined
# here so the entire script can be updated from a single location.

# ── Survey years ──────────────────────────────────────────────────────────────
YEARS     <- c(2016, 2017, 2018, 2022)   # years included in analysis
ALL_YEARS <- 2016:2022                    # full x-axis range (preserves gap)

# ── Aggression outcome label ──────────────────────────────────────────────────
AGGR_OUTCOMES <- "a"

# ── MCMC settings: frequency models (Bernoulli) ───────────────────────────────
FREQ_FAMILY      <- bernoulli(link = "logit")
FREQ_CHAINS      <- 4
FREQ_CORES       <- 4
FREQ_ITER        <- 4000
FREQ_WARMUP      <- 2000
FREQ_ADAPT_DELTA <- 0.99

# ── MCMC settings: chase distance models (lognormal) ─────────────────────────
CHASE_FAMILY      <- lognormal()
CHASE_CHAINS      <- 4
CHASE_CORES       <- 4
CHASE_ITER        <- 4000
CHASE_WARMUP      <- 2000
CHASE_ADAPT_DELTA <- 0.99

# ── Priors (weakly informative) ───────────────────────────────────────────────
#
# Frequency (Bernoulli / logit scale):
#   Intercept Normal(-2, 1) → prior mode ~12% aggression on probability scale
#   Slopes    Normal(0, 1)  → logit slope of 1 ≈ 2.7× odds ratio
#
# Chase distance (lognormal / log scale):
#   Intercept Normal(1, 1)  → prior mode ~2.7 m on response scale
#   Slopes    Normal(0, 0.5)→ conservative; multiplicative change ~1.6× per SD
#   sigma     Exponential(1)→ positive, right-skewed, mildly regularising

prior_freq <- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0, 1),  class = b)
)

prior_chase <- c(
  prior(normal(1, 1),    class = Intercept),
  prior(normal(0, 0.5),  class = b),
  prior(exponential(1),  class = sigma)
)

# ── Visual constants ──────────────────────────────────────────────────────────
BASE_SIZE <- 12

site_colours <- c("Nata"      = "#7bafd4",
                  "Unarizaki" = "#e8957a",
                  "Sonai"     = "#74b87a")

site_shapes  <- c("Nata" = 16, "Unarizaki" = 17, "Sonai" = 15)

year_colours <- c("2016" = "#4d4d4d", "2017" = "#d73027",
                  "2018" = "#f46d43", "2022" = "#1a9850")

year_shapes  <- c("2016" = 16, "2017" = 17, "2018" = 15, "2022" = 18)



# 3. PREPARE AGGRESSION DATA ==================================================
#
# One row per encounter. unique_fish identifies the focal individual within a
# focal follow for the repeated-measures random effect.

agg_clean <- agg %>%
  mutate(reef = recode(reef, "Sonnai" = "Sonai")) %>%
  mutate(year = year(as.Date(date))) %>%
  filter(year %in% YEARS, encounter >= 1) %>%
  mutate(
    encounter_type = ifelse(focal.species == encountered.species,
                            "conspecific", "heterospecific"),
    aggressive     = as.integer(outcome %in% AGGR_OUTCOMES),
    chase_m        = as.numeric(chase.distance.m),
    unique_fish    = paste(reef, year, focal.species, observation, sep = "_")
  )

cat("\nEncounters per year × encounter type:\n")
print(table(agg_clean$year, agg_clean$encounter_type, agg_clean$reef))



# 4. PREPARE CORAL DATA =======================================================
#
# Coral cover is the primary predictor. We also test whether Acropora only
# improves predictive capacity of the models.

coral_proc <- coral %>%
  mutate(
    reef_name = recode(reef, "Sonnai" = "Sonai"),
    year      = year(as.Date(date))
  ) %>%
  filter(year %in% YEARS) %>%
  mutate(
    is_acropora   = str_detect(taxa, "^Acropora"),
    is_hard_coral = coral == 1
  )

total_pts <- coral_proc %>%
  group_by(transect) %>%
  summarise(total_points = sum(points, na.rm = TRUE), .groups = "drop")

cover_pts <- coral_proc %>%
  group_by(transect, reef_name, year) %>%
  summarise(
    acro_points  = sum(points[is_acropora],   na.rm = TRUE),
    coral_points = sum(points[is_hard_coral], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(total_pts, by = "transect") %>%
  mutate(
    acropora_pct  = 100 * acro_points  / total_points,
    hardcoral_pct = 100 * coral_points / total_points
  )

site_year_coral <- cover_pts %>%
  group_by(site = reef_name, year) %>%
  summarise(
    n_transects = n(),
    acro_mean   = mean(acropora_pct,  na.rm = TRUE),
    acro_se     = sd(acropora_pct,    na.rm = TRUE) / sqrt(n_transects),
    coral_mean  = mean(hardcoral_pct, na.rm = TRUE),
    coral_se    = sd(hardcoral_pct,   na.rm = TRUE) / sqrt(n_transects),
    .groups = "drop"
  ) %>%
  mutate(
    acro_std  = scale(acro_mean)[, 1],
    coral_std = scale(coral_mean)[, 1]
  )

# Store scaling parameters for back-transforming predictions to % cover
ACRO_CENTER   <- mean(site_year_coral$acro_mean)
ACRO_SCALE    <- sd(site_year_coral$acro_mean)
CORAL_CENTER  <- mean(site_year_coral$coral_mean)
CORAL_SCALE   <- sd(site_year_coral$coral_mean)

cat("\nAcropora cover per site × year:\n")
print(site_year_coral %>% select(site, year, acro_mean, acro_se, acro_std))
cat("\nCoral cover per site × year:\n")
print(site_year_coral %>% select(site, year, coral_mean, coral_se, coral_std))



# 5. COMBINE DATASETS =========================================================
#
# Join Acropora cover to raw encounter data at the site × year level.
# Split into separate data frames for conspecific and heterospecific frequency
# models. Chase distance data are filtered to aggressive encounters only.

model_df <- agg_clean %>%
  rename(site = reef) %>%
  inner_join(site_year_coral, by = c("site", "year")) %>%
  mutate(
    year_f         = factor(year, levels = YEARS),
    site           = factor(site, levels = c("Nata", "Unarizaki", "Sonai")),
    encounter_type = factor(encounter_type,
                            levels = c("conspecific", "heterospecific"))
  )

# Separate encounter-type subsets for frequency models
model_df_con <- model_df %>% filter(encounter_type == "conspecific")
model_df_het <- model_df %>% filter(encounter_type == "heterospecific")

# Chase distance: aggressive encounters only, strictly positive distances
# (lognormal requires y > 0)
model_df_chase <- model_df %>%
  filter(aggressive == 1, !is.na(chase_m), chase_m > 0)

cat("\nConspecific encounters:", nrow(model_df_con),
    "| Heterospecific:", nrow(model_df_het), "\n")
cat("Chase encounters:", nrow(model_df_chase), "\n")
cat("Unique focal fish:", length(unique(model_df$unique_fish)), "\n")



# 6. BAYESIAN MODELS ==========================================================
#
# LOO-CV (leave-one-out cross-validation):
#   Higher ELPD = better out-of-sample predictive accuracy.
#   ELPD difference > 4 SE considered practically meaningful.
#
# Bayes R²: proportion of variance in the outcome explained by the model,
#   computed from posterior draws (Gelman et al. 2019).

## 6a. Frequency — conspecific encounters --------------------------------------
#
# Conspecific aggression = direct resource competition between territorial rivals.
# Formula: aggressive ~ coral_std + site + (1 | unique_fish)

fit_freq_con <- brm(
  formula = aggressive ~ coral_std + site + (1 | unique_fish),
  data    = model_df_con,
  family  = FREQ_FAMILY,
  prior   = prior_freq,
  control = list(adapt_delta = FREQ_ADAPT_DELTA),
  chains  = FREQ_CHAINS, cores = FREQ_CORES,
  iter    = FREQ_ITER,   warmup = FREQ_WARMUP,
  seed = 42, file = "fit_freq_con_Q1"
)

cat("\n--- Conspecific frequency model summary ---\n")
print(summary(fit_freq_con))

## 6b. Frequency — heterospecific encounters -----------------------------------
#
# Heterospecific aggression = competitor recognition and species discrimination.
# Separate model because the biological process and baseline rate differ
# fundamentally from conspecific aggression (Keith et al. 2023).

fit_freq_het <- brm(
  formula = aggressive ~ coral_std + site + (1 | unique_fish),
  data    = model_df_het,
  family  = FREQ_FAMILY,
  prior   = prior_freq,
  control = list(adapt_delta = FREQ_ADAPT_DELTA),
  chains  = FREQ_CHAINS, cores = FREQ_CORES,
  iter    = FREQ_ITER,   warmup = FREQ_WARMUP,
  seed = 42, file = "fit_freq_het_Q1"
)

cat("\n--- Heterospecific frequency model summary ---\n")
print(summary(fit_freq_het))

## 6c. Chase distance — separate conspecific and heterospecific ----------------
#
# Split by encounter type to avoid the encounter_type coefficient dominating
# the model and obscuring the coral cover relationship.
# Heterospecific chase data are sparse — wide credible intervals are expected
# and results should be interpreted cautiously.

model_df_chase_con <- model_df_chase %>% filter(encounter_type == "conspecific")
model_df_chase_het <- model_df_chase %>% filter(encounter_type == "heterospecific")

cat("\nChase encounters — conspecific:", nrow(model_df_chase_con),
    "| heterospecific:", nrow(model_df_chase_het), "\n")

fit_chase_con <- brm(
  formula = chase_m ~ coral_std + site + (1 | unique_fish),
  data    = model_df_chase_con,
  family  = CHASE_FAMILY,
  prior   = prior_chase,
  control = list(adapt_delta = CHASE_ADAPT_DELTA),
  chains  = CHASE_CHAINS, cores = CHASE_CORES,
  iter    = CHASE_ITER,   warmup = CHASE_WARMUP,
  seed = 42, file = "fit_chase_con_Q1v7"
)

cat("\n--- Chase distance model summary (conspecific) ---\n")
print(summary(fit_chase_con))

fit_chase_het <- brm(
  formula = chase_m ~ coral_std + site + (1 | unique_fish),
  data    = model_df_chase_het,
  family  = CHASE_FAMILY,
  prior   = prior_chase,
  control = list(adapt_delta = CHASE_ADAPT_DELTA),
  chains  = CHASE_CHAINS, cores = CHASE_CORES,
  iter    = CHASE_ITER,   warmup = CHASE_WARMUP,
  seed = 42, file = "fit_chase_het_Q1v7"
)

cat("\n--- Chase distance model summary (heterospecific) ---\n")
print(summary(fit_chase_het))

## 6d. Model comparison: Acropora vs total coral cover ------------------------
#
# Here we fit equivalent models using Acropora cover and compare 
# predictive accuracy via LOO-CV.

fit_freq_con_acro <- brm(
  formula = aggressive ~ acro_std + site + (1 | unique_fish),
  data    = model_df_con,
  family  = FREQ_FAMILY,
  prior   = prior_freq,
  control = list(adapt_delta = FREQ_ADAPT_DELTA),
  chains  = FREQ_CHAINS, cores = FREQ_CORES,
  iter    = FREQ_ITER,   warmup = FREQ_WARMUP,
  seed = 42, file = "fit_freq_con_acro_Q1"
)
cat("\n--- Frequency conspecific ACROPORA model summary ---\n")
print(summary(fit_freq_con_acro))

fit_freq_het_acro <- brm(
  formula = aggressive ~ acro_std + site + (1 | unique_fish),
  data    = model_df_het,
  family  = FREQ_FAMILY,
  prior   = prior_freq,
  control = list(adapt_delta = FREQ_ADAPT_DELTA),
  chains  = FREQ_CHAINS, cores = FREQ_CORES,
  iter    = FREQ_ITER,   warmup = FREQ_WARMUP,
  seed = 42, file = "fit_freq_het_acro_Q1"
)
cat("\n--- Frequency heterospecific ACROPORA model summary ---\n")
print(summary(fit_freq_het_acro))

fit_chase_con_acro <- brm(
  formula = chase_m ~ acro_std + site + (1 | unique_fish),
  data    = model_df_chase_con,
  family  = CHASE_FAMILY,
  prior   = prior_chase,
  control = list(adapt_delta = CHASE_ADAPT_DELTA),
  chains  = CHASE_CHAINS, cores = CHASE_CORES,
  iter    = CHASE_ITER,   warmup = CHASE_WARMUP,
  seed = 42, file = "fit_chase_con_acro_Q1v7"
)
cat("\n--- Chase conspecific ACROPORA model summary ---\n")
print(summary(fit_chase_con_acro))

fit_chase_het_acro <- brm(
  formula = chase_m ~ acro_std + site + (1 | unique_fish),
  data    = model_df_chase_het,
  family  = CHASE_FAMILY,
  prior   = prior_chase,
  control = list(adapt_delta = CHASE_ADAPT_DELTA),
  chains  = CHASE_CHAINS, cores = CHASE_CORES,
  iter    = CHASE_ITER,   warmup = CHASE_WARMUP,
  seed = 42, file = "fit_chase_het_acro_Q1v7"
)
cat("\n--- Chase heterospecific ACROPORA model summary ---\n")
print(summary(fit_chase_het_acro))


## 6f. Posterior predictive checks --------------------------------------------
# check the observed line (black) falls within predicted lines (light blue) 
pp_check(fit_freq_con,       type = "dens_overlay", ndraws = 100)
pp_check(fit_freq_het,       type = "dens_overlay", ndraws = 100)
pp_check(fit_chase_con,          type = "dens_overlay", ndraws = 100)
pp_check(fit_chase_het,          type = "dens_overlay", ndraws = 100)
pp_check(fit_freq_con_acro, type = "dens_overlay", ndraws = 100)
pp_check(fit_freq_het_acro, type = "dens_overlay", ndraws = 100)
pp_check(fit_chase_con_acro,    type = "dens_overlay", ndraws = 100)
pp_check(fit_chase_het_acro,    type = "dens_overlay", ndraws = 100)

## 6g. LOO-CV: Acropora vs total coral cover ----------------------------------

fit_freq_con       <- add_criterion(fit_freq_con,       "loo")
fit_freq_het       <- add_criterion(fit_freq_het,       "loo")
fit_chase_con          <- add_criterion(fit_chase_con,   "loo")
fit_chase_het          <- add_criterion(fit_chase_het,  "loo")
fit_freq_con_acro <- add_criterion(fit_freq_con_acro, "loo")
fit_freq_het_acro <- add_criterion(fit_freq_het_acro, "loo")
fit_chase_con_acro    <- add_criterion(fit_chase_con_acro,"loo")
fit_chase_het_acro    <- add_criterion(fit_chase_het_acro,"loo")

# run moment matching to reduce the effect of outliers (Vehtari et al 2024)
fit_freq_con       <- add_criterion(fit_freq_con,       "loo", moment_match = TRUE)
fit_freq_het       <- add_criterion(fit_freq_het,       "loo", moment_match = TRUE)
fit_chase_con          <- add_criterion(fit_chase_con,  "loo", moment_match = TRUE)
fit_chase_het          <- add_criterion(fit_chase_het,  "loo", moment_match = TRUE)
fit_freq_con_acro <- add_criterion(fit_freq_con_acro, "loo", moment_match = TRUE)
fit_freq_het_acro <- add_criterion(fit_freq_het_acro, "loo", moment_match = TRUE)
fit_chase_con_acro    <- add_criterion(fit_chase_con_acro,    "loo", moment_match = TRUE)
fit_chase_het_acro    <- add_criterion(fit_chase_het_acro,    "loo", moment_match = TRUE)

cat("\nPareto k diagnostics (coral models):\n")
print(pareto_k_table(fit_freq_con$criteria$loo))
print(pareto_k_table(fit_freq_het$criteria$loo))
print(pareto_k_table(fit_chase_con$criteria$loo))
print(pareto_k_table(fit_chase_het$criteria$loo))

cat("\n--- LOO-CV: Conspecific frequency (Acropora vs total coral) ---\n")
print(loo_compare(fit_freq_con, fit_freq_con_acro))

cat("\n--- LOO-CV: Heterospecific frequency (Acropora vs total coral) ---\n")
print(loo_compare(fit_freq_het, fit_freq_het_acro))

cat("\n--- LOO-CV: Chase distance (Acropora vs total coral) ---\n")
print(loo_compare(fit_chase_con, fit_chase_con_acro))

cat("\n--- LOO-CV: Chase distance (Acropora vs total coral) ---\n")
print(loo_compare(fit_chase_het, fit_chase_het_acro))

## 6h. Coefficient plot --------------------------------------------------------
#
# Posterior distributions of all fixed-effect coefficients from the three
# primary models, faceted by outcome. The dashed line at zero is the key
# reference. Coefficients on the logit scale (frequency models) and log scale
# (chase distance) — not directly comparable across facets but interpretable
# within each.

# Helper to extract all b_ coefficients from a model
extract_coefs <- function(fit, outcome_label) {
  as_draws_df(fit) %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      outcome   = outcome_label,
      parameter = str_remove(parameter, "^b_") %>%
        str_replace_all("_", " ") %>%
        str_replace("coral std",     "Coral cover") %>%
        str_replace("encounter typeheterospecific", "Encounter: heterospecific") %>%
        str_replace("siteSonai",      "Site: Sonai") %>%
        str_replace("siteUnarizaki",  "Site: Unarizaki") %>%
        str_replace("Intercept",      "Intercept")
    )
}

coef_data <- bind_rows(
  extract_coefs(fit_freq_con,  "(a) Frequency\n(conspecific)"),
  extract_coefs(fit_freq_het,  "(b) Frequency\n(heterospecific)"),
  extract_coefs(fit_chase_con, "(c) Chase distance\n(conspecific)"),
  extract_coefs(fit_chase_het, "(d) Chase distance\n(heterospecific)")
) %>%
  filter(parameter != "Intercept") %>%
  mutate(
    outcome   = factor(outcome, levels = c("(a) Frequency\n(conspecific)",
                                           "(b) Frequency\n(heterospecific)",
                                           "(c) Chase distance\n(conspecific)",
                                           "(d) Chase distance\n(heterospecific)")),
    parameter = factor(parameter, levels = c("Coral cover",
                                             "Site: Unarizaki",
                                             "Site: Sonai"))
  )

# Compute whether 95% CI excludes zero for each coefficient × outcome
ci_flag <- coef_data %>%
  group_by(outcome, parameter) %>%
  summarise(
    lo95 = quantile(value, 0.025),
    hi95 = quantile(value, 0.975),
    sig  = lo95 > 0 | hi95 < 0,
    pp   = pmax(mean(value > 0), mean(value < 0)),
    .groups = "drop"
  )

coef_data <- coef_data %>%
  left_join(ci_flag %>% select(outcome, parameter, sig, pp),
            by = c("outcome", "parameter"))
coef_summary %>%
  mutate(outcome = str_replace(outcome, "\n", " ")) %>%
  write.csv("coef_summary.csv", row.names = FALSE)

print(coef_summary)


# shared plot elements used in both halves
coef_halfeye <- list(
  stat_halfeye(
    colour              = "grey25",
    slab_alpha          = 0.80,
    .width              = c(0.66, 0.95),
    point_interval      = median_qi,
    point_size          = 2.5,
    interval_size_range = c(0.6, 1.2),
    normalize           = "panels"
  ),
  scale_fill_gradient(
    low    = "white",
    high   = "#4d9b87",
    limits = c(0.5, 1.0),
    breaks = c(0.5, 0.75, 1.0),
    labels = c("0.50", "0.75", "1.00"),
    name   = "Posterior\nprobability"
  ),
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey55", linewidth = 0.5),
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))),
  facet_wrap(~ outcome, ncol = 2, scales = "free_x"),
  theme_classic(base_size = 14),
  theme(
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_blank(),
    panel.spacing    = unit(1.5, "lines"),
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 14),
    axis.text.y      = element_text(size = 14)
  )
)

fig_coef_freq <- ggplot(
  coef_data %>% filter(str_detect(outcome, "Frequency")),
  aes(x = value, y = parameter, fill = pp)
) +
  coef_halfeye +
  geom_text(
    data = ci_flag %>% filter(str_detect(outcome, "Frequency")),
    aes(x = -Inf, y = parameter, label = sprintf("PP = %.2f", pp)),
    hjust = -0.1, vjust = -4,
    size = 4, colour = "grey20",
    inherit.aes = FALSE
  ) +
  labs(x = "Posterior coefficient (logit scale)",
       y = NULL) +
  theme(legend.position = "right")

fig_coef_chase <- ggplot(
  coef_data %>% filter(str_detect(outcome, "Chase")),
  aes(x = value, y = parameter, fill = pp)
) +
  coef_halfeye +
  geom_text(
    data = ci_flag %>% filter(str_detect(outcome, "Chase")),
    aes(x = -Inf, y = parameter, label = sprintf("PP = %.2f", pp)),
    hjust = -0.1, vjust = -4,
    size = 4, colour = "grey20",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.05)),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  labs(x = "Posterior coefficient (log scale)",
       y = NULL) +
  theme(legend.position  = "none",
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank())

fig_coef <- (fig_coef_freq | fig_coef_chase) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("Figure_coef_Q1.pdf", fig_coef, width = 12, height = 5,
       units = "in", device = "pdf")
ggsave("Figure_coef_Q1.png", fig_coef, width = 12, height = 5,
       units = "in", dpi = 300)



# 7. POSTERIOR SUMMARIES FOR FIGURES ==========================================

## 7a. Observed site × year summaries (for plotting over model predictions) ----

obs_freq <- model_df %>%
  group_by(site, year, year_f, encounter_type) %>%
  summarise(
    freq    = mean(aggressive),
    n       = n(),
    freq_se = sqrt(freq * (1 - freq) / n),
    .groups = "drop"
  )

obs_chase <- model_df_chase %>%
  group_by(site, year, year_f, encounter_type) %>%
  summarise(
    chase_mu = mean(chase_m),
    chase_se = sd(chase_m) / sqrt(n()),
    n        = n(),
    .groups  = "drop"
  )

# Join coral cover to observed summaries (needed for conditional effects fig)
obs_freq_coral  <- obs_freq  %>%
  left_join(site_year_coral %>% select(site, year, coral_mean), by = c("site","year"))

obs_chase_coral <- obs_chase %>%
  left_join(site_year_coral %>% select(site, year, coral_mean), by = c("site","year"))

# coral_mean already present from Section 5 join — assign directly
model_df_coral       <- model_df
model_df_chase_coral <- model_df_chase

## 7b. Conditional effects posteriors for Figure 2 ----------------------------
#
# Predictions across the observed coral range, separately per site.
# Site is in newdata (fixed effect), unique_fish is marginalised (re_formula=NA).

coral_seq <- seq(min(model_df$coral_mean), max(model_df$coral_mean), length.out = 80)

nd_coral_base <- expand_grid(
  site       = factor(levels(model_df$site), levels = levels(model_df$site)),
  coral_mean = coral_seq
) %>%
  mutate(coral_std = (coral_mean - CORAL_CENTER) / CORAL_SCALE)

pred_freq_con <- nd_coral_base %>%
  add_epred_draws(fit_freq_con, ndraws = 2000, re_formula = NA) %>%
  group_by(site, coral_mean) %>%
  summarise(mean_est = mean(.epred),
            lo95 = quantile(.epred, 0.025),
            hi95 = quantile(.epred, 0.975), .groups = "drop") %>%
  mutate(encounter_type = "conspecific")

pred_freq_het <- nd_coral_base %>%
  add_epred_draws(fit_freq_het, ndraws = 2000, re_formula = NA) %>%
  group_by(site, coral_mean) %>%
  summarise(mean_est = mean(.epred),
            lo95 = quantile(.epred, 0.025),
            hi95 = quantile(.epred, 0.975), .groups = "drop") %>%
  mutate(encounter_type = "heterospecific")

pred_freq_cond <- bind_rows(pred_freq_con, pred_freq_het) %>%
  mutate(encounter_type = factor(encounter_type,
                                 levels = c("conspecific","heterospecific")))

pred_chase_con <- nd_coral_base %>%
  add_epred_draws(fit_chase_con, ndraws = 2000, re_formula = NA) %>%
  group_by(site, coral_mean) %>%
  summarise(mean_est = mean(.epred),
            lo95 = quantile(.epred, 0.025),
            hi95 = quantile(.epred, 0.975), .groups = "drop") %>%
  mutate(encounter_type = "conspecific")

pred_chase_het <- nd_coral_base %>%
  add_epred_draws(fit_chase_het, ndraws = 2000, re_formula = NA) %>%
  group_by(site, coral_mean) %>%
  summarise(mean_est = mean(.epred),
            lo95 = quantile(.epred, 0.025),
            hi95 = quantile(.epred, 0.975), .groups = "drop") %>%
  mutate(encounter_type = "heterospecific")

pred_chase_cond <- bind_rows(pred_chase_con, pred_chase_het) %>%
  mutate(encounter_type = factor(encounter_type,
                                 levels = c("conspecific", "heterospecific")))



# 8. FIGURE 1: TEMPORAL TRAJECTORIES ==========================================
#
# Panel A: Acropora and total coral cover per site × year
# Panel B: Aggression frequency per site × year, faceted by encounter type
# Panel C: Chase distance per site × year, faceted by encounter type
#
# Model ribbons show population-level predictions from the Acropora model
# evaluated at each year's mean Acropora cover across sites.

## 8a. Shared elements ---------------------------------------------------------

x_scale <- scale_x_continuous(
  limits       = c(2015.7, 2022.5),
  breaks       = ALL_YEARS,
  labels       = ALL_YEARS,
  minor_breaks = NULL
)

bleach_line <- geom_vline(
  xintercept = 2016.5, linetype = "dashed",
  colour = "grey45", linewidth = 0.55
)

strip_theme <- theme(
  strip.background = element_blank(),
  strip.text       = element_text(size = BASE_SIZE + 2, face = "bold")
)

## 8b. Panel A: coral cover ----------------------------------------------------

coral_long <- site_year_coral %>%
  pivot_longer(cols = c(acro_mean, coral_mean),
               names_to = "metric", values_to = "cover_mean") %>%
  mutate(
    cover_se = ifelse(metric == "acro_mean", acro_se, coral_se),
    metric   = recode(metric,
                      "acro_mean"  = "Acropora",
                      "coral_mean" = "Total hard coral")
  )

fig_a <- ggplot(coral_long,
                aes(x = year, y = cover_mean, colour = site,
                    group = interaction(site, metric),
                    shape = metric, linetype = metric)) +
  bleach_line +
  annotate("text", x = 2016.55, y = Inf, hjust = 0, vjust = 1.5,
           label = "Bleaching", colour = "grey25", size = 5, fontface = "bold") +
  geom_line(linewidth = 0.65, alpha = 0.7) +
  geom_errorbar(aes(ymin = cover_mean - 1.96 * cover_se,
                    ymax = cover_mean + 1.96 * cover_se),
                width = 0.15, alpha = 0.5) +
  geom_point(size = 2.8, stroke = 0.8, fill = "white") +
  scale_shape_manual(values = c("Acropora" = 16, "Total hard coral" = 21),
                     name = "Coral metric") +
  scale_linetype_manual(values = c("Acropora" = "solid", "Total hard coral" = "dashed"),
                        name = "Coral metric") +
  scale_colour_manual(values = site_colours, name = "Site") +
  guides(colour = "none") +
  x_scale +
  labs(x = NULL, y = "Cover (%)", tag = "(a)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.tag = element_text(face = "bold", size = BASE_SIZE + 2),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35))

## 8c. Panel B: aggression frequency -------------------------------------------

fig_b <- ggplot() +
  bleach_line +
  geom_line(data = obs_freq,
            aes(x = year, y = freq, colour = site, group = site),
            linewidth = 0.6, alpha = 0.65) +
  geom_errorbar(data = obs_freq,
                aes(x = year,
                    ymin = pmax(freq - 1.96 * freq_se, 0),
                    ymax = pmin(freq + 1.96 * freq_se, 1),
                    colour = site),
                width = 0.15, alpha = 0.5) +
  geom_point(data = obs_freq,
             aes(x = year, y = freq, colour = site, shape = site),
             size = 2.8) +
  facet_wrap(~ encounter_type, scales = "free_y") +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(values = site_shapes, name = "Site") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  guides(colour = "none", shape = "none") +
  x_scale +
  labs(x = NULL, y = "Aggression frequency\n(prop. encounters)", tag = "(b)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.tag = element_text(face = "bold", size = BASE_SIZE + 2),
        panel.grid.major.x = element_line(colour = "grey98", linewidth = 0.35)) +
  strip_theme

## 8d. Panel C: chase distance -------------------------------------------------

fig_c <- ggplot() +
  bleach_line +
  geom_line(data = obs_chase,
            aes(x = year, y = chase_mu, colour = site, group = site),
            linewidth = 0.6, alpha = 0.65) +
  geom_errorbar(data = obs_chase,
                aes(x = year,
                    ymin = pmax(chase_mu - 1.96 * chase_se, 0),
                    ymax = chase_mu + 1.96 * chase_se,
                    colour = site),
                width = 0.15, alpha = 0.5) +
  geom_point(data = obs_chase,
             aes(x = year, y = chase_mu, colour = site, shape = site),
             size = 2.8) +
  facet_wrap(~ encounter_type, scales = "free_y") +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(values = site_shapes, name = "Site") +
  x_scale +
  labs(x = "Year", y = "Chase distance (m)", tag = "(c)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2),
        panel.grid.major.x = element_line(colour = "grey98", linewidth = 0.35)) +
  strip_theme

## 8e. Assemble and save -------------------------------------------------------

fig_1 <- (fig_a / fig_b / fig_c) +
  plot_layout(heights = c(1.1, 1, 1), guides = "collect") &
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = BASE_SIZE + 2),
        legend.title     = element_text(size = BASE_SIZE + 2, face = "bold"),
        legend.key.size  = unit(1.2, "lines"),
        axis.text        = element_text(size = BASE_SIZE + 2),
        axis.title       = element_text(size = BASE_SIZE + 2))

ggsave("Figure1_temporal_Q1.pdf",
       fig_1, width = 7.5, height = 10, units = "in", device = "pdf")
ggsave("Figure1_temporal_Q1.png",
       fig_1, width = 7.5, height = 10, units = "in", dpi = 300)



# 9. FIGURE 2: CORAL–AGGRESSION RELATIONSHIP ================================
#
# Conditional effects plots from the Acropora models.
# X-axis: coral cover (%). Y-axis: aggression outcome.
# Lines: model-fitted relationship per site; ribbon: 95% CI.
# Jittered points: raw encounter-level data (density visible).
# Larger symbols: observed site × year means (shaped by year).
#
# Panel A: aggression frequency (conspecific and heterospecific from separate models)
# Panel B: chase distance (conspecific and heterospecific from separate models)

## 9a. Panel A: frequency ------------------------------------------------------

fig_2a <- ggplot() +
  geom_ribbon(data = pred_freq_cond,
              aes(x = coral_mean, ymin = lo95, ymax = hi95, group = site),
              fill = "grey85", alpha = 0.4) +
  geom_line(data = pred_freq_cond,
            aes(x = coral_mean, y = mean_est, colour = site),
            linewidth = 1.0) +
  geom_jitter(data = model_df_coral,
              aes(x = coral_mean, y = aggressive, colour = site),
              width = 0.4, height = 0.02, size = 0.6, alpha = 0.2) +
  geom_point(data = obs_freq_coral,
             aes(x = coral_mean, y = freq, colour = site, shape = year_f),
             size = 3.5) +
  facet_wrap(~ encounter_type, scales = "free_y") +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(values = year_shapes, name = "Year") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       y = "Aggression frequency\n(prop. encounters)",
       tag = "(a)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag     = element_text(face = "bold", size = BASE_SIZE + 2),
        legend.position = "none") +
  strip_theme

## 9b. Panel B: chase distance -------------------------------------------------

fig_2b <- ggplot() +
  geom_ribbon(data = pred_chase_cond,
              aes(x = coral_mean, ymin = lo95, ymax = hi95, group = site),
              fill = "grey85", alpha = 0.4) +
  geom_line(data = pred_chase_cond,
            aes(x = coral_mean, y = mean_est, colour = site),
            linewidth = 1.0) +
  geom_jitter(data = model_df_chase_coral,
              aes(x = coral_mean, y = chase_m, colour = site),
              width = 0.4, height = 0, size = 0.6, alpha = 0.2) +
  geom_point(data = obs_chase_coral,
             aes(x = coral_mean, y = chase_mu, colour = site, shape = year_f),
             size = 3.5) +
  facet_wrap(~ encounter_type, scales = "free_y") +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(values = year_shapes, name = "Year") +
  labs(x = "Coral cover (%)",
       y = "Chase distance (m)",
       tag = "(b)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(plot.tag        = element_text(face = "bold", size = BASE_SIZE + 2),
        strip.text      = element_blank(),   # already labelled in panel A
        strip.background = element_blank(),
        legend.position = "right")

## 9c. Assemble and save -------------------------------------------------------

fig_2 <- (fig_2a / fig_2b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.text     = element_text(size = BASE_SIZE),
        legend.title    = element_text(size = BASE_SIZE, face = "bold"))

ggsave("Figure2_coral_aggression_Q1.pdf",
       fig_2, width = 8, height = 8, units = "in", device = "pdf")
ggsave("Figure2_coral_aggression_Q1.png",
       fig_2, width = 8, height = 8, units = "in", dpi = 300)

cat("\nQ1v7 complete.\n")
