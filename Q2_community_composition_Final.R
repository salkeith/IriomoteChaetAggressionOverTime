# =============================================================================
# Q2: To what extent can observed changes in aggression be explained by
#     altered butterflyfish community composition and/or abundance?
#
# ANALYTICAL APPROACH
# -------------------
# Three nested questions:
#   (a) Did total butterflyfish density change across years?
#       → Linear mixed effects model (lme4), LRT vs null
#
#   (b) Did feeding guild composition change across years?
#       → Guild abundance aggregated to transect level; LMM testing
#         year × guild interaction vs additive null, LRT.
#         Proportional guild composition plotted to show what changed.
#
#   (c) Did species-level community composition change across years?
#       → PERMANOVA (vegan::adonis2) on Bray–Curtis dissimilarities;
#         site fitted first (sequential SS) to account for between-site
#         baseline differences. Pairwise comparisons with BH correction.
#         Homogeneity of dispersion assessed via PERMDISP (betadisper).
#         Proportional species composition plotted to show what changed.
#
# =============================================================================

rm(list = ls())

# 0. PACKAGES =================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)      # adonis2, vegdist, betadisper, cmdscale
library(patchwork)
library(scales)
library(ggtext)
library(dplyr)
library(tidyr)

set.seed(42)


# 1. LOAD AND PREPARE DATA ====================================================

## 1a. Inspect fish object -----------------------------------------------------

load("SimplifiedAggFishCoral.RData")

# Rename "reef" to "site" for consistency with Q1 code
colnames(fish)[1] <- "site"
cat("Fish object columns:\n")
print(names(fish))
cat("\nFirst rows:\n")
print(head(fish))

YEARS <- sort(unique(fish$year))
cat("\nYears detected:", YEARS, "\n")
cat("\nGuilds present:\n")
print(table(fish$guild, useNA = "ifany"))

## 1b. Set factor levels -------------------------------------------------------

fish <- fish %>%
  mutate(
    site   = factor(site,  levels = c("Nata", "Unarizaki", "Sonai")),
    year_f = factor(year,  levels = YEARS),
    guild  = factor(guild, levels = c("obligate", "facultative"))
  )

# Check all species have a guild assignment
unmatched <- fish %>%
  filter(is.na(guild)) %>%
  distinct(species)

if (nrow(unmatched) > 0) {
  cat("\nWARNING — species with no guild assignment (excluded from guild analyses):\n")
  print(unmatched)
} else {
  cat("\nAll species have a guild assignment.\n")
}


## 1c. Summarise for analyses --------------------------------------------------

# Total density per transect — zero counts retained as valid observations.
fish_transect <- fish %>%
  group_by(site, year, year_f, transect) %>%
  summarise(total_count = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(density = total_count / 250)

# Guild density aggregated to transect level — zeros retained (for LMMs).
fish_guild <- fish %>%
  filter(!is.na(guild)) %>%
  group_by(site, year, year_f, transect, guild) %>%
  summarise(guild_count = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(guild_density = guild_count / 250)

# Guild proportional composition per transect (for guild LMM summary output).
# Proportion of each guild out of total fish per transect; transects with
# zero total fish excluded (proportion undefined).
guild_prop <- fish %>%
  filter(!is.na(guild)) %>%
  group_by(site, year_f, transect, guild) %>%
  summarise(count = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(site, year_f, transect) %>%
  mutate(
    total = sum(count),
    prop  = ifelse(total > 0, count / total, NA_real_)
  ) %>%
  filter(!is.na(prop)) %>%
  ungroup()

# Mean proportional guild composition per year (pooled across sites and transects)
guild_prop_summary <- guild_prop %>%
  group_by(year_f, guild) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

cat("\nGuild proportional composition per year:\n")
print(guild_prop_summary)


# 2. TOTAL ABUNDANCE ===========================================================

## 2a. LM: total abundance ~ site (site differences) ---------------------------
# Tests whether sites differ in overall abundance and whether trajectories
# diverged across years. No repeated transects so no random effect needed.
# Density is right-skewed; log-transform applied (check hist/QQ above first).

# Check for zeros before log-transforming
if (any(fish_transect$density == 0)) {
  eps <- min(fish_transect$density[fish_transect$density > 0]) / 2
  fish_transect <- fish_transect %>% mutate(log_density = log(density + eps))
  cat(sprintf("\nZero densities present — log(density + %.4f) used\n", eps))
} else {
  fish_transect <- fish_transect %>% mutate(log_density = log(density))
  cat("\nNo zero densities — log(density) used\n")
}

mod_abu_site_null <- lm(log_density ~ year_f,        data = fish_transect)
mod_abu_site_add  <- lm(log_density ~ site + year_f, data = fish_transect)
mod_abu_site_int  <- lm(log_density ~ site * year_f, data = fish_transect)

lrt_abu_site <- anova(mod_abu_site_null, mod_abu_site_add, mod_abu_site_int)

cat("\n--- Total abundance by site: LRT ---\n")
print(lrt_abu_site)
cat("\n--- Total abundance by site: additive model summary ---\n")
print(summary(mod_abu_site_add))
anova(mod_abu_site_add)

# Residual check
par(mfrow = c(1, 2))
hist(residuals(mod_abu_site_add), main = "Residuals", xlab = "")
qqnorm(residuals(mod_abu_site_add)); qqline(residuals(mod_abu_site_add))
par(mfrow = c(1, 1))



# 3. GUILD COMPOSITION =========================================================

## 3a. LMM: guild abundance ~ year × guild -------------------------------------
# Response: density per guild per transect.
# The year × guild interaction tests whether the mix of guilds shifted
# differently across years.
#
# Additive (null): guild_density ~ year_f + guild + (1|site)
# Full:            guild_density ~ year_f * guild + (1|site)

mod_guild_null <- lmer(
  guild_density ~ 1 + (1 | site),
  data = fish_guild, REML = FALSE
)

mod_guild_add <- lmer(
  guild_density ~ year_f + guild + (1 | site),
  data = fish_guild, REML = FALSE
)

mod_guild_int <- lmer(
  guild_density ~ year_f * guild + (1 | site),
  data = fish_guild, REML = FALSE
)

lrt_guild <- anova(mod_guild_null, mod_guild_add, mod_guild_int)

cat("\n--- Guild composition: LRT (year × guild interaction) ---\n")
print(lrt_guild)


# 4. SPECIES COMPOSITION (PERMANOVA) ===========================================

## 4a. Build species × transect matrix -----------------------------------------

species_wide <- fish %>%
  group_by(site, year, year_f, transect, species) %>%
  summarise(count = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = species,
    values_from = count,
    values_fill = 0
  )

meta        <- species_wide %>% select(site, year, year_f, transect)
comm_matrix <- species_wide %>%
  select(-site, -year, -year_f, -transect) %>%
  as.matrix()

# Transects with zero total abundance are excluded from the community matrix
# only — Bray–Curtis dissimilarity is undefined for all-zero samples.
# These transects remain in the LMM analyses above.
row_totals <- rowSums(comm_matrix)
keep       <- row_totals > 0

if (any(!keep)) {
  cat(sprintf(
    "\n%d transect(s) with zero fish excluded from PERMANOVA (retained in LMMs)\n",
    sum(!keep)
  ))
  meta        <- meta[keep, ]
  comm_matrix <- comm_matrix[keep, ]
}

bc_dist <- vegdist(comm_matrix, method = "bray")

cat(sprintf("\nCommunity matrix: %d transects × %d species\n",
            nrow(comm_matrix), ncol(comm_matrix)))


## 4b. PERMANOVA: full community -----------------------------------------------
# Site fitted first (sequential SS) so between-site compositional differences
# are partitioned before testing the temporal effect of year.

perm_full <- adonis2(
  bc_dist ~ site + year_f,
  data         = meta,
  permutations = 999,
  by           = "terms"
)

cat("\n--- Full community PERMANOVA (all species) ---\n")
print(perm_full)


## 4c. PERMDISP: homogeneity of multivariate dispersion ------------------------
# A significant result means the PERMANOVA F-statistic may partly reflect
# a dispersion effect and should be interpreted cautiously.

disp_year <- betadisper(bc_dist, meta$year_f)
perm_disp <- permutest(disp_year, permutations = 999)

cat("\n--- PERMDISP: multivariate dispersion by year ---\n")
print(perm_disp)
cat("\nGroup dispersions:\n")
print(disp_year$group.distances)


## 4d. Pairwise PERMANOVA (BH-corrected) ----------------------------------------

pairwise_adonis2 <- function(dist_obj, meta, group_var, permutations = 999) {
  dist_mat <- as.matrix(dist_obj)
  groups   <- sort(unique(as.character(meta[[group_var]])))
  pairs    <- combn(groups, 2, simplify = FALSE)

  results <- map_dfr(pairs, function(pair) {
    idx      <- as.character(meta[[group_var]]) %in% pair
    sub_dist <- as.dist(dist_mat[idx, idx])
    sub_meta <- meta[idx, , drop = FALSE]
    res      <- adonis2(
      sub_dist ~ sub_meta[[group_var]],
      permutations = permutations
    )
    tibble(
      comparison = paste(pair, collapse = " vs "),
      F_value    = round(res$F[1], 3),
      R2         = round(res$R2[1], 3),
      p_value    = res$`Pr(>F)`[1]
    )
  })

  results$p_adj_BH <- p.adjust(results$p_value, method = "BH")
  results
}

pairwise_perm <- pairwise_adonis2(
  dist_obj     = bc_dist,
  meta         = meta,
  group_var    = "year_f",
  permutations = 999
)

cat("\n--- Pairwise PERMANOVA by year (BH-corrected) ---\n")
print(pairwise_perm)


## 4e. Species and guild composition by site -----------------------------------

# Mean density per species per site (pooled across years and transects)
spp_by_site <- fish %>%
  group_by(site, species) %>%
  summarise(mean_density = mean(abundance / 250, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site, values_from = mean_density, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total))

cat("\nMean species density by site (fish per 250m²):\n")
print(spp_by_site, n = Inf)

# Guild composition by site
guild_by_site <- fish %>%
  filter(!is.na(guild)) %>%
  group_by(site, guild) %>%
  summarise(mean_density = mean(abundance / 250, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  mutate(prop = mean_density / sum(mean_density))

cat("\nGuild composition by site:\n")
print(guild_by_site)

# Pairwise PERMANOVA by site
pairwise_site <- pairwise_adonis2(
  dist_obj     = bc_dist,
  meta         = meta,
  group_var    = "site",
  permutations = 999
)

cat("\n--- Pairwise PERMANOVA by site (BH-corrected) ---\n")
print(pairwise_site)


## 4f. Shannon diversity by site and year ----------------------------------------
# H' (Shannon) captures richness + evenness per transect.
# No repeated transects — lm() used throughout, mirroring section 2b.

# Reuse the filtered community matrix from section 4a (zero-filled, empty rows removed)
div_meta <- species_wide[keep, ] %>%
  select(site, year_f, transect) %>%
  mutate(
    H        = diversity(comm_matrix, index = "shannon"),
    richness = specnumber(comm_matrix),
    J        = ifelse(richness > 1, H / log(richness), NA_real_)
  )

cat("\nMean Shannon H' by site:\n")
print(
  div_meta %>%
    group_by(site) %>%
    summarise(mean_H = round(mean(H), 3), sd_H = round(sd(H), 3),
              mean_J = round(mean(J, na.rm = TRUE), 3), .groups = "drop")
)

cat("\nMean Shannon H' by site and year:\n")
print(
  div_meta %>%
    group_by(site, year_f) %>%
    summarise(mean_H = round(mean(H), 3), .groups = "drop") %>%
    pivot_wider(names_from = site, values_from = mean_H)
)

# Shannon H'
mod_H_null <- lm(H ~ year_f,        data = div_meta)
mod_H_site <- lm(H ~ site + year_f, data = div_meta)
mod_H_int  <- lm(H ~ site * year_f, data = div_meta)

cat("\n--- Shannon diversity: sequential LRT ---\n")
print(anova(mod_H_null, mod_H_site, mod_H_int))
cat("\n--- Shannon diversity: additive model ---\n")
print(anova(mod_H_site))
print(summary(mod_H_site))


# 5. ORDINATION ================================================================

## 5a. PCoA on Bray–Curtis dissimilarities -------------------------------------
# PCoA produces a unique solution and variance explained per axis is
# quantifiable from eigenvalues.

pcoa_out <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_eig <- pcoa_out$eig
pcoa_var <- round(pcoa_eig / sum(abs(pcoa_eig)) * 100, 1)

pcoa_df <- as_tibble(pcoa_out$points, .name_repair = "minimal") %>%
  setNames(c("PCo1", "PCo2")) %>%
  bind_cols(meta)

cat(sprintf("\nPCoA: PCo1 = %.1f%%, PCo2 = %.1f%% of variation\n",
            pcoa_var[1], pcoa_var[2]))


# 6. FIGURES ===================================================================

# ── Shared constants and palettes — consistent with Q1 ───────────────────────

BASE_SIZE <- 12

year_colours <- c("2016" = "#4d4d4d", "2017" = "#d73027",
                  "2018" = "#f46d43", "2022" = "#1a9850")

site_colours <- c("Nata"      = "#7bafd4",
                  "Unarizaki" = "#e8957a",
                  "Sonai"     = "#74b87a")

site_shapes  <- c("Nata" = 16, "Unarizaki" = 17, "Sonai" = 15)

guild_colours <- c("obligate"    = "#4aada8",
                   "facultative" = "#7dbf8e")

# Number of top species shown in both species figures (panels b and c)
N_TOP <- 10

strip_theme <- theme(
  strip.background = element_blank(),
  strip.text       = element_text(size = BASE_SIZE + 2, face = "bold")
)


# ── Panel (a): total abundance over time ──────────────────────────────────────
# Boxplot of raw transect densities; points coloured by site show spatial spread.

fig_abu <- ggplot(fish_transect,
                  aes(x = year_f, y = density)) +
  geom_boxplot(
    width         = 0.55,
    alpha         = 0.60,
    colour        = "grey30",
    fill          = "grey90",
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = site, shape = site),
    width = 0.18, size = 1.8, alpha = 0.70
  ) +
  scale_colour_manual(values = site_colours, name = "Site") +
  scale_shape_manual(values  = site_shapes,  name = "Site") +
  labs(
    x   = "Year",
    y   = expression("Chaetodon spp. density (fish m"^{-2}*")"),
    tag = "(a)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.tag           = element_text(face = "bold", size = BASE_SIZE + 2),
    panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.35),
    legend.position    = "none"
  )


# ── Panel (b): proportional species composition over time ─────────────────────
# Proportion of each top-N species relative to the total community per
# transect, averaged across transects per year. Shows whether the dominant
# species changed in relative abundance across the disturbance–recovery cycle.

top_species <- fish %>%
  group_by(species) %>%
  summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total, n = N_TOP) %>%
  arrange(desc(total)) %>%
  pull(species)

cat("\nTop", N_TOP, "species:", paste(top_species, collapse = ", "), "\n")

# Proportion of each top species as share of the total community per transect,
# averaged across transects per year
spp_prop <- fish %>%
  group_by(site, year_f, transect) %>%
  mutate(transect_total = sum(abundance, na.rm = TRUE)) %>%
  filter(species %in% top_species, transect_total > 0) %>%
  group_by(site, year_f, transect, species) %>%
  summarise(
    prop = sum(abundance, na.rm = TRUE) / first(transect_total),
    .groups = "drop"
  ) %>%
  group_by(year_f, species) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

# Species ordered by guild then abundance within guild
obligate_spp    <- c("lunulatus", "unimaculatus", "trifascialis", "ornatissimus", "plebeius")
facultative_spp <- c("vagabundus", "citrinellus", "rafflesii", "argentatus", "ephippium")
species_order   <- c(obligate_spp, facultative_spp)

spp_prop <- spp_prop %>%
  mutate(species = factor(species, levels = species_order))

# Guild-coded colour palette: warm tones (obligate), cool tones (facultative)
species_palette <- c(
  "lunulatus"    = "#c05a40",
  "unimaculatus" = "#d47c5c",
  "trifascialis" = "#e8a07a",
  "ornatissimus" = "#f0c09a",
  "plebeius"     = "#f5d8bc",
  "vagabundus"   = "#4e9da8",
  "citrinellus"  = "#7ec0ba",
  "rafflesii"    = "#6898c0",
  "argentatus"   = "#9ab8d4",
  "ephippium"    = "#b8d0c8"
)

# Transparent header entries create bold guild-label rows in the legend
species_palette_full <- c(
  "hdr_obl" = "#ffffff00",
  species_palette[obligate_spp],
  "hdr_fac" = "#ffffff00",
  species_palette[facultative_spp]
)

legend_breaks <- names(species_palette_full)

legend_labels <- c(
  "hdr_obl"      = "**Obligate**",
  "lunulatus"    = "*C. lunulatus*",
  "unimaculatus" = "*C. unimaculatus*",
  "trifascialis" = "*C. trifascialis*",
  "ornatissimus" = "*C. ornatissimus*",
  "plebeius"     = "*C. plebeius*",
  "hdr_fac"      = "**Facultative**",
  "vagabundus"   = "*C. vagabundus*",
  "citrinellus"  = "*C. citrinellus*",
  "rafflesii"    = "*C. rafflesii*",
  "argentatus"   = "*C. argentatus*",
  "ephippium"    = "*C. ephippium*"
)

# override.aes makes header rows invisible (no colour box, just the bold text)
override_fill   <- c("transparent", species_palette[obligate_spp],
                     "transparent", species_palette[facultative_spp])
override_colour <- c("transparent", rep("grey20", length(obligate_spp)),
                     "transparent", rep("grey20", length(facultative_spp)))

# Catch any top species not in the pre-defined palette and assign fallback colours
extra_spp  <- setdiff(top_species, names(species_palette))
extra_cols <- if (length(extra_spp) > 0) {
  setNames(hue_pal()(length(extra_spp)), extra_spp)
} else {
  c()
}

fig_spp <- ggplot(spp_prop,
                  aes(x = year_f, y = mean_prop, fill = species)) +
  geom_col(position = "stack", width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values   = species_palette_full,
    breaks   = legend_breaks,
    limits   = legend_breaks,
    labels   = legend_labels,
    na.value = "transparent",
    name     = NULL,
    guide    = guide_legend(
      ncol = 1,
      override.aes = list(
        fill   = override_fill,
        colour = override_colour
      )
    )
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Year", y = "Proportional abundance", tag = "(b)") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.tag           = element_text(face = "bold", size = BASE_SIZE + 2),
    panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.35),
    legend.position    = "right",
    legend.text        = element_markdown(size = BASE_SIZE - 1)
  )


# ── Panel (c): PCoA ordination ─────────────────────────────────────────────────
# 95% confidence ellipses per year show the degree of temporal overlap in
# community composition. Points shaped by site.

fig_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, colour = year_f)) +
  stat_ellipse(
    aes(fill = year_f),
    geom   = "polygon",
    alpha  = 0.10,
    colour = NA,
    level  = 0.95
  ) +
  stat_ellipse(
    level     = 0.95,
    linewidth = 0.5,
    lty       = 2
  ) +
  geom_point(
    aes(shape = site),
    size  = 2.5,
    alpha = 0.85
  ) +
  scale_colour_manual(values = year_colours, name = "Year") +
  scale_fill_manual(values   = year_colours, name = "Year") +
  scale_shape_manual(values  = site_shapes,  name = "Site") +
  labs(
    x   = sprintf("PCo1 (%.1f%%)", pcoa_var[1]),
    y   = sprintf("PCo2 (%.1f%%)", pcoa_var[2]),
    tag = "(c)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.tag         = element_text(face = "bold", size = BASE_SIZE + 2),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = BASE_SIZE - 1),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line        = element_line(colour = "black")
  )


# ── Assemble and save — 2×2 layout ────────────────────────────────────────────
# Row 1: (a) total abundance over time  |  
# Row 2: (b) proportional composition   |  (c) PCoA ordination

fig_community <- (fig_abu | fig_spp) / fig_pcoa +
  plot_layout(guides = "keep") &
  theme(
    text       = element_text(size = BASE_SIZE),
    axis.text  = element_text(size = BASE_SIZE),
    axis.title = element_text(size = BASE_SIZE),
    plot.tag   = element_text(face = "bold", size = BASE_SIZE + 2)
  )

ggsave("Figure_Q2_community.pdf", fig_community, width = 14, height = 5,
       units = "in", device = "pdf")
ggsave("Figure_Q2_community.png", fig_community, width = 14, height = 5,
       units = "in", dpi = 300)

cat("\nFigure saved: Figure_Q2_community.pdf / .png\n")
cat("\nQ2 community composition complete.\n")
