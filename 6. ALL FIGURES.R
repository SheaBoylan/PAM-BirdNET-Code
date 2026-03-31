# =============================================================================
# DISSERTATION FIGURES & TABLES
# Requires: all models and data objects from analysis script in environment
# Run top to bottom to produce all figures and tables
#
# Figure 1  — Detection probability ~ Distance x Stand type
# Figure 2  — Detection probability ~ Community vocal activity
# Figure 3  — Vegetation coefficient plot
# Figure 4  — Classification accuracy ~ Distance x Stand type
# Figure 5  — Species-specific classification accuracy
# Figure 6  — Top 10 misclassification pairs
# Figure S1 — Basal area & detection probability (two panels)
# Figure S2 — Classification accuracy ~ Community vocal activity (null result)
# Figure S3 — Random effects variance components
# Figure S4 — Detection rate x distance by species
# Table 1   — AIC model comparisons
# Table 2   — Fixed effects for best models
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(glmmTMB)

out_dir <- "C:/Users/Boyla/OneDrive/Desktop/University work/Year 4/Dissertation/Figures"
dir.create(out_dir, showWarnings = FALSE)


# =============================================================================
# SHARED SETUP
# =============================================================================

col_dec  <- "#2E7D52"
col_pine <- "#8B6914"
col_low  <- "#a8d5b5"
col_mid  <- "#2E7D52"
col_high <- "#1a3a1a"

stand_labels  <- c(deciduous = "Deciduous", pine = "Pine")
desk_accuracy <- 0.93
raw_breaks    <- c(0, 10, 50, 100, 200, 344)
dist_seq      <- seq(0, 100, length.out = 200)

mean_dist    <- mean(full_zero_filled$distance_m,           na.rm = TRUE)
mean_ba      <- mean(full_zero_filled$total_ba_plot,         na.rm = TRUE)
mean_can_LAI <- mean(full_zero_filled$canopy_LAI_conversion, na.rm = TRUE)
mean_grd_LAI <- mean(full_zero_filled$ground_LAI_conversion, na.rm = TRUE)

ba_q     <- quantile(full_zero_filled$total_ba_plot,         probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
ba_range <- range(full_zero_filled$total_ba_plot,            na.rm = TRUE)

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold", size = 11),
    legend.text        = element_text(size = 11),
    axis.title         = element_text(face = "bold", size = 12),
    axis.text          = element_text(size = 11),
    plot.tag           = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank()
  )

bar_theme <- theme_classic(base_size = 13) +
  theme(
    axis.title.x       = element_text(face = "bold", size = 12),
    axis.text.y        = element_text(face = "italic", size = 10),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    panel.grid.major.y = element_blank()
  )

y_prob <- scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1),
                             expand = expansion(add = c(0, 0)))
x_dist <- scale_x_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100),
                             expand = expansion(add = c(0, 3)))


# =============================================================================
# HELPERS
# =============================================================================

# Population-level predictions from glmmTMB (fixed effects only).
# re.form = NA marginalises over random effects.
# Models fitted with weights = n_total require n_total = 1 as a dummy in
# newdata — does not affect predicted probabilities, only SE scaling.
predict_glmmTMB <- function(model, newdata) {
  pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE,
                  re.form = NA, allow.new.levels = TRUE)
  newdata$fit <- plogis(pred$fit)
  newdata$lwr <- plogis(pred$fit - 1.96 * pred$se.fit)
  newdata$upr <- plogis(pred$fit + 1.96 * pred$se.fit)
  return(newdata)
}

# Prepend a (distance = 0, fit = anchor_val) row to a stand_type prediction grid
add_anchor <- function(pred_df, anchor_val = 1) {
  bind_rows(
    data.frame(distance_m = 0, stand_type = c("deciduous", "pine"),
               fit = anchor_val, lwr = anchor_val, upr = anchor_val),
    pred_df
  ) %>% arrange(stand_type, distance_m)
}

# Extract fixed effects + 95% CI from a glmmTMB model as a tidy data frame
extract_fixed <- function(model, model_name) {
  coef_df      <- as.data.frame(summary(model)$coefficients$cond)
  coef_df$term <- rownames(coef_df)
  ci           <- as.data.frame(confint(model, parm = "beta_"))
  ci$term      <- rownames(ci)
  colnames(ci) <- c("lwr_95", "upr_95", "estimate_ci", "term")
  coef_df %>%
    left_join(ci %>% select(term, lwr_95, upr_95), by = "term") %>%
    mutate(model = model_name, .before = term) %>%
    select(model, term, estimate = Estimate, se = `Std. Error`,
           z = `z value`, p_value = `Pr(>|z|)`, lwr_95, upr_95) %>%
    mutate(across(where(is.numeric), ~ round(., 4)),
           sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                           p_value < 0.05  ~ "*",   p_value < 0.1  ~ ".",
                           TRUE            ~ ""))
}

# Build an AIC comparison section as a tidy data frame
make_aic_section <- function(section_name, model_names, ...) {
  AIC(...) %>% as.data.frame() %>%
    mutate(section   = section_name,
           model     = model_names,
           delta_AIC = round(AIC - min(AIC), 2),
           AIC       = round(AIC, 2)) %>%
    select(section, model, df, AIC, delta_AIC)
}


# =============================================================================
# FIGURE 1: Detection probability ~ Distance x Stand type
# Model: stand_comparison_veg_int5
# Vegetation covariates at means; anchor at (0, 100%)
# =============================================================================

pred_det <- expand.grid(distance_m = dist_seq, stand_type = c("deciduous", "pine")) %>%
  mutate(total_ba_plot = mean_ba, canopy_LAI_conversion = mean_can_LAI,
         ground_LAI_conversion = mean_grd_LAI, scientific_name = NA, site = NA) %>%
  predict_glmmTMB(model = stand_comparison_veg_int5) %>%
  add_anchor(anchor_val = 1)

fig1 <- ggplot(pred_det, aes(x = distance_m, colour = stand_type, fill = stand_type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_colour_manual(values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_fill_manual(  values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  x_dist + y_prob +
  labs(x = "Distance from speaker (m)", y = "Probability of detection") +
  base_theme

ggsave(file.path(out_dir, "Figure1_Detection.png"), fig1, width = 6, height = 5, dpi = 300, bg = "white")
cat("✓ Figure 1 saved\n")


# =============================================================================
# FIGURE 2: Detection probability ~ Community vocal activity
# Model: detection_mask_veg_log
# X-axis: log(n_nontarget + 1) labelled with raw counts; distance at mean
# =============================================================================

pred_activity <- expand.grid(n_nontarget = seq(0, 344, length.out = 200),
                             stand_type  = c("deciduous", "pine")) %>%
  mutate(distance_m = mean_dist, total_ba_plot = mean_ba,
         canopy_LAI_conversion = mean_can_LAI, ground_LAI_conversion = mean_grd_LAI,
         scientific_name = NA, site = NA) %>%
  predict_glmmTMB(model = detection_mask_veg_log)

fig2 <- ggplot(pred_activity, aes(x = log(n_nontarget + 1), colour = stand_type, fill = stand_type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_colour_manual(values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_fill_manual(  values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_x_continuous(breaks = log(raw_breaks + 1), labels = raw_breaks, expand = expansion(add = c(0, 0.1))) +
  y_prob +
  labs(x = "Non-target detections (community vocal activity)", y = "Probability of detection") +
  base_theme

ggsave(file.path(out_dir, "Figure2_Community_Vocal_Activity.png"), fig2, width = 6, height = 5, dpi = 300, bg = "white")
cat("✓ Figure 2 saved\n")


# =============================================================================
# FIGURE 3: Vegetation coefficient plot
# Model: veg_model_ba_LAI_int
# Distance x vegetation interaction terms with 95% CI and significance stars
# =============================================================================

ci_veg      <- as.data.frame(confint(veg_model_ba_LAI_int2))
ci_veg$term <- rownames(ci_veg)
colnames(ci_veg) <- c("lwr", "upr", "estimate", "term")

coef_veg      <- as.data.frame(summary(veg_model_ba_LAI_int2)$coefficients$cond)
coef_veg$term <- rownames(coef_veg)

interactions <- ci_veg %>%
  left_join(coef_veg %>% select(term, `Pr(>|z|)`), by = "term") %>%
  rename(p_value = `Pr(>|z|)`) %>%
  filter(term %in% c("distance_m:total_ba_plot",
                     "distance_m:canopy_LAI_conversion",
                     "distance_m:ground_LAI_conversion",
                     "canopy_LAI_conversion:ground_LAI_conversion")) %>%
  mutate(
    label = factor(term,
                   levels = c("distance_m:total_ba_plot",
                              "distance_m:canopy_LAI_conversion",
                              "distance_m:ground_LAI_conversion",
                              "canopy_LAI_conversion:ground_LAI_conversion"),
                   labels = c("Distance x Basal area",
                              "Distance x Canopy LAI",
                              "Distance x Ground LAI",
                              "canopy_LAI_conversion:ground_LAI_conversion")),
    sig       = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                          p_value < 0.05  ~ "*",   TRUE           ~ "ns"),
    direction = ifelse(estimate > 0, "Positive", "Negative")
  )

fig3 <- ggplot(interactions, aes(x = estimate, y = label, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(x = upr, label = sig), hjust = -0.4, size = 4.5, colour = "grey30") +
  scale_colour_manual(values = c(Positive = col_dec, Negative = col_pine), name = "Effect direction") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Interaction coefficient estimate (log-odds)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.text    = element_text(size = 11),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 11),
        legend.text  = element_text(size = 11),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.y = element_blank())

ggsave(file.path(out_dir, "Figure3_Vegetation_Coefficients.png"), fig3, width = 7, height = 4, dpi = 300, bg = "white")
cat("✓ Figure 3 saved\n")


# =============================================================================
# FIGURE 4: Classification accuracy ~ Distance x Stand type
# Model: accuracy_int
# Anchor at (0, 93%) — empirical desk-based baseline
# n_total = 1 dummy required (model fitted with weights = n_total)
# =============================================================================

pred_acc <- expand.grid(distance_m = dist_seq, stand_type = c("deciduous", "pine")) %>%
  mutate(expected_species = NA, transect_id = NA, n_total = 1) %>%
  predict_glmmTMB(model = accuracy_int) %>%
  add_anchor(anchor_val = desk_accuracy)

fig4 <- ggplot(pred_acc, aes(x = distance_m, colour = stand_type, fill = stand_type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  geom_hline(yintercept = desk_accuracy, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  annotate("text", x = 2, y = desk_accuracy + 0.02,
           label = "Desk-based baseline (93%)", hjust = 0, size = 3.5, colour = "grey40") +
  scale_colour_manual(values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_fill_manual(  values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  x_dist + y_prob +
  labs(x = "Distance from speaker (m)", y = "Probability of correct classification") +
  base_theme

ggsave(file.path(out_dir, "Figure4_Classification.png"), fig4, width = 6, height = 5, dpi = 300, bg = "white")
cat("✓ Figure 4 saved\n")


# =============================================================================
# FIGURE 5: Species-specific classification accuracy
# Data: species_accuracy
# =============================================================================

tier_colours <- c("High (>=80%)" = col_dec, "Medium (50-79%)" = col_pine, "Low (<50%)" = "#8B2020")

species_plot <- species_accuracy %>%
  filter(!is.na(expected_species)) %>%
  arrange(accuracy_pct) %>%
  mutate(
    expected_species = factor(expected_species, levels = expected_species),
    tier = factor(case_when(accuracy_pct >= 80 ~ "High (>=80%)",
                            accuracy_pct >= 50 ~ "Medium (50-79%)",
                            TRUE               ~ "Low (<50%)"),
                  levels = c("High (>=80%)", "Medium (50-79%)", "Low (<50%)"))
  )

fig5 <- ggplot(species_plot, aes(x = accuracy_pct, y = expected_species, fill = tier)) +
  geom_vline(xintercept = 93, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(accuracy_pct, 1), "%")), hjust = -0.15, size = 3.2, colour = "grey30") +
  annotate("text", x = 92, y = Inf, label = "Desk-based baseline (93%)",
           hjust = 1, vjust = 0.94, size = 2.23, colour = "grey40") +
  scale_fill_manual(values = tier_colours, name = "Performance tier") +
  scale_x_continuous(limits = c(0, 102), breaks = c(0, 25, 50, 75, 100),
                     labels = percent_format(scale = 1, accuracy = 1),
                     expand = expansion(add = c(0, 0))) +
  labs(x = "Classification precision (%)", y = NULL) +
  bar_theme +
  theme(legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 11),
        legend.text     = element_text(size = 11))

ggsave(file.path(out_dir, "Figure5_Species_Precision.png"), fig5, width = 8, height = 7, dpi = 300, bg = "white")
cat("✓ Figure 5 saved\n")


# =============================================================================
# FIGURE 6: Top 10 misclassification pairs
# Data: playback_detections
# =============================================================================

top_confusions <- playback_detections %>%
  filter(expected_species != Scientific.name) %>%
  count(expected_species, Scientific.name, sort = TRUE) %>%
  head(10) %>%
  mutate(confusion_label = factor(paste0(expected_species, " -> ", Scientific.name),
                                  levels = rev(paste0(expected_species, " -> ", Scientific.name))))

fig6 <- ggplot(top_confusions, aes(x = n, y = confusion_label)) +
  geom_col(width = 0.7, fill = col_pine) +
  geom_text(aes(label = n), hjust = -0.3, size = 3.5, colour = "grey30") +
  scale_x_continuous(limits = c(0, max(top_confusions$n) * 1.15),
                     expand = expansion(add = c(0, 0))) +
  labs(x = "Number of misclassifications", y = NULL) +
  bar_theme

ggsave(file.path(out_dir, "Figure6_Confusion_Pairs.png"), fig6, width = 9, height = 5, dpi = 300, bg = "white")
cat("✓ Figure 6 saved\n")


# =============================================================================
# FIGURE S1: Basal area & detection probability
# Model: stand_comparison_veg_int5
# Panel A: Detection vs basal area (distance at mean, lines by stand type)
# Panel B: Detection vs distance (BA at quantiles, stand type marginalised)
# =============================================================================

ba_labels <- c(Low    = paste0("Low (", round(ba_q["25%"], 3), " m\u00b2)"),
               Medium = paste0("Medium (", round(ba_q["50%"], 3), " m\u00b2)"),
               High   = paste0("High (", round(ba_q["75%"], 3), " m\u00b2)"))

pred_ba_main <- expand.grid(total_ba_plot = seq(ba_range[1], ba_range[2], length.out = 200),
                            stand_type    = c("deciduous", "pine")) %>%
  mutate(distance_m = mean_dist, canopy_LAI_conversion = mean_can_LAI,
         ground_LAI_conversion = mean_grd_LAI, scientific_name = NA, site = NA) %>%
  predict_glmmTMB(model = stand_comparison_veg_int5)

pS1a <- ggplot(pred_ba_main, aes(x = total_ba_plot, colour = stand_type, fill = stand_type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_colour_manual(values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_fill_manual(  values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_x_continuous(expand = expansion(add = c(0, 0.01))) +
  y_prob +
  labs(tag = "A", x = "Basal area (m\u00b2)", y = "Probability of detection") +
  base_theme

pred_ba_dist <- expand.grid(distance_m  = dist_seq,
                            ba_quantile = c("Low", "Medium", "High")) %>%
  mutate(total_ba_plot = case_when(ba_quantile == "Low"    ~ ba_q["25%"],
                                   ba_quantile == "Medium" ~ ba_q["50%"],
                                   ba_quantile == "High"   ~ ba_q["75%"]),
         stand_type = "deciduous", canopy_LAI_conversion = mean_can_LAI,
         ground_LAI_conversion = mean_grd_LAI, scientific_name = NA, site = NA,
         ba_quantile = factor(ba_quantile, levels = c("Low", "Medium", "High"))) %>%
  predict_glmmTMB(model = stand_comparison_veg_int5)

pred_ba_dist <- bind_rows(
  data.frame(distance_m = 0,
             ba_quantile = factor(c("Low", "Medium", "High"), levels = c("Low", "Medium", "High")),
             fit = 1, lwr = 1, upr = 1,
             total_ba_plot = c(ba_q["25%"], ba_q["50%"], ba_q["75%"]),
             stand_type = "deciduous", canopy_LAI_conversion = mean_can_LAI,
             ground_LAI_conversion = mean_grd_LAI, scientific_name = NA, site = NA),
  pred_ba_dist
) %>% arrange(ba_quantile, distance_m)

pS1b <- ggplot(pred_ba_dist, aes(x = distance_m, colour = ba_quantile, fill = ba_quantile)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  scale_colour_manual(values = c(Low = col_low, Medium = col_mid, High = col_high), labels = ba_labels, name = "Basal area") +
  scale_fill_manual(  values = c(Low = col_low, Medium = col_mid, High = col_high), labels = ba_labels, name = "Basal area") +
  x_dist + y_prob +
  labs(tag = "B", x = "Distance from speaker (m)", y = "Probability of detection") +
  base_theme

figS1 <- (pS1a | pS1b) + plot_layout(guides = "keep")
ggsave(file.path(out_dir, "FigureS1_BasalArea.png"), figS1, width = 12, height = 5, dpi = 300, bg = "white")
cat("✓ Figure S1 saved\n")


# =============================================================================
# FIGURE S2: Classification accuracy ~ Community vocal activity (null result)
# Model: accuracy_activity_log
# accuracy_int (without activity) was preferred by AIC — included as
# supplementary demonstration that vocal activity does not affect accuracy
# =============================================================================

pred_acc_activity <- expand.grid(n_nontarget = seq(0, 344, length.out = 200),
                                 stand_type  = c("deciduous", "pine")) %>%
  mutate(distance_m = mean(accuracy_summary$distance_m, na.rm = TRUE),
         expected_species = NA, transect_id = NA, n_total = 1) %>%
  predict_glmmTMB(model = accuracy_activity_log)

figS2 <- ggplot(pred_acc_activity, aes(x = log(n_nontarget + 1), colour = stand_type, fill = stand_type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, colour = NA) +
  geom_line(aes(y = fit), linewidth = 1.1) +
  geom_hline(yintercept = desk_accuracy, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  annotate("text", x = 0.1, y = desk_accuracy + 0.02,
           label = "Desk-based baseline (93%)", hjust = 0, size = 3.5, colour = "grey40") +
  scale_colour_manual(values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_fill_manual(  values = c(deciduous = col_dec, pine = col_pine), labels = stand_labels, name = "Stand type") +
  scale_x_continuous(breaks = log(raw_breaks + 1), labels = raw_breaks, expand = expansion(add = c(0, 0.1))) +
  y_prob +
  labs(x = "Non-target detections (community vocal activity)", y = "Probability of correct classification") +
  base_theme

ggsave(file.path(out_dir, "FigureS2_Classification_VocalActivity.png"), figS2, width = 6, height = 5, dpi = 300, bg = "white")
cat("✓ Figure S2 saved\n")


# =============================================================================
# FIGURE S3: Random effects variance components
# Models: stand_comparison_veg_int5 (detection), accuracy_int (classification)
# =============================================================================

variance_data <- data.frame(
  model   = c("Detection", "Detection", "Classification", "Classification"),
  group   = c("Species", "Site", "Species", "Transect"),
  std_dev = c(1.34763, 0.67194, 1.15372, 0.39615)
) %>%
  mutate(variance = std_dev^2,
         model    = factor(model, levels = c("Detection", "Classification")),
         group    = factor(group, levels = c("Species", "Site", "Transect")))

figS3 <- ggplot(variance_data, aes(x = group, y = variance, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(variance, 3)), vjust = -0.5, size = 3.5, colour = "grey30") +
  scale_fill_manual(values = c(Species = col_high, Site = col_pine, Transect = col_pine), guide = "none") +
  scale_y_continuous(limits = c(0, max(variance_data$variance) * 1.2), expand = expansion(add = c(0, 0))) +
  facet_wrap(~ model, scales = "free_x") +
  labs(x = NULL, y = "Variance (SD\u00b2)") +
  theme_classic(base_size = 13) +
  theme(axis.title.y     = element_text(face = "bold", size = 12),
        axis.text        = element_text(size = 11),
        strip.text       = element_text(face = "bold", size = 12),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank())

ggsave(file.path(out_dir, "FigureS3_Variance_Components.png"), figS3, width = 7, height = 5, dpi = 300, bg = "white")
cat("✓ Figure S3 saved\n")


# =============================================================================
# FIGURE S4: Detection rate x distance by species (spaghetti plot)
# Three species with steepest drop-off at 100m highlighted automatically
# =============================================================================

det_by_species <- full_zero_filled %>%
  group_by(scientific_name, distance_m) %>%
  summarise(det_rate = mean(detected, na.rm = TRUE), .groups = "drop")

extreme_species <- det_by_species %>%
  filter(distance_m == 100) %>%
  arrange(det_rate) %>%
  slice(1:3) %>%
  pull(scientific_name)

highlight_colours <- c("#2E7D52", "#8B6914", "#8B2020", "grey80")
names(highlight_colours) <- c(extreme_species, "Other")

det_by_species <- det_by_species %>%
  mutate(highlight = factor(ifelse(scientific_name %in% extreme_species, scientific_name, "Other"),
                            levels = c(extreme_species, "Other")))

figS4 <- ggplot(det_by_species,
                aes(x = distance_m, y = det_rate, group = scientific_name,
                    colour = highlight, alpha = highlight, linewidth = highlight)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_colour_manual(values = highlight_colours, name = NULL,
                      guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) +
  scale_alpha_manual(   values = c(rep(1, 3), 0.4), guide = "none") +
  scale_linewidth_manual(values = c(rep(1.2, 3), 0.5), guide = "none") +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100), expand = expansion(add = c(2, 3))) +
  y_prob +
  labs(x = "Distance from speaker (m)", y = "Detection rate") +
  theme_classic(base_size = 13) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(face = "italic", size = 10),
        axis.title       = element_text(face = "bold", size = 12),
        axis.text        = element_text(size = 11),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank())

ggsave(file.path(out_dir, "FigureS4_Detection_Species.png"), figS4, width = 7, height = 5, dpi = 300, bg = "white")
cat("✓ Figure S4 saved\n")


# =============================================================================
# TABLE 1: AIC model comparisons
# =============================================================================

aic_table <- bind_rows(
  make_aic_section("Detection: Stand Type",
                   c("Distance x Stand Type", "+ Vegetation (main effects)",
                     "+ Distance x Basal Area", "+ Distance x Ground LAI",
                     "+ Distance x Canopy LAI", "+ Canopy LAI x Ground LAI (best)"),
                   stand_comparison_int, stand_comparison_veg, stand_comparison_veg_int2,
                   stand_comparison_veg_int3, stand_comparison_veg_int4, stand_comparison_veg_int5),
  
  make_aic_section("Detection: Vegetation",
                   c("Distance + Vegetation (main effects)", "+ Distance x Basal Area",
                     "+ Distance x Canopy LAI + Distance x Ground LAI (best)"),
                   veg_model_only, veg_model_ba_int, veg_model_ba_LAI_int),
  
  make_aic_section("Detection: Masking",
                   c("Best model without masking", "+ Activity (linear)", "+ Activity (log)",
                     "+ Distance x Activity (linear)", "+ Distance x Activity (log)",
                     "+ Vegetation + Activity (linear)", "+ Vegetation + Activity (log) (best)",
                     "+ Vegetation + Distance x Activity (linear)", "+ Vegetation + Distance x Activity (log)"),
                   stand_comparison_veg_int2, detection_mask_linear, detection_mask_log,
                   detection_mask_int_linear, detection_mask_int_log, detection_mask_veg_linear,
                   detection_mask_veg_log, detection_mask_veg_int_linear, detection_mask_veg_int_log),
  
  make_aic_section("Classification: Stand Type",
                   c("Distance x Stand Type (best)", "+ Vegetation (main effects)",
                     "+ Distance x Basal Area", "+ Distance x Ground LAI",
                     "+ Distance x Canopy LAI", "+ Canopy LAI x Ground LAI"),
                   accuracy_int, accuracy_veg, accuracy_veg_int2,
                   accuracy_veg_int3, accuracy_veg_int4, accuracy_veg_int5),
  
  make_aic_section("Classification: Masking",
                   c("Best model without activity (best)", "+ Activity (linear)", "+ Activity (log)",
                     "+ Distance x Activity (linear)", "+ Distance x Activity (log)",
                     "+ Vegetation + Activity (linear)", "+ Vegetation + Activity (log)",
                     "+ Vegetation + Distance x Activity (linear)", "+ Vegetation + Distance x Activity (log)"),
                   accuracy_int, accuracy_activity_linear, accuracy_activity_log,
                   accuracy_activity_int_linear, accuracy_activity_int_log, accuracy_activity_veg_linear,
                   accuracy_activity_veg_log, accuracy_activity_veg_int_linear, accuracy_activity_veg_int_log)
)

write.csv(aic_table, file.path(out_dir, "Table1_AIC_Comparisons.csv"), row.names = FALSE)
cat("✓ Table 1 saved\n")


# =============================================================================
# TABLE 2: Fixed effects for best models
# =============================================================================

fixed_effects <- bind_rows(
  extract_fixed(stand_comparison_veg_int5, "Detection: Stand Type"),
  extract_fixed(veg_model_ba_LAI_int,      "Detection: Vegetation"),
  extract_fixed(detection_mask_veg_log,    "Detection: Masking"),
  extract_fixed(accuracy_int,              "Classification: Stand Type"),
  extract_fixed(accuracy_activity_log,     "Classification: Masking")
)

write.csv(fixed_effects, file.path(out_dir, "Table2_Fixed_Effects.csv"), row.names = FALSE)
cat("✓ Table 2 saved\n")
