###MASKING MODELS: PROBABILITY OF DETECTION###
# Model 1: Distance × Stand Type + Masking (linear)
detection_mask_linear <- glmmTMB(
  detected ~ distance_m * stand_type + n_nontarget + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_linear)
confint(detection_mask_linear)

# Model 2: Distance × Stand Type + Masking (log-transformed)
detection_mask_log <- glmmTMB(
  detected ~ distance_m * stand_type + log(n_nontarget + 1) + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_log)
confint(detection_mask_log)

# Model 3: With Distance × Masking interaction (linear)
detection_mask_int_linear <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * n_nontarget + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_int_linear)
confint(detection_mask_int_linear)

# Model 4: With Distance × Masking interaction (log)
detection_mask_int_log <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * log(n_nontarget + 1) + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_int_log)
confint(detection_mask_int_log)

# Model 5: Full model with vegetation + masking (linear)
detection_mask_veg_linear <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    n_nontarget + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_veg_linear)
confint(detection_mask_veg_linear)

# Model 6: Full model with vegetation + masking (log)
detection_mask_veg_log <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    log(n_nontarget + 1) + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_veg_log)
confint(detection_mask_veg_log)

# Model 7: Full model + Distance × Masking interaction (linear)
detection_mask_veg_int_linear <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    distance_m * n_nontarget + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_veg_int_linear)
confint(detection_mask_veg_int_linear)

# Model 8: Full model + Distance × Masking interaction (log)
detection_mask_veg_int_log <- glmmTMB(
  detected ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    distance_m * log(n_nontarget + 1) + 
    (1|scientific_name) + (1|site),
  data = full_zero_filled,
  family = binomial()
)
summary(detection_mask_veg_int_log)
confint(detection_mask_veg_int_log)

# Betabinomial masking model 1
accuracy_activity_log_bb <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + log(n_nontarget + 1) +
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary, family = betabinomial(), weights = n_total
)
summary(accuracy_activity_log_bb)
confint(accuracy_activity_log_bb)

#Betabinomial model 4
accuracy_activity_log_bb_rs <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + log(n_nontarget + 1) +
    (1 + distance_m | expected_species) + (1|transect_id),
  data = accuracy_summary, family = betabinomial(), weights = n_total
)
summary(accuracy_activity_log_bb_rs)
confint(accuracy_activity_log_bb_rs)


# Compare all masking models
AIC(stand_comparison_veg_int2,          # Best model without masking
    detection_mask_linear,
    detection_mask_log,
    detection_mask_int_linear,
    detection_mask_int_log,
    detection_mask_veg_linear,
    detection_mask_veg_log,
    detection_mask_veg_int_linear,
    detection_mask_veg_int_log)





###ACCURACY MASKING MODELS###
# Model 1: Distance × Stand Type + Community Activity (linear)
accuracy_activity_linear <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + n_nontarget + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_linear)
confint(accuracy_activity_linear)

# Model 2: Distance × Stand Type + Community Activity (log-transformed)
accuracy_activity_log <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + log(n_nontarget + 1) + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_log)
confint(accuracy_activity_log)

# Model 3: With Distance × Activity interaction (linear)
accuracy_activity_int_linear <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * n_nontarget + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_int_linear)
confint(accuracy_activity_int_linear)

# Model 4: With Distance × Activity interaction (log)
accuracy_activity_int_log <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * log(n_nontarget + 1) + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_int_log)
confint(accuracy_activity_int_log)

# Model 5: Full model with vegetation + activity (linear)
accuracy_activity_veg_linear <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    n_nontarget + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_veg_linear)
confint(accuracy_activity_veg_linear)

# Model 6: Full model with vegetation + activity (log)
accuracy_activity_veg_log <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    log(n_nontarget + 1) + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_veg_log)
confint(accuracy_activity_veg_log)

# Model 7: Full model + Distance × Activity interaction (linear)
accuracy_activity_veg_int_linear <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    distance_m * n_nontarget + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_veg_int_linear)
confint(accuracy_activity_veg_int_linear)

# Model 8: Full model + Distance × Activity interaction (log)
accuracy_activity_veg_int_log <- glmmTMB(
  accuracy_pct/100 ~ distance_m * stand_type + distance_m * total_ba_plot + 
    canopy_LAI_conversion + ground_LAI_conversion + 
    distance_m * log(n_nontarget + 1) + 
    (1|expected_species) + (1|transect_id),
  data = accuracy_summary,
  family = binomial(),
  weights = n_total
)
summary(accuracy_activity_veg_int_log)
confint(accuracy_activity_veg_int_log)

# Compare all models
AIC(accuracy_int,                      # Best model without activity
    accuracy_activity_linear,
    accuracy_activity_log,
    accuracy_activity_int_linear,
    accuracy_activity_int_log,
    accuracy_activity_veg_linear,
    accuracy_activity_veg_log,
    accuracy_activity_veg_int_linear,
    accuracy_activity_veg_int_log)
