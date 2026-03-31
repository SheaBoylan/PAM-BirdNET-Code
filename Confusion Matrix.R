###CONFUSION MATRIX###
library(dplyr)
library(tidyr)
library(ggplot2)

# Create confusion matrix data
confusion_data <- playback_detections %>%
  count(expected_species, Scientific.name) %>%
  # Calculate percentages by expected species (row percentages)
  group_by(expected_species) %>%
  mutate(
    total = sum(n),
    percentage = (n / total) * 100
  ) %>%
  ungroup()

# View top confusions
confusion_data %>%
  filter(expected_species != Scientific.name) %>%  # Only misclassifications
  arrange(desc(n)) %>%
  head(20)

# Create matrix format
confusion_matrix <- confusion_data %>%
  select(expected_species, Scientific.name, n) %>%
  pivot_wider(
    names_from = Scientific.name,
    values_from = n,
    values_fill = 0
  )

print(confusion_matrix)


# Heatmap of confusion matrix
ggplot(confusion_data, aes(x = Scientific.name, y = expected_species, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "darkred", 
                      name = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    x = "BirdNET Prediction",
    y = "Expected Species (Playback)",
    title = "BirdNET Classification Confusion Matrix"
  )

ggsave("Confusion_Matrix.png", width = 12, height = 10, dpi = 300)


# Overall confusion summary
confusion_summary <- playback_detections %>%
  mutate(
    correct = (expected_species == Scientific.name)
  ) %>%
  summarise(
    total_predictions = n(),
    correct_predictions = sum(correct),
    incorrect_predictions = sum(!correct),
    overall_accuracy = mean(correct) * 100
  )

print(confusion_summary)

# Species-specific accuracy
species_accuracy <- playback_detections %>%
  mutate(correct = (expected_species == Scientific.name)) %>%
  group_by(expected_species) %>%
  summarise(
    n_predictions = n(),
    n_correct = sum(correct),
    accuracy_pct = mean(correct) * 100,
    .groups = 'drop'
  ) %>%
  arrange(desc(accuracy_pct))

print(species_accuracy)

# What are the top misclassifications?
top_confusions <- playback_detections %>%
  filter(expected_species != Scientific.name) %>%  # Only errors
  count(expected_species, Scientific.name, sort = TRUE) %>%
  head(10)

print(top_confusions)

# Visualize top confusions
ggplot(top_confusions, aes(x = reorder(paste(expected_species, "→", Scientific.name), n), 
                           y = n)) +
  geom_col(fill = "coral") +
  coord_flip() +
  labs(
    x = "Confusion (Expected → Predicted)",
    y = "Count",
    title = "Top 10 BirdNET Misclassifications"
  ) +
  theme_minimal()
