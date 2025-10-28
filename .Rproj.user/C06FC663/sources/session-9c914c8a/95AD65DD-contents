library(readxl)
library(ggpubr)
library(ggrepel)

# ======================================================
#       TREE STRESS ANALYSIS PIPELINE WITH ANOVA
# ======================================================
# --- Step 1: Compute cumulative metrics per tree ---
T_test_data <- trend_summer_2025 %>%
  group_by(Dendro_number) %>%
  mutate(
    TWD_CUM = cumsum(TWD),
    DV_CUM  = cumsum(DV)
  ) %>%
  ungroup()

# --- Step 2: Identify and remove outliers ---
outliers <- boxplot.stats(T_test_data$TWD_CUM)$out
cat("Outliers detected:", length(outliers), "\n")

DV_combined_clean <- T_test_data %>%
  filter(!TWD_CUM %in% outliers)  # using TWD for outlier removal

# --- Step 3: Aggregate by tree (to ensure independence) ---
tree_summary <- DV_combined_clean %>%
  group_by(Dendro_number, SP_CODE, Stress_Level) %>%
  summarise(
    TWD_CUM_final = max(TWD_CUM, na.rm = TRUE),
    DV_CUM_final  = max(DV_CUM, na.rm = TRUE),
    .groups = "drop"
  )

# --- Step 4: Define stress bins including 0 ---
tree_summary <- tree_summary %>%
  filter(!is.na(Stress_Level)) |>
  mutate(range = cut(
    Stress_Level,
    breaks = c(-Inf, 0, 9, 10),
    labels = c("0", "1-9", "10"),
    include.lowest = TRUE
  ))

# ======================================================
#                 ANOVA TESTS
# ======================================================

# --- ANOVA for TWD ---
anova_TWD <- aov(TWD_CUM_final ~ range, data = tree_summary)
summary(anova_TWD)

# Pairwise comparisons (Tukey)
tukey_TWD <- TukeyHSD(anova_TWD)
tukey_TWD

# --- ANOVA for DV ---
anova_DV <- aov(DV_CUM_final ~ range, data = tree_summary)
summary(anova_DV)

# Pairwise comparisons (Tukey)
tukey_DV <- TukeyHSD(anova_DV)
tukey_DV

# ======================================================
#                  BOXPLOTS
# ======================================================

# TWD boxplot with ANOVA significance
ggboxplot(tree_summary, x = "range", y = "TWD_CUM_final", fill = "range",
          add = "jitter", alpha = 0.7) +
  stat_compare_means(method = "anova", label = "p.signif") +
  labs(title = "TWD across Stress Levels",
       x = "Stress Level", y = "Cumulative TWD (µm)") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    #axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )

# DV boxplot with ANOVA significance
ggboxplot(tree_summary, x = "range", y = "DV_CUM_final", fill = "range",
          add = "jitter", alpha = 0.7) +
  stat_compare_means(method = "anova", label = "p.signif") +
  labs(title = "DV across Stress Levels",
       x = "Stress Level", y = "Cumulative DV (µm)") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    #axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )


