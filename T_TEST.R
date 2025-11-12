library(readxl)
library(ggpubr)
library(ggrepel)

# ======================================================
#      BOXPLOT OF TWD VS DV FOR STRESS LEVEL 1 AND 10
# ======================================================

# Filter data for Stress Levels 1 and 10
tree_summary <- trend_summer_2025 %>%
  filter(Stress_Level %in% c(1, 10))
summary(tree_summary)

# ======================================================
#       TREE STRESS ANALYSIS PIPELINE WITH ANOVA
# ======================================================
# --- Step 1: Compute cumulative metrics per tree ---
T_test_data <- trend_summer_2025 %>%
  group_by(Dendro_number) %>%
  mutate(
    TWD_CUM = cumsum(TWD), #cunsum
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
    breaks = c(-Inf, 5, 9, 10),
    labels = c( "1-5", "6-9", "10"),
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

# Make a dataframe for stat_pvalue_manual
tukey_TWD_df <- data.frame(
  group1 = "10",
  group2 = "6-9",
  p.adj = 0.024,
  y.position = max(tree_summary$TWD_CUM_final) * 0.98  # position above boxplot
)


# TWD boxplot with ANOVA significance
ggboxplot(tree_summary, x = "range", y = "TWD_CUM_final", fill = "range",
          add = "jitter", alpha = 0.7) +
  stat_pvalue_manual(tukey_TWD_df, label = "p.adj", tip.length = 0.03) +
  labs(title = "TWD across Stress Levels",
       x = "Stress Level", y = "Cumulative TWD (µm)") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
  )

# Define the comparisons
# DV boxplot with ANOVA significance
ggboxplot(tree_summary, x = "range", y = "DV_CUM_final", fill = "range", add = "jitter", alpha = 0.7) + 
  stat_compare_means(method = "anova", label = "p.signif") + 
  labs(title = "DV across Stress Levels", x = "Stress Level", y = "Cumulative DV (µm)")+ 
  theme( legend.position = "right",
         plot.title = element_text(face = "bold", size = 18, hjust = 0.5), 
         axis.title = element_text(face = "bold", size = 14), 
         axis.text = element_text(face = "bold", size = 12, color = "black"),
         legend.title = element_text(face = "bold", size = 12), 
         legend.text = element_text(face = "bold", size = 10), 
         #axis.line = element_line(color = "black", linewidth = 0.8), 
         panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5) )
