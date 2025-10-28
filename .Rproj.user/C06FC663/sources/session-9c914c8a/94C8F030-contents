##---------
# find covarites causing trend

# Load CSV
summer_2025 <- read.csv("/scratch/ope4/BIO_599/SUMMER/SUMMER_2025.csv")
colnames(summer_2025)

summer_2025 <- summer_2025 |>
  select(Dendrometer,SP_CODE, PLOT, Stress_Level) |>
  na.omit()

##join summer 2025 dataset to our full datset
summer_2025_clean <- summer_2025 %>%
  group_by(Dendrometer, SP_CODE, PLOT) %>%
  summarise(
    Stress_Level = first(Stress_Level),
    .groups = "drop"
  )

trend_summer_2025 <- trend_all %>%
  left_join(summer_2025_clean, by = c("Dendro_number" = "Dendrometer"))

##look at summary stat
summary(trend_summer_2025)
# Look at the distribution of samples for each factor
table(trend_summer_2025[ , c("SP_CODE", "PLOT")])

##check for collinairty
# Select the relevant columns
predictors <- trend_summer_2025[, c("min_temp", "max_temp", "soil_moisture_am", "evapotranspiration", "vpd", "precipitation")]

# Check pairwise correlations
cor_matrix <- cor(predictors, use = "complete.obs")  # exclude missing values
print(cor_matrix)

# Log-transform (common for right-skewed positive variables)
trend_summer_2025$z_max_temp <- scale(trend_summer_2025$max_temp)
trend_summer_2025$z_soil_moisture_am <- scale(trend_summer_2025$soil_moisture_am + 1)  # add 1 if zeros exist
trend_summer_2025$z_evapotranspiration <- scale(trend_summer_2025$evapotranspiration + 1)
trend_summer_2025$z_precipitation <- scale(trend_summer_2025$precipitation + 1)
trend_summer_2025$z_vpd <- scale(trend_summer_2025$vpd + 1)
summary(trend_summer_2025)

# Fit LME with tree nested in species and location as random effects
# You can start simple: random intercept per tree
M1.TWD <- lme(
  TWD ~ z_max_temp + z_soil_moisture_am + 
    z_evapotranspiration + z_precipitation + z_vpd,
  random = ~1 | Dendro_number,   # random intercept per tree
  method = "REML",
  na.action = na.omit,
  data = trend_summer_2025
)
summary(M1.TWD)
TWD_png <- tab_model(M1.TWD)
TWD_png


#DV
M1.DV <- lme(
  DV ~ z_max_temp + z_soil_moisture_am + 
    z_evapotranspiration + z_precipitation+ z_vpd,
  random = ~1 | Dendro_number,   # random intercept per tree
  method = "REML",
  na.action = na.omit,
  data = trend_summer_2025
)
summary(M1.DV)

DV_png <- tab_model(M1.DV)
DV_png


##look at how the species and plot play a role in TWD
##model mixed model for plot and sp
# TWD with same structure
trend_summer_2025$SP_CODE <- as.factor(trend_summer_2025$SP_CODE)
trend_summer_2025$PLOT <- as.factor(trend_summer_2025$PLOT)

#table 2
# Refit the model
M2_TWD <- lme(TWD ~ (z_max_temp + z_soil_moisture_am + 
                       z_evapotranspiration + z_precipitation+ z_vpd) * SP_CODE,
              random = ~1 | Dendro_number,
              data = trend_summer_2025,
              method = "REML",
              na.action = na.omit)
summary(M2_TWD)
tab_model(M2_TWD)
# Prepare the table
sig_table <- summary(M2_TWD)$tTable %>%
  as_tibble(rownames = "variable") %>%
  filter(`p-value` < 0.05) %>%
  mutate(
    # Add stars for significance
    signif = case_when(
      `p-value` < 0.001 ~ "***",
      `p-value` < 0.01  ~ "**",
      `p-value` < 0.05  ~ "*"
    ),
    # Round numeric columns
    Value     = round(Value, 0),
    Std.Error = round(Std.Error, 0),
    `t-value` = round(`t-value`, 3),
    `p-value` = round(`p-value`, 3),
    # Shorten variable names
    variable = str_replace(variable, "SP_CODE", ""),
  )

# Create flextable
ft <- flextable(sig_table[, c("variable", "Value", "Std.Error", "DF", "t-value", "p-value", "signif")]) %>%
  autofit()
ft
# Export to Word (can save as PDF from Word)
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "/scratch/ope4/BIO_599/Figure/Table_2.docx")

#Table 3
# Refit the model
M3_TWD <- lme(TWD ~ (z_max_temp + z_soil_moisture_am + 
                       z_evapotranspiration + z_precipitation+ z_vpd) * PLOT,
              random = ~1 | Dendro_number,
              data = trend_summer_2025,
              method = "REML",
              na.action = na.omit)
summary(M3_TWD)
tab_model(M3_TWD)
# Prepare the table
sig_table <- summary(M3_TWD)$tTable %>%
  as_tibble(rownames = "variable") %>%
  filter(`p-value` < 0.05) %>%
  mutate(
    # Add stars for significance
    signif = case_when(
      `p-value` < 0.001 ~ "***",
      `p-value` < 0.01  ~ "**",
      `p-value` < 0.05  ~ "*"
    ),
    # Round numeric columns
    Value     = round(Value, 0),
    Std.Error = round(Std.Error, 0),
    `t-value` = round(`t-value`, 3),
    `p-value` = round(`p-value`, 3),
    # Shorten variable names
    variable = str_replace(variable, "PLOT", ""),
  )

# Create flextable
ft <- flextable(sig_table[, c("variable", "Value", "Std.Error", "DF", "t-value", "p-value", "signif")]) %>%
  autofit()
ft
# Export to Word (can save as PDF from Word)
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "/scratch/ope4/BIO_599/Figure/Table_3.docx")


