
##----------------------------------------------------------------------------------------------------------------------
#functional mechanism for tree water deficit

summer_2025 <- read.csv("/scratch/ope4/BIO_599/SUMMER/SUMMER_2025.csv")
summer_2024 <- read.csv("/scratch/ope4/BIO_599/SUMMER/SUMMER_2024.csv")

colnames(summer_2025)
colnames(summer_2024)
# Convert FW and DW to numeric if they are not
summer_2025$FW <- as.numeric(summer_2025$FW)
# Replace commas with dots (for decimal), convert empty strings to NA
summer_2025$DW <- as.numeric(ifelse(summer_2025$DW == "", NA, gsub(",", ".", summer_2025$DW)))
summer_2025$LDMC <- summer_2025$DW/summer_2025$FW
summer_2024$FW <- as.numeric(summer_2024$FW)
summer_2024$DW <- as.numeric(summer_2024$DW)
summer_2024$LDMC <- summer_2024$DW/summer_2024$FW
summer_2024$DBH <- summer_2024$DBH * 2.54


summer_2025 <- summer_2025 |>
  dplyr::select(-X.1, -X.2, -X.3, -X.4, -FW, -DW, -Species) |>
  dplyr::filter(!is.na(REPS))

summer_2024 <- summer_2024 |>
  dplyr::select(-TW, -FW, -DW, -Species) |>
  dplyr::filter(!is.na(REPS))

# Extract Tree_name and coordinates from summer_2025
coords_2025 <- summer_2025 %>%
  select(Tree_name, X, Y) %>%
  distinct()  # keep only one row per Tree_name

# Join X and Y into summer_2024
summer_2024 <- summer_2024 %>%
  left_join(coords_2025, by = "Tree_name")

# Row-bind the two datasets
merged_summer <- bind_rows(summer_2024, summer_2025)


##look at summary stat
summary(merged_summer)
# Look at the distribution of samples for each factor
table(merged_summer[ , c("SP_CODE", "PLOT")])

#EDA
# Select numeric columns that could be predictors
predictors <- merged_summer %>%
  select(DBH, Height, WP, LWC, RWC, LT, LDMC)
# Correlation matrix (pairwise, excluding NAs)
cor_matrix <- cor(predictors, use = "pairwise.complete.obs")
# Print the correlation matrix
print(round(cor_matrix, 2))

##get the absolute value of LWC
merged_summer$LWC <- abs(merged_summer$LWC)
merged_summer$WP <- -abs(merged_summer$WP)

# Select only Dendro_number and Trend from trend_all
trend_subset <- trend_all %>%
  select(Dendro_number, Trend, Slope) |>
  distinct()

# Join Trend into merged_summer by matching dendrometer numbers
merged_summer <- merged_summer %>%
  left_join(trend_subset, by = c("Dendrometer" = "Dendro_number")) 

merged_summer <- merged_summer %>%
  filter(!is.na(Trend))  %>%
  filter(Trend != "not enough data")
length(unique(merged_summer$Dendrometer))
summary(merged_summer)



# Summarise predictors by Year and Tree_name
merged_summary <- merged_summer %>%
  group_by(Tree_name, Time, Year) %>%
  summarise(
    PLOT = first(PLOT),
    SP_CODE = first(SP_CODE),
    Trend = first(Trend),
    Slope = first(Slope),
    Stress_Level = first(Stress_Level),
    DBH = mean(DBH, na.rm = TRUE),
    Height = mean(Height, na.rm = TRUE),
    LWC = mean(LWC, na.rm = TRUE),
    RWC = mean(RWC, na.rm = TRUE),
    LT = mean(LT, na.rm = TRUE),
    WP = mean(WP, na.rm = TRUE),
    LDMC = mean(LDMC, na.rm = TRUE),
    .groups = "drop"
  )

##scale all predictors
##scale the varaibes
merged_summary <- merged_summary %>%
  mutate(across(c(LWC, LDMC, DBH, Height, RWC, LT, WP),
                ~scale(.),
                .names = "z_{.col}")) 
summary(merged_summary)
##Make Trend a factor with Stable as reference
merged_summer$Trend <- factor(merged_summer$Trend)

# Remove rows with NAs in predictors (just like multinom did)
model_data <- merged_summary %>%
  filter(!is.na(z_DBH) & !is.na(z_Height) & !is.na(z_LWC) & !is.na(z_WP) & !is.na(z_LDMC) & !is.na(Trend))
summary(model_data)

numeric_data <- model_data %>%
  select(Slope, z_WP, z_DBH, z_LWC, z_LDMC, z_Height) %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%  # <- this unlists matrices
  drop_na()

ggpairs(
  numeric_data,
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.6))
)


model_data$Stress_Level <- factor(model_data$Stress_Level, ordered = TRUE)
gls_2 <- gls(Slope~z_WP + z_DBH + z_LWC + z_LDMC + z_Height + Stress_Level + SP_CODE, data = model_data )
check_model(gls_2)
#summary(gls_2)
H2 <- lme(
  Slope ~ z_WP + z_DBH + z_LWC + z_LDMC + z_Height + Stress_Level + SP_CODE,
  random = ~1 | PLOT,
  method = "REML",
  na.action = na.omit,
  data = model_data
)
summary(H2)
tab_model(H2)
check_model(H2)

#summary(M1.nested)
lme2_ml <- update(H2, method = "ML")
drop1(lme2_ml, test = "Chisq")
#Refit with REML and validate
#After selecting the optimal fixed structure, refit the model using REML for final coefficient estimates.
lme_final_2 <- update(lme2_ml, fixed = Slope ~ z_LWC + Stress_Level + SP_CODE,
                      method = "REML")
summary(lme_final_2)
tab_model(lme_final_2, show.se=TRUE)
check_model(lme_final_2)


