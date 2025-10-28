library(sf)
library(zoo)
library(Kendall)
library(trend)
library(lme4)
library(sjPlot)
library(nlme)
library(effects)
library(ggeffects)
library(RColorBrewer)
#library(flextable)
library(knitr)
library(officer)
library(webshot2)
#install.packages("systemfonts")
library(nnet)
library(MASS)
library(performance)
library(GGally)
library(tidyverse)


# 1. List all CSV files in the folder
files <- list.files(
  path = "/scratch/ope4/BIO_599/dendrometers/",
  pattern = "_data_.*\\.csv$",
  full.names = TRUE
)

# 2. Function to read a file and add dendro_number
read_dendro_file <- function(file) {
  # Read the CSV
  sub <- read.csv(file, sep = ";", header = FALSE, stringsAsFactors = FALSE)
  
  # Extract dendro number from filename
  dendro_number <- as.numeric(gsub(".*_data_([0-9]+)_.*", "\\1", file))
  
  # Add column
  sub$Dendro_number <- dendro_number
  
  return(sub)
}

# 3. Apply to all files and combine into one data frame
all_data <- lapply(files, read_dendro_file) %>% bind_rows()

# select col of interest
combined_data <- all_data |>
  select(Dendro_number, V2, V4, V7) |>
  filter(V7 != 0) 


##rename the cols
colnames(combined_data) <- c("Dendro_number", "DateTime", "Temp", "dendro")

summary(combined_data)

# Convert DateTime to POSIXct if not already
combined_data$DateTime <- as.POSIXct(combined_data$DateTime, format = "%Y.%m.%d %H:%M")
combined_data$Temp <- as.numeric(combined_data$Temp)  # safe if it's already numeric
combined_data$Dendro_number <- as.numeric(combined_data$Dendro_number)  # safe if it's already numeric
str(combined_data)
head(combined_data$DateTime, 20)


# 1. Compute the maximum per dendrometer across all time
dendro_max <- combined_data %>%
  group_by(Dendro_number) %>%
  summarise(dmax = max(dendro, na.rm = TRUE), .groups = 'drop')

# 2. Join dmax and calculate TWD for each measurement
combined_data <- combined_data %>%
  left_join(dendro_max, by = "Dendro_number") %>%
  mutate(TWD = dmax - dendro)  # TWD for each 15-min measurement

daily_variation <- combined_data %>%
  mutate(Date = as.Date(DateTime)) %>%       # add date-only column
  group_by(Dendro_number, Date) %>%         # group by dendrometer and day
  summarise(
   min_dendro = min(dendro, na.rm = TRUE),
    max_dendro = max(dendro, na.rm = TRUE),
    DV = max_dendro - min_dendro,
    TWD = mean(TWD, na.rm = TRUE),
   min_temp = min(Temp, na.rm = TRUE),
   max_temp = max(Temp, na.rm = TRUE),
   temp_range = max_temp - min_temp,
     ) %>%
  ungroup()

summary(daily_variation)


length(unique(combined_data$Dendro_number))
length(unique(daily_variation$Dendro_number))

###---------------------

#add the climate dataset
precip <- read_csv("/scratch/ope4/BIO_599/Dendro Climate Points/precipitation_time_series.csv")

summary(precip)
length(unique(precip$Dndrmtr))

# Step 2: Prepare precipitation (keep only relevant columns)
precip_daily <- precip %>%
  select(Dndrmtr, date, precipitation, vpd) %>%
  rename(Dendro_number = Dndrmtr,
         Date = date)

# Find date range of dendrometer data
date_range <- range(daily_variation$Date, na.rm = TRUE)

# Clip precip to fall within that range
precip_clipped <- precip_daily %>%
  filter(Date >= date_range[1] & Date <= date_range[2])

####evapo
#add the climate dataset
evapo <- read_csv("/scratch/ope4/BIO_599/Dendro Climate Points/evapotranspiration_time_series.csv")

summary(evapo)
length(unique(evapo$Dndrmtr))

##interpolate to daily instead of every 8 days
# Make sure date is Date type
evapo <- evapo %>%
  mutate(date = as.Date(date))

# Interpolate daily for each dendrometer
evapo_daily <- evapo %>%
  group_by(Dndrmtr) %>%
  arrange(date) %>%
  # Create full daily sequence for each dendrometer
  complete(date = seq(min(date), max(date), by = "1 day")) %>%
  # Interpolate evapotranspiration
  mutate(evapotranspiration = na.approx(evapotranspiration, x = date, maxgap = Inf)) %>%
  ungroup()

# Step 2: Prepare precipitation (keep only relevant columns)
evapo_daily <- evapo_daily %>%
  select(Dndrmtr, date, evapotranspiration) %>%
  rename(Dendro_number = Dndrmtr,
         Date = date)

# Find date range of dendrometer data
date_range <- range(daily_variation$Date, na.rm = TRUE)

# Clip precip to fall within that range
evapo_clipped <- evapo_daily %>%
  filter(Date >= date_range[1] & Date <= date_range[2])


#### soil moisture
#add the climate dataset
soil <- read_csv("/scratch/ope4/BIO_599/Dendro Climate Points/soil_moisture_time_series.csv")

summary(soil)
length(unique(soil$Dndrmtr))

# Step 2: Prepare precipitation (keep only relevant columns)
soil <- soil %>%
  select(Dndrmtr, date, soil_moisture_am) %>%
  rename(Dendro_number = Dndrmtr,
         Date = date)

# Find date range of dendrometer data
date_range <- range(daily_variation$Date, na.rm = TRUE)

# Clip precip to fall within that range
soil_clipped <- soil %>%
  filter(Date >= date_range[1] & Date <= date_range[2])

######### merge all climate dataset---------------------------
merged <- reduce(list(soil_clipped, evapo_clipped, precip_clipped), left_join, by = c("Dendro_number", "Date"))

colnames(merged)

# Step 3: Merge
merged_data <- daily_variation %>%
  left_join(merged, by = c("Dendro_number", "Date"))


##step 4: test for trend in daily variation of dendrometer
## Step 4: test for trend in daily variation of dendrometer
check_trend <- function(x, alpha = 0.05) {
  # Drop NAs
  x <- x[!is.na(x)]
  
  # Require at least 14 points for decompose (2 periods if freq=7)
  if (length(x) < 14) {
    return(data.frame(Trend = "not enough data", Slope = NA))
  }
  
  # Turn into time series
  ts <- ts(x, frequency = 7)
  
  # Perform additive decomposition
  decomp_add <- decompose(ts, type = "additive")
  
  # Use the trend component (remove NAs first)
  trend_vals <- decomp_add$trend
  trend_vals <- trend_vals[!is.na(trend_vals)]
  
  if (length(trend_vals) < 3) {
    return(data.frame(Trend = "not enough trend data", Slope = NA))
  }
  
  mk <- trend::mk.test(trend_vals)
  sen <- trend::sens.slope(trend_vals)
  slope <- sen$estimates
  
  # Determine trend direction
  if (mk$p.value < alpha) {
    if (slope > 0) {
      return(data.frame(Trend = "increasing", Slope = slope))
    } else if (slope < 0) {
      return(data.frame(Trend = "decreasing", Slope = slope))
    }
  }
  return(data.frame(Trend = "stable", Slope = slope))
}


## Apply by Dendro_number
trend_results <- merged_data %>%
  group_by(Dendro_number) %>%
  summarise(
    tmp = list(check_trend(TWD)),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(tmp)  # expands Trend and Slope into separate columns

table(trend_results$Trend)

## Merge back with main data
trend_all <- merged_data %>%
  left_join(trend_results, by = "Dendro_number") %>%
  mutate(Trend = Trend)  # keep trend label for each record
#

##---------------------------------------------------------------------------------------------------------------------------
# Load CSV
df <- read.csv("/scratch/ope4/BIO_599/SUMMER_2025.csv")

table(df$Dendrometer)

df <- df |>
  select(Dendrometer,Tree_name, X, Y) |>
  na.omit()

# Keep only unique rows
df_unique <- df[!duplicated(df), ]

# Convert to sf object using X (longitude) and Y (latitude)
gdf <- st_as_sf(df_unique, coords = c("X", "Y"), crs = 4326)  # WGS84

# Export to shapefile
st_write(gdf, "/scratch/ope4/BIO_599/SHP/SUMMER_2025.shp")
