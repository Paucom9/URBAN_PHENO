# Clean environment and set the working directory to the location of the data files
rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment
setwd("D:/URBAN TRENDS/BMS data/BMS DATA 2024")  

# Load required libraries
# ----
library(data.table)  # For efficient data handling
library(mgcv)        # For generalized additive models
library(dplyr)       # For data manipulation
library(tidyr)       # For data tidying
library(broom)       # For converting statistical analysis objects into tidy data frames
library(stringr)     # For string manipulation
library(lubridate)   # For easy and intuitive work with dates and times
library(doParallel)  # For increasing loop performance
library(changepoint) # For change point analyses
# ----


# ---- Data Import and Preparation --- #
#----
# eBMS data

# Import butterfly count data
ebms_count_df <- read.csv("ebms_count.csv", sep = ",", dec = ".")
# Import visit data
ebms_visit_df <- read.csv("ebms_visit.csv", sep = ",", dec = ".")
# Import climate region data
ebms_clim_df <- read.csv("ebms_transect_climate.csv", sep = ",", dec = ".")
# Import transect coordinates
ebms_coord_df <- read.csv("ebms_transect_coord.csv", sep = ",", dec = ".")
# Import country codes
country_codes <- read.csv("country_codes.csv", sep = ";", dec = ".")

# Extract bms_id from transect_id and select relevant columns
ebms_clim_df <- ebms_clim_df %>%
  mutate(bms_id = str_extract(transect_id, "^[^.]*")) %>%
  dplyr::select(bms_id, transect_id, genzname)

# uBMS data
# Import count data
ubms_count_df <- read.csv("output_count_table.csv", sep = ",", dec = ".")
# Assign a unique identifier to uBMS data
ubms_count_df$bms_id <- "ES-uBMS"
# Select and reorder columns to match the structure of eBMS count data
ubms_count_df <- ubms_count_df %>%
  dplyr::select(visit_id, bms_id, transect_id, visit_date, year, month, day, species_name, count)

# Import visit data and perform necessary transformations
ubms_visit_df <- read.csv("raw_ubms_ebms_visit.csv", sep = ",", dec = ".")
ubms_visit_df$bms_id <- "ES-uBMS"  # Assign the uBMS identifier

# Rename columns and calculate date components
ubms_visit_df <- ubms_visit_df %>%
  rename(visit_date = date_of_visit) %>%
  mutate(
    visit_date = ymd(visit_date),
    year = year(visit_date),
    month = month(visit_date),
    day = day(visit_date),
    week = week(visit_date),
    ebms_partner = TRUE
  )

# Adjust column order to match eBMS visit data structure
ubms_visit_df <- ubms_visit_df %>%
  dplyr::select(visit_id, bms_id, transect_id, visit_date, year, month, day, ebms_partner, week)

# Import transect cooordinates
ubms_coord <- read.csv("ubms_sites.csv", sep = ";", dec = ".")

ubms_coord <- ubms_coord %>%mutate(bms_id = "ES-uBMS")

# Filter out rows with NA in coordinates before transformation
ubms_coord_filtered <- ubms_coord %>%
  filter(!is.na(transect_longitude) & !is.na(transect_latitude))

# Convert to an sf object
ubms_coord_sf <- st_as_sf(ubms_coord_filtered, coords = c("transect_longitude", "transect_latitude"), crs = 4326, remove = FALSE)

# Transform coordinates to EPSG:3035
ubms_coord_transformed <- st_transform(ubms_coord_sf, crs = 3035)

ubms_coord_final <- data.frame(ubms_coord_transformed) %>%
  transmute(
    bms_id,
    transect_id,
    transect_length,
    transect_lon = st_coordinates(ubms_coord_transformed)[, 1],
    transect_lat = st_coordinates(ubms_coord_transformed)[, 2]
  )

## Concatenate rows

m_count_df <- rbind(ebms_count_df, ubms_count_df)
m_visit_df<- rbind(ebms_visit_df, ubms_visit_df)
m_clim_df<- ebms_clim_df
m_coord<- rbind(ebms_coord_df, ubms_coord_final)

## Transform data frames to data tables

m_count <- data.table(m_count_df)
m_visit <- data.table(m_visit_df)
m_clim <- data.table(m_clim_df)
dt_country_cod <- data.table(country_codes)

## Change column names

setnames(m_visit, c('transect_id', 'visit_date'), c('SITE_ID', 'DATE'))
setnames(m_count, c('transect_id', 'visit_date','species_name', 'count'),
         c('SITE_ID', 'DATE', 'SPECIES', 'COUNT'))
setnames(m_clim, c('transect_id', 'genzname'),
         c('SITE_ID', 'RCLIM'))

# Perform a left join to add RCLIM from m_clim to m_visit based on bms_id and SITE_ID
m_visit <- m_visit[m_clim, on = .(bms_id, SITE_ID), nomatch = 0]

## Perform a left join to merge m_clim into m_count

m_count <- left_join(m_count, m_clim, by = c("SITE_ID", "bms_id"))

# Merge m_count with dt_country_cod to include country_code
m_count <- merge(m_count, dt_country_cod, by = "bms_id", all.x = TRUE)

#NA values in RCLIM correspond to uBMS transects. All of them are part of climate region K
m_count$RCLIM <- replace(m_count$RCLIM, is.na(m_count$RCLIM), "K. Warm temperate and mesic")

# Create and ID factor
m_count$ID <- paste(m_count$SPECIES, m_count$SITE_ID, sep = "_")

# Year as factor
m_count$year <- as.factor(m_count$year)
m_visit$year <- as.factor(m_visit$year)

# Check count and visit data
head(m_count)
head(m_visit)

# --- Filter data by univoltine species --- #
#----

setwd("D:/URBAN TRENDS/Traits data")  
traits <- read.csv("Traits_table.csv", sep = ";", dec = ".")

traits <- traits %>% rename(SPECIES = Taxon) %>% mutate(SPECIES = gsub("_", " ", SPECIES))

# Identify univoltine species
uni_sp <- traits %>%
  filter(Vol_max == 1.0) %>%
  pull(SPECIES)

# Filter m_count for unvoltine species
m_count_univol <- m_count %>%
  filter(SPECIES %in% uni_sp)
#----

# --- Calculate pheno estimates for univoltines butterflies --- #

# Create an anchor argument to add zeros (0) before and after the monitoring season 
# This ensures that the flight curve starts and end at zero

anchor <- data.table(
  bms_id = rep(NA_character_, each = 151),
  visit_id = rep(NA_integer_, each = 151),
  SITE_ID = rep(NA_character_, each = 151),
  DATE = rep(as.Date(NA), each = 151),
  year = rep(NA_integer_, each = 151),
  month = rep(NA_integer_, each = 151),
  day = rep(NA_integer_, each = 151),
  SPECIES = rep(NA_character_, each = 151),
  COUNT = rep(0, each = 151),
  RCLIM = rep(NA_character_, each = 151),
  country_code = rep(NA_character_, each = 151),
  Country.Name = rep(NA_character_, each = 151),
  ID = rep(NA_character_, each = 151),
  ebms_partner = rep(NA, each = 151),
  week = rep(NA_integer_, each = 151),
  julian_day = c(rep(1:59, times = n_distinct(all_counts$bms_id)), rep(274:365, times = n_distinct(all_counts$bms_id)))
)

# Initialize an empty data frame to store phenology estimates
phenology_estimates <- data.frame(ID = numeric(), 
                                  YEAR = numeric(), 
                                  SPECIES = character(), 
                                  SITE_ID = character(), 
                                  ONSET_mean = numeric(),
                                  ONSET_var = numeric(), 
                                  OFFSET_mean = numeric(), 
                                  OFFSET_var = numeric(), 
                                  PEAKDAY = numeric(), 
                                  FLIGHT_LENGTH_mean = numeric(),
                                  FLIGHT_LENGTH_var = numeric(), 
                                  stringsAsFactors = FALSE)


# Create a progress bar with an estimated total count
pb <- progress_estimated(length(unique(m_count_univol$ID)))

for(id in unique(m_count_univol$ID)){
  
  sub_count <- m_count_univol[ID == id]
  site <- unique(sub_count$SITE_ID)
  
  # Update the progress bar
  pb$tick()$print()
  
  for(YEAR in unique(sub_count$year)){
    
    sub_count_year<-  sub_count[year == YEAR]
    sub_visit <- m_visit[year == YEAR & SITE_ID == site]
    
    if(nrow(sub_count_year) >= 3 & nrow(sub_visit) >= 10){
      
      # Identify dates in sub_visit not present in sub_count_year
      missing_dates <- sub_visit[!sub_visit$DATE %in% sub_count_year$DATE, ]
      
      # Create new rows in sub_count_year for missing dates with COUNT = 0
      missing_rows <- missing_dates %>%
        mutate(COUNT = 0)
      
      # Include zero counts in the species*site df
      all_counts <- bind_rows(sub_count_year, missing_rows)
      
      # Convert DATE column to Date class
      all_counts <- all_counts %>%
        mutate(DATE = as.Date(DATE))
      
      # Calculate Julian day
      all_counts <- all_counts %>%
        mutate(julian_day = yday(DATE))
      
      # Bind the anchor rows with the counts
      
      all_counts <- rbind(all_counts, anchor)
      
      tryCatch({
      
      # Fit a GAM model
      gam_model <- gam(COUNT ~ s(julian_day), data = all_counts, family = nb)
      
      if (inherits(gam_model, "try-error")) {
      peak_day <- NA
      onset_var <- NA
      onset_mean <- NA
      offset_var <- NA
      offset_mean <- NA
      flight_length_var <- NA
      flight_length_mean <- NA
      
      } else {
      
      # Create a sequence of Julian days from 1 to 365
      julian_days <- 1:365
      
      # Create a sequence of julian_day from 1 to 365
      new_data <- data.frame(julian_day = julian_days)
      
      # Predict COUNT using the GAM model
      predict_count <- predict(gam_model, newdata = new_data, type = "response")
      
      pheno <- data.frame(Julian_Day = julian_days, Predicted_Count = predict_count)
      
      # --- Phenology estimates --- #
      
      peak_day <- which.max(pheno$Predicted_Count) # Day of the flight curve in which the abundance is maximum
      
      # Identify onset and offset of butterfly activity using change point analysis
      cp_mean <- cpt.mean(pheno$Predicted_Count, method = "PELT", penalty = "Manual", pen.value = 0.05)
      cp_var <- cpt.var(pheno$Predicted_Count, method = "PELT", penalty = "Manual", pen.value = 0.05)
      
      # Extract the change points
      change_points_mean <- cpts(cp_mean)
      change_points_var <- cpts(cp_var)
      
      if (inherits(change_points_mean, "try-error")) {
        onset_mean <- NA
        offset_mean <- NA
        flight_length_mean <- NA
      } else {
        onset_mean <- min(change_points_mean) # First day of appearance with cp_mean method
        offset_mean <- max(change_points_mean) # Last day of appearance with cp_var method
        flight_length_mean <- length(onset_mean:offset_mean)   # Length of the flight period
      }
      
      if (inherits(change_points_var, "try-error")) {
        onset_var <- NA
        offset_var <- NA
        flight_length_var <- NA
      } else {
        onset_var <- min(change_points_var) # First day of appearance with cp_mean method
        offset_var <- max(change_points_var) # Last day of appearance with cp_var method
        flight_length_var <- length(onset_var:offset_var)   # Length of the flight period
      }
      
      
      }
      
      }, error = function(e) {

        print(paste("Error for ", id, " in ", YEAR, ":", e))
        
      })
      
      # Save the results into the data frame
      phenology_estimates <- rbind(phenology_estimates, 
                                   data.frame(ID = id, 
                                              YEAR = YEAR, 
                                              SPECIES = unique(sub_count_year$SPECIES), 
                                              SITE_ID = unique(sub_count_year$SITE_ID), 
                                              ONSET_mean = onset_mean,
                                              ONSET_var = onset_var, 
                                              OFFSET_mean = offset_mean, 
                                              OFFSET_var = offset_var, 
                                              PEAKDAY = peak_day, 
                                              FLIGHT_LENGTH_mean = flight_length_mean,
                                              FLIGHT_LENGTH_var = flight_length_var))
      
      #cat(sprintf("Phenology estimates calculated for %s in year %s\n", id, YEAR))
      
      
    } else {
      
      #cat(sprintf("Not enough data to calculate phenology estimates for %s in year %s\n", id, YEAR))
      
    }

  }
  
}


#----

head(phenology_estimates)

# Specify the file path and name
file_path <- "D:/URBAN TRENDS/Urban_pheno/pheno_estimates.csv"

# Save as a CSV file
write.csv(phenology_estimates, file = file_path, row.names = FALSE)




