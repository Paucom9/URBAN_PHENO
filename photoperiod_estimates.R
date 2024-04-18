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
library(sf)
library(suncalc)
library(sp)
library(progress)
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
#----


# --- Photoperiod calculation --- #

head(m_coord)

# Define the EPSG:3035 projection
proj <- CRS("+init=epsg:3035")

# Function to convert coordinates from EPSG:3035 to latitude and longitude
#----
convert_coordinates <- function(lon, lat) {
  # Create a SpatialPoints object
  coordinates <- SpatialPoints(matrix(c(lon, lat), ncol = 2), proj4string = proj)
  # Convert coordinates to geographical coordinate system (latitude and longitude)
  geographical_coordinates <- spTransform(coordinates, CRS("+proj=longlat"))
  # Extract converted coordinates
  converted_coordinates <- coordinates(geographical_coordinates)
  return(converted_coordinates)
}
#----

# Function to calculate the photoperiod at the beginning of each season
#----
calculate_photoperiod <- function(latitude, longitude, date) {
  # Convert coordinates to latitude and longitude
  converted_coordinates <- convert_coordinates(longitude, latitude)
  converted_latitude <- converted_coordinates[2]
  converted_longitude <- converted_coordinates[1]
  
  # Calculate the position of the sun at sunrise and sunset
  sun_position <- getSunlightTimes(date = date, lat = converted_latitude, lon = converted_longitude)
  
  # Extract sunrise and sunset times
  sunrise <- sun_position$times[1]
  sunset <- sun_position$times[3]
  
  # Calculate the duration of the day (photoperiod) in minutes
  photoperiod <- as.numeric(sunset - sunrise, units = "minutes")
  
  return(photoperiod)
}
#----

# Initialize a progress bar
pb <- progress_bar$new(format = "[:bar] :percent :eta", total = length(unique(m_coord$transect_id)))

# Create an empty dataframe to store the results
photoperiod_df <- data.frame(bms_id = character(),
                             SITE_ID = character(),
                             latitude = numeric(),
                             photo_spring = numeric(),
                             photo_summer = numeric(),
                             photo_autumn = numeric(),
                             stringsAsFactors = FALSE)

# Define dates for the spring equinox, summer solstice, and autumn equinox
spring_equinox<-as.Date("2024-03-21")
summer_equinox<-as.Date("2024-06-21")
autumn_equinox<- as.Date("2024-09-21")

# Loop over unique transect IDs
for (bms_id in unique(m_coord$transect_id)) {
  subset_data <- m_coord[m_coord$transect_id == bms_id, ]
  
  # Check if any NA values exist in transect_lon or transect_lat columns
  if (any(is.na(subset_data$transect_lon)) || any(is.na(subset_data$transect_lat))) {
    # If NA values exist, assign NA to photoperiod variables and skip calculation
    photo_spring <- NA
    photo_summer <- NA
    photo_autumn <- NA
    latitude <- NA
    
  } else {
  
  # Convert coordinates to latitude and longitude
  geographical_coordinates <- convert_coordinates(subset_data$transect_lon[1], subset_data$transect_lat[1])
  latitude <- geographical_coordinates[1, 2]
  
  # Calculate photoperiod
  sun_times_spr <- getSunlightTimes(spring_equinox, latitude, 0)  # Using longitude 0 for simplicity
  sun_times_summ <- getSunlightTimes(summer_equinox, latitude, 0)  # Using longitude 0 for simplicity
  sun_times_aut <- getSunlightTimes(autumn_equinox, latitude, 0)  # Using longitude 0 for simplicity
  
  photo_spring <- as.numeric(difftime(sun_times_spr$sunset, sun_times_spr$sunrise, units = "hours"))
  photo_summer <- as.numeric(difftime(sun_times_summ$sunset, sun_times_summ$sunrise, units = "hours"))
  photo_autumn <- as.numeric(difftime(sun_times_aut$sunset, sun_times_aut$sunrise, units = "hours"))
  
  }
  
  # Append the results to the dataframe
  photoperiod_df <- rbind(photoperiod_df, 
                          data.frame(bms_id = subset_data$bms_id,
                                     SITE_ID = subset_data$transect_id[1],
                                     latitude = latitude,
                                     stringsAsFactors = FALSE,
                                     photo_spring = photo_spring,
                                     photo_summer = photo_summer,
                                     photo_autumn = photo_autumn))
  
  
  # Increment the progress bar
  pb$tick()
  
}

head(photoperiod_df)

# Specify the file path and name
file_path <- "D:/URBAN TRENDS/Urban_pheno/photo_estimates.csv"

# Save as a CSV file
write.csv(photoperiod_df, file = file_path, row.names = FALSE)


