rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment


# Required packages

library(sf)
library(raster)
library(climateExtract)
library(terra)
library(magrittr)
library(dplyr)
library(data.table)
library(geodata)
library(tidyr)
library(progress)


# --- Charge data --- #
# ----
# eBMS data

setwd("D:/URBAN TRENDS/BMS data/BMS DATA 2024")  

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

# --- Merge data
head(m_coord)
head(country_codes)

m_coord_c <- merge(m_coord, country_codes, by = "bms_id", all.x = TRUE)

# ---

# ---

# Assuming 'm_coord' is your data frame and it's already been loaded into R

# Step 1: Convert to an sf object
# Replace EPSG:3035 with the correct EPSG code for ETRS89-extended / LAEA Europe if it's different
m_coord_sf <- st_as_sf(m_coord_c, coords = c("transect_lon", "transect_lat"), crs = 3035)

# Step 2: Transform to WGS 84
m_coord_sf_transformed <- st_transform(sf_point, crs = 4326)

# Extract the transformed coordinates back into the data frame
country_coord$lon_wgs84 <- st_coordinates(m_coord_sf_transformed)[,1]
country_coord$lat_wgs84 <- st_coordinates(m_coord_sf_transformed)[,2]

# View the results
head(m_coord)


# ---


setwd("D:/URBAN TRENDS/Climate data/Copernicus climate data")  


pb <- progress_bar$new(format = "[:bar] :percent :eta", total = length(unique(m_coord$country_code)))


for(country in unique(m_coord_c$country_code)){
  
  country_coord <- subset(m_coord_c, m_coord_c$country_code == country)
  
  sf_point <- st_as_sf(country_coord[!is.na(country_coord$transect_lon), ], coords = c("transect_lon", "transect_lat"), crs = 3035)
  
  m_coord_sf_transformed <- st_transform(sf_point, crs = 4326)
  country_coord$lon_wgs84 <- st_coordinates(m_coord_sf_transformed)[,1]
  country_coord$lat_wgs84 <- st_coordinates(m_coord_sf_transformed)[,2]
  
  grid_points<-st_make_grid(sf_point,cellsize = 100000 , square = TRUE,  what = "polygons")
  grid_points<- grid_points[sf_point]
  bms_bb<- st_bbox(st_transform(grid_points, 4326))
  bms_bbsf<- st_make_grid(bms_bb, n = 1, what = "polygons")
  bms_grid<- st_make_grid(bms_bb, cellsize = 25000 , square = TRUE,  what = "polygons")
  bms_bbsf<-st_as_sf(bms_bbsf)
  
  
  plot(sf_point$geometry)
  plot(bms_bbsf, add=T)
  plot(grid_points[sf_point], add=T, col = '#ff000088')
  
  #a character string for specific time period to be downloaded. 
  #Chunk available are "2011-2023", "1995-2010", "1980-1994", "1965-1979", and "1950-1964", but check what's available at
  
  ## minimum
  # 2000-2010
  climate_data = extract_nc_value(first_year = 2011, last_year = 2021,local_file = FALSE, file_path = NULL,sml_chunk = "2011-2023", spatial_extent = bms_bbsf,
                                  clim_variable = "min temp",statistic = "mean", grid_size = 0.25, ecad_v = NULL, write_raster = TRUE, out = "raster_min_temp1.tiff",
                                  return_data = TRUE)
  
  ## maximm
  # 2000-2010
  climate_data = extract_nc_value(first_year = 2011, last_year = 2021,local_file = FALSE,file_path = NULL,sml_chunk = "2011-2023",spatial_extent = bms_bbsf,
                                  clim_variable = "max temp",statistic = "mean", grid_size = 0.25, ecad_v = NULL,write_raster = TRUE,out = "raster_max_temp1.tiff",
                                  return_data = TRUE)
  
  rbk_min = terra::rast("raster_min_temp1.tiff")
  rbk_max = terra::rast("raster_max_temp1.tiff")
  
  be_gdd<- gdd_extract(base_temp = 5, min_temp = rbk_min, max_temp = rbk_max, gdd_method = 'be')
  
  tp_index <- get_layer_indice(x = be_gdd,
                               date_format = "%Y-%m-%d",
                               indice_level = "month")
  
  month_cumsum_gdd <- cumsum_rb(be_gdd, indices = tp_index)
  
  # get the indices for the last day of each month - calculate the number of day in each month along the time-series
  last_day_index <- as.numeric(table(as.numeric(factor(lubridate::floor_date(as.Date(names(month_cumsum_gdd), "%Y-%m-%d"), "month"))))) 
  
  last_day_gdd<- month_cumsum_gdd[[cumsum(last_day_index)]]
  
  # Initialize tibble
  gdd_df <- tibble(
    Year = character(),
    Month = character(),
    GDD5 = numeric(),
    SITE_ID = character(),
    bms_id = character()
  )
  
  for(site in unique(country_coord$transect_id)){
    
    sub_country_coord <- subset(country_coord, country_coord$transect_id == site)
    
    point <- unname(unlist(country_coord[1, 8:9]))
    
    values <- terra::extract(last_day_gdd, point) # Extract values at the specified point
    
    # Pivot the data frame to a long format
    values_long <- values %>%
      pivot_longer(
        cols = everything(), 
        names_to = "Date", 
        values_to = "GDD5"
      )%>%
      # Extract year and month from the Date
      # Assuming the date format is YYYY-MM-DD
      separate(Date, into = c("Year", "Month", "Day"), sep = "-") %>%
      # Select only the columns we want (exclude Day)
      dplyr::select(-Day) %>%
      # Arrange the data if necessary
      arrange(Year, Month)
    
    values_long$bms_id <- unique(sub_country_coord$bms_id)
    values_long$SITE_ID <- unique(sub_country_coord$transect_id)
    
    
    gdd_df<- rbind(gdd_df, values_long)
    
    
  }
  
  
  unique_country_code <- unique(country_coord$country_code)[1]
  
  # Create the filename string with the unique country code
  filename <- paste0("gdd_df_", unique_country_code, ".csv")
  
  # Use write_csv from the readr package to write the file
  write.csv(gdd_df, filename, row.names = FALSE)
  
  # Increment the progress bar
  pb$tick()
  
  
}









































