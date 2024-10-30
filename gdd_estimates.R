# --- Calculation of GDD across sites and years --- #


rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment


# Required packages
#----
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
#----

# --- Charge data --- #
#----

# eBMS data

setwd("D:/URBAN TRENDS/BMS data/BMS DATA 2024")  

ebms_coord_df <- read.csv("ebms_transect_coord.csv", sep = ",", dec = ".")
# Import country codes
country_codes <- read.csv("country_codes.csv", sep = ";", dec = ".")

# uBMS data
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
m_coord<- rbind(ebms_coord_df, ubms_coord_final)

#Merge data
head(m_coord)
head(country_codes)

m_coord_c <- merge(m_coord, country_codes, by = "bms_id", all.x = TRUE)
m_coord_c<- na.omit(m_coord_c)
#----



# --- Loop --- #


# Set a permanent directory for output files
permanent_dir <- "D:/URBAN TRENDS/Climate data/Copernicus climate data"

Period <- c("2011-2023", "1995-2010", "1980-1994", "1965-1979")
First_year<- c(2011, 1995, 1980, 1976)
Last_year <- c(2023, 2010, 1994, 1979)

subset_m_coord_c <- m_coord_c[m_coord_c$Country.Name %in% c("Spain", "Belgium"), ]
head(subset_m_coord_c)



for(country in unique(m_coord_c$country_code)){
  
  # Create and set a new temporary directory for this iteration
  temp_dir <- tempdir()
  setwd(temp_dir)
  
  country_coord <- subset(m_coord_c, m_coord_c$country_code == country)
  
  # Output a message indicating the start of processing for the current country
  cat("Processing country: ", unique(country_coord$Country.Name), "\n")
  
  sf_point <- st_as_sf(country_coord[!is.na(country_coord$transect_lon), ], coords = c("transect_lon", "transect_lat"), crs = 3035)
  
  m_coord_sf_transformed <- st_transform(sf_point, crs = 4326)
  country_coord$lon_wgs84 <- st_coordinates(m_coord_sf_transformed)[,1]
  country_coord$lat_wgs84 <- st_coordinates(m_coord_sf_transformed)[,2]
  
  grid_points<-st_make_grid(sf_point,cellsize = 100000 , square = TRUE,  what = "polygons")
  grid_points<- grid_points[sf_point]
  bms_bb<- st_bbox(st_transform(grid_points, 4326))
  bms_bbsf<- st_make_grid(bms_bb, n = 1, what = "polygons")
  bms_grid<- st_make_grid(bms_bb, cellsize = 10000 , square = TRUE,  what = "polygons") # replace cellsize
  bms_bbsf<-st_as_sf(bms_bbsf)
  
  
  plot(sf_point$geometry)
  plot(bms_bbsf, add=T)
  plot(grid_points[sf_point], add=T, col = '#ff000088')
  
  
  
  for(i in seq(1:4)){
    
    Period_chunk <- Period[i]
    First_year_chunk <- First_year[i]
    Last_year_chunk <- Last_year[i]
    
    cat("           Processing period: ", Period_chunk, "\n")
    
    
    ## minimum
    climate_data = extract_nc_value(first_year = First_year_chunk, last_year = Last_year_chunk,local_file = FALSE, file_path = NULL,sml_chunk = Period_chunk, spatial_extent = bms_bbsf,
                                    clim_variable = "min temp",statistic = "mean", grid_size = 0.25, ecad_v = NULL, write_raster = TRUE, out = "raster_min_temp1.tiff",
                                    return_data = TRUE)
    
    ## maximum
    climate_data = extract_nc_value(first_year = First_year_chunk, last_year = Last_year_chunk,local_file = FALSE,file_path = NULL,sml_chunk = Period_chunk,spatial_extent = bms_bbsf,
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
      site_location <- vect(cbind(point[1], point[2]), crs="EPSG:4326")
      
      values <- terra::extract(last_day_gdd, site_location) # Extract values at the specified point
      
      # Pivot the data frame to a long format
      values_long <- values %>%
        pivot_longer(
          cols = -ID, # Exclude the ID column
          names_to = "Date", 
          values_to = "GDD5"
        ) %>%
        # Extract year and month from the Date
        # Assuming the date format is YYYY-MM-DD
        separate(Date, into = c("Year", "Month", "Day"), sep = "-") %>%
        # Select only the columns we want (exclude Day and ID)
        dplyr::select(-Day, -ID) %>%
        # Arrange the data if necessary
        arrange(Year, Month)
      
      values_long$bms_id <- unique(sub_country_coord$bms_id)
      values_long$SITE_ID <- unique(sub_country_coord$transect_id)
      
      
      gdd_df<- rbind(gdd_df, values_long)
      
      
    }
    
    # Create the filename string with the unique country code
    filename <- paste0(permanent_dir, "/", "gdd5_025_", country, "_", Period_chunk, ".csv")  # Added "/" separator
    
    # Use write_csv from the readr package to write the file
    write.csv(gdd_df, filename, row.names = FALSE)
    
  }
    
  
}




# --- Data compilation --- #

setwd("E:/URBAN TRENDS/Climate data/Copernicus climate data")


# Get a list of R files in the working directory starting with "gdd5_"
file_list <- list.files(pattern = "^gdd5_025_", full.names = TRUE)

# Read each CSV file and bind them row-wise
combined_data <- do.call(rbind, lapply(file_list, function(file) {
  read.csv(file, header = TRUE)  # Use read.csv for CSV files
}))


head(combined_data)



# Create new column names for GDD5
new_col_names <- paste0("GDD5_", 1:12)

# Pivot the data wider to create separate columns for each Month's GDD5 value
combined_data <- combined_data %>%
  pivot_wider(names_from = Month, values_from = GDD5, names_prefix = "GDD5_") %>%
  mutate(across(starts_with("GDD5"), as.numeric))  # Convert GDD5 columns to numeric if needed

# Save the new dataframe as a CSV file in the working directory
write.csv(combined_data, "gdd5_025_compiled.csv", row.names = FALSE)



library(dplyr)

# Count NA rows for each 'bms_id' and total number of rows
na_rows_by_bms_id <- combined_data %>%
  group_by(bms_id) %>%
  summarize(total_rows = n(),
            num_na_rows = sum(is.na(GDD5)))

# Print the result
print(na_rows_by_bms_id, n = 22)








