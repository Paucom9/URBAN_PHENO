rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment

#
#----
# Load necessary libraries
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
library(tidyr)
library(lubridate)
#----

# Function to calculate GDDs
#----

# Función para extraer y sumar GDD por periodo para cada transect_id y año
extract_gdd_periods <- function(transect_id, lon, lat, be_gdd, periods) {
  # Extraer las coordenadas del sitio (lon, lat) en be_gdd
  coords <- data.frame(x = lon, y = lat)
  
  # Extraer los valores de GDD del raster be_gdd para las coordenadas dadas
  site_gdd <- terra::extract(be_gdd, coords)
  
  # Obtener las fechas del raster be_gdd y convertirlas en días julianos
  dates <- as.Date(names(be_gdd))  # Extraer los nombres de las capas (fechas)
  day_of_year <- as.numeric(format(dates, "%j"))  # Convertir a días julianos (1-365)
  
  # Inicializar un dataframe para almacenar los resultados
  gdd_sums <- data.frame()
  
  # Loop a través de los años disponibles en be_gdd
  for (year in unique(format(dates, "%Y"))) {
    # Filtrar las capas del año actual
    year_idx <- which(format(dates, "%Y") == year)
    gdd_year <- site_gdd[, year_idx]
    day_of_year_year <- day_of_year[year_idx]  # Días julianos de ese año
    
    # Sumar GDD para cada periodo definido (late winter, spring, summer, autumn)
    gdd_latewinter <- sum(gdd_year[day_of_year_year %in% periods$late_winter], na.rm = TRUE)
    gdd_spring <- sum(gdd_year[day_of_year_year %in% periods$spring], na.rm = TRUE)
    gdd_summer <- sum(gdd_year[day_of_year_year %in% periods$summer], na.rm = TRUE)
    gdd_autumn <- sum(gdd_year[day_of_year_year %in% periods$autumn], na.rm = TRUE)
    
    # Agregar una fila al dataframe con los resultados para este año y sitio
    gdd_sums <- rbind(gdd_sums, data.frame(
      transect_id = transect_id,
      year = year,
      GDD_latewinter = gdd_latewinter,
      GDD_spring = gdd_spring,
      GDD_summer = gdd_summer,
      GDD_autumn = gdd_autumn
    ))
  }
  
  return(gdd_sums)
}


#----


# --- Charge data --- #
#----

# eBMS data

setwd("E:/URBAN TRENDS/BMS data/BMS DATA 2024")  

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

# Pheno data
setwd("E:/URBAN TRENDS/Urban_pheno")  
pheno_df <- read.csv("pheno_estimates.csv", sep = ",", dec = ".")
pheno_df$YEAR <- as.factor(pheno_df$YEAR)
pheno_df$SITE_ID <- as.factor(pheno_df$SITE_ID)

head(pheno_df)


# Calculate the mean values of pheno variables for each ID, SPECIES, and SITE_ID
pheno_mean_dates <- pheno_df %>%
  group_by(ID, SPECIES, SITE_ID) %>%
  summarise(
    ONSET_mean_avg = mean(ONSET_mean, na.rm = TRUE),
    OFFSET_mean_avg = mean(OFFSET_mean, na.rm = TRUE),
    FLIGHT_LENGTH_mean_avg = mean(FLIGHT_LENGTH_mean, na.rm = TRUE),
    PEAKDAY_avg = mean(PEAKDAY, na.rm = TRUE)
  )

# View the summary
print(pheno_mean_dates)



# Set a permanent directory for output files
permanent_dir <- "E:/URBAN TRENDS/Climate data/Copernicus climate data2"


# Definir los periodos usando días julianos
periods <- list(
  late_winter = 1:59,     # Enero 1 - Febrero 28/29 (días julianos 1 - 59)
  spring = 60:151,        # Marzo 1 - Mayo 31 (días julianos 60 - 151)
  summer = 152:243,       # Junio 1 - Agosto 31 (días julianos 152 - 243)
  autumn = 244:304        # Septiembre 1 - Octubre 31 (días julianos 244 - 304)
)


# Begin looping over each country
for(country in unique(m_coord_c$country_code)){
  
  # Create and set a new temporary directory for this iteration
  temp_dir <- tempdir()
  setwd(temp_dir)
  
  country_coord <- subset(m_coord_c, m_coord_c$country_code == country)
  
  # Output a message indicating the start of processing for the current country
  cat("Processing country:", unique(country_coord$Country.Name), "\n")
  
  # Convert to sf object and set CRS to 3035
  sf_point <- st_as_sf(country_coord[!is.na(country_coord$transect_lon), ], coords = c("transect_lon", "transect_lat"), crs = 3035)
  
  # Transform to WGS84 for coordinates
  m_coord_sf_transformed <- st_transform(sf_point, crs = 4326)
  
  # Update country_coord with WGS84 coordinates
  country_coord$lon_wgs84 <- st_coordinates(m_coord_sf_transformed)[,1]
  country_coord$lat_wgs84 <- st_coordinates(m_coord_sf_transformed)[,2]
  
  # Create grid and bounding box as per your code
  grid_points <- st_make_grid(sf_point, cellsize = 100000, square = TRUE, what = "polygons")
  grid_points <- grid_points[sf_point]
  
  # Transform grid_points to WGS84
  grid_points_wgs84 <- st_transform(grid_points, crs = 4326)
  
  # Get the bounding box of the grid points
  bms_bb <- st_bbox(grid_points_wgs84)
  
  # Create a single polygon from the bounding box
  bms_bbsf <- st_as_sfc(bms_bb)
  
  # Ensure bms_bbsf is an sf object
  bms_bbsf <- st_as_sf(bms_bbsf)
  
  plot(sf_point$geometry)
  plot(bms_bbsf, add=T)
  plot(grid_points[sf_point], add=T, col = '#ff000088')
  
  # Subset phenological data for the current country
  country_pheno <- pheno_mean_dates %>%
    filter(SITE_ID %in% country_coord$transect_id)
  
  # Check if there are any sites to process
  if(nrow(country_pheno) == 0){
    cat("No phenological data to process for country:", country, "\n")
    next
  }
  
  # Merge phenological data with coordinates
  country_pheno <- country_pheno %>%
    left_join(
      country_coord %>% dplyr::select(transect_id, lon_wgs84, lat_wgs84),
      by = c("SITE_ID" = "transect_id")
    )
  
  # Identify and exclude SITE_IDs without coordinates
  missing_coords <- country_pheno %>%
    filter(is.na(lon_wgs84) | is.na(lat_wgs84)) %>%
    pull(SITE_ID) %>%
    unique()
  
  if(length(missing_coords) > 0){
    message("The following SITE_IDs do not have coordinates and will be excluded: ", paste(missing_coords, collapse = ", "))
    country_pheno <- country_pheno %>%
      filter(!(SITE_ID %in% missing_coords))
  }
  
  # Loop over each period
  for(i in seq_along(Period)){
    
    Period_chunk <- Period[i]
    First_year_chunk <- First_year[i]
    Last_year_chunk <- Last_year[i]
    
    cat("           Processing period:", Period_chunk, "\n")
    
    min_output_filename <- paste0("raster_min_temp_", country, "_", Period_chunk, ".tiff")
    climate_data_min <- extract_nc_value(
      first_year = First_year_chunk,
      last_year = Last_year_chunk,
      local_file = FALSE,
      file_path = NULL,
      sml_chunk = Period_chunk,
      spatial_extent = bms_bbsf,  # Use bms_bbsf as spatial extent
      clim_variable = "min temp",
      statistic = "mean",
      grid_size = 0.25,
      ecad_v = NULL,
      write_raster = TRUE,
      out = min_output_filename,
      return_data = TRUE
    )
    
    max_output_filename <- paste0("raster_max_temp_", country, "_", Period_chunk, ".tiff")
    climate_data_max <- extract_nc_value(
      first_year = First_year_chunk,
      last_year = Last_year_chunk,
      local_file = FALSE,
      file_path = NULL,
      sml_chunk = Period_chunk,
      spatial_extent = bms_bbsf,  # Use bms_bbsf as spatial extent
      clim_variable = "max temp",
      statistic = "mean",
      grid_size = 0.25,
      ecad_v = NULL,
      write_raster = TRUE,
      out = max_output_filename,
      return_data = TRUE
    )
    
    # Load the temperature rasters
    rbk_min <- terra::rast(min_output_filename)
    rbk_max <- terra::rast(max_output_filename)
    
    # Compute GDD using the be method
    be_gdd <- gdd_extract(base_temp = 5, min_temp = rbk_min, max_temp = rbk_max, gdd_method = 'be')
    
    # Aplicar el bucle for para recorrer cada transect_id en country_coord
    gdd_results_list <- list()
    
    for (i in 1:nrow(country_coord)) {
      # Extraer los GDD para el transecto actual
      gdd_results_list[[i]] <- extract_gdd_periods(
        transect_id = country_coord$transect_id[i],   # ID del transecto
        lon = country_coord$lon_wgs84[i],             # Longitud
        lat = country_coord$lat_wgs84[i],             # Latitud
        be_gdd = be_gdd,                              # Raster de GDD
        periods = periods                             # Periodos de GDD (late winter, spring, summer, autumn)
      )
      
      # Imprimir el progreso cada 10 transectos
      if (i %% 10 == 0) {
        cat("Procesados", i, "de", nrow(country_coord), "transectos\n")
      }
    }
    
    gdd_results <- do.call(rbind, gdd_results_list)
    
    gdd_results$Country <- country
    
    # Create the filename for saving results
    filename <- paste0(permanent_dir, "/", "gdd_results_", country, "_", Period_chunk, ".csv")
    
    # Save the results to a CSV file
    write.csv(gdd_results, filename, row.names = FALSE)
    
  }
}

# Set working directory to the permanent directory
setwd(permanent_dir)

# Get a list of CSV files starting with "gdd_results_"
file_list <- list.files(pattern = "^gdd_results_.*\\.csv$", full.names = TRUE)

# Read each CSV file and bind them row-wise
combined_data <- do.call(rbind, lapply(file_list, function(file) {
  read.csv(file, header = TRUE)
}))

# View the combined data
head(combined_data)

# Save the combined data as a CSV file
write.csv(combined_data, "gdd_seasonal_periods.csv", row.names = FALSE)

# Optional: Analyze NA rows for each 'SPECIES' and 'SITE_ID'
na_rows_by_site <- combined_data %>%
  group_by(transect_id) %>%
  summarize(total_rows = n(),
            num_na_rows = sum(is.na(GDD_latewinter ) | is.na(GDD_spring ) | is.na(GDD_summer ) |  is.na(GDD_autumn ) ))

# Print the result
print(na_rows_by_site)


unique(combined_data$Country)







