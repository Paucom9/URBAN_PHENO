rm(list = ls())  # Remove all objects from the current R session to ensure a clean working environment

#Libraries required
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(car)   # For vif() function
library(MuMIn)
library(performance)


# --- Charge data & select and transform variables --- #
#----
# Pheno data
setwd("E:/URBAN TRENDS/Urban_pheno")  
pheno_df <- read.csv("pheno_estimates.csv", sep = ",", dec = ".")
pheno_df$YEAR <- as.factor(pheno_df$YEAR)
pheno_df$SITE_ID <- as.factor(pheno_df$SITE_ID)

# Climate data
setwd("E:/URBAN TRENDS/Climate data/Copernicus climate data")  
gdd_df <- read.csv("gdd5_025_compiled.csv", sep = ",", dec = ".")
gdd_df <- gdd_df %>% select(bms_id, SITE_ID, Year, GDD5_3, GDD5_6, GDD5_9)
gdd_df$Year <- as.factor(gdd_df$Year)
gdd_df$SITE_ID <- as.factor(gdd_df$SITE_ID)
gdd_df$bms_id <- as.factor(gdd_df$bms_id)
gdd_df <- gdd_df %>% rename(YEAR = Year)
setwd("E:/URBAN TRENDS/Climate data/Copernicus climate data2")  
gdd_df2 <- read.csv("gdd_seasonal_periods.csv", sep = ",", dec = ".")
gdd_df2 <- gdd_df2 %>%  rename(SITE_ID = transect_id )
gdd_df2 <- gdd_df2 %>%  rename(YEAR = year )
gdd_df2$YEAR <- as.factor(gdd_df2$YEAR)
gdd_df2$SITE_ID <- as.factor(gdd_df2$SITE_ID)

# Photoperiod data
setwd("E:/URBAN TRENDS/Urban_pheno")
photo_df <- read.csv("photo_estimates.csv", sep = ",", dec = ".")
photo_df$SITE_ID <- as.factor(photo_df$SITE_ID)
photo_df$bms_id <- as.factor(photo_df$bms_id)

# Urbanisation data
setwd("E:/URBAN TRENDS/Urbanisation data")
urban_df <- read.csv("urban_year_data.csv", sep = ",", dec = ".")
urban_df$year <- as.factor(urban_df$year)
urban_df$SITE_ID <- as.factor(urban_df$SITE_ID)
urban_df <- urban_df %>% rename(YEAR = year)
urban_df2 <- read.csv("population_year_data.csv", sep = ",", dec = ".")
urban_df2$YEAR <- as.factor(urban_df2$YEAR)
urban_df2$SITE_ID <- as.factor(urban_df2$SITE_ID)

#----

# --- Merge data --- #

head(pheno_df)
head(gdd_df)
head(photo_df)
head(urban_df)

# Merge pheno_df, gdd_df, photo_df, and urban_df

final_df<- NULL

final_df <- merge(pheno_df, gdd_df, by = c("SITE_ID", "YEAR"), all.x = TRUE) # Merge pheno_df and gdd_df
final_df <- merge(final_df, gdd_df2, by = c("SITE_ID", "YEAR"), all.x = TRUE) # Merge pheno_df and gdd_df
final_df <- merge(final_df, photo_df, by = c("SITE_ID", "bms_id"), all.x = TRUE) # Merge with photo_df
final_df <- merge(final_df, urban_df, by = c("SITE_ID", "YEAR"), all.x = TRUE) # Merge with urban_df
final_df <- merge(final_df, urban_df2, by = c("SITE_ID", "YEAR"), all.x = TRUE) # Merge with urban_df


# View the first few rows of the merged dataframe
head(final_df)
str(final_df)

# Count rows with NA values:
# Assuming your dataframe is named df
num_na_rows <- sum(rowSums(is.na(final_df)) > 0)
print(num_na_rows)

# Remove rows with NA values
final_df <- na.omit(final_df)
str(final_df)
final_df <- final_df[is.finite(final_df$ONSET_mean), ]
final_df <- final_df[is.finite(final_df$OFFSET_mean), ]
final_df <- final_df[is.finite(final_df$PEAKDAY), ]
final_df <- final_df[is.finite(final_df$FLIGHT_LENGTH_mean), ]

str(final_df)
head(final_df)


# Load necessary package
library(dplyr)

# Define a function to remove extreme outliers based on 3 times IQR
remove_extreme_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)  # 1st Quartile
  Q3 <- quantile(x, 0.75, na.rm = TRUE)  # 3rd Quartile
  IQR <- Q3 - Q1  # Interquartile range
  lower_bound <- Q1 - 3 * IQR  # Extreme lower bound
  upper_bound <- Q3 + 3 * IQR  # Extreme upper bound
  x > lower_bound & x < upper_bound  # Return TRUE for non-extreme outliers
}

# Apply the function to each SPECIES to remove only extreme outliers in ONSET_mean and OFFSET_mean
final_df <- final_df %>%
  group_by(SPECIES) %>%
  filter(remove_extreme_outliers(ONSET_mean) & remove_extreme_outliers(OFFSET_mean)) %>%
  ungroup()


# Scale specific variables (e.g., GDD_latewinter, GDD_spring, GDD_summer, GDD_autumn)
final_df <- final_df %>%
  mutate_at(vars(GDD_latewinter, GDD_spring, GDD_summer, GDD_autumn, built_500,built_2000, pop_500, pop_1000, latitude.x), scale)



# Define a function to extract and format model summaries
extract_model_summary <- function(model) {
  if (is.null(model)) {
    return(rep("Model did not converge", 9))  # Return message for all columns if the model didn't converge
  }
  
  summary_data <- summary(model)
  coefficients <- summary_data$coefficients$cond  # Extract conditional model coefficients
  
  # Format the output as "Estimate (p-value)" or "Estimate (<0.001)" if p-value is smaller than 0.001
  formatted_output <- apply(coefficients, 1, function(row) {
    estimate <- round(row["Estimate"], 4)
    p_value <- ifelse(is.na(row["Pr(>|z|)"]), "NA", 
                      ifelse(row["Pr(>|z|)"] < 0.001, "<0.001", round(row["Pr(>|z|)"], 4)))
    return(paste0(estimate, " (", p_value, ")"))
  })
  
  return(formatted_output)
}

# Define a function to run models for a single species and extract results
run_models_for_species <- function(species_name, data) {
  
  # Subset the data for the specific species
  subset_species <- data %>% filter(SPECIES == species_name)
  
  # Define a function to fit a model with error handling
  fit_model <- function(formula, data, family) {
    tryCatch({
      model <- glmmTMB(formula, data = data, family = family)
      return(model)
    }, error = function(e) {
      message("Model did not converge for species: ", species_name)
      return(NULL)
    }, warning = function(w) {
      message("Model did not converge for species: ", species_name)
      return(NULL)
    })
  }
  
  # Fit the models with error handling
  glmm_onset <- fit_model(ONSET_mean ~ pop_500 * GDD_latewinter + pop_500 * GDD_spring + pop_500 * latitude.x + (1|SITE_ID) + (1|YEAR), subset_species, poisson)
  
  glmm_peak <- fit_model(PEAKDAY ~ pop_500 * GDD_spring + pop_500 * GDD_summer + pop_500 * latitude.x + (1|SITE_ID) + (1|YEAR), subset_species, poisson)
  
  glmm_offset <- fit_model(OFFSET_mean ~ pop_500 * GDD_summer + pop_500 * GDD_autumn + pop_500 * latitude.x + (1|SITE_ID) + (1|YEAR), subset_species, poisson)
  
  glmm_flight <- fit_model(FLIGHT_LENGTH_mean ~ pop_500 * GDD_spring + pop_500 * GDD_summer + GDD_autumn + pop_500 * latitude.x + (1|SITE_ID) + (1|YEAR), subset_species, poisson)
  
  # Extract and format the coefficients and p-values for each model
  onset_summary <- extract_model_summary(glmm_onset)
  peak_summary <- extract_model_summary(glmm_peak)
  offset_summary <- extract_model_summary(glmm_offset)
  flight_summary <- extract_model_summary(glmm_flight)
  
  # Collect results into a data frame for each species
  results <- data.frame(
    Model = c("Onset Mean", "Peak Day", "Offset Mean", "Flight Length"),
    pop_500 = c(onset_summary["pop_500"], peak_summary["pop_500"], offset_summary["pop_500"], flight_summary["pop_500"]),
    GDD_latewinter = c(onset_summary["GDD_latewinter"], NA, NA, NA),
    GDD_spring = c(onset_summary["GDD_spring"], peak_summary["GDD_spring"], NA, flight_summary["GDD_spring"]),
    GDD_summer = c(NA, peak_summary["GDD_summer"], offset_summary["GDD_summer"], flight_summary["GDD_summer"]),
    GDD_autumn = c(NA, NA, offset_summary["GDD_autumn"], flight_summary["GDD_autumn"]),
    Latitude = c(onset_summary["latitude.x"], peak_summary["latitude.x"], offset_summary["latitude.x"], flight_summary["latitude.x"]),
    pop_500_GDD_latewinter = c(onset_summary["pop_500:GDD_latewinter"], NA, NA, NA),
    pop_500_GDD_spring = c(onset_summary["pop_500:GDD_spring"], peak_summary["pop_500:GDD_spring"], NA, flight_summary["pop_500:GDD_spring"]),
    pop_500_Latitude = c(onset_summary["pop_500:latitude.x"], peak_summary["pop_500:latitude.x"], offset_summary["pop_500:latitude.x"], flight_summary["pop_500:latitude.x"])
  )
  
  return(results)
}

# Now, run the models for each species in final_df
species_list <- unique(final_df$SPECIES)
total_species <- length(species_list)  # Total number of unique species

# Initialize a list to store the results for all species
all_species_results <- list()

for (i in seq_along(species_list)) {
  species <- species_list[i]
  
  # Debugging: Print which species is being processed along with progress
  print(paste("Processing species", i, "of", total_species, ":", species))
  
  # Run the models and get the coefficients and p-values for each species
  species_results <- run_models_for_species(species, final_df)
  
  # Store the results in the list
  all_species_results[[species]] <- species_results
  
  # Debugging: Print completion of species processing
  print(paste("Completed processing species", i, "of", total_species, ":", species))
}



# Combine all the results into a single data frame
final_results_df <- bind_rows(lapply(names(all_species_results), function(species_name) {
  species_results <- all_species_results[[species_name]]
  # Add a column for the species name
  species_results <- mutate(species_results, SPECIES = species_name)
  return(species_results)
}))

# View the combined data frame
print(head(final_results_df))

# Reorder the columns to put SPECIES first
final_results_df <- final_results_df %>%
  select(SPECIES, everything())  # Moves SPECIES to the first column

# View the head of the updated data frame
print(head(final_results_df))


# Save the final_results_df as a CSV file
write.csv(final_results_df, "E:/URBAN TRENDS/Urban_pheno/Results/pop500_speciesmodel_results.csv", row.names = FALSE)



