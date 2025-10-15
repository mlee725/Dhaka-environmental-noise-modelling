##############################
#
## Source script for sound prediction, SES, and population data
#
##############################

# This script preps (cleans and transforms) the sound prediction and population and SES data so it can be used in analysis

##################
## Load in packages =====
#################
library(sf)
library(raster)
library(terra)
library(dplyr)
library(tidyterra)
library(sp)

##################
## Load in data =====
#################

# 50m sound predictions
layers <- read.csv("~/noise_lur_layers.csv")

# Thana: percent living in poverty (HCRUpper(%))
HCR2022 <- st_read("~/WB_Upazila_Wise_Poverty_Estimate_2022_Dhaka.shp")

# Ward shapefile 2022
adim2022 <- st_read("~/Dhaka_city_corporation_wards.shp")

# Population
WorldPop2020 <- raster("~/bgd_ppp_2020_UNadj.tif")
WorldPop2020 <- projectRaster(WorldPop2020, crs = "EPSG:32646")
WorldPop2020<- terra::rast(WorldPop2020)

# Suburban
suburban <- st_read("~/sub-urban_from_1km_hex_grid.shp")

#################
## Functions =====
#################

## convert sf data to utm
wgS84_to_utmz46N <- function(x){
  out <- st_make_valid(x)
  out <- st_transform(out, crs = "EPSG:32646")
  return(out)
}

## convert point data of sound predictions to rasters
points_to_raster <- function(x){
  out <- raster::rasterFromXYZ(layers[,c("x_utm", "y_utm", x)], crs = "EPSG:32646", digits=4) # UTM zone 46N
  out <- terra::rast(out)
  outline <- st_make_valid(adim2022)
  outline <- st_transform(outline, crs = "EPSG:32646") 
  out_masked <- mask(x = out, mask = outline)
  return(out_masked)
}

##################
## Modify data =====
##################

# -- sound prediction raster surfaces -- #
noise_all_masked <- points_to_raster("pred_overall_noise")
noise_day_masked <- points_to_raster("pred_day_noise")
noise_night_masked <- points_to_raster("pred_night_noise")

# -- change crs -- #
HCR2022 <- wgS84_to_utmz46N(HCR2022)
adim2022 <- wgS84_to_utmz46N(adim2022)
suburban <- wgS84_to_utmz46N(suburban)

# -- modify SES data -- #
## extract data from sound predictions to thanas and define SES groups
HCR2022 <- HCR2022 %>%
  mutate(noise_All_mean = exactextractr::exact_extract(noise_all_masked, HCR2022, 'mean')) %>%
  mutate(noise_All_med = exactextractr::exact_extract(noise_all_masked, HCR2022, 'median')) %>%
  mutate(noise_day_mean = exactextractr::exact_extract(noise_day_masked, HCR2022, 'mean')) %>%
  mutate(noise_day_med = exactextractr::exact_extract(noise_day_masked, HCR2022, 'median')) %>%
  mutate(noise_night_mean = exactextractr::exact_extract(noise_night_masked, HCR2022, 'mean')) %>%
  mutate(noise_night_med = exactextractr::exact_extract(noise_night_masked, HCR2022, 'median')) %>%
  # suburban thanas
  mutate(drop = ifelse(upazila_th %in% c("Sabujbag", "Mugda", "Shah Ali", "Bhatara", "Khilgaon", "Demra", "Turag", "Uttarkhan", "Badda",
                                         "Khilkhet"), 1, 0)) %>% 
  mutate(drop_anysub = ifelse(upazila_th %in% c("Pallabi", "Dakkhinkhan", "Jatrabari", "Rupnagar","Sabujbag", "Mugda", "Shah Ali", "Bhatara", "Khilgaon", "Demra", "Turag", "Uttarkhan", "Badda",
                                                "Khilkhet"), 1, 0)) %>% # thanas with any suburban areas
  mutate(drop_50sub = ifelse(upazila_th %in% c("Turag", "Uttarkhan", "Badda", "Khilkhet"), 1, 0)) %>% # thanas with greater than 50% of the area suburban
  mutate(HCR_UP_22 = X2022.HCR.U) %>%
  mutate(catq4 = case_when(HCR_UP_22 >= 0 & HCR_UP_22 < 3.225~ "low",
                           HCR_UP_22 >= 3.225 & HCR_UP_22 < 6.3 ~ "med-low",
                           HCR_UP_22 >= 6.3 & HCR_UP_22 < 9.275 ~ "med-high",
                           HCR_UP_22 >= 9.275 & HCR_UP_22 <= 19.2 ~ "high")) %>%
  mutate(catlmh_d = case_when(HCR_UP_22 >= 0 & HCR_UP_22 < 3.225~ "low",
                              HCR_UP_22 >= 3.225 & HCR_UP_22 < 9.275 ~ "med",
                              HCR_UP_22 >= 9.275 & HCR_UP_22 <= 19.2 ~ "high")) %>%
  mutate(catq4 = ordered(catq4, levels = c("low", "med-low", "med-high", "high"))) %>%
  mutate(catlmh_d = ordered(catlmh_d, levels = c("low", "med", "high"))) %>%
  mutate(popden_census = Population/total_area)
