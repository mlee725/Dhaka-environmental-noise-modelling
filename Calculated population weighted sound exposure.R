##############################
#
## Calculated population weighted sound exposure
#
##############################

##################
## Load in packages =====
#################
library(report)
library(sf)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(raster)
library(terra)

##################
## Load in data =====
#################
source("~/Source script for exposure analysis data.R")

##################
## Functions =====
#################
## function to match raster grid resolutions of population and sound predictions and then get median exposure with population grid cells
extractedpopsound <- function(pop_raster, mask_vector, sound_raster){
  pop_raster_masked <- mask(x = pop_raster, mask = mask_vector)
  sound_resampled <- resample(sound_raster, pop_raster_masked, method = "med")
  sound_resampled <- mask(x = sound_resampled, mask = mask_vector)
  extracted2 <- terra::zonal(sound_resampled, pop_raster_masked, fun="median", na.rm=T, touches=FALSE, as.raster=F)
  return(extracted2)
}

##################
# -- A -- #
# population level exposure

# extract with rasters
extractedpopsound_all <- extractedpopsound(WorldPop2020, adim2022, noise_all_masked)
extractedpopsound_day <- extractedpopsound(WorldPop2020, adim2022, noise_day_masked)
extractedpopsound_night <- extractedpopsound(WorldPop2020, adim2022, noise_night_masked)
extractedpopsound_den <- extractedpopsound(WorldPop2020, adim2022, noise_den_masked)

# calculate pop exposure
types <- list(df1 = extractedpopsound_all, 
              df2 = extractedpopsound_day, 
              df3 = extractedpopsound_night, 
              df4 = extractedpopsound_den)
extract_df <- data.frame()
for (i in c(1:length(types))){
  pull <- rep(types[[i]][,2], types[[i]][,1]) %>%
    as.data.frame() %>%
    summary() %>%
    as.data.frame() %>%
    dplyr::select(Var2, Freq) %>%
    rename(stat = Var2,
           vals = Freq) %>%
    mutate(type = sub(".*?_(.*?)_.*", "\\1", names(types[[i]])[2]))
  extract_df <- rbind(extract_df, pull[,2:3])
}

##################
# -- B -- #
# cumulative density plots at grid 

expanded_data_night <- rep(extractedpopsound_night$pred_night_noise, extractedpopsound_night$bgd_ppp_2020_UNadj)
expanded_n_df <- data.frame(exposure = expanded_data_night)

expanded_data_day <- rep(extractedpopsound_day$pred_day_noise, extractedpopsound_day$bgd_ppp_2020_UNadj)
expanded_d_df <- data.frame(exposure = expanded_data_day)

expanded_data_all <- rep(extractedpopsound_all$pred_overall_noise, extractedpopsound_all$bgd_ppp_2020_UNadj)
expanded_a_df <- data.frame(exposure = expanded_data_all)

expanded_data_den <- rep(extractedpopsound_den$pred_den_noise, extractedpopsound_den$bgd_ppp_2020_UNadj)
expanded_den_df <- data.frame(exposure = expanded_data_den)

plot <- ggplot() +
  stat_ecdf(data=expanded_n_df, aes(exposure), alpha=0.8, lty="solid", color = "#E2D200", linewidth=1.1) +
  stat_ecdf(data=expanded_d_df, aes(exposure), alpha=0.8, lty="solid", color = "#46ACC8", linewidth=1.1) +
  stat_ecdf(data=expanded_den_df, aes(exposure), alpha=0.8, lty="solid", color = "#E58601", linewidth=1.1) +
  stat_ecdf(data=expanded_a_df, aes(exposure), alpha=0.8, lty="solid", color = "#B40F20", linewidth=1.1) +
  geom_vline(xintercept =55, lty="dashed", color = "#E2D200", linewidth=0.7) +
  geom_vline(xintercept =45, lty="dashed", color = "#46ACC8", linewidth=0.7) +
  geom_vline(xintercept =53, lty="dashed", color = "#E58601", linewidth=0.7) +
  geom_vline(xintercept =70, lty="dashed", color = "#29211F", linewidth=0.7) +
  labs(
    title = "",
    x = "Sound in dBA",
    y = "Cumulative proportion of population exposed"
  ) +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) +
  theme_classic() #+ xlim(c(50,84))
quartz()
plot
