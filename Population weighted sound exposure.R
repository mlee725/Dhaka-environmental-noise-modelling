##############################
#
## Population weighted sound exposure
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
  pop_raster_masked <- mask(x = pop_raster, mask = mask_vector) # to reduce the size of the population raster and cut computing time
  sound_resampled <- resample(sound_raster, pop_raster_masked, method = "med")
  sound_resampled <- mask(x = sound_resampled, mask = mask_vector)
  extracted <- terra::zonal(sound_resampled, pop_raster_masked, fun="median", na.rm=T, touches=FALSE, as.raster=T)
  test <- crds(extracted, df=T)
  extracted2 <- terra::zonal(sound_resampled, pop_raster_masked, fun="median", na.rm=T, touches=FALSE, as.raster=F)
  x <- cbind(test, extracted2)
  return(x)
}

## creates cumulative density plots - I pulled some of this code from this source: https://rpubs.com/clement_hironimus/probdist
cumden_plot <- function(exposure, pop, type, pop2){ # type either "Total" or "Stratified"
  # Expand to individual exposure values
  expanded_data <- rep(exposure, pop)
  expanded_df <- data.frame(exposure = expanded_data)
  sample.range <- floor(min(expanded_df$exposure)):ceiling(max(expanded_df$exposure))
  iq.mean <- mean(expanded_df$exposure)
  iq.sd <- sd(expanded_df$exposure)
  cdf <- pnorm(sample.range, iq.mean, iq.sd)
  iq.dist <- dnorm(sample.range, mean = iq.mean, sd = iq.sd)
  iq.df <- data.frame("sound" = sample.range, "Density" = iq.dist)
  iq.df <- cbind(iq.df, "CDF_LowerTail" = cdf)
  if(type == "Total"){
    ## just one curve for total population
    plot <- ggplot() +
      stat_ecdf(data=expanded_df, aes(exposure), alpha=0.8, lty="11", color = "#E7B800") +
      geom_line(data=iq.df, aes(sound, CDF_LowerTail), lty="solid", color = "#E7B800", linewidth=0.8) +
      labs(
        title = "",
        x = expression(paste("Average ", LA[eq], " (dBA) exposure")),
        y = "Cumulative propertion of population exposed"
      ) +
      theme_classic() 
  }else if(type == "Stratified"){
    ## add in additional population
    expanded_data2 <- rep(exposure, pop2)
    expanded_df2 <- data.frame(exposure = expanded_data2)
    sample.range2 <- floor(min(expanded_df2$exposure)):ceiling(max(expanded_df2$exposure))
    iq.mean2 <- mean(expanded_df2$exposure)
    iq.sd2 <- sd(expanded_df2$exposure)
    cdf2 <- pnorm(sample.range2, iq.mean2, iq.sd2)
    iq.dist2 <- dnorm(sample.range2, mean = iq.mean2, sd = iq.sd2)
    iq.df2 <- data.frame("sound" = sample.range2, "Density" = iq.dist2)
    iq.df2 <- cbind(iq.df2, "CDF_LowerTail" = cdf2)
    
    plot <- ggplot() +
      stat_ecdf(data=expanded_df, aes(exposure), alpha=0.8, lty="11", color = "#E7B800") +
      geom_line(data=iq.df, aes(sound, CDF_LowerTail), lty="solid", color = "#E7B800", linewidth=0.8) +
      stat_ecdf(data=expanded_df2, aes(exposure), alpha=0.8, lty="11", color = "#00AFBB") +
      geom_line(data=iq.df2, aes(sound, CDF_LowerTail), lty="solid", color = "#00AFBB",  linewidth=0.8) +
      labs(
        title = "",
        x = expression(paste("Average ", LA[eq], " (dBA) exposure")),
        y = "Cumulative % of population exposed"
      ) +
      theme_classic() 
  }
  return(plot)
}

##################
## population analysis =====
##################

extractedpopsound_all <- extractedpopsound(WorldPop2020, adim2022, noise_all_masked)
extractedpopsound_day <- extractedpopsound(WorldPop2020, adim2022, noise_day_masked)
extractedpopsound_night <- extractedpopsound(WorldPop2020, adim2022, noise_night_masked)

## cumulative density plots at grid-level
total_plot <- cumden_plot(extractedpopsound$pred_overall_noise, extractedpopsound$bgd_ppp_2020_UNadj, "Total")
total_plot
