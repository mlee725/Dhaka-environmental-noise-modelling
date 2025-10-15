##############################
#
## Source script for Dhaka PM2.5 and noise data
#
##############################

# This script preps (cleans and transforms) the PM2.5 and noise data so it can be used in analysis

##################
## Load in packages =====
#################
library(dplyr)
library(lubridate)
library(sf)

##################
## Load in data =====
#################
Ldata <- read.csv("~/LUR data Oct 2024.csv")
Zefan <- read.csv("~/Zefan_Jan27_filtadjupdate.csv")
noise <- read.csv("~/Updated_noise_Nov15 2025.csv")
usem <- read.csv("~/Dhaka_PM2.5_2023_YTD.csv")
met <- read.csv("~/NOAAdhakametexport.csv")
deployment <- read.csv("~/SamplerDeployments.csv")
Ssites <- st_read("~/All_Site_Geometry.shp")
TAF <- read.csv("~/TAF_jan27 temporal adi update.csv")
era5 <- read.csv("~/martha_era5_format.csv")
cat <- read.csv("~/Site category.csv")
clusters_12_all <- read.csv("~/clusters_from_Ben.csv")

##################
## Modify data that is needed for PM2.5 and noise =====
#################

## sampling sites
Ssites <- st_transform(Ssites, crs = 4326) # change crs

## meteorology clusters
clusters_12_all <- clusters_12_all %>% # using annual data
  mutate(Time_UTC_clusters = ymd_hm(Time_UTC_clusters, tz = "UTC")) %>%
  mutate(Time_clusters = with_tz(Time_UTC_clusters, tz = "Asia/Dhaka")) # change timezone

## sensor deployment and collection sheet
deployment <- deployment %>%
  mutate(DDateTime = dmy_hm(paste(DeploymentDate, DeploymentTime), tz = "Asia/Dhaka")) %>%
  mutate(DDate = dmy(DeploymentDate, tz = "Asia/Dhaka")) %>%
  mutate(CDateTime = dmy_hm(paste(CollectionDate, CollectionTime), tz = "Asia/Dhaka")) %>%
  mutate(runningtime_check = time_length(interval(DDateTime, CDateTime), unit = "hour")) %>%
  filter(!is.na(CDateTime)) %>% # 13-Feb-2023 deployment at DU was stolen, remove
  mutate(Key = paste0(NoiseID, "_", DDate, "_", Site, ".xlsx")) 

## us embassy PM2.5 monitor
usem <- usem %>%
  mutate(DateTime = ymd_hm(Date, tz = "Asia/Dhaka")) %>%
  mutate(week= week(DateTime)) %>%
  mutate(day = day(DateTime)) %>%
  mutate(yday = yday(DateTime))
usem_week <- aggregate(data = usem, NowCastConc ~ week, mean) # weekly mean
usem_day <- aggregate(data = usem, NowCastConc ~ yday, mean) # daily mean
colnames(usem_day) <- c("yday", "NowCastConc_yday")

## meteorology (airport)
met <- met %>%
  filter(!is.na(latitude)) %>%
  mutate(TP = ymd_hms(date, tz = "Asia/Dhaka")) 

## meteorology (ERA5)
era5 <- era5 %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "Asia/Dhaka")) %>%
  mutate(DateTime_h = DateTime - dseconds(1)) 

##################
## Noise data =====
#################

## fix deployment code
noise[noise$List=="mk3-01_29-Jan-2023_01_29_DU-F3.xlsx", "List"] <- "mk3-01_29-Jan-2023_DU-F3.xlsx"
noise[noise$List=="mk3-05_18-Jan_2023_BLC-F2.xlsx", "List"] <- "mk3-05_18-Jan-2023_BLC-F2.xlsx"

## create noise data for analysis
noise2 <- noise %>%
  mutate(DateTime = ymd_hms(DateTime, tz = "Asia/Dhaka")) %>%
  mutate(TimeofDay = ifelse(hour >= 6 & hour < 21, "Day", "Night")) %>%
  mutate(Date = ymd(Date, tz = "Asia/Dhaka")) %>%
  mutate(mon = lubridate::month(Date, label=TRUE)) %>%
  mutate(min = minute(DateTime)) %>%
  mutate(deplydt = gsub(".*[_]([^.]+)[_].*", "\\1", List)) %>%
  mutate(deplydt = dmy(deplydt)) %>%
  mutate(deplydtsite = gsub(".*?[_]([^_]+_[^\\.]+)\\..*", "\\1", List)) %>%
  mutate(week = lubridate::week(Date)) %>%
  mutate(wday = wday(Date)) %>%
  mutate(weeknd = ifelse(wday == 6| wday == 7, 1, 0)) %>%
  left_join(cat, by = c("Site" = "Site")) %>%
  group_by(List) %>%
  mutate(runningtime = time_length(interval(min(DateTime), max(DateTime)), unit = "hour")) %>%
  mutate(runningtime_cur = time_length(interval(min(DateTime), DateTime), unit = "hour")) %>%
  mutate(start = min(DateTime)) %>%
  mutate(stop = max(DateTime)) %>%
  mutate(RH = mean(met$RH[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>% # 
  mutate(temp = mean(met$air_temp[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(ws = mean(met$ws[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(visibility = mean(met$visibility[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(u10_m.s = mean(era5$u10_m.s[era5$DateTime >= start & era5$DateTime <= stop], na.rm = TRUE)) %>%
  mutate(v10_m.s = mean(era5$v10_m.s[era5$DateTime >= start & era5$DateTime <= stop], na.rm = TRUE)) %>%
  mutate(DeploymentDate = min(Date)) %>%
  mutate(Key = paste0(ID, "_", deplydt, "_", Site, ".xlsx")) %>%
  as.data.frame()

noise_Hourly <- aggregate(data = noise2, LEQ ~ Site + hour, mean)
noise_site <- aggregate(data = noise_Hourly, LEQ ~ Site, mean)

##################
## PM2.5 data =====
#################

## fix deployment code
Zefan[Zefan$ID=="Z2201-017_12-Sep-2023-UTTARA-R5-NAHIN.txt", "ID"] <- "Z2201-017_12-Sep-2023_UTTARA-NAHIN-R5.txt"
Zefan[Zefan$ID=="Z1910-056_12-Sep-2023-KAOLABAZAR-R46.txt", "ID"] <- "Z1910-056_12-Sep-2023_KAOLABAZAR-R46.txt"
Zefan[Zefan$ID=="Z1910-197_12-Sep-2023-CHALABON-R49.txt", "ID"] <- "Z1910-197_12-Sep-2023_CHALABON-R49.txt"
Zefan[Zefan$ID=="Z2201-027_10-Sep-2023-BARC-GROUND-F1.txt", "ID"] <- "Z2201-027_10-Sep-2023_BARC-GROUND-F1.txt"

## create PM2.5 data for analysis
Zefan2 <- Zefan %>%
  mutate(DateTime = ymd_hms(paste(Date, Time), tz = "Asia/Dhaka")) %>%
  mutate(TimeofDay = ifelse(hour >= 6 & hour < 21, "Day", "Night")) %>%
  mutate(deplydtsite = gsub(".*?[_]([^_]+_[^\\.]+)\\..*", "\\1", ID)) %>%
  mutate(deplydt = gsub(".*[_]([^.]+)[_].*", "\\1", ID)) %>%
  mutate(deplydt = dmy(deplydt)) %>%
  mutate(minute = minute(DateTime)) %>%
  mutate(Date = ymd(Date, tz = "Asia/Dhaka")) %>%
  mutate(wday = wday(Date)) %>%
  mutate(yday = yday(Date)) %>%
  mutate(weeknd = ifelse(wday == 6| wday == 7, 1, 0)) %>%
  mutate(DateTime_h = ymd_hms(paste0(Date, " ", hour, ":00:00"), tz = "Asia/Dhaka")) %>%
  left_join(TAF, by = c("week" = "week")) %>%
  mutate(us_w_pm_all = case_when(season == "Dry" ~ 181.2457672,
                                 season == "Wet" ~ 48.27640693)) %>%
  mutate(us_d_pm_all = case_when(season == "Dry" ~ 182.380791,
                                 season == "Wet" ~ 49.32570513)) %>%
  group_by(season) %>%
  mutate(cday = yday(deplydt) - min(yday)+1) %>%
  mutate(cday_24 = yday(Date) - min(yday)+1) %>%
  mutate(us_pm_all = mean(usem$NowCastConc[usem$DateTime >= min(DateTime) & usem$DateTime <= max(DateTime)])) %>% 
  # remove the runs with less that 24 hours and then add back in data for each hour with the filter values
  mutate(pm_filtersub = case_when(ID == "Z2201-029_20-Feb-2023_SHANIR-AKHRA-R42.txt" ~ 265.8721419,
                                  TRUE ~ PM2.5.filter.adjusted2)) %>%
  group_by(ID) %>%
  mutate(runningtime = time_length(interval(min(DateTime), max(DateTime)), unit = "hour")) %>%
  mutate(runningtime_cur = time_length(interval(min(DateTime), DateTime), unit = "hour")) %>%
  mutate(all_adj = pm_filtersub/TAF_uodated) %>%
  mutate(all_adj_old = pm_filtersub/TAF) %>%
  mutate(start = case_when(ID == "Z2201-029_20-Feb-2023_SHANIR-AKHRA-R42.txt" ~ ymd_hms("2023-02-20 13:49:00", tz = "Asia/Dhaka"),
                           TRUE ~ min(DateTime))) %>%
  mutate(startsub24 = start - hours(24)) %>%
  mutate(stop = case_when(ID == "Z2201-029_20-Feb-2023_SHANIR-AKHRA-R42.txt" ~ ymd_hms("2023-02-23 13:49:00", tz = "Asia/Dhaka"),
                          TRUE ~ max(DateTime))) %>%
  mutate(stop2 = max(DateTime)) %>%
  mutate(cluster_12_all = DescTools::Mode(clusters_12_all$icluster[clusters_12_all$Time_clusters >= start & clusters_12_all$Time_clusters <= stop])[1]) %>%
  mutate(key = paste0(Site, "__", season, "__", date(start))) %>%
  mutate(key2 = paste0(Site, "__", season)) %>%
  mutate(us_pm_slice = mean(usem$NowCastConc[usem$DateTime >= start & usem$DateTime <= stop])) %>% 
  mutate(CF = us_pm_all/us_pm_slice) %>%
  mutate(adj_pm_c = pm_filtersub*CF)%>%
  left_join(noise_site, by = c("Site" = "Site"), keep = T) %>% # "JUR-RKC-R3" and "POSTOGOLA-R37" don't have noise measurements
  mutate(RH = mean(met$RH[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>% 
  mutate(temp = mean(met$air_temp[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(ws = mean(met$ws[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(visibility = mean(met$visibility[met$TP >= start & met$TP <= stop], na.rm = TRUE)) %>%
  mutate(u10_m.s = mean(era5$u10_m.s[era5$DateTime >= start & era5$DateTime <= stop], na.rm = TRUE)) %>%
  mutate(v10_m.s = mean(era5$v10_m.s[era5$DateTime >= start & era5$DateTime <= stop], na.rm = TRUE)) %>%
  mutate(usem_pm_r24 = mean(usem$NowCastConc[usem$DateTime >= startsub24 & usem$DateTime <= start], na.rm = TRUE)) %>%
  left_join(usem_week, by = c("week" = "week")) %>%
  mutate(CFw = us_w_pm_all/NowCastConc) %>%
  mutate(adj_w_pm = pm_filtersub*CFw) %>%
  dplyr::rename("Site" = "Site.x") %>%
  mutate(yday = yday(DateTime)) %>%
  left_join(usem_day, by = c("yday" = "yday")) %>%
  mutate(CFd = us_d_pm_all/NowCastConc_yday) %>%
  mutate(adj_d_pm_main = pm_filtersub*CFd) %>%
  left_join(usem[,c("NowCastConc", "DateTime")], by = c("DateTime_h" = "DateTime")) %>%
  left_join(era5[,c("u10_m.s", "v10_m.s", "ws_m.s", "t2m_degC", "rh_.", "DateTime_h")], by = c("DateTime_h" = "DateTime_h")) %>%
  left_join(clusters_12_all[,c("icluster", "Time_clusters")], by = c("DateTime_h" = "Time_clusters")) %>%
  as.data.frame()
