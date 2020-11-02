
library(sf)
library(foreach)
library(units)
library(tidyverse)
library(lubridate)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

sites_2018 <- st_read(paste0(data_folder, "/Impact_Sites/Final_2018.shp"))
sites_2018 <- st_zm(sites_2018, drop=TRUE)
sites_2018$data_year <- 2018

sites_2019 <- st_read(paste0(data_folder, "/Impact_Sites/Final_2019.shp"))
sites_2019 <- st_zm(sites_2019, drop=TRUE)
sites_2019$data_year <- 2019

sites_2020 <- st_read(paste0(data_folder, "/Impact_Sites/Final_2020.shp"))
sites_2020 <- st_zm(sites_2020, drop=TRUE)
sites_2020$data_year <- 2020

sites <- rbind(sites_2018, sites_2019)
sites <- rbind(sites, sites_2020)

sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'

sites <- st_transform(sites_cea, 4326)

table(sites_cea$area_cea < as_units(100, 'hectares'))

sites %>%
    select(CI_ID,
           Data_Year=data_year,
           Area=Area_Name_,
           #CI_Start_Date=CI_Start_D,
           #ECI_End_Date=CI_End_Dat,
           Area_ha=area_cea) -> sites
sites$CI_ID <- factor(sites$CI_ID)

sites$CI_Start_Date_clean <- as.character(sites$CI_Start_Date)
sites$CI_Start_Date_clean <- str_replace(sites$CI_Start_Date_clean, '^([0-9]{4})$', '1/1/\\1')
sites$CI_Start_Date_clean <- mdy(sites$CI_Start_Date_clean)
# Set all start dates that are missing to 2016 (the median year)
sites$CI_Start_Date_clean[is.na(sites$CI_Start_Date_clean)] <- mdy('1/1/2016')
sites$CI_Start_Year <- year(sites$CI_Start_Date_clean)

sites$CI_End_Date_clean <- as.character(sites$CI_End_Date)
sites$CI_End_Date_clean <- str_replace(sites$CI_End_Date_clean, '^([0-9]{4})$', '1/1/\\1')
sites$CI_End_Date_clean <- mdy(sites$CI_End_Date_clean)

# Set all end dates that are greater than 12/31/2019 to NA, so they are treated 
# as ongoing
sites$CI_End_Date_clean[sites$CI_End_Date_clean > mdy('12/31/2019')] <- NA
sites$CI_End_Year <- year(sites$CI_End_Date_clean)

saveRDS(sites, 'sites.RDS')

# Check for overlaps
# intersections <- foreach (year in c(2019, 2019)) %do% {
#     these_sites <- sites_cea[(!sites$area_cea_lt_100ha) & (sites$data_year == 2018), ]
#     return st_intersects(these_sites, these_sites)
# }
