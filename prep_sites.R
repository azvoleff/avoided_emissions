library(sf)
library(foreach)
library(units)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

sites_2018 <- st_read(paste0(data_folder, "/Impact_Sites/Total_2018.shp"))
sites_2018 <- st_zm(sites_2018, drop=TRUE)
sites_2018$data_year <- 2018

sites_2019 <- st_read(paste0(data_folder, "/Impact_Sites/Total_2019.shp"))
sites_2019 <- st_zm(sites_2019, drop=TRUE)
sites_2019$data_year <- 2019

sites <- rbind(sites_2018, sites_2019)

sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'

sites_cea$area_cea_lt_100ha <- sites_cea$area_cea < as_units(100, 'hectares')
sites <- st_transform(sites_cea, 4326)

table(sites_cea$area_cea_lt_100ha)

save(sites, sites_cea, file='sites.RData')

# Check for overlaps
intersections <- foreach (year in c(2019. 2019)) %do% {
    these_sites <- sites_cea[(!sites$area_cea_lt_100ha) & (sites$data_year == 2018)]
    return st_intersects(these_sites, these_sites)
}
