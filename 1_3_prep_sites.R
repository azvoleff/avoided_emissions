library(sf)
library(foreach)
library(units)
library(tidyverse)
library(lubridate)

data_folder <- 'D:/Documents and Settings/azvoleff/Desktop/Guardian Stuff'

alto_comparison <- st_read(paste0(data_folder, "/AltoMayo_ComparisonArea (1).kml"))
alto_comparison <- st_zm(alto_comparison, drop=TRUE)
alto_comparison <- transmute(alto_comparison, name='Alto', type='comparison')
alto_site <- st_read(paste0(data_folder, "/AMPF_boundary_2008.shp"))
alto_site <- st_zm(alto_site, drop=TRUE)
alto_site <- transmute(alto_site, name='Alto', type='site')
alto_site <- st_transform(alto_site, 4326)
sites <- rbind(alto_comparison, alto_site)

chyulu_comparison <- st_read(paste0(data_folder, "/ChyuluHills_ComparisonArea.kml"))
chyulu_comparison <- st_zm(chyulu_comparison, drop=TRUE)
chyulu_comparison <- transmute(chyulu_comparison, name='chyulu', type='comparison')
chyulu_site <- st_read(paste0(data_folder, "/ProjectArea.shp"))
chyulu_site <- st_zm(chyulu_site, drop=TRUE)
chyulu_site <- st_transform(chyulu_site, 4326)
chyulu_site <- transmute(chyulu_site, name='chyulu', type='site')
sites <- rbind(sites, chyulu_comparison)
sites <- rbind(sites, chyulu_site)

sites_cea <- st_transform(sites, '+proj=cea')
sites_cea$area_cea <- st_area(sites_cea)
units(sites_cea$area_cea) <- 'hectares'
sites <- st_transform(sites_cea, 4326)

st_write(sites, paste0(data_folder, '/sites_combined.shp'))
saveRDS(sites, 'sites.RDS')
