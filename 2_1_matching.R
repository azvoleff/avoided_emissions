library(raster)
library(fasterize)
library(sf)
library(tidyverse)
library(units)
library(foreach)
library(mapview)
library(exactextractr)
library(raster)
library(tictoc)

options("optmatch_max_problem_size"=Inf)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

gdal_crop <- function(r, s) {
    if (filename(r) == '') {
        stop('input raster must have a valid filename')
    }
    e <- extent(s)
    band_names <- names(r)
    out_file <- tempfile(fileext='.tif')
    sf::gdal_utils('warp', filename(r), out_file,
                   options=c('-te', e@xmin, e@ymin, e@xmax, e@ymax, '-multi', 
                             '-co', 'COMPRESS=LZW'))
    out_r <- brick(out_file)
    names(out_r) <- names(r)
    return(out_r)
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}

###############################################################################
### Final data setup

f <- treatment ~ lc_2015_agriculture + precip + temp + elev + slope + 
    dist_cities + dist_roads + crop_suitability + pop_2015 + pop_growth + 
    total_biomass
saveRDS(f, 'Output/formula.RDS')

covariates <- brick('covariates_covariates.tif')
names(covariates) <- read_csv('covariates_covariates.csv')$names
lc_2015 <- brick('covariates_lc_2015.tif')
names(lc_2015) <- read_csv('covariates_lc_2015.csv')$names
fc <- brick('covariates_fc.tif')
names(fc) <- read_csv('covariates_fc.csv')$names
fc_change <- brick('covariates_fc_change.tif')
names(fc_change) <- read_csv('covariates_fc_change.csv')$names

d <- stack(covariates, lc_2015, fc, fc_change)
# Ensure only layers in the formula are included (so extra data isn't being 
# passed around)
d <- d[[c(get_names(f),
          'region',
          'ecoregion',
          'pa',
           paste0('fc_0', seq(0, 9)),
           paste0('fc_', seq(10, 19)),
           paste0('fcc_0', seq(0, 9)),
           paste0('fcc_', seq(10, 19)))]]
write_csv(data.frame(names=names(d)), file='all_covariates_names.csv')

###############################################################################
###  Load sites and covariates

sites <- readRDS('sites.RDS')
dim(sites)

# Drop sites with no overlap with GADM (since they'd throw errors later during 
# the extraction) - these are marine sites
sites <- filter(sites, !(CI_ID %in% c('242002', '242114')))

# Filter to only sites over 100 ha
sites <- sites[!sites$Area_ha < as_units(100, 'hectares'), ]
dim(sites)

# Select only sites that are not rangeland restoration
sites$Rangeland <- FALSE
sites$Rangeland[sites$Restoration == 'Rangeland Restoration'] <- TRUE
table(sites$Rangeland)
sites <- sites[!sites$Rangeland, ]
dim(sites)

regions <- readRDS('regions.RDS')
regions_rast <- fasterize(regions, raster(d[[1]]), field='level1_ID')
names(regions_rast) <- 'region'

d <- stack(d, regions_rast)

###############################################################################
###  Load sites and covariates

# Run extractions of treatment points individually to make catching any polygon 
# errors easier
treatment_key <- foreach(n=1:nrow(sites), .combine=rbind) %do% {
    print(n)
    exact_extract(d$region, sites[n, ], include_cell=TRUE, 
                  include_cols=c('CI_ID', 'Data_Year'))[[1]]
}
treatment_key <- treatment_key[!is.na(treatment_key$region), ]
saveRDS(treatment_key, 'Output/treatment_cell_key.RDS')

# Run extraction of control and treatment data by region to make the problem 
# tractable in-memory
out <- foreach(this_region_ID=unique(treatment_key$region), 
        .packages=c('exactextractr', 'sf')) %do% {
    if (file.exists(paste0('Output/treatments_and_controls_', this_region_ID, 
                           '.RDS'))) {
        print(paste0('Skipping ', this_region_ID, '. Already processed.'))
    } else {
        this_region <- regions[regions$level1_ID == this_region_ID, ]
        vals <- exact_extract(d, this_region, include_cell=TRUE)[[1]]
        saveRDS(vals, paste0('Output/treatments_and_controls_', 
                             this_region_ID, '.RDS'))
    }
}
