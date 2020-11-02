library(sf)
library(tidyverse)
library(raster)
library(fasterize)
library(gdalUtils)
library(foreach)
    
data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

load_as_vrt <- function(folder, pattern, band=FALSE, raster=TRUE) {
    vrt_file <- tempfile(fileext='.vrt')
    files <- list.files(folder, pattern=pattern)
    if (length(files) == 0) {
        stop('No files found')
    }
    if (band) {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file, b=band)
        r <- raster(vrt_file)
    } else {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file)
        r <- raster::stack(vrt_file)
    }
    if (raster) {
        return(r)
    } else {
        return(vrt_file)
    }
}

# Function used to get IDs from a rasterized set of polygons (to determine 
# which polygons were lost due to rasterization (very small polygons drop out)
get_unique <- function(x) {
    bs <- raster::blockSize(x)
    n_blocks <- bs$n
    for (block_num in 1:n_blocks) {
        these_vals <- unique(raster::getValues(x,
                                               row=bs$row[block_num], 
                                               nrows=bs$nrows[block_num]))
        if (block_num == 1) {
            out <- these_vals
        } else {
            out <- unique(c(out, these_vals))
        }
    }
    return (na.omit(out))
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    vars <- strsplit(f, split='[+ ~]')[[1]]
    vars[vars != '']
}

###############################################################################

# Load covariates
covariates_1 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_covariates_01[-.0-9]*tif')
names(covariates_1) <- c('precip',
                         'temp',
                         'elev',
                         'slope', 
                         'biome',
                         'ecoregion')
NAvalue(covariates_1) <- -32768

covariates_2 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_covariates_02[-.0-9]*tif')
names(covariates_2) <- c('dist_cities',
                         'dist_roads',
                         'travel_cities',
                         'pa', 
                         'crop_suitability')
NAvalue(covariates_2) <- -32768

population <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_pop_2000_2015_cnt[-.0-9]*tif')
names(population) <- c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015')
NAvalue(population) <- -32768

biomass <- raster('biomass.tif')
names(biomass) <- c('total_biomass')
extent(biomass) <- extent(covariates_1)

population_growth <- raster('population_growth.tif')
names(population_growth) <- c('pop_growth')

covariates <- stack(covariates_1, covariates_2, biomass, population, 
                    population_growth)
writeRaster(covariates, filename='covariates_covariates.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(covariates)), 'covariates_covariates.csv')

fc_change <- brick('fc_change.tif')
extent(fc_change) <- extent(covariates_1)
names(fc_change) <- c(paste0('fcc_0', seq(0, 9)),
                      paste0('fcc_', seq(10, 19)))
writeRaster(fc_change, filename='covariates_fc_change.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(fc_change)), 'covariates_fc_change.csv')

# lc_2001 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_lc2001_ha[-.0-9]*tif')


# lc_2001 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_lc2001_ha[-.0-9]*tif')
# names(lc_2001) <- c('lc_2001_forest',
#                     'lc_2001_grassland',
#                     'lc_2001_agriculture',
#                     'lc_2001_wetlands',
#                     'lc_2001_artificial', 
#                     'lc_2001_other',
#                     'lc_2001_water')

lc_2015 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_lc2015_ha[-.0-9]*tif')
names(lc_2015) <- c('lc_2015_forest',
                    'lc_2015_grassland',
                    'lc_2015_agriculture',
                    'lc_2015_wetlands',
                    'lc_2015_artificial', 
                    'lc_2015_other',
                    'lc_2015_water')
writeRaster(lc_2015, filename='covariates_lc_2015.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(lc_2015)), 'covariates_lc_2015.csv')

fc00_09 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc00_09_ha[-.0-9]*tif')
NAvalue(fc00_09) <- -32768
names(fc00_09) <- paste0('fc_0', seq(0, 9))
extent(fc00_09) <- extent(covariates_1)

fc10_19 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc10_19_ha[-.0-9]*tif')
NAvalue(fc10_19) <- -32768
names(fc10_19) <- paste0('fc_', seq(10, 19))
extent(fc10_19) <- extent(covariates_1)

fc <- stack(fc00_09, fc10_19)
writeRaster(fc, filename='covariates_fc.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
write_csv(data.frame(names=names(fc)), 'covariates_fc.csv')

###############################################################################
### Load GADM boundaries
regions <- st_read(paste0(data_folder, "/gadm36_levels_gpkg/gadm36_levels.gpkg"), layer="level1")
regions$level0_ID <- as.numeric(factor(regions$GID_0))
regions$level1_ID <- as.numeric(factor(regions$GID_1))
regions_rast <- fasterize(regions, raster(covariates_1[[1]]),
                          field='level1_ID')
region_IDs_after_rasterization <- get_unique(regions_rast)
regions <- regions[regions$level1_ID %in% region_IDs_after_rasterization, ]


regions$level0_ID <- as.numeric(factor(as.character(regions$GID_0)))
regions$level1_ID <- as.numeric(factor(as.character(regions$GID_1)))
# Now re-rasterize boundaries (with ID's that will disappear dropped) to ensure
# that all IDs are sequential and that they match between the data.frame and 
# the raster.
regions_rast <- fasterize(regions, raster(covariates_1[[1]]),
                          field='level1_ID')
names(regions_rast) <- 'region'
region_IDs_after_rasterization <- get_unique(regions_rast)
stopifnot(sort(region_IDs_after_rasterization) == sort(regions$level1_ID))
saveRDS(regions, file='regions.RDS')
