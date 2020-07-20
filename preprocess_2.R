library(sf)
library(tidyverse)
library(units)
library(raster)
library(fasterize)
library(gdalUtils)
library(foreach)
library(optmatch)
library(mapview)
    
options("optmatch_max_problem_size"=Inf)

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

fc_change <- raster('fc_change.tif')
names(fc_change) <- c('fc_change')
extent(fc_change) <- extent(covariates_1)

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

fc00_09 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc00_09_ha[-.0-9]*tif')
NAvalue(fc00_09) <- -32768
names(fc00_09) <- c('fc_2000',
                    'fc_2001',
                    'fc_2002',
                    'fc_2003',
                    'fc_2004',
                    'fc_2005',
                    'fc_2006',
                    'fc_2007',
                    'fc_2008',
                    'fc_2009')
extent(fc00_09) <- extent(covariates_1)
fc10_18 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc10_18_ha[-.0-9]*tif')
NAvalue(fc10_18) <- -32768
names(fc10_18) <- c('fc_2010',
                    'fc_2011',
                    'fc_2012',
                    'fc_2013',
                    'fc_2014',
                    'fc_2015',
                    'fc_2016',
                    'fc_2017',
                    'fc_2018')
extent(fc10_18) <- extent(covariates_1)

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

###############################################################################
### Final data setup
# Match on each treatment pixel

f <- treatment ~ fc_change + fc_2015 + lc_2015_agriculture + precip + temp + 
    elev + slope + dist_cities + dist_roads + crop_suitability + pop_2015 + 
    pop_growth + total_biomass

save(f, regions, file='inputs.Rdata')

d <- stack(fc00_09, fc10_18, lc_2015, covariates_1, covariates_2, biomass, 
           population, population_growth, fc_change, regions_rast)
# Ensure only layers in the formula are included (so extra data isn't being 
# passed around)
d <- d[[c(get_names(f),
          'region',
          'ecoregion',
          'pa',
          'fc_2000',
          'fc_2001',
          'fc_2002',
          'fc_2003',
          'fc_2004',
          'fc_2005',
          'fc_2006',
          'fc_2007',
          'fc_2008',
          'fc_2009',
          'fc_2010',
          'fc_2011',
          'fc_2012',
          'fc_2013',
          'fc_2014',
          # 'fc_2015', omitted as already in stack because it is a predictor in 
          # the formula
          'fc_2016',
          'fc_2017',
          'fc_2018')]]

write.table(names(d), file='all_covariates_names.txt', row.names=FALSE, 
            col.names=FALSE)
writeRaster(d, filename='all_covariates.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")
