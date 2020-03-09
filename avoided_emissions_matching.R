library(sf)
library(raster)
library(fasterize)
library(gdalUtils)
library(foreach)
library(dplyr)
library(doParallel)
library(optmatch)
    
registerDoParallel(4)

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

# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
foreach_rbind <- function(d1, d2) {
    if (is.null(d1) & is.null(d2)) {
        return(NULL)
    } else if (!is.null(d1) & is.null(d2)) {
        return(d1)
    } else if (is.null(d1) & !is.null(d2)) {
        return(d2)
    } else  {
        return(bind_rows(d1, d2))
    }
}

match_ae <- function(d, f) {
    # Filter out sites without at least one treatment unit or without at
    # least one control unit
    d <- d %>%
        filter(complete.cases(.)) %>%
        group_by(region) %>%
        mutate(n_treatment=sum(treatment),
               n_control=sum(!treatment)) %>%
        filter(n_treatment >= 1, n_control >= 1)

    # Note custom combine to handle iterations that don't return any value
    ret <- foreach (this_region=unique(d$region),
                    .packages=c('optmatch', 'dplyr'),
                    .combine=foreach_rbind, .inorder=FALSE) %do% {
        this_d <- filter(d, region == this_region)
        d_CI <- filter(this_d, treatment)
        # Filter out climates and land covers that don't appear in the CI
        # sample, and drop these levels from the factors
        this_d <- filter(this_d,
            biome %in% unique(d_CI$biome),
            ecoregion %in% unique(d_CI$ecoregion),
            pa %in% unique(d_CI$pa))
        this_d$biome <- droplevels(this_d$biome)
        this_d$ecoregion <- droplevels(this_d$ecoregion)
        this_d$pa <- droplevels(this_d$pa)
        # Can't stratify by land cover or climate if they only have one level
        if (nlevels(this_d$biome) >= 2) {
            f <- update(f, ~ . + strata(biome))
        } else {
            f <- update(f, ~ . - biome)
        }
        if (nlevels(this_d$ecoregion) >= 2) {
            f <- update(f, ~ . + strata(ecoregion))
        } else {
            f <- update(f, ~ . - ecoregion)
        }
        if (nlevels(this_d$pa) >= 2) {
            f <- update(f, ~ . + strata(pa))
        } else {
            f <- update(f, ~ . - pa)
        }
        if (nrow(d_CI) > 2) {
            model <- glm(f, data=this_d)
            dists <- match_on(model, data=this_d)
        } else {
            # Use Mahalanobis distance if there aren't enought points to run a
            # glm
            dists <- match_on(f, data=this_d)
        }
        dists <- caliper(dists, 2)
        # If the controls are too far from the treatments (due to the caliper) 
        # then the matching may fail. Can test for this by seeing if subdim 
        # runs successfully
        subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                 error=function(e) return(FALSE))
        if (subdim_works) {
            m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
            this_d <- this_d[matched(m), ]
        } else {
            this_d <- data.frame()
        }
        # Need to handle the possibility that there were no matches for this 
        # treatment, meaning this_d will be an empty data.frame
        if (nrow(this_d) == 0) {
            return(NULL)
        } else {
            return(this_d)
        }
    }
    return(ret)
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    vars <- strsplit(f, split='[+ ~]')[[1]]
    vars[vars != '']
}

###############################################################################

# Load polygons - note they are reprojected to WGS84
# TODO: both need to be planar for st_intersects, so will need to reproject the 
# rasters to CEA
sites <- st_read("Data/impact_sites.gpkg", layer="impact_sites")
# Drop the z dimension
sites <- st_zm(sites, drop=TRUE)

# Load covariates
forest <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_forest_2000_2015_ha[-.0-9]*tif')
names(forest) <- c('forest_2000',
                   'forest_2005',
                   'forest_2010',
                   'forest_2015')
NAvalue(forest) <- -32768

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
#covariates_2 <- crop(covariates_2, regions)
NAvalue(covariates_2) <- -32768

population <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_pop_2000_2015_cnt[-.0-9]*tif')
names(population) <- c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015')
#population <- crop(population, regions)
NAvalue(population) <- -32768

lc_bands <- c('forest',
              'grassland',
              'agriculture',
              'wetlands',
              'artificial', 
              'other',
              'water')
lc_2001 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_lc2001_ha[-.0-9]*tif')
names(lc_2001) <- lc_bands
lc_2015 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_lc2015_ha[-.0-9]*tif')
names(lc_2015) <- lc_bands

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

# TODO: some places can't be matched on level 2 since they are ALL of level 2, 
# like the Galapagos for example
#for (n in 1:3)) {
#m <- foreach(n=1:nrow(sites),
ae <- foreach(row_num=1:50,
#ae <- foreach(row_num=1:200,
             .packages=c('raster', 'rgeos', 'optmatch', 'dplyr', 'foreach'),
             .combine=foreach_rbind) %do% {
    print(row_num )
    p <- sites[row_num, ]
    # TODO: st_intersects expects planar coords
    inter <- st_intersects(p, regions)[[1]]
    if (length(inter) > 0) {
        r <- crop(regions_rast, regions[inter, ])
    } else {
        # Some areas won't overlap the GADM at all (marine sites for example). 
        # Return no results for these areas.
        return(NULL)
    }
    p_rast <- rasterize(p, r, background=0)
    # Note that the below will drop any portions of a site that don't fall 
    # within a region in GADM- this means that marine areas will be dropped. 
    # Given this calculation is for avoided emissions, likely not a problem... 
    #
    # Treatment is 1, control is zero
    treat_or_control <- mask(p_rast, r)
    names(treat_or_control) <- 'treatment'

    # Crop matching vars to the area of interest (the region(s) the site falls 
    # within) and then extract values
    fc <- crop(forest, r)
    cv_1 <- crop(covariates_1, r)
    cv_2 <- crop(covariates_2, r)
    vals <- getValues(stack(r, treat_or_control, fc, cv_1, cv_2))
    vals <- data.frame(vals)
    vals <- vals[!is.na(vals$treatment), ]
    vals$treatment <- as.logical(vals$treatment)
    vals$biome <- as.factor(vals$biome)
    vals$ecoregion <- as.factor(vals$ecoregion)
    vals$pa <- as.factor(vals$pa)

    # Match on each treatment pixel
    f <- treatment ~ forest_2000 + precip + temp + elev + slope + dist_cities + dist_roads + crop_suitability
    m <- match_ae(vals, f)
    if (is.null(m)) {
        return(NULL)
    } else {
        m <- select(m, c(get_names(f), 'forest_2015'))
        m$CI_ID <- as.character(p$CI_ID)
        return(m %>% select(CI_ID, everything()))
    }
}

write.csv(ae, file='ae_raw.csv', row.names=FALSE)

emissions_details <- ae %>%
    group_by(CI_ID, treatment) %>%
    summarise(forest_initial=sum(forest_2000, na.rm=TRUE),
              forest_final=sum(forest_2015, na.rm=TRUE),
              forest_loss = forest_initial - forest_final,
              n= n())
write.csv(emissions_details, file='ae_details.csv', row.names=FALSE)
emissions_summary <- emissions_details %>% group_by(CI_ID) %>%
        summarise(forest_ha_treat_minus_control=forest_loss[treatment] - forest_loss[!treatment],
                  n_treatment = n[treatment],
                  n_control = n[!treatment])
write.csv(emissions_summary, file='ae_summary.csv', row.names=FALSE)
