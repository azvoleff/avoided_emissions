library(sf)
library(tidyverse)
library(units)
library(raster)
library(foreach)
library(optmatch)
library(mapview)
library(doParallel)
#registerDoParallel(4)


MAX_TREATMENT <- 2000
    
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
        mutate(n_treatment=sum(treatment),
               n_control=sum(!treatment)) %>%
        filter(n_treatment >= 1, n_control >= 1)

    # Note custom combine to handle iterations that don't return any value
    d_CI <- filter(d, treatment)
    # Filter out climates and land covers that don't appear in the CI
    # sample, and drop these levels from the factors
    d <- filter(d,
        region %in% unique(d_CI$region),
        biome %in% unique(d_CI$biome),
        ecoregion %in% unique(d_CI$ecoregion),
        pa %in% unique(d_CI$pa))
    d$region <- droplevels(d$region)
    d$biome <- droplevels(d$biome)
    d$ecoregion <- droplevels(d$ecoregion)
    d$pa <- droplevels(d$pa)
    # Can't stratify by land cover or climate if they only have one level
    if (nlevels(d$region) >= 2) {
        f <- update(f, ~ . + strata(region))
    } else {
        f <- update(f, ~ . - region)
    }
    if (nlevels(d$biome) >= 2) {
        f <- update(f, ~ . + strata(biome))
    } else {
        f <- update(f, ~ . - biome)
    }
    if (nlevels(d$ecoregion) >= 2) {
        f <- update(f, ~ . + strata(ecoregion))
    } else {
        f <- update(f, ~ . - ecoregion)
    }
    if (nlevels(d$pa) >= 2) {
        f <- update(f, ~ . + strata(pa))
    } else {
        f <- update(f, ~ . - pa)
    }
    if (nrow(d_CI) > 2) {
        model <- glm(f, data=d)
        dists <- match_on(model, data=d)
    } else {
        # Use Mahalanobis distance if there aren't enough points to run a
        # glm
        dists <- match_on(f, data=d)
    }
    dists <- caliper(dists, 2)
    # If the controls are too far from the treatments (due to the caliper) 
    # then the matching may fail. Can test for this by seeing if subdim 
    # runs successfully
    subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                             error=function(e) return(FALSE))
    if (subdim_works) {
        m <- fullmatch(dists, min.controls=1, max.controls=1, data=d)
        d <- d[matched(m), ]
    } else {
        d <- data.frame()
    }
    # Need to handle the possibility that there were no matches for this 
    # treatment, meaning d will be an empty data.frame
    if (nrow(d) == 0) {
        return(NULL)
    } else {
        return(d)
    }
}

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    vars <- strsplit(f, split='[+ ~]')[[1]]
    vars[vars != '']
}

###############################################################################
###  Load sites and covariates

d <- brick('all_covariates.tif')
names(d) <- read.table('all_covariates_names.txt')$V1

# Load f and regions sf data (after normalizing IDs to match rasterization)
load('inputs.Rdata')

load('sites.RData')
dim(sites)

# Filter to only sites over 100 ha
sites <- sites[!sites$Area_ha < as_units(100, 'hectares'), ]
dim(sites)

# Select only sites that are not rangeland restoration
sites$Rangeland <- FALSE
sites$Rangeland[sites$Restoration == 'Rangeland Restoration'] <- TRUE
table(sites$Rangeland)
sites <- sites[!sites$Rangeland, ]
dim(sites)

###############################################################################
###  Run matching
set.seed(31)
# message_log <- file("log_message.txt", 'w')
# output_log <- file("log_output.txt", 'w')
ae <- foreach(row_num=1:nrow(sites),
             .packages=c('sf', 'tidyverse', 'optmatch', 'raster'),
             .combine=foreach_rbind, .inorder=FALSE) %do% {
    p <- sites[row_num, ]
    CI_ID <- as.character(p$CI_ID)
    Data_Year <- p$Data_Year
    if (file.exists(paste0('Output/m_', CI_ID, '_', Data_Year, '.RData'))) {
        print(paste('Skipping', CI_ID, Data_Year))
        next
    }
    # sink(output_log, append=TRUE, type='output')
    # sink(message_log, append=TRUE, type="message")
    print(paste0(row_num, ': Cropping'))
    # st_intersects expects planar coords, but approximate match from WGS84 is 
    # fine here since the regions are all much larger than the sites - will 
    # convert to CEA later when it matters
    inter <- tryCatch(st_intersects(p, regions)[[1]],
                      error=function(e) return(FALSE))
    if ((length(inter) == 0) || !inter) {
        # Some areas won't overlap the GADM at all (marine sites for example). 
        # Return no results for these areas.
        return(NULL)
    } else {
        d_crop <- gdal_crop(d, regions[inter, ])
    }

    print(paste0(row_num, ': Rasterizing and masking'))
    p_rast <- rasterize(p, d_crop, background=0)
    # Note that the below will drop any portions of a site that don't fall 
    # within a region in GADM- this means that marine areas will be dropped. 
    # Given this calculation is for avoided emissions, likely not a problem... 
    #
    # Treatment is 1, control is zero
    #
    # TODO: Need to ensure we also mask out other areas where CI has ongoing 
    # interventions.
    treat_or_control <- mask(p_rast, d_crop$region)
    names(treat_or_control) <- 'treatment'
    d_crop <- stack(treat_or_control, d_crop)

    # Project all items to cylindrical equal area
    # d_crop <- projectRaster(d_crop, crs=CRS('+proj=cea'), method='ngb')

    print(paste0(row_num, ': Reading values'))
    # Process the pixels in blocks as some regions are large
    
    treatment_vals <- extract(d_crop, p, df=TRUE)
    treatment_vals <- treatment_vals[names(treatment_vals) != 'ID']
    print(paste0(row_num, ': ', nrow(treatment_vals), ' treatment cells'))

    ###
    # Sample the treatment cells if there are more than MAX_TREATMENT
    if (nrow(treatment_vals) > MAX_TREATMENT) {
        print(paste0(row_num, ': sampling treatment cells'))
        treatment_vals <- sample_n(treatment_vals, MAX_TREATMENT)
    }
    ###
    # Sample the control cells if there are more than 100*MAX_TREATMENT
    #
    # Note the below factor is set slightly higher than the one used with 
    # sampling to ensure that there are enough cells to sample from
    if (ncell(d_crop) < (120 * nrow(treatment_vals))) {
        print(paste0(row_num, ': using all possible control points'))
        control_vals <- data.frame(getValues(d_crop))
        control_vals <- dplyr::filter(control_vals, treatment == 0)
    } else {
        print(paste0(row_num, ': sampling control points'))
        control_vals <- as.data.frame(sampleRegular(d_crop, 100 * nrow(treatment_vals), useGDAL=TRUE))
        control_vals <- dplyr::filter(control_vals, treatment == 0)
    }
    vals <- rbind(treatment_vals, control_vals)
    vals <- vals[!is.na(vals$treatment), ]
    vals$treatment <- as.logical(vals$treatment)
    vals$region <- as.factor(vals$region)
    vals$biome <- as.factor(vals$biome)
    vals$ecoregion <- as.factor(vals$ecoregion)
    vals$pa <- as.factor(vals$pa)

    print(paste0(row_num, ': Matching'))
    m <- match_ae(vals, f)

    print(paste0(row_num, ': Formatting output'))
    if (is.null(m)) {
        print(paste0(row_num, ': no matches'))
    } else {
        m <- dplyr::select(m, c(get_names(f),
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
                         'fc_2015',
                         'fc_2016',
                         'fc_2017',
                         'fc_2018'))
        m$CI_ID <- CI_ID
        m$Data_Year <- Data_Year
        m <- m %>% dplyr::select(CI_ID, everything())
        print(paste0(row_num, ': saving output'))
    }
    save(m, file=paste0('Output/m_', CI_ID, '_', Data_Year, '.RData'))
    return(NULL)
}
