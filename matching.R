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
            # Use Mahalanobis distance if there aren't enough points to run a
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
###  Load sites and covariates

d <- stack('all_covariates.tif')
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


ae <- foreach(row_num=21:22,
#ae <- foreach(row_num=1:nrow(sites),
             .packages=c('raster', 'rgeos', 'optmatch', 'dplyr', 'foreach'),
             .combine=foreach_rbind) %do% {
    print(row_num )
    p <- sites[row_num, ]
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
        r <- crop(regions_rast, regions[inter, ])
    }

    p_rast <- rasterize(p, r, background=0)
    # Note that the below will drop any portions of a site that don't fall 
    # within a region in GADM- this means that marine areas will be dropped. 
    # Given this calculation is for avoided emissions, likely not a problem... 
    #
    # Treatment is 1, control is zero
    #
    # TODO: Need to ensure we also mask out other areas where CI has ongoing 
    # interventions.
    treat_or_control <- mask(p_rast, r)
    names(treat_or_control) <- 'treatment'

    # Crop matching vars to the area of interest (the region(s) the site falls 
    # within) and then extract values
    d <- crop(d, r)
    d <- stack(r, treat_or_control, d)

    # Project all items to cylindrical equal area
    # d<- projectRaster(d, crs=CRS('+proj=cea'), method='ngb')

    # Process the pixels in blocks as some regions are large
    bs <- blockSize(d)
    m <- foreach (i=1:bs$n) %do% {
        #vals <- data.frame(getValues(d))
        vals <- data.frame(getValues(d, row=bs$row[i], nrows=bs$nrows[i]))
        vals <- vals[!is.na(vals$treatment), ]
        vals$treatment <- as.logical(vals$treatment)
        vals$biome <- as.factor(vals$biome)
        vals$ecoregion <- as.factor(vals$ecoregion)
        vals$pa <- as.factor(vals$pa)

        m <- match_ae(vals, f)
    }

    if (is.null(m)) {
        return(NULL)
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
        m$CI_ID <- as.character(p$CI_ID)
        return(m %>% dplyr::select(CI_ID, everything()))
    }
}

dplyr::select(sites, CI_ID, CI_Start_Year, CI_End_Year, Intervention, 
              Intervention_1, Intervention_2, Restoration) %>%
    right_join(ae) %>%
    mutate(cell_id=rownames(.)) -> ae

# Select initial and final forest cover based on start and end date fields for 
# each project
dplyr::select(as.data.frame(ae), cell_id, CI_Start_Year, CI_End_Year, 
              starts_with('fc_'), -fc_change) %>%
    gather(year, fc, starts_with('fc_')) %>%
    mutate(year=as.numeric(str_replace(year, 'fc_', ''))) %>%
    group_by(cell_id) %>%
    summarise(forest_initial=fc[match(CI_Start_Year[1], year)],
              forest_final=fc[match(CI_End_Year[1], year)]) -> fc

save(ae, file='ae_raw.RData')
write.csv(ae, file='ae_raw.csv', row.names=FALSE)

# TODO: Need to ensure agb change calculations are based on initial and final 
# forest cover for correct years for intervention and year of CI sites
emissions_details <- ae %>%
    group_by(CI_ID, treatment) %>%
    summarise(forest_initial=sum(forest_initial, na.rm=TRUE),
              forest_final=sum(forest_final, na.rm=TRUE),
              forest_loss = forest_initial - forest_final,
              agb_initial = sum(agb, na.rm=TRUE),
              agb_final = agb_initial * ((forest_final - forest_initial)/forest_initial + 1),
              agb_change = agb_final - agb_initial,
              n = n())
write.csv(emissions_details, file='ae_details.csv', row.names=FALSE)

emissions_summary <- emissions_details %>% group_by(CI_ID) %>%
        summarise(forest_loss_ha_treat_minus_control=forest_loss[treatment] - forest_loss[!treatment],
                  agb_loss_treat_minus_control=agb_change[treatment] - agb_change[!treatment],
                  n_treatment = n[treatment],
                  n_control = n[!treatment])
write.csv(emissions_summary, file='ae_summary.csv', row.names=FALSE)
