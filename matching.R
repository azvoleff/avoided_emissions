library(sf)
library(tidyverse)
library(units)
library(raster)
library(foreach)
library(optmatch)
library(mapview)
    
options("optmatch_max_problem_size"=Inf)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

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
save(ae, file='output_raw_matches.RData')

dplyr::select(sites, CI_ID, CI_Start_Year, CI_End_Year, Intervention, 
              Intervention_1, Intervention_2, Restoration) %>%
    right_join(ae) %>%
    mutate(cell_id=rownames(.)) -> ae

####TEMPORARY
load('ae_raw.RData')
####TEMPORARY

# Select initial and final forest cover based on start and end date fields for 
# each project
# TODO: need to add in bgb carbon calculation
dplyr::select(as.data.frame(ae), cell_id, CI_Start_Year, CI_End_Year, agb,
              starts_with('fc_'), -fc_change) %>%
    gather(year, forest_at_year_end, starts_with('fc_')) %>%
    mutate(year=as.numeric(str_replace(year, 'fc_', '')),
           CI_End_Year=ifelse(is.na(CI_End_Year), 2099, CI_End_Year)) %>%
    group_by(cell_id) %>%
    filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
    arrange(cell_id, year) %>%
    mutate(forest_loss_during_year=c(NA, diff(forest_at_year_end)),
           forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
           agb_at_year_end = agb * forest_frac_remaining,
           C_change=c(NA, diff(agb_at_year_end)) * .5,
           C_emissions_MgCO2e=C_change * -3.67) -> emissions
save(emissions, file='output_emissions_raw.Rdata')

emissions %>%
    filter(year >= CI_Start_Year) %>% # filter out the year prior to project start as no longer needed
    dplyr::select(cell_id, year, emissions) -> emissions_filtered
write.csv(emissions, file='output_emissions_filtered.csv', row.names=FALSE)

# TODO:
# 1) double-check calculations - need a ratio?
#  to convert biomass to carbon * .5
#  to convert change in C to CO2e * 3.67
# 2) what are units of avoided emissions - the agb is biomass above/below in 
#    tons


    group_by(cell_id) %>%
    summarise(forest_initial=fc[match(CI_Start_Year[1], year)],
              forest_final=fc[match(CI_End_Year[1], year)]) -> fc

write.csv(ae, file='ae_raw.csv', row.names=FALSE)

write.csv(emissions_details, file='ae_details.csv', row.names=FALSE)

emissions_summary <- emissions_details %>% group_by(CI_ID) %>%
        summarise(forest_loss_ha_treat_minus_control=forest_loss[treatment] - forest_loss[!treatment],
                  agb_loss_treat_minus_control=agb_change[treatment] - agb_change[!treatment],
                  n_treatment = n[treatment],
                  n_control = n[!treatment])
write.csv(emissions_summary, file='ae_summary.csv', row.names=FALSE)
