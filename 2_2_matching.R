library(tidyverse)
library(optmatch)
library(foreach)
library(lubridate)
library(biglm)
library(tictoc)

options("optmatch_max_problem_size"=Inf)

MAX_TREATMENT <- 2000
CONTROL_MULTIPLIER <- 50
    
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

# Basic function to extract variable names from a formula object
get_names <- function(f) {
    f <- paste0(as.character(f), collapse=' ')
    v <- strsplit(f, split='[+ ~]')[[1]]
    v <- v[v != '']
    gsub('strata\\(([a-zA-Z_]*)\\)', '\\1', v)
}

match_ae <- function(d, f) {
    if (sum(d$treatment) > 30) {
        if (nlevels(d$group) >= 2) {
            f <- update(f, ~ . + strata(group))
        } else {
            f <- update(f, ~ . - strata(group))
        }
        model <- glm(f, data=d, family=binomial())
        dists <- match_on(model, data=d)
    } else {
        # Use Mahalanobis distance if there aren't enough points to run a
        # glm
        dists <- match_on(f, data=d)
    }

    # If the controls are too far from the treatments (due to a caliper) then 
    # the matching may fail. Can test for this by seeing if subdim runs 
    # successfully
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

###############################################################################
###  Load sites and covariates

treatment_key <- readRDS('Output/treatment_cell_key.RDS')
readRDS('sites.RDS') %>%
    select(-geometry) %>%
    as_tibble() -> sites

# Filter to include only values of group that appear in the treatment pixels, 
# and to not include values that appear only in the treatment pixels
filter_groups <- function(vals) {
    vals$group <- interaction(vals$region, vals$ecoregion, vals$pa)
    vals <- filter(vals, group %in% unique(filter(vals, treatment)$group))
    treatment_groups <- unique(filter(vals, treatment)$group)
    control_groups <- unique(filter(vals, !treatment)$group)
    vals <- filter(vals, group %in% treatment_groups[treatment_groups %in% control_groups])
    # Filter out values of group that appear ONLY in the treatment pixels
    vals$group <- droplevels(vals$group)
    return(vals)
}

###############################################################################
###  Run matching

set.seed(31)
# TODO: Need to account for the data year
year <- 2018

###############################################################################
# TODO: Remove this, for testing only
ids <- unique(treatment_key$CI_ID)[c(12, 13, 30)]
treatment_key <- treatment_key %>% filter(CI_ID %in% ids)
###############################################################################

ae <- foreach(this_CI_ID=unique(treatment_key$CI_ID),
              .combine=foreach_rbind, .inorder=FALSE) %do% {

    tic()
    site <- filter(sites,
                   CI_ID == this_CI_ID,
                   Data_Year == year)
    if (file.exists(paste0('Output/m_', this_CI_ID, '_', year, '.RDS'))) {
        print(paste0('Skipping ', this_CI_ID, ' for year ', year, '. Already processed.'))
        return(NULL)
    } else {
        print(paste0('Processing ', this_CI_ID, ' for year ', year, '.'))
    }

    treatment_cell_IDs <- filter(treatment_key,
                                 CI_ID == this_CI_ID,
                                 !is.na(region),
                                 Data_Year == year)

    n_treatment_cells_total <- nrow(treatment_cell_IDs)
    vals <- foreach(this_region = unique(treatment_cell_IDs$region),
                    .combine=rbind) %do% {
        v <- readRDS(paste0('Output/treatments_and_controls_', this_region, 
                            '.RDS'))
        filter(v, region == this_region)
    }

    # TODO: need to knockout other CI sites from control sample by setting 
    # those areas to NA in treatment flag
    vals %>% full_join(
            treatment_cell_IDs %>%
                select(cell, Data_Year) %>%
                mutate(treatment=TRUE)
            , by='cell') -> vals
    vals$treatment <- as.logical(vals$treatment)
    vals$treatment[is.na(vals$treatment)] <- FALSE
    vals$Data_Year <- year

    # Eliminate any groups that are only in the control pixels, or only in the 
    # treatment pixels
    vals <- filter_groups(vals)

    sample_sizes <- vals %>%
        count(treatment, group)
    # Sample the treatment cells if there are more than MAX_TREATMENT pixels, 
    # and the control cells if there are more than CONTROL_MULTIPLIER * 
    # MAX_TREATMENT pixels
    bind_rows(
            filter(vals, treatment)  %>%
                group_by(group) %>%
                sample_n(min(MAX_TREATMENT, n())),
            filter(vals, !treatment)  %>%
                group_by(this_group=group) %>%
                sample_n(min(CONTROL_MULTIPLIER * filter(sample_sizes,
                                                         treatment == TRUE,
                                                         group == this_group[1])$n,
                             n()))
            ) %>% 
        ungroup() -> vals
    # Refilter in case any groups were lost due to the sampling
    vals <- filter_groups(vals)

    # Project all items to cylindrical equal area
    # d_crop <- projectRaster(d_crop, crs=CRS('+proj=cea'), method='ngb')
    
    f <- readRDS('Output/formula.RDS')
    f <- update(f, ~ . + group)

    # For sites that were established in or after 2005, match on the five years 
    # of deforestation data preceding the year of establishment. For sites 
    # estab prior to 2005, don't match on defor rate
    estab_year <- year(site$CI_Start_Date_clean)
    if (estab_year >= 2005) {
        init <- vals[, grepl(paste0('fc_', substr(estab_year - 5, 3, 4)), names(vals))]
        final <- vals[,grepl(paste0('fc_', substr(estab_year, 3, 4)), names(vals))]
        defor_pre_intervention <- ((final - init) / init) * 100
        # Correct for division by zero in places that had no forest cover in 
        # year 0
        defor_pre_intervention[init == 0] <- 0
        names(defor_pre_intervention) <- 'defor_pre_intervention'
        vals <- cbind(vals, defor_pre_intervention)
        f <- update(f, ~ . + defor_pre_intervention)
    }
    #vals <- vals %>% select(-starts_with('fc_'), -starts_with('fcc_'))
    
    sample_sizes <- vals %>%
        count(treatment, group)
    print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, treatment)$n, collapse=', '), ' treatment pixels'))
    print(paste0(this_CI_ID, ': ', paste(filter(sample_sizes, !treatment)$n, collapse=', '), ' control pixels'))

    if (nrow(filter(vals, treatment)) == 0) {
        print(paste0(this_CI_ID, ': No treatment values remaining after filtering'))
        return(NULL)
    } else {
        print(paste0(this_CI_ID, ': Matching'))
        m <- match_ae(vals, f)
        print(paste0(this_CI_ID, ': Formatting output'))
        if (is.null(m)) {
            print(paste0(this_CI_ID, ': no matches'))
        } else {
            m$CI_ID <- this_CI_ID
            m$Data_Year <- year
            m <- m %>% dplyr::select(CI_ID, everything())
            print(paste0(this_CI_ID, ': saving output'))
        m$sampled_fraction <- n_treatment_cells_sample / n_treatment_cells_total
        saveRDS(m, paste0('Output/m_', this_CI_ID, '_', year, '.RDS'))
        }
    }
    toc()
    return(m)
}
