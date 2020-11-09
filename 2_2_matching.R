library(tidyverse)
library(optmatch)
library(foreach)

MAX_TREATMENT <- 2000
    
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
    # Filter out sites without at least one treatment unit or without at
    # least one control unit
    d <- d %>%
        filter(complete.cases(.)) %>%
        mutate(n_treatment=sum(treatment),
               n_control=sum(!treatment)) %>%
        filter(n_treatment >= 1, n_control >= 1)

    # Note custom combine to handle iterations that don't return any value
    d_treatment <- filter(d, treatment)
    # Filter out climates and land covers that don't appear in the CI
    # sample, and drop these levels from the factors
    d <- filter(d,
        region %in% unique(d_treatment$region),
        ecoregion %in% unique(d_treatment$ecoregion),
        pa %in% unique(d_treatment$pa))
    d$region <- droplevels(d$region)
    d$ecoregion <- droplevels(d$ecoregion)
    d$pa <- droplevels(d$pa)
    # Can't stratify by land cover or climate if they only have one level
    if (nlevels(d$region) >= 2) {
        f <- update(f, ~ . + strata(region))
    } else {
        f <- update(f, ~ . - region)
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
    if (nrow(d_treatment) > 2) {
        model <- glm(f, data=d)
        dists <- match_on(model, data=d)
    } else {
        # Use Mahalanobis distance if there aren't enough points to run a
        # glm
        dists <- match_on(f, data=d)
    }
    # If the controls are too far from the treatments (due to the caliper) then 
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

f <- readRDS('Output/formula.RDS')
treatment_key <- readRDS('Output/treatment_cell_key.RDS')

###############################################################################
###  Run matching
set.seed(31)
this_CI_ID <- unique(treatment_key$CI_ID)[1]

year <- 2018
this_CID_ID <- unique(treatment_key$CI_ID)[12]
# TODO: Need to account for the data year
ae <- foreach(this_CI_ID in unique(treatment_key$CI_ID),
              .combine=foreach_rbind, .inorder=FALSE) %do% {

    if (file.exists(paste0('Output/m_', this_CI_ID, '_', year, '.RDS'))) {
        print(paste('Skipping', this_CI_ID, year, '. Already processed.'))
        return(NULL)
    }

    # TODO: Need to account for the data year
    treatment_cell_IDs <- filter(treatment_key,
                                 CI_ID == this_CI_ID,
                                 Data_Year == year)
    n_treatment_cells_total <- nrow(treatment_cell_IDs)
    print(paste0(this_CI_ID, ': ', n_treatment_cells_total, ' total treatment cells'))
    ###
    # Sample the treatment cells if there are more than MAX_TREATMENT
    if (n_treatment_cells_total <= MAX_TREATMENT) {
        n_treatment_cells_sample <- n_treatment_cells_total
    } else {
        n_treatment_cells_sample <- MAX_TREATMENT
        treatment_cell_IDs <- sample_n(treatment_cell_IDs, MAX_TREATMENT)
        print(paste0(this_CI_ID, ': sampled ', nrow(treatment_cell_IDs), ' treatment cells'))
    }

    vals <- foreach(this_region = unique(treatment_cell_IDs$region),
                    .combine=rbind) %do% {
        v <- readRDS(paste0('Output/treatments_and_controls_', this_region, 
                            '.RDS'))
        filter(v, region == this_region)
    }

    # Sample the control cells if there are more than 100*MAX_TREATMENT
    #
    # Note the below factor is set slightly higher than the one used with 
    # sampling to ensure that there are enough cells to sample from
    if (nrow(vals) < (120 * n_treatment_cells_sample)) {
        print(paste0(this_CI_ID, ': using all possible control points'))
    } else {
        vals <- sample_n(vals, 100 * n_treatment_cells_sample)
        print(paste0(this_CI_ID, ': sampled ', nrow(vals), ' control cells'))
    }

    # TODO: need to knockout other CI sites from control sample by setting 
    # those areas to NA in treatment flag
    vals %>% full_join(
            treatment_cell_IDs %>%
                select(cell, Data_Year) %>%
                mutate(treatment=TRUE)
            ) -> vals
    vals$treatment[is.na(vals$treatment)] <- FALSE

    # Project all items to cylindrical equal area
    # d_crop <- projectRaster(d_crop, crs=CRS('+proj=cea'), method='ngb')

    # Remove any points with missing data for treatment indicator
    vals <- vals[!is.na(vals$treatment), ]
    print(paste0(this_CI_ID, ': ', nrow(vals), ' total cells in matching'))
    vals$treatment <- as.logical(vals$treatment)
    vals$region <- as.factor(vals$region)
    vals$ecoregion <- as.factor(vals$ecoregion)
    vals$pa <- as.factor(vals$pa)

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
            m <- dplyr::select(m, c(get_names(f),
                             'region',
                             'ecoregion',
                             'pa',
                             'fc_00',
                             'fc_01',
                             'fc_02',
                             'fc_03',
                             'fc_04',
                             'fc_05',
                             'fc_06',
                             'fc_07',
                             'fc_08',
                             'fc_09',
                             'fc_10',
                             'fc_11',
                             'fc_12',
                             'fc_13',
                             'fc_14',
                             'fc_15',
                             'fc_16',
                             'fc_17',
                             'fc_18',
                             'fc_19'))
            m$CI_ID <- CI_ID
            m$Data_Year <- Data_Year
            m <- m %>% dplyr::select(CI_ID, everything())
            print(paste0(this_CI_ID, ': saving output'))
        m$sampled_fraction <- n_treatment_cells_sample / n_treatment_cells_total
        saveRDS(m, paste0('Output/m_', CI_ID, '_', Data_Year, '.RDS')
        }
    }

    return(m)
}

d_crop <- as.data.table(d)
