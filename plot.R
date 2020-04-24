library(tidyverse)
library(foreach)
library(sf)

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

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

load('sites.RData')

ae <- foreach(f=list.files('Output'), .combine=foreach_rbind) %do% {
    load(paste0('Output/', f))
    return(m)
}
dplyr::select(sites, CI_ID, Data_Year, CI_Start_Year, CI_End_Year, Intervention, 
              Intervention_1, Intervention_2, Restoration) %>%
    right_join(ae) %>%
    mutate(cell_id=rownames(.)) -> ae
save(ae, file='output_raw_matches.RData')

# Select initial and final forest cover based on start and end date fields for 
# each project
dplyr::select(ae, cell_id, treatment, CI_ID, Data_Year, CI_Start_Year, 
              CI_End_Year, total_biomass,
              starts_with('fc_'), -fc_change) %>%
    gather(year, forest_at_year_end, starts_with('fc_')) %>%
    mutate(year=as.numeric(str_replace(year, 'fc_', '')),
           CI_End_Year=ifelse(is.na(CI_End_Year), 2099, CI_End_Year)) %>%
    group_by(cell_id) %>%
    filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
    arrange(cell_id, year) %>%
    mutate(forest_loss_during_year=c(NA, diff(forest_at_year_end)),
           forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
           biomass_at_year_end = total_biomass * forest_frac_remaining,
           #  to convert biomass to carbon * .5
           C_change=c(NA, diff(biomass_at_year_end)) * .5,
           #  to convert change in C to CO2e * 3.67
           C_emissions_MgCO2e=C_change * -3.67) -> e
save(e, file='output_emissions_raw.Rdata')

as.data.frame(e) %>%
    group_by(CI_ID, Data_Year, year) %>%
    filter(between(year, CI_Start_Year, CI_End_Year)) %>% # filter out the year prior to project start as no longer needed
    summarise(forest_loss_ha_treat_minus_control=sum(forest_loss_during_year[treatment], na.rm=TRUE) - sum(forest_loss_during_year[!treatment], na.rm=TRUE),
              C_emissions_MgCO2e_treat_minus_control=sum(C_emissions_MgCO2e[treatment], na.rm=TRUE) - sum(C_emissions_MgCO2e[!treatment], na.rm=TRUE)) -> e_avoided
save(e_avoided, file='output_emissions_avoided.RData')
write.csv(e_avoided, file='output_emissions_avoided.csv', row.names=FALSE)

# Size of intervention groups
ae %>%
    group_by(Intervention) %>%
    summarise(n=length(unique(CI_ID))) %>%
    ggplot() +
    geom_bar(aes(Intervention, n), stat='identity')


# Plot of results by pixel
ae %>%
    group_by(Intervention, CI_ID) %>%
    summarise(treat_minus_control=agb_change[treatment] - agb_change[!treatment]) %>%
    group_by(Intervention) %>%
    summarise(sum_treat_minus_control=sum(treat_minus_control, na.rm=TRUE),
              good=sum(treat_minus_control, na.rm=TRUE) < 0) %>%
    ggplot() +
    geom_bar(aes(Intervention, sum_treat_minus_control, fill=good), stat='identity') +
    ylab('Total emissions in treatment relative to control (units??)')

# Plot of results by pixel
ae %>%
    group_by(Intervention, CI_ID) %>%
    summarise(treat_minus_control=agb_change[treatment] - agb_change[!treatment]) %>%
    group_by(Intervention) %>%
    summarise(mean_treat_minus_control=mean(treat_minus_control, na.rm=TRUE),
              good=sum(treat_minus_control, na.rm=TRUE) < 0) %>%
    ggplot() +
    geom_bar(aes(Intervention, mean_treat_minus_control, fill=good), stat='identity') +
    ylab('Mean emissions in treatment relative to control (units??)')

# Plot of avoided deforestation
ae_summary %>%
    ggplot() +
    geom_point(aes(CI_ID, forest_loss_ha_treat_minus_control, colour=good))

# Plot of avoided deforestation
    
# Plot of sample sizes
ae_details %>%
    ggplot() +
    geom_point(aes(CI_ID, n, colour=treatment), 
               position=position_dodge(width=.5)) +
    theme(legend.title = element_blank())
    
