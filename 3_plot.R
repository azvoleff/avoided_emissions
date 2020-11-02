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

load('sites.RData')

ae <- foreach(f=list.files('Output'), .combine=foreach_rbind) %do% {
    load(paste0('Output/', f))
    return(m)
}

ae %>%
    group_by(CI_ID) %>%
    summarise(sampled_fraction=sampled_fraction[1]) %>%
    ggplot(aes(sampled_fraction)) +
    geom_histogram() +
    xlab('Fraction of site used in sample') + ylab('Number of sites') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_fraction_sampled.png', width=10, height=6)

# Join the site data to the tibble
dplyr::select(sites, CI_ID, Data_Year, CI_Start_Year, CI_End_Year, Intervention, 
              Intervention_1, Intervention_2, Restoration) %>%
    right_join(ae) %>%
    mutate(cell_id=rownames(.)) -> ae

# Calculate emissions for each year
dplyr::select(ae, cell_id, treatment, sampled_fraction, CI_ID, Data_Year, 
              CI_Start_Year, CI_End_Year, total_biomass, sampled_fraction, 
              starts_with('fc_'), -fc_change) %>%
    gather(year, forest_at_year_end, starts_with('fc_')) %>%
    mutate(year=as.numeric(str_replace(year, 'fc_', '')),
           CI_End_Year=ifelse(is.na(CI_End_Year), 2099, CI_End_Year)) %>%
    group_by(cell_id) %>%
    filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
    arrange(cell_id, year) %>%
    mutate(forest_at_year_end=forest_at_year_end  * (1 / sampled_fraction), # Correct for only having a sample
           total_biomass=total_biomass * (1 / sampled_fraction), # Correct for only having a sample
           forest_loss_during_year=c(NA, diff(forest_at_year_end)),
           forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
           biomass_at_year_end = total_biomass * forest_frac_remaining,
           #  to convert biomass to carbon * .5
           C_change=c(NA, diff(biomass_at_year_end)) * .5,
           #  to convert change in C to CO2e * 3.67
           Emissions_MgCO2e=C_change * -3.67) -> ae

# Calculate avoided emissions for each year for each site
as.data.frame(ae) %>%
    group_by(CI_ID, Data_Year, year, treatment) %>%
    filter(between(year, CI_Start_Year[1], CI_End_Year[1])) %>% # filter out the year prior to project start as no longer needed
    summarise(CI_Start_Year=CI_Start_Year[1],
              CI_End_Year=CI_End_Year[1],
              # correct totals for areas where only a partial sample was used 
              # by taking into account the fraction sampled
              forest_loss_ha=sum(forest_loss_during_year, na.rm=TRUE) * (1 / sampled_fraction[1]),
              Emissions_MgCO2e=sum(Emissions_MgCO2e, na.rm=TRUE) * (1 / sampled_fraction[1])) %>%
    group_by(CI_ID, Data_Year, year) %>%
    summarise(CI_Start_Year=CI_Start_Year[1],
              CI_End_Year=CI_End_Year[1],
              forest_loss_ha_CI_minus_control=forest_loss_ha[treatment] - forest_loss_ha[!treatment],
              Emissions_MgCO2e_CI_minus_control=Emissions_MgCO2e[treatment] - Emissions_MgCO2e[!treatment]) -> ae_site
save(ae_site, file='output_emissions_avoided.Rdata')
write_csv(ae_site, path='output_emissions_avoided.csv')

ae_site %>%
    select(-forest_loss_ha_CI_minus_control)%>%
    pivot_wider(names_from=year,
                values_from=Emissions_MgCO2e_CI_minus_control,
                names_prefix='c_oavoid_') %>%
    mutate(c_oavoid=c_oavoid_2018) %>% # Duplicate the c_oavoid column as requested
    select(order(colnames(.))) -> ae_site_wide
save(ae_site_wide, file='output_emissions_avoided_wide.Rdata')
write_csv(ae_site_wide, path='output_emissions_avoided_wide.csv')

table(ae_site_wide$c_oavoid_2018 > 0)
summary(ae_site_wide$c_oavoid_2018)
sum(ae_site_wide$c_oavoid_2018)

ae_site_wide %>%
    ggplot() +
    geom_histogram(aes(c_oavoid_2018)) +
    ylab('Number of sites') + xlab('Avoided emissions (negative is good)') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_all.png', width=10, height=6)

ae_site_wide %>%
    filter(c_oavoid > -1e5) %>%
    filter(c_oavoid < 1e5) %>%
    ggplot() +
    geom_histogram(aes(c_oavoid_2018)) +
    ylab('Number of sites') + xlab('Avoided emissions (negative is good)') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_all_filter.png', width=10, height=6)


ae_site %>%
    group_by(CI_ID, Data_Year) %>%
    summarise(Emissions_MgCO2e_CI_minus_control=sum(Emissions_MgCO2e_CI_minus_control)) -> ae_site_alltime
save(ae_site_alltime, file='output_emissions_avoided_life_of_project.Rdata')
write_csv(ae_site_alltime, path='output_emissions_avoided_life_of_project.csv')

# Size of intervention groups
ae_site %>%
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
    