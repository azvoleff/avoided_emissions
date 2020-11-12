library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyverse)
library(foreach)
library(tictoc)

readRDS('sites.RDS') %>%
    select(CI_ID,
           Data_Year,
           CI_Start_Date_clean, 
           CI_End_Date_clean) %>%
    mutate(CI_Start_Year=year(CI_Start_Date_clean),
           CI_End_Year=ifelse(is.na(year(CI_End_Date_clean)), 2099, year(CI_End_Date_clean))) %>%
    select(-CI_Start_Date_clean, -CI_End_Date_clean) %>%
    lazy_dt() %>%
    full_join(
        readRDS(paste0('Output/m_ALL.RDS')) %>%
        select(cell,
               CI_ID,
               Data_Year, 
               treatment,
               sampled_fraction,
               total_biomass,
               starts_with('fc_'))
        ) %>%
    as_tibble() -> m

m %>%
    group_by(CI_ID) %>%
    summarise(sampled_fraction=sampled_fraction[1]) %>%
    ggplot(aes(sampled_fraction)) +
    geom_histogram() +
    xlab('Fraction of site used in sample') + ylab('Number of sites') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_fraction_sampled.png', width=10, height=6)

foreach (this_data_year=c(2018, 2019, 2020)) %do% {
    tic()
    m %>%
        filter(Data_Year == this_data_year) %>%
        gather(year, forest_at_year_end, starts_with('fc_')) %>%
        lazy_dt() %>%
        mutate(year=2000 + as.numeric(str_replace(year, 'fc_', ''))) %>%
        group_by(cell) %>%
        filter(between(year, CI_Start_Year[1] - 1, CI_End_Year[1])) %>% # include one year prior to project start to get initial forest cover
        arrange(cell, year) %>%
        mutate(forest_at_year_end=forest_at_year_end  * (1 / sampled_fraction), # Correct for only having a sample
               total_biomass=total_biomass * (1 / sampled_fraction), # Correct for only having a sample
               forest_loss_during_year=c(NA, diff(forest_at_year_end)),
               forest_frac_remaining = forest_at_year_end / forest_at_year_end[1],
               biomass_at_year_end = total_biomass * forest_frac_remaining,
               #  to convert biomass to carbon * .5
               C_change=c(NA, diff(biomass_at_year_end)) * .5,
               #  to convert change in C to CO2e * 3.67
               Emissions_MgCO2e=C_change * -3.67) %>%
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
                  Emissions_MgCO2e_CI_minus_control=Emissions_MgCO2e[treatment] - Emissions_MgCO2e[!treatment]) %>%
        as_tibble() -> m_site
    toc()
    saveRDS(m_site, file=paste0('output_emissions_avoided_', this_data_year, '.RDS'))
    write_csv(m_site, paste0('output_emissions_avoided_', this_data_year, '.csv'))

    m_site %>%
        select(-forest_loss_ha_CI_minus_control)%>%
        pivot_wider(names_from=year,
                    values_from=Emissions_MgCO2e_CI_minus_control,
                    names_prefix='c_oavoid_') %>%
        select(order(colnames(.))) -> m_site_wide
    write_csv(m_site_wide, paste0('output_emissions_avoided_', this_data_year, '_wide.csv'))
    saveRDS(m_site_wide, file=paste0('output_emissions_avoided_', this_data_year, '_wide.RDS'))
}


# Plot how our sites due as a percentage of the emissions of the control sites


table(m_site_wide$c_oavoid_2018 > 0)
summary(m_site_wide$c_oavoid_2018)
sum(m_site_wide$c_oavoid_2018, na.rm=TRUE)

m_site_wide %>%
    ggplot() +
    geom_histogram(aes(c_oavoid_2018)) +
    ylab('Number of sites') + xlab('Avoided emissions (negative is good)') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_all.png', width=10, height=6)

m_site_wide %>%
    filter(c_oavoid > -1e5) %>%
    filter(c_oavoid < 1e5) %>%
    ggplot() +
    geom_histogram(aes(c_oavoid_2018)) +
    ylab('Number of sites') + xlab('Avoided emissions (negative is good)') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_all_filter.png', width=10, height=6)


m_site %>%
    group_by(CI_ID, Data_Year) %>%
    summarise(Emissions_MgCO2e_CI_minus_control=sum(Emissions_MgCO2e_CI_minus_control)) -> m_site_alltime
save(m_site_alltime, file='output_emissions_avoided_life_of_project.Rdata')
write_csv(m_site_alltime, 'output_emissions_avoided_life_of_project.csv')

# Size of intervention groups
m_site %>%
    group_by(Intervention) %>%
    summarise(n=length(unique(CI_ID))) %>%
    ggplot() +
    geom_bar(aes(Intervention, n), stat='identity')


# Plot of results by pixel
m %>%
    group_by(Intervention, CI_ID) %>%
    summarise(treat_minus_control=agb_change[treatment] - agb_change[!treatment]) %>%
    group_by(Intervention) %>%
    summarise(sum_treat_minus_control=sum(treat_minus_control, na.rm=TRUE),
              good=sum(treat_minus_control, na.rm=TRUE) < 0) %>%
    ggplot() +
    geom_bar(aes(Intervention, sum_treat_minus_control, fill=good), stat='identity') +
    ylab('Total emissions in treatment relative to control (units??)')

# Plot of results by pixel
m %>%
    group_by(Intervention, CI_ID) %>%
    summarise(treat_minus_control=agb_change[treatment] - agb_change[!treatment]) %>%
    group_by(Intervention) %>%
    summarise(mean_treat_minus_control=mean(treat_minus_control, na.rm=TRUE),
              good=sum(treat_minus_control, na.rm=TRUE) < 0) %>%
    ggplot() +
    geom_bar(aes(Intervention, mean_treat_minus_control, fill=good), stat='identity') +
    ylab('Mean emissions in treatment relative to control (units??)')

# Plot of avoided deforestation
m_summary %>%
    ggplot() +
    geom_point(aes(CI_ID, forest_loss_ha_treat_minus_control, colour=good))

# Plot of avoided deforestation
    
# Plot of sample sizes
m_details %>%
    ggplot() +
    geom_point(aes(CI_ID, n, colour=treatment), 
               position=position_dodge(width=.5)) +
    theme(legend.title = element_blank())
    
