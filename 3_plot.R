library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tidyverse)
library(foreach)
library(tictoc)

m_site <- readRDS('output_raw_by_site.RDS')

m_site %>%
    group_by(name, year) %>%
    mutate(forest_loss_ha=abs(forest_loss_ha)) %>%  # Make forest loss positive since it is a measure in ha
    summarise(CI_Start_Year=CI_Start_Year[1],
              CI_End_Year=CI_End_Year[1],
              forest_loss_ha_CI_minus_control=forest_loss_ha[treatment] - forest_loss_ha[!treatment],
              Emissions_MgCO2e_CI_minus_control=Emissions_MgCO2e[treatment] - Emissions_MgCO2e[!treatment],
              fc_loss_pct_of_control=(sum(forest_loss_ha[treatment], na.rm=TRUE) /
                                      sum(forest_loss_ha[!treatment], na.rm=TRUE))*100) %>%
    relocate(year, .after=CI_End_Year) -> m_by_year
saveRDS(m_by_year, file=paste0('output_emissions_avoided.RDS'))
write_csv(m_by_year, paste0('output_emissions_avoided.csv'))
m_by_year %>%
    select(-forest_loss_ha_CI_minus_control, -fc_loss_pct_of_control) %>%
    pivot_wider(names_from=year,
                values_from=Emissions_MgCO2e_CI_minus_control,
                names_prefix='c_oavoid_') -> m_by_year_wide
write_csv(m_by_year_wide, paste0('output_emissions_avoided_wide.csv'))
saveRDS(m_by_year_wide, file=paste0('output_emissions_avoided_wide.RDS'))

# Plot how our sites due as a percentage of the emissions of the control sites
m_site %>%
    filter(Data_Year == year + 1) %>%
    select(-Data_Year) %>%
    group_by(year, CI_ID) %>%
    summarise(fc_loss_pct_of_control=(sum(forest_loss_ha[treatment], na.rm=TRUE) /
                                      sum(forest_loss_ha[!treatment], na.rm=TRUE))*100,
              lt100=factor(sum(fc_loss_pct_of_control < 100, na.rm=TRUE))) %>%
    group_by(year) %>%
    arrange(fc_loss_pct_of_control) %>%
    filter(fc_loss_pct_of_control <= 1000) %>%
    mutate(order=1:n()) %>%
    ggplot() +
    geom_bar(aes(x=order, fc_loss_pct_of_control, colour=lt100, fill=lt100), stat='identity') +
    facet_grid(year ~ ., scales='free')
ggsave('sites_filtered.png')

m_by_year %>%
    group_by(year) %>%
    summarise(gt100=sum(fc_loss_pct_of_control > 100, na.rm=TRUE),
              lt100=sum(fc_loss_pct_of_control < 100, na.rm=TRUE),
              lt0=sum(fc_loss_pct_of_control < 0, na.rm=TRUE),
              avoided_em=sum(Emissions_MgCO2e_CI_minus_control, na.rm=TRUE))

m %>%
    group_by(CI_ID) %>%
    summarise(sampled_fraction=sampled_fraction[1]) %>%
    ggplot(aes(sampled_fraction)) +
    geom_histogram() +
    xlab('Fraction of site used in sample') + ylab('Number of sites') +
    theme_bw(base_size=18)
ggsave('output_emissions_avoided_fraction_sampled.png', width=10, height=6)

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
    
