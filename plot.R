library(ggplot2)
library(dplyr)
library(sf)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

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
    
