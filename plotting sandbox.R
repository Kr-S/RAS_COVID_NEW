librarian::shelf(foreign, ggarrange, tidyverse, dplyr, lubridate, stringr, readxl, ggthemes, lme4, purrr, broom, magrittr, ggpubr, ggbeeswarm, arsenal, knitr, lmerTest)



smooth_plot <- function(data, variable, ylab = ''){  
  p <- data %>%
    filter(days_since_test %>% between(0,21)) %>%
    ggplot(aes_string(x = 'days_since_test', y = variable, color = 'class', group = 'pat_id')) +
    geom_line(aes(linetype = class), alpha = .3, size = .8) + 
    geom_smooth(aes(group = as.factor(class)), alpha= .2) +
    scale_y_log10() +
    theme_pubr() +
    scale_color_ptol() +
    labs(
      color = '',
      linetype = '',
      y = ylab,
      x = 'days since first test'
    )
  ggsave(paste0('plots/figure2B_', variable, '.png'), width = 5, height = 4)
  p
}


p_b1 <- ras_df %>% smooth_plot('Ang_1_7', 'Angiotensin 1-7')
p_b2 <- ras_df %>% smooth_plot('ACE2', 'ACE2 concentration')
p_b3 <- ras_df %>% smooth_plot('ACE_S', 'ACE activity')
p_b4 <- ras_df %>% smooth_plot('AngII', 'Angiotensin II')
p_b5 <- ras_dfpot %>% smooth_plot('potassium', 'potassium')
p_b6 <- ras_dfpot %>% smooth_plot('Aldosterone', 'aldosterone')
p_b7 <- ras_dfpot %>% smooth_plot('potassium_urine', 'potassium_urine')

ggarrange(p_b5, p_b6, p_b7)