trajectory_df <- lab_df_w_klotho_bapi %>%
  filter(
    !is.na(LBVALUE),
    !is.na(time_since_baseline),
    as.numeric(time_since_baseline) %>% between(0,365)
  ) %>% 
  mutate(
    facet_var = paste0(LBTEST, ' [', LBUNIT, ']'),
    week = as.numeric(time_since_baseline)/7,
    week = ifelse(week > 0 & week <1, 1, week),
    week_int = cut2(week, week_bounds),
    week_int_pos = str_sub(week_int, 5,6) %>% as.numeric(),
    week_int_pos = if_else(week_int_pos == 1, 0, week_int_pos, 0),
    LBVALUE = ifelse(LBTEST == 'Klotho' & LBVALUE == 0,
                     .05,
                     LBVALUE),
    DRUG = case_when(
      DRUG == 'Etel' ~ 'ETL',
      DRUG == 'Alfa' ~ 'ALFA'
    ),
    DRUG = factor(DRUG, levels = c('ETL', 'ALFA')),
    med = ifelse(DRUG == 'ETL', paste0('A', DRUG), paste0('B', DRUG)) #just to fix the plot
  )


trajectory_df %>%
  mutate(facet_var = factor(facet_var, 
                            levels = c("FGF23 [pg/ml]",
                                       "PTH [ng/l]",
                                       "Calcium [mmol/l]",
                                       "Phosphate [mmol/l]",
                                       "BAP [ng/ml]",
                                       "Klotho [ng/ml]",
                                       "1,25-(OH)2-Vit-D [pg/ml]",
                                       "25-OH-Vit-D [nmol/l]"))) %>%
  ggplot(aes(x = week, y = LBVALUE, color = DRUG)) +
  facet_wrap(~factor(facet_var), scales = 'free', strip.position = 'left', ncol = 4) +
  geom_smooth(aes(group = DRUG), alpha = .2, color = NA, span = 5) +
  geom_line(stat = "smooth", alpha = .8, span = 5) +
  geom_boxplot(
    data = trajectory_df %>%
      filter(week_int_pos != 52),
    aes(x = week_int_pos, y = LBVALUE, group = paste(week_int_pos, LBTEST, med)),
    position = position_dodge(width = 7),
    width = 4,
    alpha = 1,
    fill = NA
  ) +
  theme_pubr() +
  scale_color_ptol() +
  theme(
    text = element_text(family = "Rockwell"),
    strip.placement = "outside"
  ) +
  labs(
    x = 'weeks after baseline',
    y = '',
    color = ''
  ) +
  scale_x_continuous(breaks = seq(0,48,12)) +
  # scale_y_continuous(trans = 'log2', breaks = pretty_breaks(n = 5))
  scale_y_continuous(trans = 'log2', breaks = trans_breaks('log2', function(x) 2^x), 
                     labels = function(x) round(x, 1))