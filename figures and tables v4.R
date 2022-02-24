# figures and tables for SARS FP paper, as laid out in analysis plan 3
# author: Sebastian Hödlmoser
# date: 04.08.20


# libs and data ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stringr)
library(readxl)
library(ggthemes)
library(lme4)
library(purrr)
library(broom)
library(magrittr)
library(ggpubr)
library(ggbeeswarm)
library(arsenal)
library(knitr)
library(lmerTest)

ras_df <- readRDS('data/RAS_df.rds') %>%
  filter(days_since_test >= 0)

controls <- readRDS('data/controls_df.rds')

print_p <- function(p){
  case_when(
    p < 0.001 ~ '<0.001*',
    p < 0.01  ~ '<0.01*',
    p < 0.05  ~ '<0.05*',
    TRUE            ~ paste0(round(p,3))
  )}

display <- function(x, digits = 2){
  sprintf(paste0('%.', digits, 'f'), x)
}

# end libs and data -----------------------------------------------------------

################
#### TABLES ####
################

# medians of angiotensins -----------------------------------------------------
# mean per patienet, median over all patients of all angiotensins

ras_fingerprints <- ras_df %>%
  group_by(pat_id) %>%
  summarise_at(
    vars(
      ACE2,
      ACE2_ang_based,
      AngI,
      AngII, 
      AngIII,
      AngIV,
      Ang_1_5,
      Ang_1_7,
      Ang_1_7_plus_1_5,
      Aldosterone,
      AA2_ratio,
      PRA_S,
      ACE_S
    ),
    mean, na.rm = T
  ) %>%
  summarise_at(
    vars(
      ACE2,
      ACE2_ang_based,
      AngI,
      AngII, 
      AngIII,
      AngIV,
      Ang_1_5,
      Ang_1_7,
      Ang_1_7_plus_1_5,
      Aldosterone,
      AA2_ratio,
      PRA_S,
      ACE_S
    ),
    median, na.rm = T
  )

write_csv(ras_fingerprints, 'tables/ras_fingerprints_covid19_patients.csv')  


# end medians of angiotensins -------------------------------------------------


# n observations summary ------------------------------------------------------

tmp <- ras_df %>%
  bind_rows(controls %>% filter(type %in% c('healthy', 'influenza'))) %>%
  select(
    pat_id,
    ACE2,
    ACE2_ang_based,
    AngI,
    AngII, 
    AngIII,
    AngIV,
    Ang_1_5,
    Ang_1_7,
    ACE_S
  ) %>%
  pivot_longer(
    -pat_id,
    names_to = 'angiotensin',
    values_to = 'value'
  ) %>%
  na.omit()

tmp %>%
  group_by(angiotensin, pat_id) %>%
  summarise(
    nobs = n()
  ) %>%
  group_by(angiotensin) %>%
  summarise(
    nobs_all = sum(nobs),
    nobs_median = round(median(nobs), 2),
    `Q1/Q3` = paste0(round(quantile(nobs, c(.25, .75)), 2), collapse = ' / '),
    `min/max` = paste0(c(min(nobs), max(nobs)), collapse = ' / ')
  ) %>%
  write_csv('tables/n_observations_summary_with_controls.csv')

tmp <- ras_df %>%
  select(
    pat_id,
    ACE2,
    ACE2_ang_based,
    AngI,
    AngII, 
    AngIII,
    AngIV,
    Ang_1_5,
    Ang_1_7,
    ACE_S
  ) %>%
  pivot_longer(
    -pat_id,
    names_to = 'angiotensin',
    values_to = 'value'
  ) %>%
  na.omit()

tmp %>%
  group_by(angiotensin, pat_id) %>%
  summarise(
    nobs = n()
  ) %>%
  group_by(angiotensin) %>%
  summarise(
    nobs_all = sum(nobs),
    nobs_median = round(median(nobs), 2),
    `Q1/Q3` = paste0(round(quantile(nobs, c(.25, .75)), 2), collapse = ' / '),
    `min/max` = paste0(c(min(nobs), max(nobs)), collapse = ' / ')
  ) %>%
  write_csv('tables/n_observations_summary_without_controls.csv')

# end n observations summary --------------------------------------------------


# table1 ----------------------------------------------------------------------
# patient characteristics

table_1a_df <- controls %>%
  filter(
    type %in% c('healthy', 'influenza')
  ) %>%
  mutate(
    survival_28days = death,
    hydroxychloroquin = NA, 
    Kaletra = NA, 
    Remdesivir = NA,
    Plasma = NA, 
    Tocilizumab = NA, 
    Camostat = NA, 
    n_obs = NA, 
    length_of_hospital_stay = NA,
    steriods = NA,
    bmi = NA
  )

table_1a <- tableby(
  type ~ age + gender + bmi + ace_inhibitor + arb + diabetes + hypertension + copd + 
    hydroxychloroquin + Kaletra + Remdesivir + Plasma + Tocilizumab + Camostat + steriods +
    n_obs + length_of_hospital_stay + death,
  data = table_1a_df,
  numeric.stats=c('mean', 'sd', "median","q1q3", 'range'),
  total = FALSE,
  test = FALSE
)


table_1b_df <- ras_df %>%
  distinct(pat_id, gender, age, bmi, class, ace_inhibitor, arb, diabetes, hypertension, copd, 
           hydroxychloroquin, Kaletra, Remdesivir, Plasma, Tocilizumab, Camostat, steroids, 
           n_obs, length_of_hospital_stay, death)

table_1b <- tableby(
  class ~ age + gender + bmi + ace_inhibitor + arb + diabetes + hypertension + copd + 
    hydroxychloroquin + Kaletra + Remdesivir + Plasma + Tocilizumab + Camostat + steroids +
    n_obs + length_of_hospital_stay + death,
  data = table_1b_df,
  numeric.stats=c('mean', 'sd', "median","q1q3", 'range')
) 

write2word(table_1a %>% summary(text = TRUE), paste0(getwd(),'/tables/table_1a.docx'))
write2word(table_1b %>% summary(text = TRUE), paste0(getwd(),'/tables/table_1b.docx'))


test_df_healty <- table_1b_df %>% 
    mutate(type = 'COVID19') %>%
    bind_rows(
      controls %>% filter(type == 'healthy')
    )

test_df_influenza <- table_1b_df %>% 
  mutate(type = 'COVID19') %>%
  bind_rows(
    controls %>% filter(type == 'influenza')
  )


table1_tests <- tribble(
  ~term, ~`COVID19 vs. healthy / p.value`, ~`COVID19 vs. influenza / p.value`,
  'age', t.test(age ~ type, data = test_df_healty)$p.value, t.test(age ~ type, data = test_df_influenza)$p.value,
  'gender', fisher.test(test_df_healty$gender, test_df_healty$type)$p.value, chisq.test(test_df_influenza$gender, test_df_influenza$type)$p.value,
  'ace_inhibitor', NA, fisher.test(test_df_influenza$ace_inhibitor, test_df_influenza$type)$p.value,
  'arb', NA, fisher.test(test_df_influenza$arb, test_df_influenza$type)$p.value,
  'diabetes', NA, chisq.test(test_df_influenza$diabetes, test_df_influenza$type)$p.value,
  'hypertension', NA, chisq.test(test_df_influenza$hypertension, test_df_influenza$type)$p.value,
  'copd', NA, fisher.test(test_df_influenza$copd, test_df_influenza$type)$p.value
) %>%
  mutate_at(
    2:3,
    ~case_when(
      is.na(.) ~ '',
      .<0.001  ~ '<0.001*',
      .<0.01   ~ '<0.01*',
      .<0.05   ~ '<0.05*',
      TRUE     ~ paste(round(.,3))
    )
  ) 

write2word(table1_tests %>% kable(), paste0(getwd(),'/tables/table_1_tests.docx'))
  

# end table1 ------------------------------------------------------------------


# table 2 ---------------------------------------------------------------------
# model output of parismonious models for AngII and ACE2

# ACE2
fit_parsim_ace2 <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class +
                     gender + diabetes + hypertension + (1|pat_id),
                   data = ras_df %>% filter(days_since_test < 22))

ci <- confint(fit_parsim_ace2)[-1:-2, ]


fit_parsim_ace2 %>%
  anova() %>%
  tidy() %>%
  select(term, p.value) %>%
  mutate(
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.001 ~ '<0.001*',
      p.value < 0.01  ~ '<0.01*',
      p.value < 0.05  ~ '<0.05*',
      TRUE            ~ ''
    )
  ) 

tmp <- fit_parsim_ace2 %>% 
  fixef 


ACE2_effects <- bind_cols(
  term = rownames(summary(fit_parsim_ace2)$coefficients),
  summary(fit_parsim_ace2)$coefficients[,c('Estimate', 'Pr(>|t|)')]
) %>%
  rename(
    p_value = `Pr(>|t|)`
  ) %>%
  bind_cols(
    est_low = ci[,1],
    est_high = ci[,2]
  ) %>%
  mutate(
    est_trans = exp(Estimate),
    est_low_trans = exp(est_low),
    est_high_trans = exp(est_high),
    p_value = case_when(
      p_value < 0.001 ~ '<0.001*',
      p_value < 0.01  ~ '<0.01*',
      p_value < 0.05  ~ '<0.05*',
      TRUE            ~ paste(round(p_value,3))
    )
  ) %>%
  select(
    term, Estimate, est_low, est_high, est_trans, est_low_trans, est_high_trans, p_value
  ) %>%
  mutate_if(
    is.numeric, 
    round, 4
  )

write_csv(ACE2_effects, 'tables/table_2A_ACE2_effects.csv')
write2word(ACE2_effects %>%
             kable(), 
           paste0(getwd(),'/tables/table_2A_ACE2_effects.docx'))

# AngII

fit_parsim_angII <- lmer(log(AngII) ~ days_since_test * class + ace_inhibitor + (1|pat_id),
                   data = ras_df %>% filter(days_since_test < 22))

ci <- confint(fit_parsim_angII)[-1:-2, ]

fit_parsim_angII %>%
  anova() %>%
  tidy() %>%
  select(term, p.value) %>%
  mutate(
    p.value = round(p.value, 4),
    significance = case_when(
      p.value < 0.001 ~ '<0.001*',
      p.value < 0.01  ~ '<0.01*',
      p.value < 0.05  ~ '<0.05*',
      TRUE            ~ ''
    )
  ) %>%
  kable()



AngII_effects <- bind_cols(
  term = rownames(summary(fit_parsim_angII)$coefficients),
  summary(fit_parsim_angII)$coefficients[,c('Estimate', 'Pr(>|t|)')]
) %>%
  rename(
    p_value = `Pr(>|t|)`
  ) %>%
  bind_cols(
    est_low = ci[,1],
    est_high = ci[,2]
  ) %>%
  mutate(
    est_trans = exp(Estimate),
    est_low_trans = exp(est_low),
    est_high_trans = exp(est_high),
    p_value = case_when(
      p_value < 0.001 ~ '<0.001*',
      p_value < 0.01  ~ '<0.01*',
      p_value < 0.05  ~ '<0.05*',
      TRUE            ~ paste(round(p_value,3))
    )
  ) %>%
  select(
    term, Estimate, est_low, est_high, est_trans, est_low_trans, est_high_trans, p_value
  ) %>%
  mutate_if(
    is.numeric, 
    round, 4
  )


write_csv(AngII_effects, 'tables/table_2A_AngII_effects.csv')
write2word(AngII_effects %>%
             kable(), 
           paste0(getwd(),'/tables/table_2A_AngII_effects.docx'))

# end table 2 -----------------------------------------------------------------


# table S1 --------------------------------------------------------------------

tmp <- ras_df %>%
  mutate(type = 'SARS-CoV2 - all') %>%
  bind_rows(
    controls
  ) %>%
  group_by(pat_id, type) %>%
  summarise_at(
    vars(
      ACE2,
      ACE2_ang_based,
      AngI,
      AngII, 
      AngIII,
      AngIV,
      Ang_1_5,
      Ang_1_7,
      ACE_S
    ),
    ~mean(.)
  ) %>%
  bind_rows(
    ras_df %>%
      filter(class == 'critical') %>%
      mutate(type = 'SARS-CoV2 - critical') %>%
      group_by(pat_id, type) %>%
      summarise_at(
        vars(
          ACE2,
          ACE2_ang_based,
          AngI,
          AngII, 
          AngIII,
          AngIV,
          Ang_1_5,
          Ang_1_7,
          ACE_S
        ),
        ~mean(.)
      )
  ) %>%
  pivot_longer(
    c(-pat_id, -type),
    names_to = 'angiotensin',
    values_to = 'value'
  )
  
table_S1 <- tmp %>% 
  group_by(angiotensin, type) %>%
  summarise(
    avg = paste0(median(value, na.rm = T) %>% display(), ' (', 
                 paste0(quantile(value, c(.25,.75), na.rm = T) %>% display(), collapse = ','), ')'),
    avg = ifelse(str_detect(avg, 'NA'), '', avg)
  ) %>%
  ungroup %>%
  pivot_wider(
    id_cols = 1,
    names_from = type,
    values_from = avg
  ) %>%
  rename(
    'angiotension' = angiotensin
  ) %>%
  select(1,3,6,7,5,2,4)

write_csv(table_S1, 'tables/table_S1.csv')


# end table S1 ----------------------------------------------------------------


# table S2 --------------------------------------------------------------------

tmp <- ras_df %>%
  group_by(pat_id, class) %>%
  summarise_at(
    vars(
      ACE2,
      ACE2_ang_based,
      AngI,
      AngII, 
      AngIII,
      AngIV,
      Ang_1_5,
      Ang_1_7,
      ACE_S
    ),
  ~mean(.)
  ) %>%
  pivot_longer(
    c(-pat_id, -class),
    names_to = 'angiotensin',
    values_to = 'value'
  )

table_S2 <- tmp %>% 
  group_by(angiotensin, class) %>%
  summarise(
    avg = paste0(median(value, na.rm = T) %>% display(), ' (', 
                 paste0(quantile(value, c(.25,.75), na.rm = T) %>% display(), collapse = ','), ')')
  ) %>%
  ungroup %>%
  pivot_wider(
    id_cols = 1,
    names_from = class,
    values_from = avg
  ) %>% 
  left_join(
    tmp %>%
      split(.$angiotensin) %>%
      map_df(
        ~tibble(angiotensin = unique(.$angiotensin),
                p_value = print_p(kruskal.test(value ~ class, .)$p.value)
        )
      )
  )  %>%
  rename(
    'angiotension, Mean / Median (Q1,Q3)' = angiotensin
  )

write_csv(table_S2, 'tables/table_S2.csv')

# end table S2 ----------------------------------------------------------------


# table S3 --------------------------------------------------------------------

length_interval <- 5

ras_time_int_df <- ras_df %>%
  filter(days_since_test %>% between(1,21)) %>% 
  mutate(
    time_int = cut(days_since_test, seq(1,21,length_interval), include.lowest = T),
    class = case_when(
      class == 'mild'     ~ 'asymptomatic',
      class == 'critical' ~ 'critical',
      TRUE                ~ 'mild+severe'
    ),
    class = factor(class, levels = c('asymptomatic', 'mild+severe', 'critical')),
    x_fct = case_when(
      class == 'asymptomatic' & time_int == '[1,6]'   ~ 1, 
      class == 'mild+severe'  & time_int == '[1,6]'   ~ 2, 
      class == 'critical'     & time_int == '[1,6]'   ~ 3,
      #
      class == 'asymptomatic' & time_int == '(6,11]'  ~ 5, 
      class == 'mild+severe'  & time_int == '(6,11]'  ~ 6, 
      class == 'critical'     & time_int == '(6,11]'  ~ 7,
      #
      class == 'asymptomatic' & time_int == '(11,16]' ~ 9, 
      class == 'mild+severe'  & time_int == '(11,16]' ~ 10, 
      class == 'critical'     & time_int == '(11,16]' ~ 11, 
      #
      class == 'asymptomatic' & time_int == '(16,21]' ~ 13, 
      class == 'mild+severe'  & time_int == '(16,21]' ~ 14, 
      class == 'critical'     & time_int == '(16,21]' ~ 15, 
    )
  ) %>% 
  group_by(class, time_int, pat_id, x_fct) %>%
  summarise_at(
    vars(ACE2_ang_based, ACE2, ACE_S, AngI, AngII, AngIII, AngIV, Ang_1_5, Ang_1_7),
    ~ mean(., n.rm = TRUE)
  ) 

interval_tests <- function(marker){
  tibble(
    marker = rep(marker,4),
    ras_time_int_df %>%
      rename_('marker' = marker) %>%
      filter(class != 'asymptomatic') %>%
      split(.$time_int) %>%
      map_df(
        ~ cbind(time_int = unique(.$time_int), tidy(wilcox.test(marker ~ class, .)))[,c('time_int', 'p.value')]
      )
  ) %>%
    mutate(p.value = print_p(p.value)) %>%
    rename(
      angiotensin = marker
    ) %>%
    left_join(
      ras_time_int_df %>%
        rename_('marker' = marker) %>%
        filter(class != 'asymptomatic') %>%
        group_by(time_int, class) %>%
        summarise(
          text = paste0(display(mean(marker, na.rm=T), 2), ' (', 
                        display(quantile(marker, .25, na.rm=T), 2), ',',
                        display(quantile(marker, .75, na.rm=T), 2), ')'
                 )
        ) %>%
        pivot_wider(
          names_from = 'class',
          values_from = 'text'
        )
    ) %>%
    select(1,2,4,5,3)
}

table_S3 <- bind_rows(
  interval_tests('ACE_S'),
  interval_tests('ACE2'),
  interval_tests('ACE2_ang_based'),
  interval_tests('Ang_1_5'),
  interval_tests('Ang_1_7'),
  interval_tests('AngI'),
  interval_tests('AngII'),
  interval_tests('AngIII'),
  interval_tests('AngIV'),
)

write_csv(table_S3, 'tables/table_S3.csv')

# end table S3 ----------------------------------------------------------------


# table S4 --------------------------------------------------------------------
# likelihood ratio tests from minimal model to full model

# ACE2

model_data <- ras_df %>% filter(days_since_test %>% between(0,21),
                                !pat_id %in% c('SARS-FP-112', 'SARS-FP-123'), 
                                !is.na(steroids))

fit_null    <- lmer(log(ACE2) ~ days_since_test + (1|pat_id),
                    data = model_data)

fit_minimal <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + (1|pat_id),
                    data = model_data)

fit_age <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                  + age + (1|pat_id),
                data = model_data)

fit_gender <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                     + age + gender + (1|pat_id),
                   data = model_data)

fit_ace_inh <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                      + age + gender + ace_inhibitor + (1|pat_id),
                    data = model_data)

fit_arb <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                  + age + gender + ace_inhibitor + arb + (1|pat_id),
                data = model_data)

fit_steroids <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                       + age + gender + ace_inhibitor + arb + steroids + (1|pat_id),
                     data = model_data)

fit_diabetes <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                       + age + gender + ace_inhibitor + arb + steroids + diabetes + (1|pat_id),
                     data = model_data)

fit_hypertension <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + 
                           + age + gender + ace_inhibitor + arb + steroids + diabetes + hypertension + (1|pat_id),
                         data = model_data)

fit_copd <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class +
                   + age + gender + ace_inhibitor + arb + steroids + diabetes + hypertension + copd + (1|pat_id),
                 data = model_data)


lrt_ace2 <- tibble(Angiotensin = 'ACE2') %>%
  bind_cols(
    anova(
      fit_null,
      fit_minimal,
      fit_age,
      fit_gender,
      fit_ace_inh,
      fit_arb,
      fit_steroids,
      fit_diabetes,
      fit_hypertension,
      fit_copd          #fit_maximal
    ) %>%
      tidy %>%
      select(term, statistic, df, p.value) %>%
      rename(model = term,
             Chisq = statistic) %>%
      mutate(
        p.value = print_p(p.value)
      ) %>%
      bind_cols(
        formula = c( 
          paste(summary(fit_null)$call[2]),
          paste(summary(fit_minimal)$call[2]),
          paste(summary(fit_age)$call[2]),
          paste(summary(fit_gender)$call[2]),
          paste(summary(fit_ace_inh)$call[2]),
          paste(summary(fit_arb)$call[2]),
          paste(summary(fit_steroids)$call[2]),
          paste(summary(fit_diabetes)$call[2]),
          paste(summary(fit_hypertension)$call[2]),
          paste(summary(fit_copd)$call[2])
        )
      )
  )

# AngII

fit_null    <- lmer(log(AngII) ~ days_since_test + (1|pat_id),
                    data = model_data)

fit_minimal <- lmer(log(AngII) ~ days_since_test * class + (1|pat_id),
                    data = model_data)

fit_age <- lmer(log(AngII) ~ days_since_test * class +
                  age + (1|pat_id),
                data = model_data)

fit_gender <- lmer(log(AngII) ~ days_since_test * class +
                     age + gender + (1|pat_id),
                   data = model_data)

fit_ace_inh <- lmer(log(AngII) ~ days_since_test * class +
                      age + gender + ace_inhibitor + (1|pat_id),
                    data = model_data)

fit_arb <- lmer(log(AngII) ~ days_since_test * class +
                  age + gender + ace_inhibitor + arb + (1|pat_id),
                data = model_data)

fit_steroids <- lmer(log(AngII) ~ days_since_test * class +
                       age + gender + ace_inhibitor + arb + steroids + (1|pat_id),
                     data = model_data)

fit_diabetes <- lmer(log(AngII) ~ days_since_test * class + 
                       age + gender + ace_inhibitor + arb + steroids + diabetes + (1|pat_id),
                     data = model_data)

fit_hypertension <- lmer(log(AngII) ~ days_since_test * class +
                           age + gender + ace_inhibitor + arb + steroids + diabetes + hypertension + (1|pat_id),
                         data = model_data)

fit_copd <- lmer(log(AngII) ~ days_since_test * class + 
                   age + gender + ace_inhibitor + arb + steroids + diabetes + hypertension + copd + (1|pat_id),
                 data = model_data)

lrt_angII <- tibble(Angiotensin = 'AngII') %>%
  bind_cols(
    anova(
      fit_null,
      fit_minimal,
      fit_age,
      fit_gender,
      fit_ace_inh,
      fit_arb,
      fit_steroids,
      fit_diabetes,
      fit_hypertension,
      fit_copd          #fit_maximal
    ) %>%
      tidy %>%
      select(term, statistic, df, p.value) %>%
      rename(model = term,
             Chisq = statistic) %>%
      mutate(
        p.value = print_p(p.value)
      ) %>%
      bind_cols(
        formula = c( 
          paste(summary(fit_null)$call[2]),
          paste(summary(fit_minimal)$call[2]),
          paste(summary(fit_age)$call[2]),
          paste(summary(fit_gender)$call[2]),
          paste(summary(fit_ace_inh)$call[2]),
          paste(summary(fit_arb)$call[2]),
          paste(summary(fit_steroids)$call[2]),
          paste(summary(fit_diabetes)$call[2]),
          paste(summary(fit_hypertension)$call[2]),
          paste(summary(fit_copd)$call[2])
        )
      )
  )

bind_rows(
  lrt_ace2,
  lrt_angII
) %>% write_csv('tables/table_s4_LRT.csv')

# end table s4 ----------------------------------------------------------------

#################
#### FIGURES ####
#################

# figure 2B -------------------------------------------------------------------
# smooth curves of Ace2, AngII, Ang1-7, ACE-S over time, truncated after 21 days after first test

  
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


# end figure 2B ---------------------------------------------------------------


# figure 2C -------------------------------------------------------------------

fig_w <- 6
fig_h <- 5

fig2c_df <- ras_df %>%
  select(
    pat_id,
    AngII,
    ACE2_ang_based,
    ACE2,
    Ang_1_7,
    ACE_S,
    phase,
    class
  ) %>%
  bind_rows(
    controls %>%
      mutate(phase = ifelse(days < 6, 'early', 'late')) %>%
      select(
        pat_id,
        AngII,
        ACE2_ang_based,
        ACE2,
        Ang_1_7,
        ACE_S,
        phase,
        class = type
      )
  ) %>%
  group_by(class, phase, pat_id) %>%
  summarise_all(
    ~mean(.,na.rm = T)
  ) %>%
  ungroup %>%
  select(-pat_id) %>%
  mutate(
    grp = ifelse(class %in% c( 'asymptomatic', 'mild', 'severe', 'critical'), 'SARS-CoV2', 'Control'),
    class = case_when(
      class == 'healthy'    ~ 'healthy',
      class == 'asymptomatic' ~ 'asymptomatic',
      class == 'mild'       ~ 'mild',
      class == 'severe'     ~ 'severe',
      class == 'critical' & phase == 'early'  ~ 'critical, early',
      class == 'critical' & phase == 'late'  ~ 'critical, late',
      class == 'influenza' & phase == 'early' ~ 'influenza, early',
      class == 'influenza' & phase == 'late' ~ 'influenza, late',
      TRUE ~ 'discard'
    )
  ) %>% 
  filter(class != 'discard') %>%
  mutate(
    class = factor(class, levels = c('healthy', 'asymptomatic', 'mild', 'severe', 'critical, early',
                                     'critical, late', 'influenza, early', 'influenza, late'))
  ) %>%
  replace_na(list(phase = 'n.a.'))

fig2c_plot <- function(data, variable, ylab = NULL){
  p <- data %>%
    ggplot(aes_string(x = 'class', y = variable, fill = 'grp')) +
    geom_beeswarm(cex = .5, alpha = .2) +
    geom_boxplot(width = .5, position=position_dodge(.7), alpha = .5) +
    scale_fill_ptol() +
    scale_y_log10() +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = .5, hjust = .5)
    ) +
    labs(x = '', fill = '', y = ifelse(is.null(ylab), variable, ylab))
  p
}

# Ang 1-7
p_c1 <- fig2c_df %>% fig2c_plot('Ang_1_7', 'Angiotensin 1-7')

p1 <- wilcox.test(Ang_1_7 ~ grp, 
       data = fig2c_df %>% 
         filter(grp == 'SARS-CoV2' | class == 'healthy'))$p.value                              #healthy vs SARS-CoV2
p2 <- wilcox.test(Ang_1_7 ~ class, 
       data = fig2c_df %>% 
         filter(class %in% c('severe', 'critical, early', 'critical, late')) %>%
         mutate(class = ifelse(class == 'severe', 'severe', 'critical')))$p.value                #severe vs critical

h1 <- 3000
h2 <- 4000
h3 <- 5000
h4 <- 6000

p_c1 + geom_segment(aes(x = 'asymptomatic', xend = 'critical, late', y = h1, yend = h1)) +
  geom_segment(aes(x = 1, xend = 'severe', y = h2, yend = h2)) +
  geom_segment(aes(x = 1, xend = 1, y = h2, yend = h1)) +
  geom_segment(aes(x = 4, xend = 4, y = h2, yend = h1)) +
  geom_text(aes(x = 2.5, y = h2*1.3, label = print_p(p1)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h3, yend = h3)) +
  geom_segment(aes(x = 'severe', xend = 5.5, y = h4, yend = h4)) +
  geom_segment(aes(x = 'severe', xend = 'severe', y = h4, yend = h3)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h4, yend = h3)) +
  geom_text(aes(x = 4.8, y = h4*1.3, label = print_p(p2)), size = 3) 

ggsave('plots/figure2C_Ang_1_7.png', width = fig_w, height = fig_h)
ggsave('plots/figure2C_Ang_1_7_clean.png', p_c1, width = fig_w, height = fig_h)


# ACE2_ang_based
p_c2 <- fig2c_df %>% fig2c_plot('ACE2_ang_based', 'ACE2 activity')

p1 <- wilcox.test(ACE2_ang_based ~ grp, 
       data = fig2c_df %>% 
         filter(grp == 'SARS-CoV2' | class == 'healthy'))$p.value                               #healthy vs SARS-CoV2
p2 <- wilcox.test(ACE2_ang_based ~ class, 
       data = fig2c_df %>% filter(class %in% c('critical, early', 'critical, late')))$p.value     #critical, early vs late
p3 <- wilcox.test(ACE2_ang_based ~ class, 
       data = fig2c_df %>% 
         filter(str_detect(class, 'critical') | str_detect(class, 'influenza')) %>%
         mutate(class = ifelse(str_detect(class, 'critical'), 'critical', 'influenza')))$p.value   #critical vs influenza

h1 <- 1.2
h2 <- 1.4
h3 <- 1.65
h4 <- 2.5
h5 <- 2.8

p_c2 + geom_segment(aes(x = 'asymptomatic', xend = 'critical, late', y = h1, yend = h1)) +
  geom_segment(aes(x = 1, xend = 'severe', y = h2, yend = h2)) +
  geom_segment(aes(x = 1, xend = 1, y = h2, yend = h1)) +
  geom_segment(aes(x = 4, xend = 4, y = h2, yend = h1)) +
  geom_text(aes(x = 2.5, y = h2*1.2, label = print_p(p1)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h3, yend = h3)) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, early', y = h3*.9, yend = h3)) +
  geom_segment(aes(x = 'critical, late', xend = 'critical, late', y = h3*.9, yend = h3)) +
  geom_text(aes(x = 5.5, y = h3*1.2, label = print_p(p2)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h4, yend = h4)) +
  geom_segment(aes(x = 'influenza, early', xend = 'influenza, late', y = h4, yend = h4)) +
  geom_segment(aes(x = 5.5, xend = 7.5, y = h5, yend = h5)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h4, yend = h5)) +
  geom_segment(aes(x = 7.5, xend = 7.5, y = h4, yend = h5)) +
  geom_text(aes(x = 6.5, y = h5*1.2, label = print_p(p3)), size = 3)

ggsave('plots/figure2C_ACE2_ang_based.png', width = fig_w, height = fig_h)
ggsave('plots/figure2C_ACE2_ang_based_clean.png', p_c2, width = fig_w, height = fig_h)
  
# ACE2
p_c2_alt <- fig2c_df %>% fig2c_plot('ACE2', 'ACE2 concentration')

p1 <- wilcox.test(ACE2 ~ grp, 
             data = fig2c_df %>% 
               filter(grp == 'SARS-CoV2' | class == 'healthy'))$p.value                               #healthy vs SARS-CoV2
p2 <- wilcox.test(ACE2 ~ class, 
             data = fig2c_df %>% filter(class %in% c('critical, early', 'critical, late')))$p.value     #critical, early vs late
p3 <- wilcox.test(ACE2 ~ class, 
             data = fig2c_df %>% 
               filter(str_detect(class, 'critical') | str_detect(class, 'influenza')) %>%
               mutate(class = ifelse(str_detect(class, 'critical'), 'critical', 'influenza')))$p.value   #critical vs influenza

h1 <- 60
h2 <- 70
h3 <- 90
h4 <- 140
h5 <- 160

p_c2_alt + geom_segment(aes(x = 'asymptomatic', xend = 'critical, late', y = h1, yend = h1)) +
  geom_segment(aes(x = 1, xend = 'severe', y = h2, yend = h2)) +
  geom_segment(aes(x = 1, xend = 1, y = h2, yend = h1)) +
  geom_segment(aes(x = 4, xend = 4, y = h2, yend = h1)) +
  geom_text(aes(x = 2.5, y = h2*1.2, label = print_p(p1)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h3, yend = h3)) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, early', y = h3*.9, yend = h3)) +
  geom_segment(aes(x = 'critical, late', xend = 'critical, late', y = h3*.9, yend = h3)) +
  geom_text(aes(x = 5.5, y = h3*1.2, label = print_p(p2)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h4, yend = h4)) +
  geom_segment(aes(x = 'influenza, early', xend = 'influenza, late', y = h4, yend = h4)) +
  geom_segment(aes(x = 5.5, xend = 7.5, y = h5, yend = h5)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h4, yend = h5)) +
  geom_segment(aes(x = 7.5, xend = 7.5, y = h4, yend = h5)) +
  geom_text(aes(x = 6.5, y = h5*1.21, label = print_p(p3)), size = 3)

ggsave('plots/figure2C_ACE2_measured.png', width = fig_w, height = fig_h)
ggsave('plots/figure2C_ACE2_measured_clean.png', p_c2_alt, width = fig_w, height = fig_h)

# ACE-S
p_c3 <- fig2c_df %>% fig2c_plot('ACE_S', 'ACE activty')

p1 <- wilcox.test(ACE_S ~ class, 
       data = fig2c_df %>% filter(class %in% c('healthy', 'asymptomatic')))$p.value              #healthy vs asymptomatic
p2 <- wilcox.test(ACE_S ~ class, 
       data = fig2c_df %>% filter(class %in% c('asymptomatic', 'mild')))$p.value                 #asymptomatic vs mild
p3 <- wilcox.test(ACE_S ~ class, 
       data = fig2c_df %>% filter(class %in% c('mild', 'severe')))$p.value                       #mild vs severe
p4 <- wilcox.test(ACE_S ~ class, 
       data = fig2c_df %>% 
         filter(class %in% c('severe', 'critical, early', 'critical, late')) %>%
         mutate(class = ifelse(class == 'severe', 'severe', 'critical')))$p.value               #severe vs critical
p5 <- wilcox.test(ACE_S ~ class, 
       data = fig2c_df %>% 
         filter(str_detect(class, 'critical') | str_detect(class, 'influenza')) %>%
         mutate(class = ifelse(str_detect(class, 'critical'), 'critical', 'influenza')))$p.value #critical vs influenza

h1 <- 20
h2 <- 15
h3 <- 12
h4 <- 9
h5 <- 8
h6 <- 14
h7 <- 16

p_c3 + geom_segment(aes(x = 'healthy', xend = 'asymptomatic', y = h1, yend = h1)) +
  geom_segment(aes(x = 'healthy', xend = 'healthy', y = h1, yend = h1*.9)) +
  geom_segment(aes(x = 'asymptomatic', xend = 'asymptomatic', y = h1, yend = h1*.9)) +
  geom_text(aes(x = 1.5, y = h1*1.2, label = print_p(p1)), size = 3) +
  geom_segment(aes(x = 'asymptomatic', xend = 'mild', y = h2, yend = h2)) +
  geom_segment(aes(x = 'asymptomatic', xend = 'asymptomatic', y = h2, yend = h2*.9)) +
  geom_segment(aes(x = 'mild', xend = 'mild', y = h2, yend = h2*.9)) +
  geom_text(aes(x = 2.5, y = h2*1.2, label = print_p(p2)), size = 3) +
  geom_segment(aes(x = 'mild', xend = 'severe', y = h3, yend = h3)) +
  geom_segment(aes(x = 'severe', xend = 'severe', y = h3, yend = h3*.9)) +
  geom_segment(aes(x = 'mild', xend = 'mild', y = h3, yend = h3*.9)) +
  geom_text(aes(x = 3.5, y = h3*1.2, label = print_p(p3)), size = 3) +
  geom_segment(aes(x = 'severe' , xend = 5.5, y = h4, yend = h4)) +
  geom_segment(aes(x = 'critical, early' , xend = 'critical, late', y = h5, yend = h5)) +
  geom_segment(aes(x = 'severe', xend = 'severe', y = h4, yend = h5)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h4, yend = h5)) +
  geom_text(aes(x = 4.8, y = h4*1.2, label = print_p(p4)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h6, yend = h6)) +
  geom_segment(aes(x = 'influenza, early', xend = 'influenza, late', y = h6, yend = h6)) +
  geom_segment(aes(x = 5.5, xend = 7.5, y = h7, yend = h7)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h7, yend = h6)) +
  geom_segment(aes(x = 7.5, xend = 7.5, y = h7, yend = h6)) +
  geom_text(aes(x = 6.5, y = h7*1.2, label = print_p(p5)), size = 3)

ggsave('plots/figure2C_ACE_S.png', width = fig_w, height = fig_h)
ggsave('plots/figure2C_ACE_S_clean.png', p_c3, width = fig_w, height = fig_h)

# AngII
p_c4 <- fig2c_df %>% fig2c_plot('AngII', 'Angiotensin II')

p1 <- wilcox.test(AngII ~ grp, data = fig2c_df %>% filter(grp == 'SARS-CoV2' | class == 'healthy'))$p.value #healthy vs SARS-CoV2
p2 <- wilcox.test(AngII ~ class, 
       data = fig2c_df %>% 
         filter(str_detect(class, 'critical') | class == 'healthy') %>%
         mutate(class = ifelse(class == 'healthy', 'healthy', 'critical')))$p.value              #healthy vs critical
p3 <- wilcox.test(AngII ~ class, 
       data = fig2c_df %>% 
         filter(str_detect(class, 'critical') | str_detect(class, 'influenza')) %>%
         mutate(class = ifelse(str_detect(class, 'critical'), 'critical', 'influenza')))$p.value #critical vs influenza

h1 <- 8000
h2 <- 10000
h3 <- 13000
h4 <- 16000
h5 <- 19000
h6 <- 23000

p_c4 + geom_segment(aes(x = 'asymptomatic', xend = 'critical, late', y = h1, yend = h1)) +
  geom_segment(aes(x = 1, xend = 'severe', y = h2, yend = h2)) +
  geom_segment(aes(x = 1, xend = 1, y = h2, yend = h1)) +
  geom_segment(aes(x = 4, xend = 4, y = h2, yend = h1)) +
  geom_text(aes(x = 2.5, y = h2*1.3, label = print_p(p1)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h3, yend = h3)) +
  geom_segment(aes(x = 1, xend = 5.5, y = h4, yend = h4)) +
  geom_segment(aes(x = 1, xend = 1, y = h3, yend = h4)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h3, yend = h4)) +
  geom_text(aes(x = 3.3, y = h4*1.3, label = print_p(p2)), size = 3) +
  geom_segment(aes(x = 'critical, early', xend = 'critical, late', y = h5, yend = h5)) +
  geom_segment(aes(x = 'influenza, early', xend = 'influenza, late', y = h5, yend = h5)) +
  geom_segment(aes(x = 5.5, xend = 7.5, y = h6, yend = h6)) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = h5, yend = h6)) +
  geom_segment(aes(x = 7.5, xend = 7.5, y = h5, yend = h6)) +
  geom_text(aes(x = 6.5, y = h6*1.3, label = print_p(p3)), size = 3)
  
ggsave('plots/figure2C_AngII.png', width = fig_w, height = fig_h)
ggsave('plots/figure2C_AngII_clean.png', p_c4, width = fig_w, height = fig_h)

# end figure 2C ---------------------------------------------------------------


# figure 3 --------------------------------------------------------------------

fig3_plot <- function(data, marker, ylab = NULL){
  
  data %>%
    filter(
      days_since_test %>% between(0,21),
      class == 'critical'
    ) %>%
    mutate(
      death = ifelse(death, 'deceased', 'recovered'),
      death = factor(death, levels = c('recovered', 'deceased'))
    ) %>%
    ggplot(aes_string(x = 'days_since_test', y = marker, color = 'death', group = 'pat_id')) +
    geom_line(aes(linetype = death), alpha = .3, size = .8) + 
    geom_smooth(aes(group = as.factor(death)), alpha= .2) +
    scale_y_log10() +
    theme_pubr() +
    scale_color_ptol() +
    labs(
      color = '',
      linetype = '',
      y = ifelse(is.null(ylab), marker, ylab),
      x = 'days since first test'
    )
  
}

ggarrange(
  ras_df %>% fig3_plot('Ang_1_7', 'Angiotensin 1-7'),
  ras_df %>% fig3_plot('ACE2', 'ACE2 concentration'),
  ras_df %>% fig3_plot('ACE_S', 'ACE activity'),
  ras_df %>% fig3_plot('AngII', 'Angiotensin II'),
  common.legend = TRUE,
  align = 'hv'
  )

ggsave('plots/figure3_alternative.png')



# end figure 3 ----------------------------------------------------------------


# figure s1 -------------------------------------------------------------------
#correlation ACE2 and ACE2 ang based 

#v1

r_sq <- summary(lm(log(ACE2) ~ log(ACE2_ang_based), 
                   data = ras_df))$r.squared %>% round(3)

ras_df %>%
  ggplot(aes(x = ACE2, y = ACE2_ang_based)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'lm', se = FALSE, color = 'darkgreen', size = 1.5, alpha = .6) +
  geom_text(aes(x = 1, y = .5, 
                label = paste0('R^2 == ', r_sq)), parse = TRUE)+
  theme_pubr() +
  labs(
    y = 'ACE2 activity',
    x = 'ACE2 concentration'
  )

ggsave('plots/figure_S1.png', width = 5, height = 5)

#v2

r_sq <- summary(lm(log(ACE2) ~ log(ACE2_activity_v2), 
                   data = ras_df))$r.squared %>% round(3)
ras_df %>%
  ggplot(aes(x = ACE2, y = ACE2_activity_v2)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'lm', se = FALSE, color = 'darkgreen', size = 1.5, alpha = .6) +
  geom_text(aes(x = .5, y = .5, 
                label = paste0('R^2 == ', r_sq)), parse = TRUE)+
  theme_pubr() +
  labs(
    y = 'ACE2 activity',
    x = 'ACE2 concentration'
  )

ggsave('plots/figure_S1_aternative 2 Ang 1_5 over AngII.png', width = 5, height = 5)

#v3

r_sq <- summary(lm(log(ACE2) ~ log(ACE2_activity_v3), 
                   data = ras_df))$r.squared %>% round(3)
ras_df %>%
  ggplot(aes(x = ACE2, y = ACE2_activity_v3)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = 'lm', se = FALSE, color = 'darkgreen', size = 1.5, alpha = .6) +
  geom_text(aes(x = .5, y = .5, 
                label = paste0('R^2 == ', r_sq)), parse = TRUE)+
  theme_pubr() +
  labs(
    y = 'ACE2 activity',
    x = 'ACE2 concentration'
  )

ggsave('plots/figure_S1_aternative 3 Ang 1_5 plus Ang 1_7 over sum of AngI AngII Ang 1_5 Ang 1_7.png', width = 5, height = 5)

# end figure S1 ---------------------------------------------------------------


# figure S2 -------------------------------------------------------------------

ras_time_int_df <- ras_df %>%
  filter(days_since_test %>% between(1,21)) %>% 
  mutate(
    time_int = cut(days_since_test, seq(1,21,length_interval), include.lowest = T),
    class = case_when(
      class == 'mild'     ~ 'asymptomatic',
      class == 'critical' ~ 'critical',
      TRUE                ~ 'mild+severe'
    ),
    class = factor(class, levels = c('asymptomatic', 'mild+severe', 'critical')),
    x_fct = case_when(
      class == 'asymptomatic' & time_int == '[1,6]'   ~ 1, 
      class == 'mild+severe'  & time_int == '[1,6]'   ~ 2, 
      class == 'critical'     & time_int == '[1,6]'   ~ 3,
      #
      class == 'asymptomatic' & time_int == '(6,11]'  ~ 5, 
      class == 'mild+severe'  & time_int == '(6,11]'  ~ 6, 
      class == 'critical'     & time_int == '(6,11]'  ~ 7,
      #
      class == 'asymptomatic' & time_int == '(11,16]' ~ 9, 
      class == 'mild+severe'  & time_int == '(11,16]' ~ 10, 
      class == 'critical'     & time_int == '(11,16]' ~ 11, 
      #
      class == 'asymptomatic' & time_int == '(16,21]' ~ 13, 
      class == 'mild+severe'  & time_int == '(16,21]' ~ 14, 
      class == 'critical'     & time_int == '(16,21]' ~ 15, 
    )
  ) %>% 
  group_by(class, time_int, pat_id, x_fct) %>%
  summarise_at(
    vars(ACE2, ACE_S, AngII, Ang_1_7),
    ~ mean(log(.), n.rm = TRUE)
  ) 


fig_S2_plot <- function(data, variable, ylab){

  data %>% 
    mutate_at(
      variable,
      exp
    ) %>%
    mutate(
      grp = paste(class, x_fct)
    ) %>%
    ggplot(aes_string(x = 'x_fct', y = variable, fill = 'class', group = 'grp')) +
    geom_beeswarm(aes(color = class), alpha = .5, size = 1, position = position_dodge()) +
    geom_boxplot(width = .7, alpha = .5) +
    theme_pubr() +
    scale_color_ptol() +
    scale_fill_ptol() +
    scale_y_log10() +
    labs(
      x = '5-day intervals since test',
      y = ylab,
      fill = '',
      color = ''
    ) +
    scale_x_continuous(breaks = c(2,6,10,14), labels = c('1-6', '7-11', '12-16', '17-21'))
  
}

ggarrange(
  ras_time_int_df %>% fig_S2_plot('Ang_1_7', 'Angiotensin 1-7'),
  ras_time_int_df %>% fig_S2_plot('ACE2', 'ACE2 concentration'),
  ras_time_int_df %>% fig_S2_plot('ACE_S', 'ACE activity'),
  ras_time_int_df %>% fig_S2_plot('AngII', 'Angiotensin II'),
  common.legend = TRUE, 
  align = 'hv',
  ncol = 2, nrow = 2
)

ggsave('plots/figure_S2.png', width = 9, height = 6)

# end figure S2 ---------------------------------------------------------------

# figure S3 trajectories by comedication --------------------------------------

figS3_plot <- function(data, marker, filter_arg, ylab = NULL){
  
  if(marker == 'ACE2'){
    limits <- c(.1,150)
    } else {
      limits <- c(1,7000)
    }
  
  data %>%
    filter(
      days_since_test %>% between(0,21),
      ace_inhibitor == filter_arg[1],
      arb == filter_arg[2]
    ) %>%
    ggplot(aes_string(x = 'days_since_test', y = marker, color = 'class', group = 'pat_id')) +
    geom_line(aes(linetype = class), alpha = .3, size = .8) + 
    geom_smooth(aes(group = as.factor(class)), alpha= .2) +
    scale_y_log10() +
    theme_pubr() +
    scale_color_ptol() +
    labs(
      color = '',
      linetype = '',
      y = ifelse(is.null(ylab), marker, ylab),
      x = 'days since first test',
      title = paste0('ACEi: ', filter_arg[1], ' / ARB: ', filter_arg[2])
    ) +
    coord_cartesian(ylim = limits)
  
}

ggarrange(
  ras_df %>% figS3_plot('ACE2', c('no', 'no'), 'ACE2 concentration'),
  ras_df %>% figS3_plot('ACE2', c('yes', 'no'), 'ACE2 Concentration') + labs(y = ''),
  ras_df %>% figS3_plot('ACE2', c('no', 'yes'), 'ACE2 Concentration') + labs(y = ''),
  ras_df %>% figS3_plot('AngII', c('no', 'no'), 'Angiotensin II'),
  ras_df %>% figS3_plot('AngII', c('yes', 'no')) + labs(y = ''),
  ras_df %>% figS3_plot('AngII', c('no', 'yes')) + labs(y = ''),
  common.legend = TRUE,
  align = 'hv',
  ncol = 3,
  nrow = 2
)

ggsave('plots/figure_S3.png', height = 8, width = 12)

# end figure SX trajectories by comedication ----------------------------------


# figure S4 minimal mixed model plot ------------------------------------------

# ACE2

fit_minimal <- lmer(log(ACE2) ~ days_since_test * class + I(days_since_test^2) * class + (1|pat_id),
                    data = ras_df %>% filter(days_since_test < 22))

new_dat <- expand.grid(
  days_since_test = 0:21,
  class = c('asymptomatic', 'mild', 'severe', 'critical')
)
new_dat$ACE2 <- exp(predict(fit_minimal, new_dat, re.form=NA))

mm <- model.matrix(terms(fit_minimal), new_dat)
predFun<-function(.) mm %*% fixef(.)
n_bootstraps <- 5000
bb<-bootMer(fit_minimal,FUN=predFun,nsim=n_bootstraps) #do this 500 times
bb_se<-apply(bb$t,2,function(x) quantile(x, c(.025, .975)))
#bb_se<-apply(bb$t,2,function(x)  x[order(x)][c(5,195)])
new_dat$LC<-bb_se[1,]
new_dat$UC<-bb_se[2,] 
new_dat$pred<-predict(fit_minimal,newdata=new_dat,re.form=NA)

p_ace2 <- new_dat %>%
  ggplot(aes(x = days_since_test, y = pred, group = class, color = class)) +
  geom_line(data = ras_df %>% filter(gender == 'm', days_since_test < 22),
            aes(x = days_since_test, y = log(ACE2), group = pat_id, linetype = class),
            alpha = .3,
            size = 1) +
  geom_ribbon(aes(ymin = LC, ymax = UC), alpha = .1) +
  geom_line(size = 1.5) +
  scale_color_ptol() +
  theme_pubr() +
  labs(
    x = 'days since test',
    y = 'log(ACE2 concentration)',
    shape = 'controls',
    # title = 'Regression Curves (i.e. fixed effects) of minimal model',
    # subtitle = 'with 95% bootstrap CI',
    color = '',
    linetype= ''
  ) +
  theme(legend.position = 'right')


# AngII

fit_minimal <- lmer(log(AngII) ~ days_since_test * class + (1|pat_id),
                    data = ras_df %>% filter(days_since_test < 22))

new_dat <-expand.grid(
  days_since_test = 0:21,
  class = c('asymptomatic', 'mild', 'severe', 'critical')
)
new_dat$AngII <- exp(predict(fit_minimal, new_dat, re.form=NA))

mm <- model.matrix(terms(fit_minimal), new_dat)
predFun<-function(.) mm %*% fixef(.)
n_bootstraps <- 500
bb<-bootMer(fit_minimal,FUN=predFun,nsim=n_bootstraps) #do this 200 times
bb_se<-apply(bb$t,2,function(x) quantile(x, c(.025, .975)))
#bb_se<-apply(bb$t,2,function(x)  x[order(x)][c(5,195)])
new_dat$LC<-bb_se[1,]
new_dat$UC<-bb_se[2,] 
new_dat$pred<-predict(fit_minimal,newdata=new_dat,re.form=NA)

p_angII <- new_dat %>%
  ggplot(aes(x = days_since_test, y = pred, group = class, color = class)) +
  geom_line(data = ras_df %>% filter(gender == 'm', days_since_test < 22),
            aes(x = days_since_test, y = log(AngII), group = pat_id, linetype = class),
            alpha = .3,
            size = 1) +
  geom_ribbon(aes(ymin = LC, ymax = UC), alpha = .1) +
  geom_line(size = 1.5) +
  scale_color_ptol() +
  theme_pubr() +
  labs(
    x = 'days since test',
    y = 'log(Angiotensin II)',
    shape = 'controls',
    # title = 'Regression Curves (i.e. fixed effects) of minimal model',
    # subtitle = 'with 95% bootstrap CI',
    color = '',
    linetype= ''
  ) +
  theme(legend.position = 'right')

ggarrange(
  p_ace2,
  p_angII,
  align = 'v',
  common.legend = TRUE,
  ncol = 1
)

ggsave('plots/figure_S4.png', width = 5, height = 8)

# end figure sx minimal mixed model plot --------------------------------------









