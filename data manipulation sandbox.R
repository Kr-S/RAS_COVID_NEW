librarian::shelf(dplyr, plotly, shiny, tidyverse, lubridate ,stringr , readxl, ggthemes)





ACOVACT_lab_all_ID_no_names_2020_09_22 <- read.delim("C:/Users/SKrenn/RAS-COVID/SARS_FP_data_repository/ACOVACT_lab_all_ID_no_names_2020_09_22.txt")


renamevars <- c(potassium = 'Kalium', potassium_urine ='Kalium /U')

laball <- ACOVACT_lab_all_ID_no_names_2020_09_22 %>%
  distinct() %>%
  pivot_wider(names_from = Type, values_from = Lab.values) %>%
  select(
    renamevars,
    pat_id = 'ID',
    date_measurement = 'X.2'
  ) %>%
  pivot_longer(cols=names(renamevars),names_to="names", values_to="values") %>%
  filter(!is.na(values)) %>%
  pivot_wider(names_from = names, values_from = values) %>%
  mutate(
    date_measurement = as.Date(strptime(date_measurement, "%d.%m.%y %H:%M")),
    potassium = as.numeric(as.character(potassium)),
    potassium_urine = as.numeric(as.character(potassium_urine))
  )




ras_dfpot <- measurement_df %>% 
  left_join(lab_values) %>%
  left_join(lab_values_2 %>% rename(il_6 = `Interleukin-6`)) %>%
  left_join(resp_values) %>%
  left_join(patients) %>%
  full_join(laball) %>%
  left_join(steroids_bmi) %>%
  group_by(pat_id) %>%
  mutate(
    days_since_test = ifelse(!is.na(date_first_test_asdate),
                             as.numeric(date_measurement - date_first_test_asdate) ,
                             ifelse(
                               !is.na(date_hospital_asdate),
                               as.numeric(date_measurement - date_hospital_asdate),
                               as.numeric(date_measurement - date_icu_asdate)
                             )
                             
    ),
    phase = ifelse(days_since_test < 6, 'early', 'late'),
    days_test_to_death = ifelse(!is.na(date_first_test_asdate),
                                as.numeric(date_death_asdate - date_first_test_asdate),
                                ifelse(
                                  !is.na(date_hospital_asdate),
                                  as.numeric(date_death_asdate - date_hospital_asdate),
                                  as.numeric(date_death_asdate - date_icu_asdate)
                                )),
    survival_28days = ifelse(
      is.na(date_death_asdate) | days_test_to_death > 28, 'survival', 'death'
    ),
    death = !is.na(date_death_asdate),
    n_obs = n(),
    obs_length = max(days_since_test),
    length_of_hospital_stay = ifelse(any(death),
                                     ifelse(!is.na(date_first_test_asdate),
                                            as.numeric(date_death_asdate - date_first_test_asdate),
                                            ifelse(
                                              !is.na(date_hospital_asdate),
                                              as.numeric(date_death_asdate - date_hospital_asdate),
                                              as.numeric(date_death_asdate - date_icu_asdate)
                                            )),
                                     ifelse(!is.na(date_first_test_asdate),
                                            as.numeric(date_discharge_asdate - date_first_test_asdate),
                                            ifelse(
                                              !is.na(date_hospital_asdate),
                                              as.numeric(date_discharge_asdate - date_hospital_asdate),
                                              as.numeric(date_discharge_asdate - date_icu_asdate)
                                            ))
    )
  ) %>%
  ungroup



