
ras_df_raw <- readRDS('data/RAS_df.rds')

ras_df_cox <- ras_df_raw %>%
      filter(
        !is.na(Aldosterone)
      ) %>%
      select(
        pat_id,
        potassium,
        Aldosterone,
        date_measurement,
        date_hospital_asdate,
        date_hospital_discharge_asdate,
        date_icu_asdate,
        date_icu_discharge_asdate,
        date_death_asdate,
        age,
        gender
      )

ras_df_cox_event <- ras_df_raw %>%
  filter(!is.na(potassium)) %>%
  arrange(date_measurement) %>%
  mutate(
    hypokalemia = potassium < 3.5
  ) %>%
  group_by(pat_id) %>%
  summarise(
    state_first_hypokalemia = as.numeric(any(hypokalemia)),
    date_first_hypokalemia = if_else(any(hypokalemia), date_measurement[which(hypokalemia)[1]], tail(date_measurement, 1)),
    hk_at_begin = date_first_hypokalemia == head(date_measurement, 1),
    n_meas = n()
  )

ras_df_cox_event %>%
  group_by(hk_at_begin) %>%
  summarise(n())

ras_df_cox %>%
  left_join(ras_df_cox_event) %>% 
  filter(date_measurement < date_first_hypokalemia) %>%
  group_by(pat_id) %>%
  arrange(pat_id, date_measurement) %>%
  mutate(
    date1 = date_measurement,
    date2 = lead(date_measurement),
    date2 = if_else(is.na(date2), date_first_hypokalemia, date2),
    state = if_else(date_first_hypokalemia > date2, 0, state_first_hypokalemia)
  )

ras_df_raw %>%
  filter(pat_id == 'SARS-FP-4') %>%
  select(date_measurement, potassium, Aldosterone) %>% 
  pivot_longer(
    -1
  ) %>%
  ggplot(aes(x= date_measurement, y = value, group = name)) +
  facet_grid(name~., scales = 'free_y') +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 3.5))
  
  
  

  
