librarian::shelf(tidyverse, lubridate ,stringr , readxl, ggthemes)


# read excel files ------------------------------------------------------------

clinical_data      <- read_excel('data/SARS-FP_clinical_data_2020_08_01_NEW_NEW.xlsx')
# lab_values_KFJ     <- read_excel('data/SARS-FP_KFJ_lab_values_all.xlsx',
#                                  col_names = c('pat_id', 'source', 'variable', 'range_normal',
#                                                'date_measurement', 'value', 'value_star', 'foo'))    # brauch ich wohl eher nicht
lab_values_KFJ     <- read_tsv('data/ACOVACT_lab_all_ID_no_names_2020_09_22.txt', 
                                 col_names = c('source', 'variable', 'range_normal', 
                                               'date_measurement', 'value', 'value_star', 'foo', 'pat_id')) 
lab_values_RAS     <- read_excel('data/SARS-FP_RAS_measurments_all_NEW_2.0_NEW_new actvity marker.xlsx')
ards_controls      <- read_excel('data/ARDS_controls_NEW.xlsx')
healthy_controls   <- read_excel('data/healthy controls_age_sexNEW.xlsx')
heart_controls     <- read_excel('data/heart_failure controls_NEW.xlsx')
hl63_conrols       <- read_excel('data/HL63_serum.xlsx') 
influenza_controls <- read_excel('data/Influenza_serum_co-morb_NEW_NEW.xlsx')
pneumonia_controls <- read_excel('data/pneumonia_serum.xlsx')
steroids_bmi       <- read_excel('data/SARS-FP steriods and bmi.xlsx',
                                 col_names = c('pat_id', 'steroids', 'bmi')) %>%
  mutate(
    bmi = as.numeric(bmi),
    steroids = ifelse(steroids == 1, 'yes', 'no')
  )
resp_values        <- read_excel('data/resp_samples_Covid_19_ACE2_BUN_normalized_and_PCR.xlsx')  %>%
  select(
    pat_id = "Patient ID",
    date_measurement = "Date-Sampling",
    ACE2_lung = "ACE-2 lung (directly measured)",
    ACE2_resp_normalized = "ACE2 resp normalized" ,
    exclude = "BUN resp LLQ (1=YES) -> möglicherweise zu sehr verdünnt"
  ) %>%
  mutate(
    date_measurement = as.Date(date_measurement)
  )

# data prep

# patient data ----------------------------------------------------------------

patients <- clinical_data %>% 
  select(pat_id = 'Pat ID', 
         gender, 
         age, 
         add_viro = `additional - viro`,
         category = `category (ICU, hospitalized, mild, asymptomatic)`,
         vent_status = `vent status`,
         date_first_symptoms = `Date first symptoms`,
         date_first_test = `Date first test`,
         date_hospital = `date hospital`,
         date_icu = `Date ICU`,
         date_intub = `date intub`,
         date_extub = `date extub`,
         date_icu_excharge = `Date ICU discharge`,
         date_hospital_discharge = `date hosptial discharge`,
         date_death = death,
         diabetes = `diabetes (y/n)`,
         hypertension = `art hypertension (y/n)`,
         copd = `COPD (y/n)`,
         ace_inhibitor = `ACEi (y/n)`,
         arb = `ARB (y/n)`,
         comedication = `1=hydroxychloroquin, 2= Kaletra, 3= Remdesivir, 4 = Plasma, 5= Tocilizumab, 6 = Camostat`
  ) %>% 
  mutate(
    date_first_symptoms_asdate = as.Date(as.numeric(date_first_symptoms), origin = "1899-12-30"),
    date_first_test_asdate = as.Date(as.numeric(date_first_test), origin = "1899-12-30"),
    date_hospital_asdate = as.Date(as.numeric(date_hospital), origin = "1899-12-30"),
    date_icu_asdate = as.Date(as.numeric(date_icu), origin = "1899-12-30"),
    date_discharge_asdate = as.Date(as.numeric(date_hospital_discharge), origin = "1899-12-30"),
    date_death_asdate = as.Date(as.numeric(date_death), origin = "1899-12-30"),
    vent_status = as.factor(vent_status),
    class = case_when(
      vent_status == 0        ~ 'asymptomatic',
      vent_status == 1        ~ 'mild',
      vent_status %in% c(2,3) ~ 'severe',
      vent_status == 4        ~ 'critical'
    ),
    class = factor(class, levels = c('asymptomatic', 'mild', 'severe', 'critical')),
    hydroxychloroquin = if_else(str_detect(comedication, '1'), 'yes', 'no', 'no'),
    Kaletra = if_else(str_detect(comedication, '2'), 'yes', 'no', 'no'),
    Remdesivir = if_else(str_detect(comedication, '3'), 'yes', 'no', 'no'),
    Plasma = if_else(str_detect(comedication, '4'), 'yes', 'no', 'no'),
    Tocilizumab = if_else(str_detect(comedication, '5'), 'yes', 'no', 'no'),
    Camostat = if_else(str_detect(comedication, '6'), 'yes', 'no', 'no'),
  )
         
# measurements ----------------------------------------------------------------

measurement_df <-
  clinical_data %>%
    select(pat_id = 'Pat ID', contains('Date')) %>%
    select(-c(14:21)) %>%
  pivot_longer(-1,
               values_to = 'date_measurement',
               names_to = 'measurement_nr') %>%
  bind_cols(
    clinical_data %>%
      select(pat_id_foo = 'Pat ID', contains('ACO-ID')) %>%
      pivot_longer(-1,
                   names_to = 'foo',
                   values_to = 'measurement_id')
  ) %>%
  select(-foo, -pat_id_foo) %>%
  filter(!is.na(measurement_id)) %>%
  mutate(
    date_measurement = as.Date(date_measurement)
  ) 


# lab data --------------------------------------------------------------------

lab_values <- lab_values_RAS %>%
  rename('measurement_id' = '...1') %>%
  select(
    measurement_id, 
    AngII = `Ang II (1-8)`,
    ACE2,
    ACE2_ang_based = `ACE2 (Ang based))`,
    Ang_1_7 = `Ang 1-7`,
    AngI = `Ang I (1-10)`,
    AngIII = `Ang III (2-8)`,
    Ang_1_5 = `Ang 1-5`,
    AngIV = `Ang IV (3-8)`,
    Aldosterone,
    Ang_1_7_plus_1_5 = `Ang 1-7 + Ang 1-5`,
    AA2_ratio = `AA2-Ratio`,
    PRA_S = `PRA-S`,
    ACE_S = `ACE-S`,
    ACE2_activity_v2 = `Ang 1-5 / AngII`,
    ACE2_activity_v3 = `(Ang1-5+Ang1-7) / (AngI+AngII+An1-5+Ang1-7)`
  ) %>%
  mutate(
    AngII = as.numeric(AngII),
    ACE2 = as.numeric(ACE2),
    ACE2_ang_based = str_replace_all(ACE2_ang_based, ',', '.'),
    ACE2_ang_based = ifelse(
      str_sub(ACE2_ang_based,1,1) == '>',
      as.numeric(str_sub(ACE2_ang_based, 3, -1L)),
      ifelse(
        str_sub(ACE2_ang_based,1,1) == '<',
        0.04,
        as.numeric(ACE2_ang_based)
      )
    ),
    AA2_ratio = as.numeric(AA2_ratio),
    ACE_S = as.numeric(ACE_S)
  )

# lab values ------------------------------------------------------------------

lab_values_2 <- lab_values_KFJ %>%
  select(pat_id, variable, date_measurement, value) %>%
  filter(
    variable %in% c('Interleukin-6')
  ) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(pat_id, date_measurement),
    names_from = variable,
    values_from = value
  ) %>%
  mutate(
    date_measurement = as.Date(date_measurement, format = '%d.%m.%y')
  ) %>%
  filter(
    !is.na(pat_id)
  ) %>% group_by(pat_id, date_measurement) %>%
  summarise_all(
    ~mean(as.numeric(str_replace(., '<', '')), na.rm = TRUE)
  )
  



# merge of patients, measurements and lab data --------------------------------

ras_df <- measurement_df %>% 
  left_join(lab_values) %>%
  left_join(lab_values_2 %>% rename(il_6 = `Interleukin-6`)) %>%
  left_join(resp_values) %>%
#  left_join(potassium_df)  hier lab values einfügen
  left_join(patients) %>%
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

saveRDS(ras_df, 'data/RAS_df.rds')


# controls --------------------------------------------------------------------

controls <- 
  bind_rows(
    ards_controls %>%
      rename(
        day = `sample ID`,
        AngII = `Ang II (1-8)`,
        ACE2_ang_based = `ACE2 (Ang based))`
      ) %>%
      select(
        day, AngII, ACE2_ang_based,
        ACE2,
        Ang_1_7 = `Ang 1-7`,
        AngI = `Ang I (1-10)`,
        AngIII = `Ang III (2-8)`,
        Ang_1_5 = `Ang 1-5`,
        AngIV = `Ang IV (3-8)`,
        Aldosterone,
        Ang_1_7_plus_1_5 = `Ang 1-7 + Ang 1-5`,
        AA2_ratio = `AA2-Ratio`,
        PRA_S = `PRA-S`,
        ACE_S = `ACE-S`
      ) %>%
      mutate(
        pat_id = str_sub(day, 1,2),
        days = as.numeric(str_sub(day, -1L, -1L)),
        AngII = as.numeric(AngII),
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        type = 'ARDS'
      ) %>%
      select(-day),
    healthy_controls %>%
      rename(
        AngII = `Ang II (1-8)`,
        ACE2_ang_based = `ACE-2 (Ang based)`
      ) %>%
      select(AngII, ACE2_ang_based, ACE2,
             Ang_1_7 = `Ang 1-7`,
             AngI = `Ang I (1-10)`,
             AngIII = `Ang III (2-8)`,
             Ang_1_5 = `Ang 1-5`,
             AngIV = `Ang IV (3-8)`,
             Aldosterone,
             Ang_1_7_plus_1_5 = `Ang 1-7 + Ang 1-5`,
             AA2_ratio = `AA2-Ratio`,
             PRA_S = `PRA-S`,
             ACE_S = `ACE-S`,
             age = Age,
             gender = Sex
      ) %>%
      mutate(
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        type = 'healthy',
        pat_id = paste0(type, 1:n())
      ),
    heart_controls %>%
      rename(
        AngII = `Ang II (1-8)`,
        ACE2_ang_based = `ACE-2 (Ang based)`
      ) %>%
      select(
        AngII, ACE2_ang_based, ACE2,
        Ang_1_7 = `Ang 1-7`,
        AngI = `Ang I (1-10)`,
        Ang_1_5 = `Ang 1-5`,
        Aldosterone,
        PRA_S = `PRA-S`
      ) %>%
      mutate(
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        type = 'heart failure',
        pat_id = paste0(type, 1:n())
      ),
    influenza_controls %>%
      rename(
        pat_id = `Patient-ID`,
        AngII = `Ang II (1-8)`,
        ACE2_ang_based = `ACE2 (Ang based))`,
        days = `Time after intub`
      ) %>%
      select(pat_id, days, AngII, ACE2_ang_based, ACE2, 
             Ang_1_7 = `Ang 1-7`,
             AngI = `Ang I (1-10)`,
             AngIII = `Ang III (2-8)`,
             Ang_1_5 = `Ang 1-5`,
             AngIV = `Ang IV (3-8)`,
             Aldosterone,
             Ang_1_7_plus_1_5 = `Ang 1-7 + Ang 1-5`,
             AA2_ratio = `AA2-Ratio`,
             PRA_S = `PRA-S`,
             ACE_S = 'ACE-S',
             age = 'Age',
             gender = 'sex',
             death,
             copd = 'COPD',
             diabetes = 'DM',
             hypertension = `art HAT`,
             ace_inhibitor = 'ACEi',
             arb = 'ARB'
             ) %>%
      mutate(
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        type = 'influenza',
        pat_id = as.character(pat_id),
        AA2_ratio = as.numeric(AA2_ratio),
        death = !is.na(death)
      ),
    # bind_cols(
    #   pneumonia_controls %>%
    #     select(pat_id = 'Pat ID', contains('ACO')) %>%
    #     pivot_longer(
    #       -1,
    #       values_to = 'measurement_id'),
    #   pneumonia_controls %>%
    #     select(pat_id = 'Pat ID', contains('Days since')) %>%
    #     pivot_longer(
    #       -1,
    #       values_to = 'days'
    #     ) %>%
    #     select(days)
    #   ) %>%
    #   left_join(lab_values) %>%
    #   select(-measurement_id, -name) %>%
    #   mutate(
    #     type = 'pneumonia',
    #     ACE2_ang_based = ifelse(
    #       str_sub(ACE2_ang_based,1,1) == '<',
    #       0.04,
    #       as.numeric(ACE2_ang_based)
    #     ),
    #     ACE2 = as.numeric(ACE2),
    #     AA2_ratio = as.numeric(AA2_ratio),
    #     ACE_S = as.numeric(ACE_S)
    #   )
  )

controls$copd[controls$copd == 'yey'] <- 'yes'

saveRDS(controls, 'data/controls_df.rds')











