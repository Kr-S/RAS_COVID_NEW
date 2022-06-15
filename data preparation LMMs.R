
library(tidyverse)
library(lubridate)
library(stringr)
library(readxl)
library(ggthemes)
library(janitor)
library(geepack)
library(geeM)
library(geeasy)
library(splitstackshape)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(naniar)
library(Hmisc)
library(corrplot)
library(plotmo)
library(splines)
library(sjlabelled)
library(stargazer)
library(table1)
library(tidyr)
library(broom)
library(scales)
library(flextable)
library(ggrepel)
library(lme4)
library(MuMIn)
library(DT)


# read excel files ------------------------------------------------------------

clinical_data      <- read_excel('data/data_update_20211202/2022-05-15_SARS-FP_clinical-data.xlsx')
group_ids          <- read_excel('data/new_category.xlsx') %>% rename(status = status3)
lab_values_KFJ     <- read_tsv('data/ACOVACT_lab_all_ID_no_names_2020_09_22.txt', 
                                 col_names = c('source', 'variable', 'range_normal', 
                                               'date_measurement', 'value', 'value_star', 'foo', 'pat_id')) 
ACOVACT_lab_all_ID_no_names_2020_09_22 <- read.delim("data/ACOVACT_lab_all_ID_no_names_2020_09_22.txt")
lab_values_RAS     <- read_excel('data/data_update_20211202/SARS-FP_RAS_measurments_all_NEW_2.0_NEW_new actvity marker_additional.xlsx') #new data
ards_controls      <- read_excel('data/ARDS_controls_NEW.xlsx')
healthy_controls   <- read_excel('data/healthy controls_age_sexNEW.xlsx')
heart_controls     <- read_excel('data/heart_failure controls_NEW.xlsx')
hl63_conrols       <- read_excel('data/HL63_serum.xlsx') 
influenza_controls <- read_excel('data/Influenza_serum_co-morb_NEW_NEW.xlsx')
pneumonia_controls <- read_excel('data/pneumonia_serum.xlsx')
insulin            <- read_excel("data/Insulin.xlsx")
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
    exclude = "BUN resp LLQ (1=YES) -> möglicherweise zu sehr verdünnt",
    viral_load = 'viral load',
    days_after_intub = 'days after intub'
  ) %>%
  mutate(
    date_measurement = as.Date(date_measurement)
  )
saveRDS(resp_values, 'data/resp_values_prep.rds')

ace_values <- read_excel('data/data_update_20211202/ACE_all_2021_12_02.xlsx') %>% #new data
  rename(
    measurement_id = `...1`
  ) %>%
  filter(measurement_id != '346') %>%
  mutate(ACE = as.numeric(ACE))



# data prep

# patient data ----------------------------------------------------------------

patients <- clinical_data %>% 
  select(pat_id = 'Pat ID', 
         sex = gender, 
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
         date_icu_discharge = `Date ICU discharge`,
         date_hospital_discharge = `date hosptial discharge`,
         date_death = death,
         Diabetes = `diabetes (y/n)`,
         hypertension = `art hypertension (y/n)`,
         copd = `COPD (y/n)`,
         ace_inhibitor = `ACEi (y/n)`,
         arb = `ARB (y/n)`,
         comedication = `1=hydroxychloroquin, 2= Kaletra, 3= Remdesivir, 4 = Plasma, 5= Tocilizumab, 6 = Camostat`
  ) %>% 
  mutate(
    arb = as.numeric(arb),
    date_first_symptoms_asdate = as.Date(as.numeric(date_first_symptoms), origin = "1899-12-30"),
    date_first_test_asdate = as.Date(as.numeric(date_first_test), origin = "1899-12-30"),
    date_hospital_asdate = as.Date(as.numeric(date_hospital), origin = "1899-12-30"),
    date_hospital_discharge_asdate = as.Date(as.numeric(date_hospital_discharge), origin = "1899-12-30"),
    date_icu_asdate = as.Date(as.numeric(date_icu), origin = "1899-12-30"),
    date_icu_discharge_asdate = as.Date(as.numeric(date_icu_discharge), origin = "1899-12-30"),
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
    RASi = factor(case_when(ace_inhibitor == "yes" | arb == "yes" ~1, TRUE ~ 0), levels=c(0,1),labels=c("No RAS-Inhibitor","RAS-Inhibitor"))
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
    pras = `PRA-S`,
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
    ACE_S = as.numeric(ACE_S),
    alt_RAS = Ang_1_7 + Ang_1_5,
    pras = AngI + AngII + AngIII + AngIV + Ang_1_7 + Ang_1_5,
    alt_RAS_ratio = (Ang_1_7 + Ang_1_5) / pras,
    AngII_to_Ang_1_7 = AngII / Ang_1_7,
    AngII_to_Ang_1_5 = AngII / Ang_1_5,
    Ang_1_7_to_AngII = Ang_1_7 / AngII
  )

# lab values ------------------------------------------------------------------

renamevars <- c(potassium = 'Kalium', 
                pot_urine ='Kalium /U', 
                il_6 = 'Interleukin-6',
                sod = 'Natrium',
                CRP = 'CRP',
                crea = "Kreatinin",               #which one?
                crea_urine = "Kreatinin /U",
                crea_serum = "Kreatinin /SM",
                d_dimer = "D-Dimer"
                )

laball <- ACOVACT_lab_all_ID_no_names_2020_09_22 %>%
  distinct() %>%
  mutate(Type= replace(Type ,Type=='Kalium art.', 'Kalium')) %>% 
  pivot_wider(names_from = Type, values_from = Lab.values) %>%
  select(
    renamevars,
    pat_id = 'ID',
    date_measurement_time = 'X.2'
  ) %>%
  pivot_longer(cols=names(renamevars),names_to="names", values_to="values") %>%
  mutate(date_measurement = as.Date(strptime(date_measurement_time, "%d.%m.%y %H:%M"))) %>%
  filter(!is.na(values)) %>%
  arrange(date_measurement_time) %>%
  group_by(pat_id, date_measurement, names) %>%
  summarise(
    values = as.numeric(head(values, 1))  #if multiple measures per day, take first value
  ) %>%
  pivot_wider(names_from = names, values_from = values)

# merge of patients, measurements and lab data --------------------------------

ras_df <- measurement_df %>% 
  full_join(lab_values) %>%
  full_join(meds) %>% 
  full_join(insulin) %>% 
  left_join(ace_values) %>%
  left_join(patients_pot) %>%
  left_join(group_ids %>% mutate(status = factor(status))) %>%
  left_join(steroids_bmi) %>%
  mutate(ras_measurement = TRUE) %>%
  full_join(resp_values) %>%
  full_join(laball) %>%
  tidyr::replace_na(list(ras_measurement = FALSE)) %>%
  mutate(ald = Aldosterone,
         ang2 = AngII,
         ang1 = AngI,
         days_since_test = ifelse(!is.na(date_first_test_asdate), as.numeric(date_measurement - date_first_test_asdate),
                                  ifelse(!is.na(date_hospital_asdate),
                                         as.numeric(date_measurement - date_hospital_asdate),
                                         as.numeric(date_measurement - date_icu_asdate))),
         days_since_hospital = ifelse(!is.na(date_hospital_asdate),as.numeric(date_measurement - date_hospital_asdate) ,
                                      ifelse(!is.na(date_first_test_asdate),as.numeric(date_measurement - date_first_test_asdate),
                                             as.numeric(date_measurement - date_icu_asdate))),
         end_of_icu = ifelse(is.na(date_icu_asdate),NA,
                             ifelse(!is.na(date_icu_discharge_asdate), date_icu_discharge_asdate, 
                                    ifelse(!is.na(date_hospital_discharge_asdate), date_hospital_discharge_asdate, date_death_asdate))),
         in_icu = ifelse(date_measurement>=date_icu_asdate, ifelse(date_measurement<=end_of_icu, 1,0),0),
         never_icu = factor(is.na(date_icu_asdate)*1, levels=c(0,1),labels=c("ICU", "Never ICU")),
         icu_or_death = ifelse(!is.na(date_death), 1, ifelse(never_icu == 0, 1,0)),
         in_icu = ifelse(is.na(in_icu),0,in_icu),
         days_from_icu_admission = as.numeric(date_measurement - date_icu_asdate),
         five_days_since_hospital = ceiling((days_since_hospital+1)/ 5),
         seven_days_since_hospital = ceiling((days_since_hospital+1)/ 7),
         hypoK = factor(ifelse(!is.na(potassium),ifelse(potassium<3.5,1,0),NA), levels=c(0,1), labels=c("No Hypokalemia", "Hypokalemia")),
         aa2r_LOQ = factor(ifelse(ald<20, 1, ifelse(ang2<2, 1, 0)), levels = c(0,1), labels = c("Regular", "LOQ")),
         ald = ifelse(ald<20, 20, ald),
         pras = ifelse(pras < 10 , 10, pras),
         ang2 = ifelse(ang2 < 2 , 2, ang2),
         aa2r = ald/ang2
  ) %>% 
  group_by(pat_id) %>% 
  #filter(ald>10) %>%
  mutate(ever_hypoK = factor(min(potassium,na.rm = TRUE)<3.5,levels=c(TRUE,FALSE),labels=c("Hypokalemia", "No Hypokalemia")),
         severe_hypoK = factor(ifelse(!is.na(potassium),ifelse(potassium<3,1,0),NA), levels=c(0,1), labels=c("No Hypokalemia", "Hypokalemia")),
         ever_severe_hypoK = factor(min(potassium,na.rm = TRUE)<3,levels=c(TRUE,FALSE),labels=c("Hypokalemia", "No Hypokalemia"))
         , first_week_hypoK = factor(case_when(any(hypoK == "Hypokalemia" & seven_days_since_hospital == 1) ~ 1, TRUE ~0) , levels=c(0,1), labels = c("No Hypokalemia during First Week","Hypokalemia During First Week"))
         , second_week_hypoK = factor(case_when(any(hypoK == "Hypokalemia" & seven_days_since_hospital == 2) ~ 1, TRUE ~0) , levels=c(0,1), labels = c("No Hypokalemia during Second Week","Hypokalemia During Second Week"))
         , third_week_hypoK = factor(case_when(any(hypoK == "Hypokalemia" & seven_days_since_hospital == 3) ~ 1, TRUE ~0) , levels=c(0,1), labels = c("No Hypokalemia during Third Week","Hypokalemia During Third Week"))
  ) %>% 
  mutate_at(vars(latest_pot = potassium, # getting last values
                 latest_ald = ald,
                 latest_ang2 = ang2,
                 latest_pras = pras,
                 latest_aa2r= aa2r,
                 latest_arb = arb,
                 latest_acei = acei,
                 latest_mra = mra,
                 latest_loop_diuretic = loop_diuretic,
                 latest_thiazid = thiazid_diuretic,
                 latest_pot_flush = pot_flush,
                 latest_pot_supp = pot_supp,
                 latest_cat = catecholamine
  ), lag) %>%
  mutate_at(vars(next_ald = ald,   # getting next values 
                 next_aa2r = aa2r,
                 next_ang2 = ang2,
                 next_pras = pras,
                 next_pot = potassium
  ), lead)%>% 
  mutate_at(vars(latest_arb, # setting 0 of binary vars to na, so the latest time always only shows time of latest given medication
                 latest_acei,
                 latest_mra,
                 latest_loop_diuretic,
                 latest_thiazid,
                 latest_pot_flush,
                 latest_pot_supp,
                 latest_cat),
            ~na_if(.,y=0)) %>% 
  mutate(latest_pot_date = as.Date(ifelse(!is.na(latest_pot),lag(date_measurement),NA), origin = "1970-01-01"), # get the times of latest and next values
         latest_ald_date = as.Date(ifelse(!is.na(latest_ald),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_arb_date = as.Date(ifelse(!is.na(latest_arb),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_acei_date = as.Date(ifelse(!is.na(latest_acei),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_mra_date = as.Date(ifelse(!is.na(latest_mra),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_loop_diuretic_date = as.Date(ifelse(!is.na(latest_loop_diuretic),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_thiazid_date = as.Date(ifelse(!is.na(latest_thiazid),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_pot_flush_date = as.Date(ifelse(!is.na(latest_pot_flush),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_pot_supp_date = as.Date(ifelse(!is.na(latest_pot_supp),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_cat_date =as.Date(ifelse(!is.na(latest_cat),lag(date_measurement),NA), origin = "1970-01-01"),
         next_ald_date = as.Date(ifelse(!is.na(next_ald),lead(date_measurement),NA), origin = "1970-01-01"),
         latest_aa2r_date = as.Date(ifelse(!is.na(latest_aa2r),lag(date_measurement),NA), origin = "1970-01-01"),
         next_aa2r_date = as.Date(ifelse(!is.na(next_aa2r),lead(date_measurement),NA), origin = "1970-01-01"),
         latest_ang2_date = as.Date(ifelse(!is.na(latest_ang2),lag(date_measurement),NA), origin = "1970-01-01"),
         next_ang2_date = as.Date(ifelse(!is.na(next_ang2),lead(date_measurement),NA), origin = "1970-01-01"),
         latest_pras_date = as.Date(ifelse(!is.na(latest_pras),lag(date_measurement),NA), origin = "1970-01-01"),
         next_pras_date = as.Date(ifelse(!is.na(next_pras),lead(date_measurement),NA), origin = "1970-01-01"),
         next_pot_date = as.Date(ifelse(!is.na(next_pot),lead(date_measurement),NA), origin = "1970-01-01"),
  ) %>% 
  tidyr::fill(latest_ald_date, latest_ald, # fill forward empty slots of the latest values and dates
              latest_pot, latest_pot_date,
              latest_aa2r_date, latest_aa2r,
              latest_arb,latest_arb_date,
              latest_acei,latest_acei_date,
              latest_mra,latest_mra_date,
              latest_loop_diuretic,latest_loop_diuretic_date,
              latest_thiazid,latest_thiazid_date,
              latest_pot_flush,latest_pot_flush_date,
              latest_pot_supp,latest_pot_supp_date,
              latest_cat, latest_cat_date,
              latest_ang2, latest_ang2_date,
              latest_pras, latest_pras_date,
              .direction = "down") %>% 
  tidyr::fill(next_ald, next_ald_date, # fill backward empty slots of next values and dates
              next_aa2r, next_aa2r_date,
              next_ang2, next_ang2_date,
              next_pras, next_pras_date,
              next_pot, next_pot_date,
              .direction = "up") %>% 
  mutate(days_since_ald = as.numeric(date_measurement - latest_ald_date), # calculate differences in days from dates and changes in vars if needed
         days_since_pot = as.numeric(date_measurement - latest_pot_date),
         days_since_arb = as.numeric(date_measurement - latest_arb_date),
         days_since_acei = as.numeric(date_measurement - latest_acei_date),
         days_since_mra = as.numeric(date_measurement - latest_mra_date),
         days_since_loop_diuretic = as.numeric(date_measurement - latest_loop_diuretic_date),
         days_since_thiazid = as.numeric(date_measurement - latest_thiazid_date),
         days_since_pot_flush = as.numeric(date_measurement - latest_pot_flush_date),
         days_since_pot_supp = as.numeric(date_measurement - latest_pot_supp_date),
         days_since_cat = as.numeric(date_measurement - latest_cat_date),
         days_to_next_ald = as.numeric(next_ald_date - date_measurement),
         days_since_aa2r = as.numeric(date_measurement - latest_aa2r_date),
         days_to_next_aa2r = as.numeric(next_aa2r_date - date_measurement),
         days_since_ang2 = as.numeric(date_measurement - latest_ang2_date),
         days_to_next_ang2 = as.numeric(next_ang2_date - date_measurement),
         days_since_pras = as.numeric(date_measurement - latest_pras_date),
         days_to_next_pras = as.numeric(next_pras_date - date_measurement),
         days_to_next_pot = as.numeric(next_pot_date - date_measurement),
         pot_change = potassium - latest_pot,
         Median_Potassium = median(potassium,na.rm=TRUE), 
         Median_Aldosterone = median(ald,na.rm=TRUE) ,
         median_ald_tert = factor(ntile(Median_Aldosterone,3), levels=c(1,2,3), labels=c("First", "Second", "Third")),
         median_pot_tert = factor(ntile(Median_Potassium,3), levels=c(1,2,3), labels=c("First", "Second", "Third"))) %>% 
  ungroup %>% 
  tidyr::replace_na(replace=list(days_since_arb = 0,
                                 days_since_acei = 0,
                                 days_since_mra = 0,
                                 days_since_loop_diuretic = 0,
                                 days_since_thiazid = 0,
                                 days_since_pot_flush = 0,
                                 days_since_pot_supp = 0,
                                 days_since_cat = 0,
                                 latest_arb = 0,
                                 latest_acei = 0,
                                 latest_mra = 0,
                                 latest_loop_diuretic = 0,
                                 latest_thiazid = 0,
                                 latest_pot_flush = 0,
                                 latest_pot_supp = 0,
                                 latest_cat = 0,
                                 arb = 0,
                                 acei = 0,
                                 mra = 0,
                                 loop_diuretic = 0,
                                 thiazid = 0,
                                 pot_flush = 0,
                                 pot_supp = 0)) %>% 
  mutate(insulin = if_else(insulin=="yes",1,0,missing = 0),
         Diabetes= factor(diabetes, levels=c("yes","no"), labels=c("Yes","No")),
         pot_tert = factor(ntile(potassium,3), levels=c(1,2,3), labels=c("First", "Second", "Third")),
         ald_tert = factor(ntile(ald,3), levels=c(1,2,3), labels=c("First", "Second", "Third")),
         pot_crea_urine = pot_urine/crea_urine
  )%>% 
  #  mutate_at(vars(latest_arb,latest_acei,latest_mra,latest_thiazid,latest_loop_diuretic,latest_pot_flush,latest_pot_supp),as.factor) %>% 
  filter(days_since_hospital %in% 0:20)%>% 
  distinct() %>% 
  relocate(any_of(c("pat_id", "sex","date_days", "days_since_hospital", "date_measurement", "date_hospital_asdate","seven_days_since_hospital", "potassium","first_week_hypoK", "second_week_hypoK", "third_week_hypoK","hypoK","pot_tert","hypoK","ever_hypoK","ald","ald_tert", "cat","latest_cat", "aa2r", "latest_aa2r", "next_aa2r", "never_icu","in_icu","date_icu_asdate","date_icu_discharge_asdate","date_hospital_discharge_asdate","date_death_asdate","num_date_death","date_death","end_of_icu","pat_id","days_since_ald", "days_since_pot", "latest_pot_date", "latest_pot", "ald", "latest_ald", "latest_ald_date" )), .after=potassium)


  

saveRDS(ras_df, 'data/data_update_20211202/RAS_df.rds')


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
        alt_RAS = `Ang 1-7 + Ang 1-5`,
        AA2_ratio = `AA2-Ratio`,
        pras = `PRA-S`,
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
             alt_RAS = `Ang 1-7 + Ang 1-5`,
             AA2_ratio = `AA2-Ratio`,
             pras = `PRA-S`,
             ACE_S = `ACE-S`,
             age = Age,
             sex = Sex
      ) %>%
      mutate(
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        alt_RAS = Ang_1_7 + Ang_1_5,
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
        pras = `PRA-S`
      ) %>%
      mutate(
        ACE2_ang_based = ifelse(
          str_sub(ACE2_ang_based,1,1) == '<',
          0.04,
          as.numeric(ACE2_ang_based)
        ),
        alt_RAS = Ang_1_7 + Ang_1_5,
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
             alt_RAS = `Ang 1-7 + Ang 1-5`,
             AA2_ratio = `AA2-Ratio`,
             pras = `PRA-S`,
             ACE_S = 'ACE-S',
             age = 'Age',
             sex = 'sex',
             death,
             copd = 'COPD',
             Diabetes = 'DM',
             hypertension = `art HAT`,
             ace_inhibitor = 'ACEi',
             arb = 'ARB'
             ) %>%
      mutate(
        arb=as.numeric(arb),
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
  ) %>%
  mutate(
    status = 'influenza',
    AngII_to_Ang_1_7 = AngII / Ang_1_7,
    AngII_to_Ang_1_5 = AngII / Ang_1_5,
    Ang_1_7_to_AngII = Ang_1_7 / AngII,
    pras = AngI + AngII + AngIII + AngIV + Ang_1_7 + Ang_1_5,
    alt_RAS_ratio = (Ang_1_7 + Ang_1_5) / pras
  )

controls$copd[controls$copd == 'yey'] <- 'yes'

saveRDS(controls, 'data/controls_df.rds')









