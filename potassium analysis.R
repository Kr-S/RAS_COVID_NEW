

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



add_columns <- function(df, columns){
  new <- rep(NA_character_, length(columns))
  names(new) <- columns
  mutate(df, !!!new)
}

# load files ----

clinical_data      <- read_excel('data/SARS-FP_clinical_data_2020_08_01_NEW_NEW.xlsx')
group_ids          <- read_excel('data/new_category.xlsx') %>% rename(status = status3)
lab_values_KFJ     <- read_tsv('data/ACOVACT_lab_all_ID_no_names_2020_09_22.txt', 
                               col_names = c('source', 'variable', 'range_normal', 
                                             'date_measurement', 'value', 'value_star', 'foo', 'pat_id')) 
ACOVACT_lab_all_ID_no_names_2020_09_22 <- read.delim("data/ACOVACT_lab_all_ID_no_names_2020_09_22.txt")
lab_values_RAS     <- read_excel('data/SARS-FP_RAS_measurments_all_NEW_2.0_NEW_new actvity marker.xlsx')
ards_controls      <- read_excel('data/ARDS_controls_NEW.xlsx')
healthy_controls   <- read_excel('data/healthy controls_age_sexNEW.xlsx')
heart_controls     <- read_excel('data/heart_failure controls_NEW.xlsx')
hl63_conrols       <- read_excel('data/HL63_serum.xlsx') 
influenza_controls <- read_excel('data/Influenza_serum_co-morb_NEW_NEW.xlsx')
pneumonia_controls <- read_excel('data/pneumonia_serum.xlsx')
raw_meds            <- read_excel('data/meds_long.xlsx')
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

ace_values <- read_excel('data/ACE_new.xlsx') %>%
  rename(
    measurement_id = `...1`
  ) %>%
  filter(measurement_id != '346') %>%
  mutate(ACE = as.numeric(ACE))

# patient data ----


patients_pot <- clinical_data %>% 
  select(pat_id = 'Pat ID', 
         sex = gender, 
         age, 
         vent = `vent status`,
         date_first_symptoms = `Date first symptoms`,
         date_first_test = `Date first test`,
         date_hospital = `date hospital`,
         date_icu = `Date ICU`,
         date_icu_discharge = `Date ICU discharge`,
         date_hospital_discharge = `date hosptial discharge`,
         date_death = death,
         diabetes = `diabetes (y/n)`,
         hypertension = `art hypertension (y/n)`,
         copd = `COPD (y/n)`
         ) %>% 
  mutate(date_first_symptoms_asdate = as.Date(as.numeric(date_first_symptoms), origin = "1899-12-30"),
         date_first_test_asdate = as.Date(as.numeric(date_first_test), origin = "1899-12-30"),
          date_hospital_asdate = as.Date(as.numeric(date_hospital), origin = "1899-12-30"),
          date_hospital_discharge_asdate = as.Date(as.numeric(date_hospital_discharge), origin = "1899-12-30"),
          date_icu_asdate = as.Date(as.numeric(date_icu), origin = "1899-12-30"),
          date_icu_discharge_asdate = as.Date(as.numeric(date_icu_discharge), origin = "1899-12-30"),
          date_discharge_asdate = as.Date(as.numeric(date_hospital_discharge), origin = "1899-12-30"),
          date_death_asdate = as.Date(as.numeric(date_death), origin = "1899-12-30"),
         age = age/10,
         severity = dplyr::case_when(
           vent %in% 0        ~ 0,
           vent %in% 1        ~ 1,
           vent %in% c(2,3)   ~ 2,
           vent %in% 4        ~ 3
         ),
         vent = factor(vent, levels = c(0,1,2,3,4), labels = c("None","Nasal O2","High Flow","Pressurized Mask","Intubation")),
         severity = factor(severity, levels = c(0,1,2,3), labels = c('Asymptomatic', 'Mild', 'Severe', 'Critical')),
         intubation = factor(ifelse(vent=="Intubation",1,0),levels=c(0,1), labels= c("Not Intubated", "Intubated"))) %>% 
  var_labels(vent="Breathing Assistance",
             intubation = "Ever Intubated",
             severity = "Severity by Respiratory Assistance")

# measurements ----

measurement_ids <- clinical_data %>%
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

# ras lab values ----

ras_lab <- lab_values_RAS %>%
  rename('measurement_id' = '...1') %>%
  select(ald = Aldosterone, pras = `PRA-S`, ang2 = `Ang II (1-8)`,
         measurement_id)%>%
  full_join(measurement_ids) %>% 
  select(pras, ald, date_measurement, pat_id,ang2)


# reformat meds ----


meds <- raw_meds %>% 
  mutate(start_date = replace(start_date, start_date=="ongoing", 
                              as.numeric(min(raw_meds$start_date[raw_meds$start_date!="ongoing"], na.rm = TRUE), na.rm = TRUE)),
         end_date = replace(end_date, end_date=="ongoing", 
                            as.numeric(max(raw_meds$end_date[raw_meds$end_date!="ongoing"], na.rm = TRUE), na.rm = TRUE)),
         start_date=excel_numeric_to_date(as.numeric(start_date)),
         end_date = excel_numeric_to_date(as.numeric(end_date)))%>%
  filter(!is.na(start_date)|!is.na(end_date)) %>% 
  mutate(drug_class = tolower(drug_class)) %>%
  filter(!is.na(start_date)|!is.na(end_date)) %>% 
  mutate(drug_class = tolower(drug_class),
         drug= tolower(drug),
         daycount = end_date - start_date) %>% 
  transmute(pat_id, drug_class,drug, date_measurement = map2(start_date, end_date, seq, by = "1 day")) %>%
  unnest(cols=c(date_measurement)) %>% 
  mutate(classtreated=1) %>% 
  pivot_wider(id_cols=c(pat_id, date_measurement),names_from=drug_class, values_from=classtreated, values_fn = function(values) first(values),values_fill = 0) %>% 
  select(-"NA", thiazid_diuretic = "thiazid diuretic", loop_diuretic = "loop diuretic")


# potassium analysis time frame ----


renamevars_pot <- c(potassium = 'Kalium')

potassium_df <- ACOVACT_lab_all_ID_no_names_2020_09_22 %>%
  distinct() %>%
  pivot_wider(names_from = Type, values_from = Lab.values) %>%
  select(
    renamevars_pot,
    pat_id = 'ID',
    date_measurement_time = 'X.2') %>%
  pivot_longer(cols=names(renamevars_pot),names_to="names", values_to="values") %>%
  mutate(date_measurement = as.Date(strptime(date_measurement_time, "%d.%m.%y %H:%M"))) %>%
  filter(!is.na(as.numeric(values))) %>%
  arrange(date_measurement_time) %>%
  group_by(pat_id, date_measurement, names) %>%
  summarise(values = as.numeric(head(values, 1))) %>%  #if multiple measures per day, take first value
  pivot_wider(names_from = names, values_from = values) %>% 
  left_join(patients_pot) %>% 
  left_join(steroids_bmi) %>% 
  mutate(days_since_test = ifelse(!is.na(date_first_test_asdate), as.numeric(date_measurement - date_first_test_asdate),
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
         never_icu = factor(is.na(in_icu)*1, levels=c(0,1),labels=c("ICU", "Never ICU")),
         icu_or_death = ifelse(!is.na(date_death), 1, ifelse(never_icu == 0, 1,0)),
         in_icu = ifelse(is.na(in_icu),0,in_icu),
         days_from_icu_admission = as.numeric(date_measurement - date_icu_asdate)
         ) %>% 
  full_join(insulin) %>%
  full_join(meds) %>% 
  full_join(ras_lab) %>% 
  group_by(pat_id) %>% 
  #filter(ald>10) %>%
  mutate(aa2r_LOQ = factor(ifelse(ald<20, 1, ifelse(ang2<2, 1, 0)), levels = c(0,1), labels = c("Regular", "LOQ")),
         ald = ifelse(ald<20, 20, ald),
         pras = ifelse(pras < 10 , 10, pras),
         ang2 = ifelse(ang2 < 2 , 2, ang2),
         aa2r = ald/ang2) %>% 
  mutate_at(vars(latest_pot = potassium, # creating last values
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
  mutate_at(vars(next_ald = ald,   # creating next values 
            next_aa2r = aa2r,
            next_ang2 = ang2,
            next_pras = pras
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
         ) %>% 
  fill(latest_ald_date, latest_ald, # fill backwards empty slots of the latest values and dates
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
  fill(next_ald, next_ald_date, # fill forward empty slots of next values and dates
       next_aa2r, next_aa2r_date,
       next_ang2, next_ang2_date,
       next_pras, next_pras_date,
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
         pot_change = potassium - latest_pot) %>% 
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
         five_days_since_hospital = ceiling((days_since_hospital+1)/ 5),
         seven_days_since_hospital = ceiling((days_since_hospital+1)/ 7),
         Diabetes= factor(diabetes, levels=c("yes","no"), labels=c("Yes","No"))) %>% 
#  mutate_at(vars(latest_arb,latest_acei,latest_mra,latest_thiazid,latest_loop_diuretic,latest_pot_flush,latest_pot_supp),as.factor) %>% 
  filter(days_since_hospital %in% 0:20)%>% 
  relocate(any_of(c("pat_id", "date_days", "potassium", "cat","latest_cat", "aa2r", "latest_aa2r", "next_aa2r", "never_icu","in_icu","date_icu_asdate","date_icu_discharge_asdate","date_hospital_discharge_asdate","date_death_asdate","num_date_death","date_death","end_of_icu","pat_id", "days_since_hospital", "seven_days_since_hospital","days_since_ald", "days_since_pot", "latest_pot_date", "latest_pot", "ald", "latest_ald", "latest_ald_date" )), .after=potassium)


potaldvars <- c("pat_id","in_icu","potassium","latest_ald","latest_pot","days_since_ald","days_since_pot","sex","age","days_since_hospital","days_since_test", "pot_supp", "days_since_arb", "days_since_acei","days_since_mra", "days_since_loop_diuretic", "days_since_thiazid", "days_since_pot_flush", "days_since_pot_supp", "latest_arb", "latest_acei", "latest_mra", "latest_loop_diuretic", "latest_thiazid", "latest_pot_flush", "latest_pot_supp", "latest_cat", "days_since_cat")



names(potassium_df)
""
## vars "pat_id","potassium","latest_ald","latest_pot","latest_ald_date","latest_pot_date","sex","age","days_since_hospital"

## meds "insulin","beta-blocker","alpha-blocker","calcium channel blocker","alpha-2-agonist","kalzium antagonist",
##      "potassium channel opener","alpha1-adrenoceptor-antagonist" "arb","acei"m"loop diuretic","thiazid diuretic","carbonic anhydrase inhibitor"
##      "sulfonamide-based diuretic","mra","catecholamine","pot_supp","pot_flush","calc","phos","mag"   




length(unique(potassium_df$pat_id))

potald <- potassium_df %>% 
  select(all_of(potaldvars)) %>% 
  drop_na() %>% 
  mutate(ID = as.numeric(factor(pat_id)))

hist(log10(potald$latest_ald))
qqnorm(log10(potald$latest_ald))
hist(log10(potald$latest_pot))
qqnorm(log10(potald$latest_pot))
hist(potald$pot_change)
qqnorm(potald$pot_change)
hist(log10(potald$potassium))
qqnorm(log10(potald$potassium))
hist(log10(potald$latest_ald*potald$days_since_ald))
qqnorm(log10(potald$latest_ald*potald$days_since_ald))

hist(potald$days_since_hospital)
qqnorm(potald$days_since_hospital)

phist<- potald %>% ggplot(aes(x=log10(potassium))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
pqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log10(potassium))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(phist, pqq, labels = "AUTO")


lahist<- potald %>% ggplot(aes(x=log10(latest_ald+5))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
laqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log10(latest_ald+5))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(lahist, laqq, labels = "AUTO")
lphist <- potald %>% ggplot(aes(x=log10(latest_pot))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
lpqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log10(latest_pot))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(lphist, lpqq, labels = "AUTO")

pchist<- potald %>% ggplot(aes(x=pot_change)) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
pcqq<-ageqq<-potald %>% ggplot(aes(sample=log10(pot_change))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(pchist, pcqq, labels = "AUTO")

dsahist<- potald %>% ggplot(aes(x=log10(days_since_ald))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
dsaqq<-ageqq<-potald %>% ggplot(aes(sample=log10(days_since_ald))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(dsahist, dsaqq, labels = "AUTO")

dshhist<-potald %>% ggplot(aes(x=days_since_hospital)) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
dshqq<-ageqq<-potald %>% ggplot(aes(sample=days_since_hospital)) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(dshhist, dshqq, labels = "AUTO")



agehist<-potald %>% ggplot(aes(x=age)) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
ageqq<-potald %>% ggplot(aes(sample=age)) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(agehist, ageqq, labels = "AUTO")


patient_sex <- potald %>% group_by(pat_id) %>% filter(row_number()==1) %>% ungroup() %>%  ggplot()+
  geom_bar(aes(x=sex))+
  labs(title="Patients")
complete_observation_sex <- potald %>%ggplot()+
  geom_bar(aes(x=sex))+
  labs(title="Observations")
plot_grid(patient_sex, complete_observation_sex, labels = "AUTO")


potald %>%ggplot()+
  geom_bar(aes(x=pot_supp))+
  labs(title="Observations")

  
dim(potald)


##IDA

IDA = potald %>% select(age, in_icu, days_since_hospital,
                                    latest_ald, days_since_ald,
                                    latest_pot, days_since_pot, 
                                    latest_arb,days_since_arb,
                                    latest_acei, days_since_acei,
                                    latest_mra,days_since_mra,
                                    latest_loop_diuretic,days_since_loop_diuretic,
                                    latest_thiazid,days_since_thiazid,
                                    latest_pot_flush,days_since_pot_flush,
                                    latest_pot_supp,days_since_pot_supp) %>% 
  mutate(ald_time_interaction =latest_ald*days_since_ald,
         pot_time_interaction = latest_pot*days_since_pot,
         arb_time_interaction = latest_arb*days_since_arb,
         acei_time_interaction = latest_acei*days_since_acei,
         mra_time_interaction = latest_mra*days_since_mra,
         loop_diuretic_time_interaction = latest_loop_diuretic*days_since_loop_diuretic,
         thiazid_time_interaction = latest_thiazid*days_since_thiazid,
         pot_flush_time_interaction = latest_pot_flush*days_since_pot_flush,
         pot_supp_time_interaction = latest_pot_supp*days_since_pot_supp
         )



IDA = data.frame(IDA)
logIDA = data.frame(lapply(IDA, log10))
logIDA2 = data.frame(lapply(logIDA, log10))
dev.off()
pdf(file="output/Predictor Distributions.pdf")
for (x in c(1:length(colnames(IDA)))){
  if (all(IDA[,x] %in% c(0,1))){
    print(ggplot(data=data.frame(IDA), aes(x = IDA[,x])) +
            geom_bar(colour="black", fill="white") +
            labs(x=colnames(IDA[x]))
          )
  } else {
  print(ggplot(data=data.frame(IDA), aes(x = IDA[,x])) +
          geom_histogram(data=data.frame(IDA)[x], aes(y=..density..), colour="black", fill="white")+
          geom_density(data=data.frame(IDA)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
          labs(title=colnames(IDA)[x],x=colnames(IDA)[x], y = "Density"))
  print(ggqqplot(IDA[,x], main=colnames(IDA)[x]))
  print(ggplot(data=data.frame(logIDA), aes(x = logIDA[,x])) +
          geom_histogram(data=data.frame(logIDA)[x], aes(y=..density..), colour="black", fill="white")+
          geom_density(data=data.frame(logIDA)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
          labs(title=paste("",colnames(logIDA)[x]),x=colnames(logIDA)[x], y = "Density"))
  print(ggqqplot(logIDA[,x], main=paste("logged",colnames(logIDA)[x])))
  print(ggplot(data=data.frame(logIDA2), aes(x = logIDA2[,x])) +
          geom_histogram(data=data.frame(logIDA2)[x], aes(y=..density..), colour="black", fill="white")+
          geom_density(data=data.frame(logIDA2)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
          labs(title=paste("double logged",colnames(logIDA2)[x]),x=colnames(logIDA2)[x], y = "Density"))
  print(ggqqplot(logIDA2[,x], main=paste("double logged",colnames(logIDA2)[x])))
  }
}
dev.off()

#### Descriptive Analysis ####


#### Table One ####

t1df <- potassium_df %>% group_by(pat_id) %>% mutate(Mean_Potassium = median(potassium,na.rm=TRUE), Mean_Aldosterone = median(ald,na.rm=TRUE))%>% 
  filter(row_number()==1) %>% 
  ungroup %>%  
  mutate(Age = age, Sex = factor(sex, levels=c("f","m"), labels=c("Female","Male")), BMI = bmi,
            Hypertension = factor(hypertension, levels=c("yes","no"), labels=c("Yes","No")),
            COPD = factor(copd, levels=c("yes","no"),labels=c("Yes","No")), 
            pat_id, ICU = never_icu, Potassium = potassium, Aldosterone = ald) %>% 
  var_labels(Age = "Age (years)",BMI = "BMI (kg/m²)", Mean_Potassium ="Serum Potassium (mmol/L)", Mean_Aldosterone= "Serum Aldosterone (pmol/L)", ICU = "ICU stay")


my_plots <- lapply(names(t1df), function(var_x){
  p <- 
    ggplot(t1df) +
    aes_string(var_x)
  if(is.numeric(t1df[[var_x]])) {
    p <- p +geom_histogram(aes(y=..density..), colour="black", fill="white")+ geom_density(alpha=.2, fill="#FF6666")
  } 
  else {
    p <- p + geom_bar(colour="black", fill="white")
  } 
})

plot_grid(my_plots)

(t1 <- table1(~ Age + BMI + Diabetes + Hypertension + COPD + ICU + Mean_Potassium + Mean_Aldosterone  | Sex, data=t1df, render.continuous="Median [Q1,Q3]"))

write.table(t1,file="output/T1.png")


#### Smooth Plots ####



total_aldplot <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2)+
  ylab("Aldosterone (pmol/L,Log Scale)")+
  xlab("Days since hospitalization")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(9,100), n.breaks=10)

total_potplot <- potassium_df%>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  theme(legend.title = element_blank())+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

total_aa2rplot <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since hospitalization")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

test_total_aldplot <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2) +
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days since test")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(9,100), n.breaks=10)

test_total_potplot <- potassium_df%>% 
  ggplot(aes(x=days_since_test,y=potassium)) +  
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since test")+  
  theme(legend.title = element_blank())+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

test_total_aa2rplot <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2) +
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since test")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

aldICUplot <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2) +
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days from ICU admission")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(9,100), n.breaks=10)

potICUplot <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=potassium)) +  
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days from ICU admission")+
  theme(legend.title = element_blank())+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aa2rICUplot <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=never_icu, linetype=never_icu),col="grey0", alpha = 0.2) +
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days from ICU admission")+
  theme(legend.title = element_blank())+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)


smoothplots_nolegend <- cowplot::plot_grid(total_potplot + theme(legend.position="none")
                                  ,test_total_potplot + theme(legend.position="none")
                                  ,potICUplot + theme(legend.position="none")
                                  ,total_aldplot + theme(legend.position="none")
                                  ,test_total_aldplot + theme(legend.position="none")
                                  ,aldICUplot + theme(legend.position="none")
                                  ,total_aa2rplot + theme(legend.position="none")
                                  ,test_total_aa2rplot + theme(legend.position="none")
                                  ,aa2rICUplot + theme(legend.position="none")
                                  , labels="AUTO", label_size = 12, ncol=3, nrow=3)


legend <- cowplot::get_legend(total_aldplot+ labs(linetype = "ICU Admission",col = "ICU Admission"))

smoothplots <- cowplot::plot_grid(smoothplots_nolegend,legend, ncol=2, rel_widths = c(2,.2))



cowplot::save_plot('output/Trend of Ald, AA2-R and Pot by ICU Admission.png', smoothplots, ncol=3, nrow=3)



#### Smooth Plots by Ventilation####



total_aldplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("Aldosterone (pmol/L,Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

total_potplot_vent <- potassium_df%>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

total_aa2rplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

test_total_aldplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

test_total_potplot_vent <- potassium_df%>% 
  ggplot(aes(x=days_since_test,y=potassium)) +  
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since test")+  
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

test_total_aa2rplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

aldICUplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

potICUplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=potassium)) +  
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days from ICU admission")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aa2rICUplot_vent <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=vent, linetype=vent, col=vent), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

legend <- cowplot::get_legend(total_aldplot_vent+ labs(linetype = "Breathing Assistance",col = "Breathing Assistance"))

smoothplots_vent_nolegend <- cowplot::plot_grid(total_potplot_vent + theme(legend.position="none"),
                                                test_total_potplot_vent + theme(legend.position="none"),
                                                potICUplot_vent + theme(legend.position="none"),
                                                total_aldplot_vent + theme(legend.position="none"),
                                                test_total_aldplot_vent + theme(legend.position="none"),
                                                aldICUplot_vent + theme(legend.position="none"),
                                                total_aa2rplot_vent + theme(legend.position="none"),
                                                test_total_aa2rplot_vent +theme(legend.position="none"),
                                                aa2rICUplot_vent + theme(legend.position="none"), labels="AUTO", label_size = 12, ncol=3, nrow=3)
smoothplots_vent <- cowplot::plot_grid(smoothplots_vent_nolegend,legend, ncol=2, rel_widths = c(2,.2))


cowplot::save_plot('output/Trend of Ald, AA2-R and Pot by Maximum Ventilation Type.png', smoothplots_vent, ncol=3, nrow=3)


#### Smooth Plots by Severity####



total_aldplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("Aldosterone (pmol/L,Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

total_potplot_severity <- potassium_df%>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

total_aa2rplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

test_total_aldplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

test_total_potplot_severity <- potassium_df%>% 
  ggplot(aes(x=days_since_test,y=potassium)) +  
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since test")+  
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

test_total_aa2rplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

aldICUplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

potICUplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=potassium)) +  
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days from ICU admission")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aa2rICUplot_severity <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=severity, linetype=severity, col=severity), alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

legend <- cowplot::get_legend(total_aldplot_severity+ labs(linetype = "Severity",col = "Severity"))

smoothplots_severity_nolegend <- cowplot::plot_grid(total_potplot_severity + theme(legend.position="none"),
                                                test_total_potplot_severity + theme(legend.position="none"),
                                                potICUplot_severity + theme(legend.position="none"),
                                                total_aldplot_severity + theme(legend.position="none"),
                                                test_total_aldplot_severity + theme(legend.position="none"),
                                                aldICUplot_severity + theme(legend.position="none"),
                                                total_aa2rplot_severity + theme(legend.position="none"),
                                                test_total_aa2rplot_severity +theme(legend.position="none"),
                                                aa2rICUplot_severity + theme(legend.position="none"), labels="AUTO", label_size = 12, ncol=3, nrow=3)
smoothplots_severity <- cowplot::plot_grid(smoothplots_severity_nolegend,legend, ncol=2, rel_widths = c(2,.2))


cowplot::save_plot('output/Trend of Ald, AA2-R and Pot by Maximum Severity.png', smoothplots_severity, ncol=3, nrow=3)


#### Smooth Plots by Intubation####

total_aldplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("Aldosterone (pmol/L,Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

total_potplot_intubation <- potassium_df%>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

total_aa2rplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since hospitalization")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

test_total_aldplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

test_total_potplot_intubation <- potassium_df%>% 
  ggplot(aes(x=days_since_test,y=potassium)) +  
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since test")+  
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

test_total_aa2rplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days since test")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

aldICUplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=ald)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("Aldosterone (pmol/L, Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(9,100), n.breaks=10)

potICUplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=potassium)) +  
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days from ICU admission")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aa2rICUplot_intubation <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=aa2r)) +
  theme_cowplot(12)+
  geom_smooth(aes(group=intubation, linetype=intubation), col="grey0", alpha = 0.2)+
  ylab("AA2-Ratio (Log Scale)")+
  xlab("Days from ICU admission")+
  scale_y_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,10), n.breaks=12)

legend <- cowplot::get_legend(total_aldplot_intubation+ labs(linetype = "Ever Intubated"))

smoothplots_intubation_nolegend <- cowplot::plot_grid(total_potplot_intubation + theme(legend.position="none"),
                                                    test_total_potplot_intubation + theme(legend.position="none"),
                                                    potICUplot_intubation + theme(legend.position="none"),
                                                    total_aldplot_intubation + theme(legend.position="none"),
                                                    test_total_aldplot_intubation + theme(legend.position="none"),
                                                    aldICUplot_intubation + theme(legend.position="none"),
                                                    total_aa2rplot_intubation + theme(legend.position="none"),
                                                    test_total_aa2rplot_intubation +theme(legend.position="none"),
                                                    aa2rICUplot_intubation + theme(legend.position="none"), labels="AUTO", label_size = 12, ncol=3, nrow=3)

smoothplots_intubation <- cowplot::plot_grid(smoothplots_intubation_nolegend,legend, ncol=2, rel_widths = c(2,.2))

cowplot::save_plot('output/Trend of Ald, AA2-R and Pot by Intubation.png', smoothplots_intubation, ncol=3, nrow=3)



#### Potassium by Aldosterone per 5 day Interval ####



scx = scale_x_log10(limits=c(19,400), n.breaks=10)
scy = scale_y_continuous(breaks=seq(2.5,5.6,by=0.5), limits=c(2.5,5.6))
theme = theme_cowplot(12)


fullscatter <-potassium_df %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

fullscatter_latest <- potassium_df%>% 
  filter(days_since_ald <2) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

fullscatter_next <- potassium_df%>% 
  filter(days_to_next_ald <2) %>% 
  ggplot(aes(y=potassium, x=next_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter1 <-potassium_df %>% 
  filter(seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest1 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_next1 <- potassium_df %>% 
  filter(days_to_next_ald <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=next_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter2 <-potassium_df %>% 
  filter(seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest2 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_next2 <- potassium_df %>% 
  filter(days_to_next_ald <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=next_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter3 <-potassium_df %>% 
  filter(seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme +geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest3 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_next3 <- potassium_df %>% 
  filter(days_to_next_ald <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=next_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=20),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

(aldpotscat <- cowplot::plot_grid( fullscatter,fullscatter_latest,fullscatter_next, scatter1,scatter_latest1,scatter_next1, scatter2,scatter_latest2,scatter_next2,scatter3,scatter_latest3,scatter_next3, label_size = 12, ncol=3, labels = "AUTO"))

cowplot::save_plot('output/Lack of Correlation Ald and Pot.png', aldpotscat, ncol=3, nrow=4)

#### Potassium by AA2R per 5 day interval ####

scx = scale_x_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,max(potassium_df$aa2r, na.rm=TRUE)+1), n.breaks=12)

scy = scale_y_continuous(breaks=seq(2.5,5.6,by=0.5), limits=c(2.5,5.6))
theme = theme_cowplot(12)


aa2r_fullscatter <- potassium_df %>% 
  ggplot(aes(y=potassium, x=aa2r, color=aa2r_LOQ))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme(legend.title = element_blank())+
  theme  + annotation_logticks(sides="b") 

aa2r_fullscatter_latest <- potassium_df%>% 
  filter(days_since_aa2r <2) %>% 
  ggplot(aes(y=potassium, x=latest_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme  + annotation_logticks(sides="b") 

aa2r_fullscatter_next <- potassium_df%>% 
  filter(days_to_next_aa2r <2) %>% 
  ggplot(aes(y=potassium, x=next_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx + 
  theme  + annotation_logticks(sides="b") 

aa2r_scatter1 <-potassium_df %>% 
  filter(seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_latest1 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=latest_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_next1 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=next_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx +
  theme + annotation_logticks(sides="b") 

aa2r_scatter2 <-potassium_df %>% 
  filter(seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_latest2 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=latest_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_next2 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=next_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter3 <-potassium_df %>% 
  filter(seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_latest3 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=latest_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

aa2r_scatter_next3 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=next_aa2r))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme  + annotation_logticks(sides="b") 

(aa2rpotscat <- cowplot::plot_grid( aa2r_fullscatter,aa2r_fullscatter_latest,aa2r_fullscatter_next, aa2r_scatter1,aa2r_scatter_latest1,aa2r_scatter_next1, aa2r_scatter2,aa2r_scatter_latest2,aa2r_scatter_next2,aa2r_scatter3,aa2r_scatter_latest3,aa2r_scatter_next3, label_size = 12, ncol=3, labels = "AUTO"))

cowplot::save_plot('output/Lack of Correlation AA2-R and Pot.png', aa2rpotscat, ncol=3, nrow=4)


#### Ang2 by AA2R per 5 day interval ####

scx = scale_x_log10(limits=c(min(potassium_df$aa2r,na.rm = TRUE)-min(potassium_df$aa2r)/100,max(potassium_df$aa2r, na.rm=TRUE)+1), n.breaks=12)

min_ang2 <- min(potassium_df$ang2, na.rm=TRUE)-min(potassium_df$ang2, na.rm=TRUE)/100
max_ang2 <- max(potassium_df$ang2, na.rm=TRUE)+min(potassium_df$ang2, na.rm=TRUE)/100
scy = scale_y_continuous(breaks=seq(min_ang2,max_ang2,by=0.5), limits=c(min_ang2,max_ang2))
theme = theme_cowplot(12)


aa2r_fullscatter <- potassium_df %>% 
  ggplot(aes(y=ang2, x=aa2r, color=aa2r_LOQ))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme(legend.title = element_blank())+
  theme  + annotation_logticks(sides="b") 

aa2r_fullscatter_latest <- potassium_df%>% 
  filter(days_since_aa2r <2) %>% 
  ggplot(aes(y=ang2, x=latest_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme  + annotation_logticks(sides="b") 

aa2r_fullscatter_next <- potassium_df%>% 
  filter(days_to_next_aa2r <2) %>% 
  ggplot(aes(y=ang2, x=next_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx + 
  theme  + annotation_logticks(sides="b") 

ang2_aa2r_scatter1 <-potassium_df %>% 
  filter(seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=ang2, x=aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_latest1 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=ang2, x=latest_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_next1 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=ang2, x=next_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx +
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter2 <-potassium_df %>% 
  filter(seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=ang2, x=aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_latest2 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=ang2, x=latest_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_next2 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=ang2, x=next_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter3 <-potassium_df %>% 
  filter(seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=ang2, x=aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Same Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_latest3 <- potassium_df %>% 
  filter(days_since_aa2r <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=ang2, x=latest_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Previous Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme + annotation_logticks(sides="b") 

ang2_aa2r_scatter_next3 <- potassium_df %>% 
  filter(days_to_next_aa2r <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=ang2, x=next_aa2r))+
  geom_point()+
  ylab("Angiotensin 2 (mmol/L)")+
  xlab("Next Day Aldosterone-Angiotensin-2-Ratio (Logarithmic Scale)")+
  scy +
  scx+
  theme  + annotation_logticks(sides="b") 

(aa2rpotscat <- cowplot::plot_grid( aa2r_fullscatter,aa2r_fullscatter_latest,aa2r_fullscatter_next, ang2_aa2r_scatter1,ang2_aa2r_scatter_latest1,ang2_aa2r_scatter_next1, ang2_aa2r_scatter2,ang2_aa2r_scatter_latest2,ang2_aa2r_scatter_next2,ang2_aa2r_scatter3,ang2_aa2r_scatter_latest3,ang2_aa2r_scatter_next3, label_size = 12, ncol=3, labels = "AUTO"))

cowplot::save_plot('output/Scatterplots AA2-R and Angiotensin 2.png', aa2rpotscat, ncol=3, nrow=4)






#### GEE without splines ####


formule1 <- potassium ~ 
  age + sex + in_icu + days_since_hospital+
  latest_ald + latest_ald:days_since_ald +
  latest_pot + latest_pot:days_since_pot+
  latest_arb+latest_arb:days_since_arb+
  latest_acei+latest_acei:days_since_acei+
  latest_mra+latest_mra:days_since_mra+
  latest_loop_diuretic+latest_loop_diuretic:days_since_loop_diuretic+
  latest_thiazid+latest_thiazid:days_since_thiazid+
  latest_pot_flush+latest_pot_flush:days_since_pot_flush+
  latest_pot_supp+latest_pot_supp:days_since_pot_supp+
  latest_cat + latest_cat:days_since_cat



potald <- potald %>% 
  mutate(hypoK= if_else(potassium <3.5, 1, 0))


geelm.control(std.err ="san.se")
fit_potald <-  geeglm(formula= formule1, 
                      data=potald, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)

summary(fit_potald)


potald$residuals_pot <-fit_potald$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()


potald_min1 <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~.-1)

summary(potald_min1)

fit_potald_min1 <-  geeglm(formula= formule1, 
                      data=potald_min1, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_min1)

potald_min1$residuals_potmin1 <-fit_potald_min1$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_min1_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_min1, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()


potald_1div <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~1-(1/.))
summary(potald_1div)

potald_1div[potald_1div == -Inf] <- 1

summary(potald_1div)

fit_potald_1div<-  geeglm(formula= formule1, 
                      data=potald_1div, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_1div)

potald_1div$residuals_pot1div <-fit_potald_1div$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_1div_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_1div, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()

plot_scatter(t1,age,potassium,fit.line ="loess")



potald_2exp <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~0.5-2**(-.))
summary(potald_2exp)

summary(potald_2exp)

fit_potald_2exp<-  geeglm(formula= formule1, 
                          data=potald_2exp, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_2exp)

potald_2exp$residuals_pot1div <-fit_potald_2exp$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_2exp_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_2exp, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()

plot_scatter(t1,age,potassium,fit.line ="loess")

#### GEE without spline & log predictors ####



potald <- potald %>% 
  mutate(latest_ald_log = log10(latest_ald))


formule2 <- potassium ~ 
  age + sex + in_icu + days_since_hospital+
  latest_ald_log + latest_ald_log:days_since_ald +
  latest_pot + latest_pot:days_since_pot+
  latest_arb+latest_arb:days_since_arb+
  latest_acei+latest_acei:days_since_acei+
  latest_mra+latest_mra:days_since_mra+
  latest_loop_diuretic+latest_loop_diuretic:days_since_loop_diuretic+
  latest_thiazid+latest_thiazid:days_since_thiazid+
  latest_pot_flush+latest_pot_flush:days_since_pot_flush+
  latest_pot_supp+latest_pot_supp:days_since_pot_supp+
  latest_cat + latest_cat:days_since_cat

    


geelm.control(std.err ="san.se")
fit_potald_log <-  geeglm(formula= formule2, 
                      data=potald, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)

summary(fit_potald_log)


potald$residuals_pot <-fit_potald_log$residuals

x=1

predictors <- attr(terms(fit_potald_log),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_log_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()


potald_min1_log <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~.-1)

summary(potald_min1_log)

fit_potald_min1_log <-  geeglm(formula= formule2, 
                           data=potald_min1_log, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_min1_log)

potald_min1_log$residuals_potmin1 <-fit_potald_min1_log$residuals

x=1

predictors <- attr(terms(fit_potald_min1_log),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_min1_log_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_min1_log, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()


potald_1div_log <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~1-(1/.))
summary(potald_1div_log)

potald_1div_log[potald_1div_log == -Inf] <- 1

summary(potald_1div_log)

fit_potald_1div_log<-  geeglm(formula= formule2, 
                          data=potald_1div_log, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_1div_log)

potald_1div_log$residuals_pot1div <-fit_potald_1div_log$residuals

x=1

predictors <- attr(terms(fit_potald_1div_log),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_1div_log_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_1div_log, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()




potald_2exp_log <- potald %>% 
  mutate_at(c("days_since_ald", "days_since_pot", "days_since_arb", "days_since_acei",
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid",
              "days_since_pot_flush","days_since_pot_supp"),~0.5-2**(-.))


summary(potald_2exp_log)


summary(potald_2exp_log)

fit_potald_2exp_log<-  geeglm(formula= formule2, 
                              data=potald_2exp_log, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_2exp_log)

potald_2exp_log$residuals_pot1div <-fit_potald_2exp_log$residuals

x=1

predictors <- attr(terms(fit_potald_2exp_log),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/potald_2exp_log_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=potald_2exp_log, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()

library(sjPlot)
library(sjmisc)
library(sjlabelled)

QIC(fit_potald)
QIC(fit_potald_log)
QIC(fit_potald_min1)
QIC(fit_potald_min1_log)
QIC(fit_potald_1div)
QIC(fit_potald_1div_log)
QIC(fit_potald_2exp)
QIC(fit_potald_2exp_log)

tab_model(fit_potald,fit_potald_log,fit_potald_min1,fit_potald_min1_log,fit_potald_1div,fit_potald_1div_log,fit_potald_2exp,fit_potald_2exp_log,
          dv.labels = c( "Regular Elapsed Days",
                         "Regular Elapsed Days and Log of Alsdosterone","Elapsed Days Minus One", "Elapsed Days Minus One and Log of Aldosterone", 
                        "One Over Elapsed Days", "One Over Elapsed Days and Log of Aldosterone", "Exponentiated Elapsed Days", "Exponentiated Elapsed Days and Log of Aldosterone"), digits = 5,
          file = "output/Elapsed Time Transformation and Log of Aldosterone GEE Results.html", linebreak = TRUE)

#### combined testing of covariate effects ####

vc.pot <- as.matrix(vcov(fit_potald_2exp_log)[6:7, 6:7])
coeff.pot <- coef(fit_potald_2exp_log)[6:7]
stat <- coeff.pot %*% solve(vc.pot) %*% coeff.pot
(pval.sim.pot <- pchisq(stat, 2, lower.tail=FALSE))


#### predictor distributions ####


IDA = potald %>% mutate(pot_int = latest_pot*days_since_pot,
                        ald_int = latest_ald*days_since_ald,
                        ald_int_log = latest_ald_log*days_since_ald,
                        arb_int = latest_arb*days_since_arb,
                        acei_int = latest_acei*days_since_acei,
                        mra_int = latest_mra*days_since_mra,
                        loop_diuretic_int = latest_loop_diuretic*days_since_loop_diuretic,
                        thiazid_int = latest_thiazid*days_since_thiazid,
                        pot_flush_int = latest_pot_flush*days_since_pot_flush,
                        pot_supp_int = latest_pot_supp*days_since_pot_supp) %>% 
  select(age,days_since_hospital,potassium,latest_pot,days_since_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log,arb_int,
         acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int)


logIDA = IDA %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
logIDA2 = logIDA %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                    ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
dev.off()

pdf(file="output/Predictor Distributions.pdf")
for (x in c(1:length(colnames(IDA)))){
  xa = colnames(IDA)[x]
  try(print(ggplot(data=data.frame(IDA), aes_string(x = xa)) +
              geom_histogram(data=data.frame(IDA)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(IDA)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=colnames(IDA)[x],x=colnames(IDA)[x], y = "Density")))
  try(print(ggplot(data=data.frame(IDA), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=colnames(IDA)[x],x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Log of ",colnames(IDA)[x]),x=paste0("Log of ",colnames(IDA)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Log of ",colnames(IDA)[x]),x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA2), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA2)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA2)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Double Log of ",colnames(IDA)[x]),x=paste0("Double Log of ",colnames(IDA)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA2), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Double Log of ",colnames(IDA)[x]),x="Theoretical", y = "Sample")))
}
dev.off()


IDA_1div = potald_1div_log %>%mutate(pot_int = latest_pot*days_since_pot,
                                 ald_int = latest_ald*days_since_ald,
                                 ald_int_log = latest_ald_log*days_since_ald,
                                 arb_int = latest_arb*days_since_arb,
                                 acei_int = latest_acei*days_since_acei,
                                 mra_int = latest_mra*days_since_mra,
                                 loop_diuretic_int = latest_loop_diuretic*days_since_loop_diuretic,
                                 thiazid_int = latest_thiazid*days_since_thiazid,
                                 pot_flush_int = latest_pot_flush*days_since_pot_flush,
                                 pot_supp_int = latest_pot_supp*days_since_pot_supp) %>% 
  select(age,days_since_hospital,potassium,latest_pot,days_since_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log,arb_int,
         acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int)


logIDA_1div = IDA_1div %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                          ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
logIDA2_1div = logIDA_1div %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                              ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
dev.off()

pdf(file="output/Predictor Distributions 1div.pdf")
for (x in c(1:length(colnames(IDA_1div)))){
  xa = colnames(IDA_1div)[x]
  try(print(ggplot(data=data.frame(IDA_1div), aes_string(x = xa)) +
              geom_histogram(data=data.frame(IDA_1div)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(IDA_1div)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=colnames(IDA_1div)[x],x=colnames(IDA_1div)[x], y = "Density")))
  try(print(ggplot(data=data.frame(IDA_1div), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=colnames(IDA_1div)[x],x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA_1div), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA_1div)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA_1div)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Log of ",colnames(IDA_1div)[x]),x=paste0("Log of ",colnames(IDA_1div)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA_1div), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Log of ",colnames(IDA)[x]),x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA2_1div), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA2_1div)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA2_1div)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Double Log of ",colnames(IDA_1div)[x]),x=paste0("Double Log of ",colnames(IDA_1div)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA_1div2), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Double Log of ",colnames(IDA_1div)[x]),x="Theoretical", y = "Sample")))
}
dev.off()


IDA_min1 = potald_min1_log %>% mutate(pot_int = latest_pot*days_since_pot,
                                  ald_int = latest_ald*days_since_ald,
                                  ald_int_log = latest_ald_log*days_since_ald,
                                  arb_int = latest_arb*days_since_arb,
                                  acei_int = latest_acei*days_since_acei,
                                  mra_int = latest_mra*days_since_mra,
                                  loop_diuretic_int = latest_loop_diuretic*days_since_loop_diuretic,
                                  thiazid_int = latest_thiazid*days_since_thiazid,
                                  pot_flush_int = latest_pot_flush*days_since_pot_flush,
                                  pot_supp_int = latest_pot_supp*days_since_pot_supp) %>% 
  select(age,days_since_hospital,potassium,latest_pot,days_since_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log,arb_int,
         acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int)

logIDA_min1 = IDA_min1 %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                          ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
logIDA2_min1 = logIDA_min1 %>% mutate_at(vars(age,days_since_hospital,potassium,days_since_pot,latest_pot,pot_int,latest_ald,days_since_ald,ald_int,ald_int_log
                                              ,arb_int,acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log10)
dev.off()

pdf(file="output/Predictor Distributions min1.pdf")
for (x in c(1:length(colnames(IDA_min1)))){
  xa = colnames(IDA_min1)[x]
  try(print(ggplot(data=data.frame(IDA_min1), aes_string(x = xa)) +
              geom_histogram(data=data.frame(IDA_min1)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(IDA_min1)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=colnames(IDA_min1)[x],x=colnames(IDA_min1)[x], y = "Density")))
  try(print(ggplot(data=data.frame(IDA_min1), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=colnames(IDA_min1)[x],x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA_min1), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA_min1)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA_min1)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Log of ",colnames(IDA_min1)[x]),x=paste0("Log of ",colnames(IDA_min1)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA_min1), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Log of ",colnames(IDA)[x]),x="Theoretical", y = "Sample")))
  try(print(ggplot(data=data.frame(logIDA2_min1), aes_string(x = xa)) +
              geom_histogram(data=data.frame(logIDA2_min1)[x], aes(y=..density..), colour="black", fill="white")+
              geom_density(data=data.frame(logIDA2_min1)[x], aes(y=..density..), alpha=.2, fill="#FF6666")+
              labs(title=paste0("Double Log of ",colnames(IDA_min1)[x]),x=paste0("Double Log of ",colnames(IDA_min1)[x]), y = "Density")))
  try(print(ggplot(data=data.frame(logIDA2_min1), aes_string(sample = xa)) +
              geom_qq()+ geom_qq_line()+
              labs(title=paste0("Double Log of ",colnames(IDA_min1)[x]),x="Theoretical", y = "Sample")))
}
dev.off()

geelm.control(std.err ="san.se")
fit_potald <-  geeglm(formula= formule1, 
                      data=potald, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)

summary(fit_potald)
plot(fit_potald)

plotres(fit_potald)
plotmo(fit_potald)




geelm.control(std.err ="san.se")
ar1_fit_potald <-  geeglm(formula= formule1, 
                          data=regpotald, family= gaussian, id=ID, corstr="ar1", waves=days_since_hospital)
ind_fit_potald <-  geeglm(formula= formule1, 
                          data=regpotald, family= gaussian, id=ID, corstr="independence", waves=days_since_hospital)
exc_fit_potald <-  geeglm(formula= formule1, 
                          data=regpotald, family= gaussian, id=ID, corstr="exchangeable", waves=days_since_hospital)



QIC(ar1_fit_potald)
QIC(ind_fit_potald)
QIC(exc_fit_potald)


#### Outlier Removal ####

outliers20 <- potald[-c(43,315,423),][-c(1,418),][-c(126,363),][-c(417),][-c(437,196,278),][-c(422,424),][-c(108,43,201),][-c(72,338),][-c(53,92),]
geelm.control(std.err ="san.se")
fit_potald <-  geeglm(formula= formule1, 
                      data=outliers20, family= gaussian, id=ID, corstr="ar1", waves=days_since_hospital)


summary(fit_potald)
plot(fit_potald)

plotres(fit_potald)


hist(fit_potald$residuals)


potald$residuals_pot <-fit_potald$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/20_outliers removed_pot_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=outliers20, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()






#### GEE before ICU without splines ####


before_ICU <- potald %>%
  mutate(after_ICU = ifelse(in_icu == 1, 1, NA)) %>% 
  group_by(ID) %>% 
  fill(after_ICU,.direction = "down") %>% 
  ungroup() %>% 
  mutate(after_ICU = replace_na(after_ICU,0)) %>% 
  filter(after_ICU==0)


geelm.control(std.err ="san.se")
fit_before_ICU <-  geeglm(formula= formule1, 
                      data=before_ICU, family= gaussian, id=ID, corstr="ar1", waves=days_since_hospital)


summary(fit_before_ICU)
plot(fit_before_ICU)

plotres(fit_before_ICU)


hist(fit_before_ICU$residuals)


before_ICU$residuals_pot <-fit_before_ICU$residuals

x=1

predictors <- attr(terms(fit_before_ICU),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="output/pot_residuals_without_splines.pdf")
for (x in 1:length(predictors)){
  xa = predictors[x]
  print(
    ggplot(data=fit_before_ICU, aes_string(x=xa,y="residuals_pot")) +
      geom_point(colour="black", fill="white") +
      geom_smooth() +
      labs(y="Residuals Potassium",
           x=predictors[x])
  )
}
dev.off()

ggplot(data=fit_before_ICU, aes(x=in_icu, y=pot_supp))+
  geom_point()+
  geom_smooth()




# 
# ## GEE with splines
# 
# formule2 <- potassium ~  ns(latest_ald,knots=c(10,50,75), Boundary.knots = c(5,100)) + 
#   ns(latest_ald,knots=c(10,50,75), Boundary.knots = c(5,100)):ns(days_since_ald,Boundary.knots= c(1,6), knots=c(3,4)) + 
#   ns(latest_pot, Boundary.knots = c(3,4.5), knots=c(3.5,4)) +
#   ns(latest_pot, Boundary.knots = c(3,4.5), knots=c(3.5,4)):ns(days_since_pot, Boundary.knots = c(1,3),knots=2) + 
#   days_since_hospital + ns(age, Boundary.knots = c(36, 85), knots=c(50,66,78)) + sex + pot_supp + in_icu
# 
# geelm.control(std.err ="san.se")
# fit_potald_spl <-  geeglm(formula= formule2, 
#                       data=potald, family= gaussian, id=ID, corstr="ar1", waves=days_since_hospital)
# 
# summary(fit_potald_spl)
# plot(fit_potald_spl)
# 
# 
# 
# plotres(fit_potald_spl)
# 
# 
# hist(fit_potald_spl$residuals)
# 
# 
# potald$residuals_pot <-fit_potald_spl$residuals
# 
# dev.off()
# pdf(file="output/pot_residuals_with_splines.pdf")
# for (x in 1:length(colnames(potald))){
#   if (!colnames(potald)[x] %in% c("pat_id")){
#     print(
#       ggplot(data=potald, aes(x=unlist(potald[,x]),y=residuals_pot)) +
#         geom_point(colour="black", fill="white") +
#         geom_smooth(aes(y=residuals_pot)) +
#         labs(y="Residuals Potassium",
#              x=ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x]),
#              title = ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x])
#         )
#     )
#   }
# }
# dev.off()
# 
# 
# ## GEE Hypokalemia without splines
# 
# potald <- potald %>% 
#   mutate(hypoK= if_else(potassium <3.5, 1, 0))
# 
# 
# 
# formule3 <- hypoK ~  latest_ald + latest_ald:days_since_ald + latest_pot + latest_pot:days_since_pot + days_since_hospital + age + sex + pot_supp
# 
# fit_hypoK <- geeglm(formula= formule3, 
#                     data=potald, family= binomial(link='logit'), id=ID, corstr="ar1", waves=days_since_hospital)
# 
# summary(fit_hypoK)
# plot(fit_hypoK)
# 
# plotres(fit_hypoK)
# 
# potald$residuals_hypoK <-fit_hypoK$residuals
# 
# dev.off()
# pdf(file="output/hypoK_residuals.pdf")
# for (x in 1:length(colnames(potald))){
#   if (!colnames(potald)[x] %in% c("pat_id")){
#     print(
#       ggplot(data=potald, aes(x=unlist(potald[,x]),y=residuals_hypoK)) +
#         geom_point(colour="black", fill="white") +
#         geom_smooth(aes(y=residuals_hypoK)) +
#         labs(y="Residuals Hypokalemia",
#              x=ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x]),
#              title = ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x])
#         )
#     )
#   }
# }
# dev.off()
# 
# 
# ## Hypokalemia GEE with splines
# 
# potald <- potald %>% 
#   mutate(hypoK= if_else(potassium <3.5, 1, 0))
# 
# 
# 
# formule4 <- hypoK ~  ns(latest_ald,knots=c(10,50,75), Boundary.knots = c(5,100)) + ns(latest_ald,knots=c(10,50,75), Boundary.knots = c(5,100)):ns(days_since_ald,Boundary.knots= c(1,6), knots=c(3,4)) + ns(latest_pot, Boundary.knots = c(3,4.5), knots=c(3.5,4)) +
#   ns(latest_pot, Boundary.knots = c(3,4.5), knots=c(3.5,4)):ns(days_since_pot, Boundary.knots = c(1,3),knots=2) + days_since_hospital + ns(age, Boundary.knots = c(36, 85), knots=c(50,66,78)) + sex + pot_supp
# 
# fit_hypoK_spl <- geeglm(formula= formule4, 
#                     data=potald, family= binomial(link='logit'), id=ID, corstr="ar1", waves=days_since_hospital)
# 
# summary(fit_hypoK_spl)
# plot(fit_hypoK_spl)
# 
# plotres(fit_hypoK_spl)
# 
# potald$residuals_hypoK <-fit_hypoK_spl$residuals
# 
# dev.off()
# pdf(file="output/hypoK_residuals_with_splines.pdf")
# for (x in 1:length(colnames(potald))){
#   if (!colnames(potald)[x] %in% c("pat_id")){
#     print(
#       ggplot(data=potald, aes(x=unlist(potald[,x]),y=residuals_hypoK)) +
#         geom_point(colour="black", fill="white") +
#         geom_smooth(aes(y=residuals_hypoK)) +
#         labs(y="Residuals Hypokalemia",
#              x=ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x]),
#              title = ifelse(label(potald)[x]=="", colnames(potald)[x],label(potald)[x])
#         )
#     )
#   }
# }
# dev.off()
# 
# 
# 
# 
# scatter.smooth(log10(potald$latest_ald), fit_potald_spl$residuals)
# qqplot(log10(potald$latest_ald), fit_potald_spl$residual)
# scatter.smooth(log10(potald$latest_ald*potald$days_since_ald), fit_potald_spl$residuals)
# qqplot(log10(potald$latest_ald*potald$days_since_ald), fit_potald_spl$residuals)
# 
# scatter.smooth(log10(potald$latest_pot), fit_potald_spl$residuals)
# qqplot(log10(potald$latest_pot), fit_potald_spl$residual)
# scatter.smooth(log10(potald$latest_pot*potald$days_since_ald), fit_potald_spl$residuals)
# qqplot(log10(potald$latest_pot*potald$days_since_ald), fit_potald_spl$residuals)
# 
# 
# 
# scatter.smooth(potald$latest_pot, fit_potald_spl$residuals)
# 
##### construction site - code ideas not for direct use
# 
# a=0
# medslong = tibble(colnames=unique(meds$drug_class))
# medslong = read_csv("\n", col_names = unique(meds$drug_class[!is.na(meds$drug_class)])) %>% 
#   add_column(pat_id=0,date_measurement=0)
# 
# for (row in 1:nrow(meds)){
#   medslong
# }
# 
# 
# do.call(rbind, 
#         lapply(seq(nrow(meds)), function(x){
#           data.frame(pat_id=meds[x,"pat_id"],
#                      drug_class=meds[x,"drug_class"],
#                      date_measurement=seq(meds[x,"start_date"], 
#                              meds[x,"end_date"],by="day")
#           )}))
# 
# dates <- meds2 %>%
#   rowwise() %>%
#   transmute(date_measurement = list(seq(as.Date(as.character(start_date), origin="1899-12-30"), as.Date(as.character(end_date), origin="1899-12-30"), by = "day"))) %>%
#   unnest(date_measurement)
# 
# unique(dates)
# 
# 
# test %>%
#   group_by(idnum) %>%
#   summarize(start=min(start),end=max(end)) %>%
#   do(data.frame(idnum=.$idnum, month=seq(.$start,.$end,by="1 day")))
# 
  