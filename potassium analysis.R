

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
         gender, 
         age, 
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
          date_death_asdate = as.Date(as.numeric(date_death), origin = "1899-12-30"))

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
  select(ald = Aldosterone, pras = `PRA-S`,
         measurement_id)%>%
  full_join(measurement_ids) %>% 
  select(pras, ald, date_measurement, pat_id) %>% 
  filter(complete.cases(.))


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
         never_icu = is.na(in_icu)*1,
         icu_or_death = ifelse(!is.na(date_death), 1, ifelse(never_icu == 0, 1,0)),
         in_icu = ifelse(is.na(in_icu),0,in_icu),
         days_from_icu_admission = as.numeric(date_measurement - date_icu_asdate)
         ) %>% 
  full_join(insulin) %>%
  full_join(meds) %>% 
  full_join(ras_lab) %>% 
  group_by(pat_id) %>% 
  #filter(ald>10) %>% 
  mutate_at(vars(latest_pot_date = date_measurement, 
                 latest_pot = potassium, 
                 latest_ald = ald,
                 latest_arb = arb,
                 latest_acei = acei,
                 latest_mra = mra,
                 latest_loop_diuretic = loop_diuretic,
                 latest_thiazid = thiazid_diuretic,
                 latest_pot_flush = pot_flush,
                 latest_pot_supp = pot_supp
                 ), lag) %>%
  mutate_at(vars(latest_arb,
                 latest_acei,
                 latest_mra,
                 latest_loop_diuretic,
                 latest_thiazid,
                 latest_pot_flush,
                 latest_pot_supp),
            ~na_if(.,y=0)) %>% 
  mutate(latest_ald_date = as.Date(ifelse(!is.na(latest_ald),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_arb_date = as.Date(ifelse(!is.na(latest_arb),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_acei_date = as.Date(ifelse(!is.na(latest_acei),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_mra_date = as.Date(ifelse(!is.na(latest_mra),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_loop_diuretic_date = as.Date(ifelse(!is.na(latest_loop_diuretic),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_thiazid_date = as.Date(ifelse(!is.na(latest_thiazid),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_pot_flush_date = as.Date(ifelse(!is.na(latest_pot_flush),lag(date_measurement),NA), origin = "1970-01-01"),
         latest_pot_supp_date = as.Date(ifelse(!is.na(latest_pot_supp),lag(date_measurement),NA), origin = "1970-01-01"),
         ) %>% 
  fill(latest_ald_date, latest_ald,
       latest_arb,latest_arb_date,
       latest_acei,latest_acei_date,
       latest_mra,latest_mra_date,
       latest_loop_diuretic,latest_loop_diuretic_date,
       latest_thiazid,latest_thiazid_date,
       latest_pot_flush,latest_pot_flush_date,
       latest_pot_supp,latest_pot_supp_date,
       .direction = "down") %>% 
  mutate(days_since_ald = as.numeric(date_measurement- latest_ald_date),
         days_since_pot = as.numeric(date_measurement - latest_pot_date),
         days_since_arb = as.numeric(date_measurement - latest_arb_date),
         days_since_acei = as.numeric(date_measurement - latest_acei_date),
         days_since_mra = as.numeric(date_measurement - latest_mra_date),
         days_since_loop_diuretic = as.numeric(date_measurement - latest_loop_diuretic_date),
         days_since_thiazid = as.numeric(date_measurement - latest_thiazid_date),
         days_since_pot_flush = as.numeric(date_measurement - latest_pot_flush_date),
         days_since_pot_supp = as.numeric(date_measurement - latest_pot_supp_date),
         pot_change = potassium - latest_pot) %>% 
  ungroup %>% 
  replace_na(replace=list(days_since_arb = 0,
             days_since_acei = 0,
             days_since_mra = 0,
             days_since_loop_diuretic = 0,
             days_since_thiazid = 0,
             days_since_pot_flush = 0,
             days_since_pot_supp = 0,
             latest_arb = 0,
             latest_acei = 0,
             latest_mra = 0,
             latest_loop_diuretic = 0,
             latest_thiazid = 0,
             latest_pot_flush = 0,
             latest_pot_supp = 0,
             arb = 0,
             acei = 0,
             mra = 0,
             loop_diuretic = 0,
             thiazid = 0,
             pot_flush = 0,
             pot_supp = 0)) %>% 
  mutate(insulin = if_else(insulin=="yes",1,0,missing = 0),
         five_days_since_hospital = ceiling((days_since_hospital+1)/ 5),
         seven_days_since_hospital = ceiling((days_since_hospital+1)/ 7)) %>% 
#  mutate_at(vars(latest_arb,latest_acei,latest_mra,latest_thiazid,latest_loop_diuretic,latest_pot_flush,latest_pot_supp),as.factor) %>% 
  filter(days_since_hospital %in% 0:20)%>% 
  relocate(any_of(c("pat_id", "date_days","never_icu","in_icu","date_icu_asdate","date_icu_discharge_asdate","date_hospital_discharge_asdate","date_death_asdate","num_date_death","date_death","end_of_icu","pat_id", "days_since_hospital", "seven_days_since_hospital","days_since_ald", "days_since_pot", "latest_pot_date", "latest_pot", "ald", "latest_ald", "latest_ald_date" )), .after=potassium) %>% 
  add_labels(never_icu, labels=c("No ICU"=1, "ICU"=0))
  
potaldvars <- c("pat_id","PRA-S","in_icu","never_icu","pot_change","potassium","latest_ald","latest_pot","days_since_ald","days_since_pot","gender","age","days_since_hospital", "pot_supp", "days_since_arb", "days_since_acei","days_since_mra", "days_since_loop_diuretic", "days_since_thiazid", "days_since_pot_flush", "days_since_pot_supp", "latest_arb", "latest_acei", "latest_mra", "latest_loop_diuretic", "latest_thiazid", "latest_pot_flush", "latest_pot_supp")



names(potassium_df)

## vars "pat_id","potassium","latest_ald","latest_pot","latest_ald_date","latest_pot_date","gender","age","days_since_hospital"

## meds "insulin","beta-blocker","alpha-blocker","calcium channel blocker","alpha-2-agonist","kalzium antagonist",
##      "potassium channel opener","alpha1-adrenoceptor-antagonist" "arb","acei"m"loop diuretic","thiazid diuretic","carbonic anhydrase inhibitor"
##      "sulfonamide-based diuretic","mra","catecholamine","pot_supp","pot_flush","calc","phos","mag"   




length(unique(potassium_df$pat_id))

potald <- potassium_df %>% 
  select(all_of(potaldvars)) %>% 
  drop_na() %>% 
  mutate(ID = as.numeric(factor(pat_id)))
#potald$latest_ald[potald$latest_ald == 10] = 5

hist(log(potald$latest_ald))
qqnorm(log(potald$latest_ald))
hist(log(potald$latest_pot))
qqnorm(log(potald$latest_pot))
hist(potald$pot_change)
qqnorm(potald$pot_change)
hist(log(potald$potassium))
qqnorm(log(potald$potassium))
hist(log(potald$latest_ald*potald$days_since_ald))
qqnorm(log(potald$latest_ald*potald$days_since_ald))

hist(potald$days_since_hospital)
qqnorm(potald$days_since_hospital)

phist<- potald %>% ggplot(aes(x=log(potassium))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
pqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log(potassium))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(phist, pqq, labels = "AUTO")


lahist<- potald %>% ggplot(aes(x=log(latest_ald+5))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
laqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log(latest_ald+5))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(lahist, laqq, labels = "AUTO")
lphist <- potald %>% ggplot(aes(x=log(latest_pot))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
lpqq<- pcqq<-ageqq<-potald %>% ggplot(aes(sample=log(latest_pot))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(lphist, lpqq, labels = "AUTO")

pchist<- potald %>% ggplot(aes(x=pot_change)) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
pcqq<-ageqq<-potald %>% ggplot(aes(sample=log(pot_change))) +
  stat_qq() + 
  stat_qq_line(col="red")
plot_grid(pchist, pcqq, labels = "AUTO")

dsahist<- potald %>% ggplot(aes(x=log(days_since_ald))) +
  geom_histogram(aes(y=..density..),color=1,fill="grey") +
  geom_density()
dsaqq<-ageqq<-potald %>% ggplot(aes(sample=log(days_since_ald))) +
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
  geom_bar(aes(x=gender))+
  labs(title="Patients")
complete_observation_sex <- potald %>%ggplot()+
  geom_bar(aes(x=gender))+
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
logIDA = data.frame(lapply(IDA, log2))
logIDA2 = data.frame(lapply(logIDA, log2))
dev.off()
pdf(file="log and nonlog.pdf")
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

t1df <- potassium_df %>% group_by(pat_id) %>% mutate(Mean_Potassium = mean(potassium,na.rm=TRUE), Mean_Aldosterone = mean(ald,na.rm=TRUE))%>% filter(row_number()==1) %>% ungroup %>%  
  transmute(Age = age, Sex = factor(gender, levels=c("f","m"), labels=c("Female","Male")), BMI = bmi, Diabetes= factor(diabetes, levels=c("yes","no"), labels=c("Yes","No")),
            Hypertension = factor(hypertension, levels=c("yes","no"), labels=c("Yes","No")),COPD = factor(copd, levels=c("yes","no"), labels=c("Yes","No")), pat_id, never_icu = factor(never_icu, levels=c(0, 1), labels=c("ICU", "No ICU")), Potassium, Aldosterone) 


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

t1 <- table1(~ Sex + Age + BMI + Diabetes + Hypertension + COPD + Aldosterone + Potassium| never_icu, data=t1df, render.continuous="Median [Min, Max]")

write.table(t1,file="tables/T1.docx")

#### Smooth Plots ####


total_aldplot <- potassium_df %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth() +
  ylab("Aldosterone")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(-100,250,by=50), limits=c(-50,100))

total_potplot <- potassium_df%>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(color="red")+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

test_total_aldplot <- potassium_df %>% 
  ggplot(aes(x=days_since_test,y=ald)) +
  theme_cowplot(12)+
  geom_smooth() +
  ylab("Aldosterone")+
  xlab("Days since test")+
  scale_y_continuous(breaks=seq(-100,250,by=50), limits=c(-50,100))

test_total_potplot <- potassium_df%>% 
  ggplot(aes(x=days_since_test,y=potassium)) +  
  geom_smooth(color="red")+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since test")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aldplot <- potassium_df[potassium_df$never_icu==1,] %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth() +
  ylab("Aldosterone")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(-100,250,by=50), limits=c(-50,100))

potplot <- potassium_df[potassium_df$never_icu==1,] %>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(color="red")+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

aldICUplot <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=ald)) +
  theme_cowplot(12)+
  geom_smooth() +
  ylab("Aldosterone")+
  xlab("Days from ICU admission")+
  scale_y_continuous(breaks=seq(-100,250,by=50), limits=c(-50,100))

potICUplot <- potassium_df %>% 
  ggplot(aes(x=days_from_icu_admission,y=potassium)) +  
  geom_smooth(color="red")+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days from ICU admission")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))


ald_everICU_plot <- potassium_df[potassium_df$never_icu==0,] %>% 
  ggplot(aes(x=days_since_hospital,y=ald)) +
  theme_cowplot(12)+
  geom_smooth() +
  ylab("Aldosterone")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(-100,250,by=50), limits=c(-50,100))

pot_everICU_plot <- potassium_df[potassium_df$never_icu==0,] %>% 
  ggplot(aes(x=days_since_hospital,y=potassium)) +  
  geom_smooth(color="red")+
  theme_cowplot(12) +
  ylab("Potassium (mmol/L)")+
  xlab("Days since hospitalization")+
  scale_y_continuous(breaks=seq(3.4,4,by=.2), limits = c(3.4, 4.2))

plot_grid(total_potplot,total_aldplot,test_total_potplot,test_total_aldplot,potplot,aldplot,pot_everICU_plot,ald_everICU_plot,potICUplot,aldICUplot, labes= c("A","B"), label_size = 12, ncol=2)


#### Potassium by Aldosterone per 5 day Interval ####

scx = scale_x_log10(limits=c(9,400), n.breaks=10)
scy = scale_y_continuous(breaks=seq(2.5,5.6,by=0.5), limits=c(2.5,5.6))
theme = theme_cowplot(12)

fullscatter <-potassium_df %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

fullscatter_latest <- potassium_df%>% 
  filter(days_since_ald <2) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter1 <-potassium_df %>% 
  filter(seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest1 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 1) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter2 <-potassium_df %>% 
  filter(seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest2 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 2)%>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx + 
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter3 <-potassium_df %>% 
  filter(seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme +geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latest3 <- potassium_df %>% 
  filter(days_since_ald <2,
         seven_days_since_hospital == 3) %>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatterICU <-potassium_df %>% 
  filter(in_icu == 1) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latestICU <- potassium_df %>% 
  filter(days_since_ald <2,
         in_icu == 1)%>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatterICU <-potassium_df %>% 
  filter(in_icu == 0) %>% 
  ggplot(aes(y=potassium, x=ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Same Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 

scatter_latestnoICU <- potassium_df %>% 
  filter(days_since_ald <2,
         in_icu == 0)%>% 
  ggplot(aes(y=potassium, x=latest_ald))+
  geom_point()+
  ylab("Potassium (mmol/L)")+
  xlab("Previous Day Aldosterone (pg/mL, Logarithmic Scale)")+
  scy +
  scx +
  theme + geom_vline(aes(xintercept=10),color="red", linetype = "dashed") + annotation_logticks(sides="b") 


(aldpotscat <- plot_grid(fullscatter,fullscatter_latest,scatter1,scatter_latest1,scatter2,scatter_latest2,scatter3,scatter_latest3, label_size = 12, ncol=2, labels = "AUTO"))

save_plot('figures/aldpotscat.png', aldpotscat, ncol=2, nrow=4)



#### GEE without splines ####


formule1 <- potassium ~ 
  age + gender + in_icu + days_since_hospital+
  latest_ald + latest_ald:days_since_ald +
  latest_pot + latest_pot:days_since_pot+ 
  latest_arb+latest_arb:days_since_arb+
  latest_acei+latest_acei:days_since_acei+
  latest_mra+latest_mra:days_since_mra+
  latest_loop_diuretic+latest_loop_diuretic:days_since_loop_diuretic+
  latest_thiazid+latest_thiazid:days_since_thiazid+
  latest_pot_flush+latest_pot_flush:days_since_pot_flush+
  latest_pot_supp+latest_pot_supp:days_since_pot_supp



potald <- potald %>% 
  mutate(hypoK= if_else(potassium <3.5, 1, 0))

hist(potald$latest_ald)
plot(potald)




geelm.control(std.err ="san.se")
fit_potald <-  geeglm(formula= formule1, 
                      data=potald, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)

summary(fit_potald)


potald$residuals_pot <-fit_potald$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="potald_residuals_without_splines.pdf")
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
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid","days_since_pot_flush","days_since_pot_supp"),~.-1)

summary(potald_min1)

fit_potald_min1 <-  geeglm(formula= formule1, 
                      data=potald_min1, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_min1)

potald_min1$residuals_pot <-fit_potald_min1$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="potald_min1_residuals_without_splines.pdf")
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
              "days_since_mra","days_since_loop_diuretic","days_since_thiazid","days_since_pot_flush","days_since_pot_supp"),~1/.)
summary(potald_1div)

potald_1div[potald_1div == Inf] <- 1
summary(potald_1div)

fit_potald_1div<-  geeglm(formula= formule1, 
                      data=potald_1div, family="gaussian", id=ID, corstr="ar1", waves=days_since_hospital)
summary(fit_potald_1div)

potald_1div$residuals_pot <-fit_potald_1div$residuals

x=1

predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")

dev.off()
pdf(file="potald_1div_residuals_without_splines.pdf")
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


library(sjPlot)
library(sjmisc)
library(sjlabelled)

QIC(fit_potald_min1)
QIC(fit_potald_1div)
QIC(fit_potald)

tab_model(fit_potald_min1,fit_potald_1div,fit_potald,
          dv.labels = c("Elapsed Days Minus One", "One Over Elapsed Days", "Regular Elapsed Days"), digits = 5,
          file = "Elapsed Time Transformation GEE Results.html")


predictors <- attr(terms(fit_potald),"term.labels") %>% str_replace_all(":","*")


IDA = potald %>% mutate(pot_int = latest_pot*days_since_pot,
                        ald_int = latest_ald*days_since_ald,
                        arb_int = latest_arb*days_since_arb,
                        acei_int = latest_acei*days_since_acei,
                        mra_int = latest_mra*days_since_mra,
                        loop_diuretic_int = latest_loop_diuretic*days_since_loop_diuretic,
                        thiazid_int = latest_thiazid*days_since_thiazid,
                        pot_flush_int = latest_pot_flush*days_since_pot_flush,
                        pot_supp_int = latest_pot_supp*days_since_pot_supp) %>% 
  select(potassium,latest_pot,latest_ald,age,days_since_hospital,pot_int,ald_int,arb_int,
         acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int)


logIDA = IDA %>% mutate_at(vars(potassium,latest_pot,latest_ald,age,days_since_hospital,pot_int,ald_int,arb_int,
                                     acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log2)
logIDA2 = logIDA %>% mutate_at(vars(potassium,latest_pot,latest_ald,age,days_since_hospital,pot_int,ald_int,arb_int,
                                 acei_int,mra_int,loop_diuretic_int,thiazid_int,pot_flush_int,pot_supp_int),log2)
pdf(file="log and nonlog.pdf")
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

## Outlier Removal

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
pdf(file="20_outliers removed_pot_residuals_without_splines.pdf")
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






## GEE before ICU without splines


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
pdf(file="pot_residuals_without_splines.pdf")
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
#   days_since_hospital + ns(age, Boundary.knots = c(36, 85), knots=c(50,66,78)) + gender + pot_supp + in_icu
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
# pdf(file="pot_residuals_with_splines.pdf")
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
# formule3 <- hypoK ~  latest_ald + latest_ald:days_since_ald + latest_pot + latest_pot:days_since_pot + days_since_hospital + age + gender + pot_supp
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
# pdf(file="hypoK_residuals.pdf")
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
#   ns(latest_pot, Boundary.knots = c(3,4.5), knots=c(3.5,4)):ns(days_since_pot, Boundary.knots = c(1,3),knots=2) + days_since_hospital + ns(age, Boundary.knots = c(36, 85), knots=c(50,66,78)) + gender + pot_supp
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
# pdf(file="hypoK_residuals_with_splines.pdf")
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
# scatter.smooth(log(potald$latest_ald), fit_potald_spl$residuals)
# qqplot(log(potald$latest_ald), fit_potald_spl$residual)
# scatter.smooth(log(potald$latest_ald*potald$days_since_ald), fit_potald_spl$residuals)
# qqplot(log(potald$latest_ald*potald$days_since_ald), fit_potald_spl$residuals)
# 
# scatter.smooth(log(potald$latest_pot), fit_potald_spl$residuals)
# qqplot(log(potald$latest_pot), fit_potald_spl$residual)
# scatter.smooth(log(potald$latest_pot*potald$days_since_ald), fit_potald_spl$residuals)
# qqplot(log(potald$latest_pot*potald$days_since_ald), fit_potald_spl$residuals)
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
  