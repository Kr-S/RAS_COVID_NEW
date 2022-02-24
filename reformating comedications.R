
library(tidyverse)
library(lubridate)
library(stringr)
library(readxl)
library(ggthemes)


meds = read_xlsx("data/meds potassium.xlsx")

x=0
for (name in colnames(meds)){
  x = x+1
  print(x)
  print(colnames(meds)[x])
  if (any(list(grepl('^Daily', colnames(meds)[x]), grepl('^Total', colnames(meds)[x]), grepl('^Form', colnames(meds)[x]),grepl('^Medication', colnames(meds)[x])))){
    colnames(meds)[x] = paste0(unlist(str_split(colnames(meds)[x-1],"\\_",n=2))[1], '_', colnames(meds)[x])
  }
}
colnames(meds)



dailynames <- colnames(meds[,grepl('Daily', colnames(meds))])


newmeds = gsub("[^\\.0-9]", "", meds[,grepl('Daily', colnames(meds))])

newmeds = meds %>%
  mutate_at(dailynames, ~ str_replace_all(., "[^,\\.\\d]", ""))


dailynames


newmeds[,grepl('Daily', colnames(newmeds))]

lapply(newmeds[,grepl('Daily', colnames(newmeds))], function(x) {x[!is.na(x)]})


lapply(meds[,grepl('Daily', colnames(meds))], function(x) {x[!is.na(x)]})


unique(unlist(lapply(meds[,grepl('Class', colnames(meds), ignore.case=TRUE)], function(x) {x[!is.na(x)]})))

