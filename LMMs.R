


print_p <- function(p){
  case_when(
    p < 0.001 ~ '<0.001*',
    p < 0.01  ~ '<0.01*',
    p < 0.05  ~ '<0.05*',
    TRUE            ~ paste0('=', round(p,3))
  )}

display <- function(x, digits = 2){
  sprintf(paste0('%.', digits, 'f'), x)
}



model_df <- ras_df %>% 
  filter(!is.na(AngII), !is.na(potassium), !is.na(in_icu), !is.na(Aldosterone)) %>% #, RASi == "No RAS-Inhibitor") %>% #, !pat_id %in% c("SARS-FP-75","SARS-FP-94")) %>% 
  mutate(LOQ = ifelse(Aldosterone<=20, 1, 0))
cor(model_df %>% select_if(is.numeric), method="spearman")
describe(model_df$LOQ)

indep <- c('Aldosterone')

for(lhs in indep){
  
  if(lhs == indep[1]){out <- tibble()}
  
  lmm_null <- lmer(as.formula(paste0('log(',lhs, ') ~ 1 + (1 | pat_id)')),
                   data = model_df) 
  
  lmm_AngII <- update(lmm_null, .~.+  AngII)
  lmm_potassium  <- update(lmm_null, .~.+  potassium)
  lmm_icu <- update(lmm_null, .~.+  in_icu)
  lmm_AngII_potassium  <- update(lmm_null, .~. + AngII + potassium)
  lmm_AngII_icu <- update(lmm_null, .~.+ AngII + in_icu)
  lmm_potassium_icu  <- update(lmm_null, .~. + potassium  + in_icu)
  lmm_full      <- update(lmm_null, .~. + AngII + potassium + in_icu)
  
  tmp <- bind_rows(
    anova(lmm_null, lmm_AngII)  %>% tidy %>% mutate(test = 'null vs. AngII'),
    anova(lmm_null, lmm_potassium )  %>% tidy %>% mutate(test = 'null vs. Potassium'),
    anova(lmm_null, lmm_icu)  %>% tidy %>% mutate(test = 'null vs. ICU'),
    anova(lmm_AngII, lmm_AngII_potassium)   %>% tidy %>% mutate(test = 'AngII vs. AngII+Potassium'),
    anova(lmm_AngII, lmm_AngII_icu)  %>% tidy %>% mutate(test = 'AngII vs. AngII+ICU'),
    anova(lmm_AngII_potassium, lmm_full)  %>% tidy %>% mutate(test = 'AngII+Potassium vs. AngII+Potassium+ICU')
  ) %>%
    mutate(
      dependent_var = lhs,
      p.value = print_p(p.value)
    ) %>%
    relocate(dependent_var, test)
  
  r_squared <- bind_rows(
    tibble(term = 'lmm_null', R2m = r.squaredGLMM(lmm_null)[1], R2c = r.squaredGLMM(lmm_null)[2]),
    tibble(term = 'lmm_AngII', R2m = r.squaredGLMM(lmm_AngII)[1], R2c = r.squaredGLMM(lmm_AngII)[2]),
    tibble(term = 'lmm_potassium',  R2m = r.squaredGLMM(lmm_potassium)[1],  R2c = r.squaredGLMM(lmm_potassium)[2]),
    tibble(term = 'lmm_ICU', R2m = r.squaredGLMM(lmm_icu)[1], R2c = r.squaredGLMM(lmm_icu)[2]),
    tibble(term = 'lmm_AngII_potassium', R2m = r.squaredGLMM(lmm_AngII_potassium)[1], R2c = r.squaredGLMM(lmm_AngII_potassium)[2]),
    tibble(term = 'lmm_AngII_icu', R2m = r.squaredGLMM(lmm_AngII_icu)[1], R2c = r.squaredGLMM(lmm_AngII_icu)[2]),
    tibble(term = 'lmm_full', R2m = r.squaredGLMM(lmm_full)[1], R2c = r.squaredGLMM(lmm_full)[2])
  )
  
  model_vars <- tibble(
    #term = c('lmm_null', 'lmm_AngII', 'lmm_potassium', 'lmm_ICU', 'lmm_full'),
    term = c('lmm_null', 'lmm_AngII', 'lmm_potassium', 'lmm_icu', 'lmm_AngII_potassium', 'lmm_AngII_icu', 'lmm_full'),
    model = c(
      'null',
      'AngII',
      'Potassium',
      'ICU',
      'AngII + Potassium',
      'AngII + ICU',
      'AngII + Potassium + ICU'
    )
  )
  
  models <- list(
    'null' = lmm_null,
    'AngII' = lmm_AngII,
    'Potassium'  = lmm_potassium,
    'ICU' = lmm_icu,
    'AngII + Potassium' = lmm_AngII_potassium,
    'AngII + ICU' = lmm_AngII_icu,
    'AngII + Potassium + ICU' = lmm_full
  )
  
  out <- bind_rows(out, tmp %>% left_join(r_squared) %>% left_join(model_vars) %>% 
                     mutate(model = factor(model, levels = rev(c(
                       'null',
                       'AngII',
                       'Potassium',
                       'ICU',
                       'AngII + Potassium',
                       'AngII + ICU',
                       'AngII + Potassium + ICU')
                     ))))
  
}

datatable(out)


pbDat <- data.frame(Aldosterone=log(model_df$Aldosterone),
                    AngII=log(model_df$AngII),
                    potassium=model_df$potassium,
                    in_icu=model_df$in_icu, 
                    pat_id=model_df$pat_id,
                    pearson=residuals(lmm_full,type="pearson"),
                    predicted=predict(lmm_full),
                    leverin_icu=hatvalues(lmm_full),
                    LOQ = factor(ifelse(model_df$Aldosterone<=20, ifelse(model_df$AngII<=2, 3, 1), ifelse(model_df$AngII<=2, 2, 0)), levels = c(0,1,2,3), labels = c("Regular", "LOQ Aldosterone", "LOQ AngII", "LOQ combined")))

lmm_full_PQL <- glmmPQL(Aldosterone ~ AngII + potassium, ~1|pat_id, family = gaussian(link = "log"), data = model_df, verbose = FALSE)


model_df %>%  select(LOQ,pat_id) %>% group_by(pat_id) %>% table

plot(lmm_full_PQL)
plot(lmm_full)
summary(lmm_full)

qqnorm(residuals(lmm_full))

ggplot(pbDat, aes(x=predicted,y=leverin_icu)) +
  geom_point(size=2)+ 
  geom_smooth(method="loess") + 
  geom_text(aes(label=pat_id))+
  theme_bw()

ggplot(pbDat, aes(x=predicted,y=log(Aldosterone))) +
  geom_point(size=2)+ 
  geom_smooth(method="loess") +
  theme_bw()

ggplot(pbDat, aes(x=predicted,y=pearson)) +
  geom_point(size=2)+ 
  geom_smooth(method="loess") +
  theme_bw()
  
ggplot(pbDat,
       aes(x=AngII,y=pearson)) +
  geom_point()+ geom_smooth(method="loess")+
  theme_bw()

ggplot(pbDat,
       aes(x=potassium,y=pearson)) +
  geom_point() + geom_smooth(method="loess") +
  theme_bw()

ggplot(pbDat,
       aes(x=in_icu,y=pearson)) +
  geom_point() + geom_smooth(method="loess") +
  theme_bw()

ggplot(pbDat,
       aes(x=log(Aldosterone),y=pearson)) +
  geom_point()+ geom_smooth(method="loess") +
  theme_bw()


means <- aggregate(pbDat[,c("AngII","potassium","in_icu")],by=list(pbDat$pat_id),FUN=mean)
lmcoefs <- summary(lm(log(Aldosterone) ~ AngII + potassium + in_icu + pat_id, data=pbDat))$coefficients[,"Estimate"]
means$effects <- c(0,lmcoefs[substr(names(lmcoefs),1,6) == "pat_id"])
means$effects <- means$effects - mean(means$effects)

cor(means[,c("AngII","potassium","in_icu","effects")],method="spearman")

ggplot(means, aes(x=AngII,y=effects)) +
  geom_point() +
  theme_bw()

ggplot(means, aes(x=potassium,y=effects)) +
  geom_point() +
  theme_bw()

ggplot(means, aes(x=in_icu,y=effects)) +
  geom_point() +
  theme_bw()



fixef(lmm_full)
lmcoefs[1:4]
#


qqnorm(residuals(lmm_full))

