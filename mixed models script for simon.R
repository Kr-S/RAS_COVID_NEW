
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



model_df <- potassium_df %>% 
  filter(!is.na(ang2), !is.na(potassium), !is.na(in_icu), !is.na(ald))

indep <- c('ald')

for(lhs in indep){
  
  if(lhs == indep[1]){out <- tibble()}
  
  lmm_null <- lmer(as.formula(paste0('log(',lhs, ') ~ 1 + (1 | pat_id)')),
                   data = model_df) 
  
  lmm_ang2 <- update(lmm_null, .~.+  log(ang2))
  lmm_potassium  <- update(lmm_null, .~.+  potassium)
  lmm_icu <- update(lmm_null, .~.+  in_icu)
  #lmm_ang2_potassium  <- update(lmm_null, .~.+ (log(AngII) + log(Potassium)))
  #lmm_ang2_icu <- update(lmm_null, .~.+ (log(AngII) + log(ICU)))
  #lmm_potassium_icu  <- update(lmm_null, .~.+ (log(Potassium)  + log(ICU)))
  lmm_full      <- update(lmm_null, .~.+ (log(ang2) + potassium + in_icu))
  
  tmp <- bind_rows(
    anova(lmm_null, lmm_ang2)  %>% tidy %>% mutate(test = 'null vs. AngII'),
    anova(lmm_null, lmm_potassium )  %>% tidy %>% mutate(test = 'null vs. Potassium'),
    anova(lmm_null, lmm_icu)  %>% tidy %>% mutate(test = 'null vs. ICU'),
    # anova(lmm_ang2, lmm_ang2_potassium)   %>% tidy %>% mutate(test = 'AngII vs. AngII+Potassium'),
    #  anova(lmm_ang2, lmm_ang2_icu)  %>% tidy %>% mutate(test = 'AngII vs. AngII+ICU'),
    anova(lmm_ang2, lmm_full)  %>% tidy %>% mutate(test = 'AngII vs. AngII+Potassium+ICU')
  ) %>%
    mutate(
      dependent_var = lhs,
      #p.value = print_p(p.value)
    ) %>%
    relocate(dependent_var, test)
  
  r_squared <- bind_rows(
    tibble(term = 'lmm_null', R2m = r.squaredGLMM(lmm_null)[1], R2c = r.squaredGLMM(lmm_null)[2]),
    tibble(term = 'lmm_ang2', R2m = r.squaredGLMM(lmm_ang2)[1], R2c = r.squaredGLMM(lmm_ang2)[2]),
    tibble(term = 'lmm_potassium',  R2m = r.squaredGLMM(lmm_potassium)[1],  R2c = r.squaredGLMM(lmm_potassium)[2]),
    tibble(term = 'lmm_ICU', R2m = r.squaredGLMM(lmm_icu)[1], R2c = r.squaredGLMM(lmm_icu)[2]),
    # tibble(term = 'lmm_ang2_potassium', R2m = r.squaredGLMM(lmm_ang2_potassium)[1], R2c = r.squaredGLMM(lmm_ang2_potassium)[2]),
    # tibble(term = 'lmm_ang2_icu', R2m = r.squaredGLMM(lmm_ang2_icu)[1], R2c = r.squaredGLMM(lmm_ang2_icu)[2]),
    tibble(term = 'lmm_full', R2m = r.squaredGLMM(lmm_full)[1], R2c = r.squaredGLMM(lmm_full)[2])
  )
  
  model_vars <- tibble(
    term = c('lmm_null', 'lmm_ang2', 'lmm_potassium', 'lmm_ICU', 'lmm_full'),
    # term = c('lmm_null', 'lmm_ang2', 'lmm_potassium', 'lmm_icu', 'lmm_ang2_potassium', 'lmm_ang2_icu', 'lmm_full'),
    model = c(
      'null',
      'AngII',
      'Potassium',
      'ICU',
      # 'AngII + Potassium',
      # 'AngII + ICU',
      'AngII + Potassium + ICU'
    )
  )
  
  models <- list(
    'null' = lmm_null,
    'AngII' = lmm_ang2,
    'Potassium'  = lmm_potassium,
    'ICU' = lmm_icu,
    #  'AngII + Potassium' = lmm_ang2_potassium,
    #  'AngII + ICU' = lmm_ang2_icu,
    'AngII + Potassium + ICU' = lmm_full
  )
  
  out <- bind_rows(out, tmp %>% left_join(r_squared) %>% left_join(model_vars) %>% 
                     mutate(model = factor(model, levels = rev(c(
                       'null',
                       'AngII',
                       'Potassium',
                       'ICU',
                       #  'AngII + Potassium',
                       #  'AngII + ICU',
                       'AngII + Potassium + ICU')
                     ))))
  
}


pbDat <- data.frame(y=log(model_df$ald),x1=log(model_df$ang2),x2=model_df$potassium,x3=model_df$in_icu, g1=model_df$pat_id,pearson=residuals(lmm_full,type="pearson"))

ggplot(pbDat,
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(pbDat,
       aes(x=x2,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(pbDat,
       aes(x=x3,y=pearson)) +
  geom_point() +
  theme_bw()


means <- aggregate(pbDat[,c("x1","x2","x3")],by=list(pbDat$g1),FUN=mean)
lmcoefs <- summary(lm(y ~ x1 + x2 + x3 + g1, data=pbDat))$coefficients[,"Estimate"]
means$effects <- c(0,lmcoefs[substr(names(lmcoefs),1,2) == "g1"])
means$effects <- means$effects - mean(means$effects)

cor(means[,c("x1","x2","x3","effects")])

ggplot(means, aes(x=x1,y=effects)) +
  geom_point() +
  theme_bw()

ggplot(means, aes(x=x2,y=effects)) +
  geom_point() +
  theme_bw()

ggplot(means, aes(x=x3,y=effects)) +
  geom_point() +
  theme_bw()



fixef(lmm_full)
lmcoefs[1:4]
#


qqnorm(residuals(lmm_full))
