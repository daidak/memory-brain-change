args = commandArgs(TRUE)

feature=as.character(args[2])
wd=as.character(args[1])


run_apoe_brain= function(wd, feature) {
  library(gratia)
  library(tidyverse)
  library(here)
  library(lmerTest)
  library(mgcv)
  
  #wd = here("data_memory_long/all")
  load(file.path(wd, "merged_reliability.rda"))
  options(bitmapType = "cairo")
  rel_parameters = "/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/reliability/reliability_parameters.rda"
  load(rel_parameters)
  
  df = wrapper_prepare_data(df.merged, feature, df.features)
  
  mod.lme.apoe.brain = lmer(delta_brain ~ apoe_status + (1|dataset), data = df)
  mod.lme.w.apoe.brain = lmer(delta_brain ~ apoe_status + (1|dataset), data = df, weights = weight.sub2)
  mod.gam.apoe.brain = 
    gamm(delta_brain ~ s(age, bs = "cr") + s(age, by = apoe_statusO) + apoe_statusO, 
         data = df, random=list(dataset=~1), weights =weight.sub2)
  mod.gam.w.apoe.brain = 
    gamm(delta_brain ~ s(age, bs = "cr") + s(age, by = apoe_statusO) + apoe_statusO, 
         data = df, random=list(dataset=~1), weights =weight.sub2)
  
  
  null.mod.gam.w.apoe.brain = 
    gamm(delta_brain ~ s(age, bs = "cr") + apoe_statusO, 
         data = df, random=list(dataset=~1), weights =weight.sub2)
  
  #### bootstrap ####
  set.seed(1234)
  nreps = 5000
  
  ### Bootstrap BRAIN x APOE ###
  print("bootstrap brainx apoe *********************************")
  g.true =  mod.gam.w.apoe.brain
  g.null = null.mod.gam.w.apoe.brain
  y.null = predict(null.mod.gam.w.apoe.brain$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap_brain_apoe(y.null,g.null,df)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[2,4]]
  (p.boot.brain.apoe.age = sum(p.orig > grot) / length(grot))
  
  
  ##### PREDICTIONS #####
  print("predictions *********************************")
  # estimate density
  d = density(mod.gam.w.apoe.brain$gam$model$age)
  
  # main estimates apoe 0 and 1 (neg and pos). estimates and derivatives with confidence intervals
  sm0 <- smooth_estimates(mod.gam.w.apoe.brain) |>
    add_confint() 
  sm0.derivatives <- derivatives(mod.gam.w.apoe.brain$gam) 
  
  
  
  #df.predict.tensory.interaction by apoe_group
  grid <- expand.grid(age = seq(25, 85, by = 5), 
                      apoe_statusO = c(0,1))
  pred.apoe <- predict(mod.gam.w.apoe.brain$gam, newdata = grid, se.fit = TRUE)
  
  sm1 = 
    data.frame(
      age = grid$age, 
      apoe_statusO = grid$apoe_statusO, 
      estimate = pred.apoe$fit, 
      se = pred.apoe$se.fit)
  sm1 = 
    sm1 %>% 
    group_by(apoe_statusO) %>% 
    nest() %>% 
    mutate(data = map(data, ~ .x %>%
                        mutate(derivative = c(NA, diff(.$estimate) / diff(.$age))))) %>% 
    unnest() 
  
  sm1$w = approx(d$x, d$y, xout = sm1$age)$y
  
  
  output <- data.frame(features = "brain")
  output <- output %>%
    mutate(lme.apoe.brain_beta = summary(mod.lme.w.apoe.brain)$coefficients[2, 1],
           lme.apoe.brain_SE = summary(mod.lme.w.apoe.brain)$coefficients[2, 2],
           lme.apoe.brain_pval = summary(mod.lme.w.apoe.brain)$coefficients[2, 5],
           gam.apoe.brain_ageF = summary(mod.gam.w.apoe.brain$gam)$s.table[2, 3],
           gam.apoe.brain_agep = summary(mod.gam.w.apoe.brain$gam)$s.table[2, 4],
           gam.apoe.brain_agep.boot = p.boot.brain.apoe.age)
  
  
  
  models = list(mod.lme.apoe.brain = mod.lme.apoe.brain,
                mod.lme.w.apoe.brain = mod.lme.w.apoe.brain,
                mod.gam.w.apoe.brain = mod.gam.w.apoe.brain,
                mod.gam.apoe.brain = mod.gam.apoe.brain)
  
  
  predictions = list(
    sm0,
    sm1,
    sm0.derivatives
    )
  
  
  ###### SAVE #######
  save(output, 
       file = file.path(wd,"mega/apoe/output/temp_models_mem",
                        paste0(feature, ".Rda")))
  
  save(models, 
       predictions,
       file = file.path(wd,"mega/apoe/output/temp_models_mem",
                        paste0("models.", feature,".Rda")))
  
}

wrapper_prepare_data = function(df, feature, df.features) {
  
  df = df%>% 
    filter(features == feature) %>% 
    unnest(data) %>% 
    filter(abs(delta_brain) < 4.5)
  df = inner_join(df, df.features)
  df= df%>% filter(!is.na(apoe_status))
  df$apoe_statusF <- as.factor(df$apoe_status)
  df$apoe_statusO <- factor(df$apoe_status, ordered = T, levels = c(0,1))
  
  #weighting variable
  df = 
    df %>% 
    mutate(bs = seD^2,
           ws = pct.err.mean^2*((time^2*n*(n+1))/ (12*(n-1)))^-1,
           icc.sub = bs /(ws + bs),
           weight.sub = 1/(1-icc.sub), 
           weight.sub2 = icc.sub^2) %>% 
    mutate(weight.sub = if_else(weight.sub > 10*min(weight.sub), 10*min(weight.sub), weight.sub), 
           weight.sub2 = if_else(weight.sub2 < .09, .09, weight.sub2))
  return(df)
}


gamm_bootstrap_brain_apoe = function(y.null,g.null, df) { 
  n = length(y.null)
  wild = sample(c(-1 ,1), size=n, replace=TRUE)
  y.resample = y.null + residuals(g.null$gam)*wild           
  df$y_resample = y.resample  %>% as.numeric() 
  
  boostrap1 = 
    gamm(y_resample ~ s(age, bs = "cr") + s(age, by = apoe_statusO) + apoe_statusO, 
         data = df, random=list(dataset=~1), weights =weight.sub2)
  bs_p = summary(boostrap1$gam)$s.table[[2,4]]
  return(bs_p)
}   
  
  
run_apoe_brain(wd, feature)


