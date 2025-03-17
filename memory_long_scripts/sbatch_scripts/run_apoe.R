args = commandArgs(TRUE)

feature=as.character(args[1])
wd=as.character(args[2])


run_apoe= function(feature, wd) {
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
  
  df.megagam = wrapper_prepare_data(df.merged, feature, df.features)
  df.megagam = df.megagam%>% filter(!is.na(apoe_status))
  df.megagam$apoe_statusF <- as.factor(df.megagam$apoe_status)
  df.megagam$apoe_statusO <- factor(df.megagam$apoe_status, ordered = T, levels = c(0,1))
  
  
  ## LME part
  print("lme models *********************************")
  mod.lme.apoe.memory = lmer(memory_slope ~ apoe_status + (1|dataset), data = df.megagam, weights = weight.sub2)
  mod.lme.apoe.brain = lmer(delta_brain ~ apoe_status + (1|dataset), data = df.megagam, weights = weight.sub2)
  
  mod.lme.apoe = lmer(memory_slope ~ delta_brain * apoe_status + (1|dataset), data = df.megagam, weights = weight.sub2)
  # Model for carriers
  mod.lme.apoe1 = lmer(memory_slope ~ delta_brain + (1|dataset), data = df.megagam[df.megagam$apoe_status == 1,], weights = weight.sub2)
  # Model for non-carriers
  mod.lme.apoe0 = lmer(memory_slope ~ delta_brain + (1|dataset), data = df.megagam[df.megagam$apoe_status == 0,], weights = weight.sub2)
  
  
  ## GAMM Part
  # Base model - apoe specific
  print("gam simple effects models *********************************")
  df.apoe0 = df.megagam %>% filter(apoe_status == 0)
  df.apoe1 = df.megagam %>% filter(apoe_status == 1)
  mod.gam_main.w.apoe0 = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr"), 
         data = df.apoe0, random=list(dataset=~1), weights =weight.sub2))
  mod.gam_main.w.apoe1 = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr"), 
         data = df.apoe1, random=list(dataset=~1), weights =weight.sub2))
  mod.gam_int.w.apoe0 = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage), 
         data = df.apoe0, random=list(dataset=~1), weights =weight.sub2))
  mod.gam_int.w.apoe1 = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage), 
         data = df.apoe1, random=list(dataset=~1), weights =weight.sub2))
  
  print("nulls for basic models *********************************")
    #nulls for basic models 
    null.mod.gam_main.w.apoe0 = suppressWarnings(
      gamm(memory_slope ~ 1, data = df.apoe0, random=list(dataset=~1), weights =weight.sub2))
    null.mod.gam_main.w.apoe1 = suppressWarnings(
      gamm(memory_slope ~ 1, data = df.apoe1, random=list(dataset=~1), weights =weight.sub2))
    null.mod.gam_int.w.apoe0 = suppressWarnings(
      gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr"), 
           data = df.apoe0, random=list(dataset=~1), weights =weight.sub2))
    null.mod.gam_int.w.apoe1 = suppressWarnings(
      gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr"), 
           data = df.apoe1, random=list(dataset=~1), weights =weight.sub2))
    
  #### Main model by APOE #####
  print("main model for apoe *********************************")
  mod.gam_main.w.apoeO = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(delta_brain, by = apoe_statusO) + apoe_statusO, 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2))
  
  # null model
  print("main null model for apoe *********************************")
  null.mod.gam_main.w.apoeO =suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + apoe_statusO, 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2))
  
  
  ##### Age x brain x APOE interaction #####
  print("apoe interaction *********************************")
  mod.gam_int.w.apoeO = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage) + ti(delta_brain, xage, by = apoe_statusO) + apoe_statusO, 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2))
  
  print("apoe interaction null model *********************************")
  null.mod.gam_int.w.apoeO = suppressWarnings(
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage) + apoe_statusO, 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2))
  
  
  #### bootstrap ####
  set.seed(1234)
  nreps = 1000
  
  ### Bootstrap BRAIN x APOE ###
  print("bootstrap brainx apoe *********************************")
  g.true =  mod.gam_main.w.apoeO 
  g.null = null.mod.gam_main.w.apoeO
  y.null = predict(null.mod.gam_main.w.apoeO$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.megagam, model = "main")}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[2,4]]
  (p.boot.main = sum(p.orig > grot) / length(grot))
  
  
  ### Bootstrap BRAIN x APOE X AGE ###
  print("bootstrap brain x apoe x age *********************************")
  g.true =  mod.gam_int.w.apoeO
  g.null =  null.mod.gam_int.w.apoeO   
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.megagam)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[4,4]]
  p.boot.int = sum(p.orig > grot) / length(grot) 
  
  
  ### Bootstrap APOE filtered groups ### 
      ## p-values are not very good for comparisons, so use with caution
  # apoe0 - main effect
  print("bootstrap apoe simple effects models *********************************")
  g.true =  mod.gam_main.w.apoe0
  g.null =  null.mod.gam_main.w.apoe0
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.apoe0,model = "main", apoe.filtered = T)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[1,4]]
  p.boot.main.apoe0 = sum(p.orig > grot) / length(grot) 
  
  # apoe1 - main effect
  g.true =  mod.gam_main.w.apoe1
  g.null =  null.mod.gam_main.w.apoe1
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.apoe1,model = "main", apoe.filtered = T)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[1,4]]
  p.boot.main.apoe1 = sum(p.orig > grot) / length(grot)
  
  # apoe0 - int effect
  g.true =  mod.gam_int.w.apoe0
  g.null =  null.mod.gam_int.w.apoe0
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.apoe0, apoe.filtered = T)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[3,4]]
  p.boot.int.apoe0 = sum(p.orig > grot) / length(grot) 
  
  g.true =  mod.gam_int.w.apoe1
  g.null =  null.mod.gam_int.w.apoe1
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.apoe1, apoe.filtered = T)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[3,4]]
  p.boot.int.apoe1 = sum(p.orig > grot) / length(grot) 
  
  
  
  ##### PREDICTIONS #####
  print("predictions *********************************")
  # estimate density
  d = density(mod.gam_main.w.apoeO$gam$model$delta_brain)
  
  # main estimates apoe 0 and 1 (neg and pos). estimates and derivatives with confidence intervals
  sm0 <- smooth_estimates(mod.gam_main.w.apoe0) |>
    add_confint() %>% 
    mutate(apoe_status = 0)
  sm1 <- smooth_estimates(mod.gam_main.w.apoe1) |>
    add_confint()  %>% 
    mutate(apoe_status = 1)
  predict.apoe.bygroup = rbind(sm0, sm1) %>% 
    mutate(apoe_statusF = factor(apoe_status), 
           features = feature)
  dm0 <- derivatives(mod.gam_main.w.apoe0) %>% 
    mutate(apoe_status = 0)
  dm1 <- derivatives(mod.gam_main.w.apoe1) %>% 
    mutate(apoe_status = 1)
  derivative.apoe.bygroup = rbind(dm0, dm1) %>% 
    mutate(apoe_statusF = factor(apoe_status))
  derivative.apoe.bygroup$w = approx(d$x, d$y, xout = derivative.apoe.bygroup$delta_brain)$y
  derivative.apoe.bygroup$features = feature 
  
  # main estimates by apoe
  predict.apoeO = smooth_estimates(mod.gam_main.w.apoeO)
  predict.apoeO$features = feature
  derivatives.apoeO =derivatives(mod.gam_main.w.apoeO)
  derivatives.apoeO$w = approx(d$x, d$y, xout = derivatives.apoeO$delta_brain)$y
  derivatives.apoeO$features = feature 
  
  
  #df.predict.tensory.interaction by apoe_group
  grid <- expand.grid(xage = seq(25, 85, by = 5), 
                      delta_brain = seq(-4.5,4.5, length.out = 100))
  pred.apoe0 <- predict(mod.gam_int.w.apoe0$gam, newdata = grid, se.fit = TRUE)
  pred.apoe1 <- predict(mod.gam_int.w.apoe1$gam, newdata = grid, se.fit = TRUE)
  
  predict.interaction.bygroup = 
    rbind(data.frame(
      xage = grid$xage, 
      delta_brain = grid$delta_brain, 
      estimate = pred.apoe0$fit, 
      se = pred.apoe0$se.fit, 
      apoe_status = 0), 
      data.frame(
        xage = grid$xage, 
        delta_brain = grid$delta_brain, 
        estimate = pred.apoe1$fit, 
        se = pred.apoe1$se.fit, 
        apoe_status = 1))
  
  
  grid2 <- expand.grid(xage = seq(25, 85, by = 5), 
                       delta_brain = seq(-4.5,4.5, length.out = 100), 
                       apoe_statusO = c(0,1))
  
  pred.apoeO <- predict(mod.gam_int.w.apoeO$gam, newdata = grid2, se.fit = TRUE)
  predict.interaction.apoeO = 
    data.frame(
      xage = grid2$xage, 
      delta_brain = grid2$delta_brain, 
      apoe_statusO = grid2$apoe_statusO,
      estimate = pred.apoeO$fit, 
      se = pred.apoeO$se.fit)
  
  
  # compute derivatives for interaction effects
  predict.interaction.apoeO <- predict.interaction.apoeO %>%
    group_by(xage, apoe_statusO) %>%
    nest() %>%
    mutate(data = map(data, ~ .x %>%
                        mutate(derivative = c(NA, diff(.$estimate) / diff(.$delta_brain))))) %>% 
    unnest() 
  
  predict.interaction.apoeO$w = approx(d$x, d$y, xout = predict.interaction.apoeO$delta_brain)$y
  predict.interaction.apoeO$features = feature
  
  
  
  predict.interaction.bygroup = 
    predict.interaction.bygroup %>% 
    group_by(xage, apoe_status) %>%
    nest() %>%
    mutate(data = map(data, ~ .x %>%
                        mutate(derivative = c(NA, diff(.$estimate) / diff(.$delta_brain))))) %>% 
    unnest() 
  predict.interaction.bygroup$w = approx(d$x, d$y, xout = predict.interaction.bygroup$delta_brain)$y
  predict.interaction.bygroup$features = feature 
  
  
  
  
  
  
  output <- data.frame(features = feature)
  output <- output %>%
    mutate(lme.apoe.memory_beta = summary(mod.lme.apoe.memory)$coefficients[2, 1],
           lme.apoe.memory_SE = summary(mod.lme.apoe.memory)$coefficients[2, 2],
           lme.apoe.memory_pval = summary(mod.lme.apoe.memory)$coefficients[2, 5],
           lme.apoe.brain_beta = summary(mod.lme.apoe.brain)$coefficients[2, 1],
           lme.apoe.brain_SE = summary(mod.lme.apoe.brain)$coefficients[2, 2],
           lme.apoe.brain_pval = summary(mod.lme.apoe.brain)$coefficients[2, 5],
           lme.apoe_beta = summary(mod.lme.apoe)$coefficients[4, 1],
           lme.apoe_SE = summary(mod.lme.apoe)$coefficients[4, 2],
           lme.apoe_pval = summary(mod.lme.apoe)$coefficients[4, 5],
           gam.main.w.apoe0F = summary(mod.gam_main.w.apoe0$gam)$s.table[1,3], 
           gam.main.w.apoe0p = summary(mod.gam_main.w.apoe0$gam)$s.table[1,4],
           gam.main.w.apoe1F = summary(mod.gam_main.w.apoe1$gam)$s.table[1,3], 
           gam.main.w.apoe1p = summary(mod.gam_main.w.apoe1$gam)$s.table[1,4],
           gam.int.w.apoe0F = summary(mod.gam_int.w.apoe0$gam)$s.table[3,3], 
           gam.int.w.apoe0p = summary(mod.gam_int.w.apoe0$gam)$s.table[3,4], 
           gam.int.w.apoe1F = summary(mod.gam_int.w.apoe1$gam)$s.table[3,3], 
           gam.int.w.apoe1p = summary(mod.gam_int.w.apoe1$gam)$s.table[3,4], 
           gam.main.w.apoeOF = summary(mod.gam_main.w.apoeO$gam)$s.table[2,3],
           gam.main.w.apoeOp = summary(mod.gam_main.w.apoeO$gam)$s.table[2,4],
           gam.int.w.apoeOF = summary(mod.gam_int.w.apoeO$gam)$s.table[4,3],
           gam.int.w.apoeOp = summary(mod.gam_int.w.apoeO$gam)$s.table[4,4],
           gam.int.w.apoeOp.bootstrap = p.boot.int,
           gam.main.w.apoeOp.bootstrap = p.boot.main,
           gam.main.w.apoe0p.bootstrap = p.boot.main.apoe0,
           gam.main.w.apoe1p.bootstrap = p.boot.main.apoe1, 
           gam.int.w.apoe0p.bootstrap = p.boot.int.apoe0,
           gam.int.w.apoe1p.bootstrap = p.boot.int.apoe1)
    
  
  
  models <- list(
    mod.lme.apoe.brain = mod.lme.apoe.brain,
    mod.lme.apoe.memory = mod.lme.apoe.memory,
    mod.lme.apoe = mod.lme.apoe,
    mod.lme.apoe0 = mod.lme.apoe0, 
    mod.lme.apoe1 = mod.lme.apoe1, 
    mod.gam_int.w.apoe0 = mod.gam_int.w.apoe0, 
    mod.gam_int.w.apoe1 = mod.gam_int.w.apoe1, 
    mod.gam_int.w.apoeO = mod.gam_int.w.apoeO, 
    mod.gam_main.w.apoe0 = mod.gam_main.w.apoe0, 
    mod.gam_main.w.apoe1 = mod.gam_main.w.apoe1, 
    mod.gam_main.w.apoeO = mod.gam_main.w.apoeO)
  
  predictions = list(
    predict.interaction.bygroup,
    predict.interaction.apoeO,
    predict.apoe.bygroup, 
    derivative.apoe.bygroup,
    predict.apoeO, 
    derivatives.apoeO
    )
  
  
  ###### SAVE #######
  save(output, 
       file = file.path(wd,"mega/apoe/output/temp_models",
                        paste0(feature, ".Rda")))
  
  save(models, 
       predictions,
       file = file.path(wd,"mega/apoe/output/temp_models",
                        paste0("models_",feature, ".Rda")))
  
}

wrapper_prepare_data = function(df, feature, df.features) {
  df =
    df %>% 
    filter(features == feature) %>% 
    unnest(data) %>% 
    filter(abs(memory_slope) < 4.5, abs(delta_brain) < 4.5)  #z score, same her
  df = inner_join(df, df.features)
  
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



gamm_bootstrap = function(y.null,g.null, df, model = "tensor_interaction", apoe.filtered = F) { 
  n = length(y.null)
  wild = sample(c(-1 ,1), size=n, replace=TRUE)
  y.resample = y.null + residuals(g.null$gam)*wild           
  df$y_resample = y.resample  %>% as.numeric() 
  
  if (apoe.filtered == F) {
    if (model == "tensor_interaction" ) {
      boostrap1 = 
        gamm(y_resample ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage) + ti(delta_brain, xage, by = apoe_statusO) + apoe_statusO, 
             data = df, random=list(dataset=~1), weights =weight.sub2)
      bs_p = summary(boostrap1$gam)$s.table[[4,4]]
    } else if (model == "main") {
      boostrap1 = 
        gamm(y_resample ~ s(delta_brain, bs = "cr") + s(delta_brain, by = apoe_statusO) + apoe_statusO, data = df, random=list(dataset=~1), weights =weight.sub2)
      bs_p = summary(boostrap1$gam)$s.table[[2,4]]
      
    }
  } else {
    if (model == "tensor_interaction" ) {
      boostrap1 = 
        gamm(y_resample ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage), 
             data = df, random=list(dataset=~1), weights =weight.sub2)
      bs_p = summary(boostrap1$gam)$s.table[[3,4]]
    } else if (model == "main") {
      boostrap1 = 
        gamm(y_resample ~ s(delta_brain, bs = "cr"), data = df, random=list(dataset=~1), weights =weight.sub2)
      bs_p = summary(boostrap1$gam)$s.table[[1,4]]
      
    }
  }
  return(bs_p)
}  
  
  
run_apoe(feature, wd)