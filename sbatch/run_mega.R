args = commandArgs(TRUE)

feature=as.character(args[1])
wd=as.character(args[2])
#wd= here::here("data_memory_long/all")
#feature = "Left-Hippocampus"

run_mega= function(feature, wd) {
  library(tidyverse)
  library(dplyr) 
  library(lmerTest)
  library(mgcv)
  nreps = 5000
  
  load(file.path(wd, "merged_reliability.rda"))
  options(bitmapType = "cairo")
  rel_parameters = "/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/reliability/reliability_parameters.rda"
  load(rel_parameters)
  
  # prepare data, merge reliability, get weights, remove outliers, etc. 
  df.megagam = wrapper_prepare_data(df.merged, feature, df.features)
  
  
  #LME 
  mod.lme = lmer(memory_slope ~ delta_brain  + (1|dataset), data = df.megagam)
  #mod.lme.slopes = lmer(memory_slope ~ delta_brain  + (1 + xage |dataset), data = df.megagam)
  
  #LME with weights
  mod.lme.w = lmer(memory_slope ~ delta_brain  + (1|dataset), data = df.megagam, weights = weight.sub2)
  mod.lme.slopes.w = lmer(memory_slope ~ delta_brain  + (1 + xage |dataset), data = df.megagam, weights = weight.sub2)
  
  #LME with weights control analysese 
  mod.lme.w.ctr.icv = lmer(memory_slope ~ delta_brain + icv + (1|dataset), data = df.megagam, weights = weight.sub2)
  mod.lme.w.ctr.icv.xage = lmer(memory_slope ~ delta_brain + icv + xage + + (1|dataset), data = df.megagam, weights = weight.sub2)
  
  # GAM
  #mod.megagam =
  #  gamm(memory_slope ~ s(delta_brain, bs = "cr"), data = df.megagam, random=list(dataset=~1))
  
  # GAM with weights
  mod.megagam.w =
    gamm(memory_slope ~ s(delta_brain, bs = "cr"), data = df.megagam, random=list(dataset=~1), weights = weight.sub2)
  # null model
  null.mod.megagam.w =
    gamm(memory_slope ~ 1, data = df.megagam, random=list(dataset=~1), weights = weight.sub2)
  
  
  
  #gam age interaction
  #mod.gam_int = 
  #  gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage), data = df.megagam, random=list(dataset=~1))
  
  #gam with weights age interaction
  mod.gam_int.w = 
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr") + ti(delta_brain, xage), 
       data = df.megagam, random=list(dataset=~1), weights =weight.sub2)
  
  ## RUN bootstrap ##
  #null gam withouth age interaction
  null.mod.gam_int.w = 
    gamm(memory_slope ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr"), 
         data = df.megagam, random=list(dataset=~1), weights =weight.sub2)
  
  
  
  ## predict matrices
  grid1 <- expand.grid(delta_brain = seq(-3,3, length.out = 20))
  
  grid2 <- expand.grid(xage = seq(25, 85, by = 5), 
                      delta_brain = seq(-3,3, length.out = 21))
  
  pred.main <- predict(mod.megagam.w$gam, newdata = grid1, se.fit = TRUE)
  pred.int <- predict(mod.gam_int.w$gam, newdata = grid2, se.fit = TRUE)
  
  
  df.predict.main <- data.frame(
    delta_brain = grid1$delta_brain,
    se = pred.main$se.fit,
    estimate = pred.main$fit)
  
  df.predict.tensor.interaction <- data.frame(
    xage = grid2$xage, 
    delta_brain = grid2$delta_brain, 
    estimate = pred.int$fit, 
    se = pred.int$se.fit)
  
  
  #gam age interaction bootstrap
 
  g.true = mod.gam_int.w 
  g.null =  null.mod.gam_int.w   
  y.null = predict(g.null$lme)
 
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.megagam)}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  
  p.orig = summary(g.true$gam)$s.table[[3,4]]
  p.boot = sum(p.orig > grot) / length(grot)
 
  
  #gam main bootstrap
  g.true = mod.megagam.w 
  g.null = null.mod.megagam.w
  y.null = predict(g.null$lme)
  grot = replicate(nreps,  tryCatch({gamm_bootstrap(y.null,g.null,df.megagam, model = "main")}, error = function(e){ NaN}))
  grot = grot[!is.nan(grot)]
  p.orig = summary(g.true$gam)$s.table[[1,4]]
  (p.boot.main = sum(p.orig > grot) / length(grot))
  
  
  
  xx <- data.frame(features = feature)
  
  xx <- xx %>%
    mutate(lme_beta = summary(mod.lme)$coefficients[2, 1],
           lme_SE = summary(mod.lme)$coefficients[2, 2],
           lme_pval = summary(mod.lme)$coefficients[2, 5],
           lme.w_beta = summary(mod.lme.w)$coefficients[2, 1],
           lme.w_SE = summary(mod.lme.w)$coefficients[2, 2],
           lme.w_pval = summary(mod.lme.w)$coefficients[2, 5],
           lme.w_beta.ctrlicv = summary(mod.lme.w.ctr.icv)$coefficients[2, 1],
           lme.w_SE.ctrlicv = summary(mod.lme.w.ctr.icv)$coefficients[2, 2],
           lme.w_pval.ctrlicv = summary(mod.lme.w.ctr.icv)$coefficients[2, 5],
           lme.w_beta.ctrlicvage = summary(mod.lme.w.ctr.icv.xage)$coefficients[2, 1],
           lme.w_SE.ctrlicvage = summary(mod.lme.w.ctr.icv.xage)$coefficients[2, 2],
           lme.w_pval.ctrlicvage = summary(mod.lme.w.ctr.icv.xage)$coefficients[2, 5],
           lme.w_beta = summary(mod.lme.w)$coefficients[2, 1],
           lme.w_SE = summary(mod.lme.w)$coefficients[2, 2],
           lme.w_pval = summary(mod.lme.w)$coefficients[2, 5],
           lme.slopes.w_beta = summary(mod.lme.slopes.w)$coefficients[2, 1],
           lme.slopes.w_SE = summary(mod.lme.slopes.w)$coefficients[2, 2],
           lme.slopes.w_pval = summary(mod.lme.slopes.w)$coefficients[2, 5],
           gamm.w.pval = summary(mod.megagam.w$gam)$s.pv[1],
           gamm.w.pval.bootstrap = p.boot.main,
           gamm.int.w.pval = summary(mod.gam_int.w$gam)$s.pv[3],
           gamm.int.w.pval.bootstrap = p.boot)
         
  
  models <- list(
    lme = mod.lme,
    lme.w = mod.lme.w,
    lme.w.ctr.icv = mod.lme.w.ctr.icv,
    lme.w.ctr.icv.age = mod.lme.w.ctr.icv.xage,
    gam.w = mod.megagam.w,
    gam_int.w = mod.gam_int.w)
  
  save(xx, 
     file = file.path(wd,"mega",
                      "temp_models",
                      paste0(feature, ".Rda")))
  
  save(models,
       df.predict.main, 
       df.predict.tensor.interaction, 
     file = file.path(wd,"mega",
                      "temp_models",
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


gamm_bootstrap = function(y.null,g.null, df, model = "tensor_interaction") { 
  n = length(y.null)
  wild = sample(c(-1 ,1), size=n, replace=TRUE)
  y.resample = y.null + residuals(g.null$gam)*wild           
  df$y_resample = y.resample  %>% as.numeric() 
  
  if (model == "tensor_interaction" ) {
  boostrap1 = 
    gamm(y_resample ~ s(delta_brain, bs = "cr") + s(xage, bs = "cr")  + ti(delta_brain, xage), data = df, random=list(dataset=~1), weights =weight.sub2)
  bs_p = summary(boostrap1$gam)$s.table[[3,4]]
  } else if (model == "main") {
    boostrap1 = 
      gamm(y_resample ~ s(delta_brain, bs = "cr"), data = df, random=list(dataset=~1), weights =weight.sub2)
    bs_p = summary(boostrap1$gam)$s.table[[1,4]]
    
  }
  return(bs_p)
}

run_mega(feature, wd)

