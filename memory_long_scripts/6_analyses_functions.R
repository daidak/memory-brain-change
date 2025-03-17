remove.vars = function(dataset){
  vars = c(
    "CSF",
    "Left-Cerebellum-Cortex",
    "Right-Cerebellum-Cortex",
    "Left-Cerebellum-White-Matter",
    "Right-Cerebellum-White-Matter",
    "Left-choroid-plexus",
    "Right-choroid-plexus",
    "Left-vessel",
    "Right-vessel",
    "SupraTentorialVol",
    "SupraTentorialVolNotVent",
    "TotalGrayVol",
    "x3rd-ventricle",
    "x4th-ventricle",
    "xRight-Cerebellum-Cortex",
    "xRight-Cerebellum-White-Matter",
    "Brain-Stem",
    "EstimatedTotalIntraCranialVol",
    "Right-VentralDC", 
    "Left-VentralDC" 
    )
  dataset <- dataset %>% filter(!features %in% vars)
  return(dataset)
  
}

set_links = function(folder) {
  tabulateddata <<- here('data-raw/tabulated')
  datasetdir <<- file.path(tabulateddata, folder)
  outdir <<- here('data_memory_long', folder)
  if (!dir.exists(outdir)) {dir.create(outdir)}
  dataset <<- folder
}


rename_var = function(df, dataset) { 
  g1 = df %>% 
    select(df.harmonize.ap[[dataset]])
  colnms <- names(g1) 
  
  names(g1) = df.harmonize.ap$harmonize
  
  return(g1)
}




wrapper_merge_mri_data = function(wd){
  load("/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/df.merged.lifespan_57K_82sites.Rda")
  load( "/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/df.all.filt.Rda")
  
  grot.ageblmri = 
    df.merge.long %>%
    group_by(rid) %>% 
    summarise(ageblmri = min(age))
  
  x = df$df$imputed.delta[[1]]
  features  = df$misc$phenotypes
  colnames(x) <- features
  
  #remove outliers
  idx = df$outliers$delta$rmsubs.miss.data
  s = df$df$base$df[idx,]
  
  df.mri.change = cbind(s,x) 
  df.mri.change = inner_join(df.mri.change,grot.ageblmri)
  save(df.mri.change,
       file = file.path(wd, "df.mri.rda"))
  return(df.mri.change)
}


wrapper_merge_mem_with_mri_data = function(wd) {
  #import memory data and merge with mri data
  load(file.path(wd, "df.memory_merged_full.rda"))
  load(file.path(wd, "df.mri.rda"))
  
  memory_out = 
    memory.merge.full %>% 
    group_by(rid) %>% 
    summarise(subject = first(subject), 
              dataset = first(dataset), 
              age = mean(age), 
              sex = first(sex), 
              total_tp_memory = max(total_tp_memory), 
              total_followup_memory = max(total_followup_memory), 
              agebl = first(agebl), 
              Int.ranef = first(Int.ranef), 
              memory_endtercept = first(memory_endtercept),
              memory_slope = first(memory_slope), 
              intercept = first(intercept)) %>% 
    filter(!is.na(memory_slope), total_followup_memory > 1.5) 
  
  merged <- inner_join(memory_out, df.mri.change) %>% 
    mutate(agediff = agebl - ageblmri,
           ageF = agebl + total_followup_memory, 
           ageFmri = ageblmri + time, 
           overlap = ifelse(agebl <= ageFmri & ageF >= ageblmri, TRUE, FALSE)) %>% 
    filter(!is.na(n)) %>% 
    filter(!(agebl < 20 | ageblmri < 20)) %>% 
    filter(abs(agediff) < 10 & abs(ageF - ageFmri) < 10) %>% 
    filter(overlap == T)
  
  save(merged,
       file = file.path(wd, "merged.rda"))
  return(merged)
}


wrapper_merge_reliability = function(wd, reliabilitydir) {
  # reliability is grouped by dataset. Not relevant for weights. 
  #load memory data
  load(file.path(wd, "merged.rda"))
  change_sign= c("Left-Inf-Lat-Vent","Left-Lateral-Ventricle","Right-Inf-Lat-Vent","Right-Lateral-Ventricle")
  #load reliability data
  load(file.path(reliabilitydir, "reliability_parameters.rda"))
  
  df = 
    df %>% select(features,
                  dataset,
                  error_variance,
                  icc)
  
  df.merged = 
    merged %>% 
    mutate(icv = scale(EstimatedTotalIntraCranialVol) %>% .[,1]) %>% 
    pivot_longer(unique(df$features), 
                 names_to = "features", 
                 values_to = "delta_brain") %>% 
    mutate(delta_brain = if_else(features %in% change_sign, -delta_brain, delta_brain))
  
  
  df.merged = 
    inner_join(df.merged, df) %>% 
    group_by(dataset, 
             features) %>% 
    mutate(delta_brain = scale(delta_brain) %>% .[,1], #z scores, select only first column since 'scale' returns a matrix
           memory_slope = scale(memory_slope) %>% .[,1], 
           Int.ranef = scale(Int.ranef) %>% .[,1],
           memory_endtercept = scale(memory_endtercept) %>% .[,1],
           intercept = scale(intercept) %>% .[,1]) %>% 
    nest()
  
  save(df.merged,
       file = file.path(wd, "merged_reliability.rda"))
  return(df.merged)
}

merge_apoe = function(wd) {
  load(file.path(wd, "merged.rda"))
  load(file.path(wd,"mega","apoe", "output","APOE_merged.rda"))
  merged_apoe_sub = merged_apoe %>% select(rid, apoe4, apoe_status)
  
  merged = inner_join(merged, merged_apoe_sub)
  save(merged,
       file = file.path(wd, "merged.rda"))
  return(merged)
  
}

prepare_df2 = function(df, p.correct, p.direction, vars_sign) {
  #change_sign= c("Left-Inf-Lat-Vent","Left-Lateral-Ventricle","Right-Inf-Lat-Vent","Right-Lateral-Ventricle")
  #df = df %>% mutate_at(vars(vars_sign), ~if_else(features %in% change_sign, -., .))
  
  #exclude to not be included in FDR correction
  df <- df %>% 
    filter(!features %in% c("SubCortGrayVol", "lh_MeanThickness_thickness","rh_MeanThickness_thickness"))
  
  for (i in p.correct) {
    #print(i)
    df[[sym(paste0(i,".fdr"))]] = p.adjust(df[[sym(i)]], method = "fdr")
    df[[sym(paste0(i,".log.fdr"))]] = -log10(df[[sym(paste0(i,".fdr"))]])
    df[[sym(paste0(i,".log"))]] = -log10(df[[sym(i)]])
  }
  
  for (i in 1:nrow(p.direction)) {
    estimate = p.direction$estimate[i]
    pvalue = p.direction$p.value[i]
    
    df <-df %>%
      mutate(!!sym(paste0(pvalue,".log.fdr")) := 
               if_else(!!sym(estimate) >0,!!sym(paste0(pvalue,".log.fdr")), -!!sym(paste0(pvalue,".log.fdr"))),
             !!sym(paste0(pvalue,".log")) := 
               if_else(!!sym(estimate) >0,!!sym(paste0(pvalue,".log")), -!!sym(paste0(pvalue,".log"))))
  }
  
  return(df)
}



change_infinite_values = function(df, vars.infinite, nreps) {
  for (i in vars.infinite) {
    if (grepl("boot", i)) {
      df = 
        df  %>% 
        mutate(!!sym(i) := if_else(is.infinite(!!sym(i)), -log10(1/nreps), !!sym(i)))
    } else {
      df = 
        df  %>% 
        mutate(!!sym(i) := if_else(is.infinite(!!sym(i)), 16, !!sym(i)))
    }
  }
  return(df)
}

mega_plots_and_derivatives = function(wd) {
  library(gratia)
  library(patchwork)
  temp_folder = here("data_memory_long/all/mega/temp_models")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    fileee <- file.path(temp_folder, mylist[i])
    load(fileee)
    
    # derivatives
    X = derivatives(models$gam.w$gam)
    
    # main gamm
    gs1 = draw(models$gam.w$gam) 
    gs2 = draw(X)
    gs = gs1[[1]] + gs2[[1]]
    fout = paste("gam_w.main", gsub("models_", "", mylist[i]), "png",sep = ".")
    ggsave(filename = file.path(wd, "mega", "plots", gsub(".Rda", "",fout)), plot = gs)
    
    
    gs3 = draw(models$gam_int.w$gam) 
    # full interaction gamm
    fout = paste("gam_w.int", gsub("models_", "", mylist[i]), "png",sep = ".")
    ggsave(filename = file.path(wd, "mega", "plots", gsub(".Rda", "",fout)), plot = gs3)
    
    dt <- df.predict.tensor.interaction %>%
      group_by(xage) %>%
      nest() %>%
      mutate(data = map(data, ~ .x %>%
                          mutate(derivative = c(NA, diff(.$estimate) / diff(.$delta_brain))))) %>% 
      unnest()
    
    dt.plot = 
      dt %>% 
      filter(xage %in% seq(40,80,by = 10)) %>% 
      mutate(xage = as.factor(xage))
    
    gs4 = 
      ggplot(dt.plot, aes(delta_brain, estimate, color = xage)) +
      geom_line(linewidth = 2) +
      geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), color = "ivory", alpha = .5) + 
      geom_line(mapping = aes(y = derivative*2), color = "lightblue", linewidth = 1.5, alpha = .5) +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      facet_grid(cols = vars(xage)) + 
      scale_y_continuous(
        name = "estimate", 
        sec.axis = sec_axis(~ . / 2, name = "derivative")
      ) +
      theme_classic() + 
      theme(legend.position = 'none')
     fout = paste("gam_w.byage", gsub("models_", "", mylist[i]), "png",sep = ".")
      ggsave(filename = file.path(wd, "mega", "plots", gsub(".Rda", "",fout)), plot = gs4)
  }
}

wrapper_get_predictions_derivative = function(wd) {
  temp_folder = here("data_memory_long/all/mega/temp_models")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    print(i)
    #i = 4
    fileee <- file.path(temp_folder, mylist[i])
    load(fileee)
    features = gsub(".Rda", "", mylist[i]) %>% gsub("models_", "", .)
    d = density(models$gam.w$gam$model$delta_brain)
    
    # predict derivative and cumulative distribution in main model
    X = derivatives(models$gam.w$gam)
    X$w = approx(d$x, d$y, xout = X$delta_brain)$y
    X$features = features
    X = X %>% 
      group_by(features) %>% 
      nest()
    
    # predict derivative from gamm interaction along delta brain by specific ages
    dt <- df.predict.tensor.interaction %>%
      group_by(xage) %>%
      nest() %>%
      mutate(data = map(data, ~ .x %>%
                          mutate(derivative = c(NA, diff(.$estimate) / diff(.$delta_brain))))) %>% 
      unnest() 
    
    dt$w = approx(d$x, d$y, xout = dt$delta_brain)$y    
    dt$features = features 
    dt = dt %>% 
      group_by(features) %>% 
      nest()
    
    predictions = left_join(X, dt, by = "features", suffix = c("main", "interaction"))
    grot[[i]] =predictions
  }
  return(grot)
}



wrapper_draw_smooth = function(wd) {
  temp_folder = here("data_memory_long/all/mega/temp_models")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    print(i)
    #i = 4
    fileee <- file.path(temp_folder, mylist[i])
    load(fileee)
    features = gsub(".Rda", "", mylist[i]) %>% gsub("models_", "", .)
    
    predict.main = smooth_estimates(models$gam.w)    
    predict.main$features = features
    df.predict.tensor.interaction$features = features
    predict.main = predict.main %>% group_by(features) %>% nest()
    predict.int = df.predict.tensor.interaction %>% group_by(features) %>% nest()
    
    predictions = left_join(predict.main, predict.int, by = "features", suffix = c("mainsmooth", "interactionsmooth"))
    grot[[i]] =predictions
  }
  return(grot)
}



wrapper_get_estimates_apoe = function(mega) {     
  mega = 
    mega %>% 
    mutate(gamm.w.apoe0.beta.m5.0.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,0), apoe_statusO == 0)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}), 
           gamm.w.apoe0.beta.m5.5.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,4.5), apoe_statusO == 0)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}),
           gamm.w.apoe0.beta.m0.5.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -0,4.5), apoe_statusO == 0)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}),
           gamm.w.apoe1.beta.m5.0.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,0), apoe_statusO == 1)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}), 
           gamm.w.apoe1.beta.m5.5.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,4.5), apoe_statusO == 1)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}),
           gamm.w.apoe1.beta.m0.5.weighted = 
             map_dbl(data.apoeO.twogroups, ~ { 
               grot = .x %>% filter(between(delta_brain, -0,4.5), apoe_statusO == 1)
               weighted.mean(grot$derivative, grot$w, na.rm = T)}),
           gamm.w.apoeInt.beta.m5.0.weighted = 
             map_dbl(data.apoeO.int, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,0), .by == "apoe_statusO")
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}), 
           gamm.w.apoeInt.beta.m5.5.weighted = 
             map_dbl(data.apoeO.int, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,4.5), .by == "apoe_statusO")
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}),
           gamm.w.apoeInt.beta.m0.5.weighted = 
             map_dbl(data.apoeO.int, ~ { 
               grot = .x %>% filter(between(delta_brain, -0,4.5), .by == "apoe_statusO")
               weighted.mean(grot$.derivative, grot$w, na.rm = T)})
    )
  return(mega)
}





wrapper_get_estimates = function(mega) {     
  mega = 
    mega %>% 
    mutate(gamm.w.beta.m2.0.weighted = 
             map_dbl(datamain, ~ { 
               grot = .x %>% filter(between(delta_brain, -2,0))
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}), 
           gamm.w.beta.m5.0.weighted = 
             map_dbl(datamain, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,0))
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}), 
           gamm.w.beta.m5.5.weighted = 
             map_dbl(datamain, ~ { 
               grot = .x %>% filter(between(delta_brain, -4.5,4.5))
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}),
           gamm.w.beta.m0.5.weighted = 
             map_dbl(datamain, ~ { 
               grot = .x %>% filter(between(delta_brain, -0,4.5))
               weighted.mean(grot$.derivative, grot$w, na.rm = T)}),
           gamm.w.int.beta.m2.0.weighted = 
             map(datainteraction, ~ {
               .x %>% filter(between(delta_brain, -2,0)) %>% 
                 group_by(xage) %>% 
                 summarise(estimate = weighted.mean(derivative, w, na.rm = T))}), 
           gamm.w.int.beta.m5.0.weighted = 
             map(datainteraction, ~ {
               .x %>% filter(between(delta_brain, -4.5,0)) %>% 
                 group_by(xage) %>% 
                 summarise(estimate = weighted.mean(derivative, w, na.rm = T))}),
           gamm.w.int.beta.m5.5.weighted = 
             map(datainteraction, ~ {
               .x %>% filter(between(delta_brain, -4.5,4.5)) %>% 
                 group_by(xage) %>% 
                 summarise(estimate = weighted.mean(derivative, w, na.rm = T))}),
           gamm.w.int.beta.m0.5.weighted = 
             map(datainteraction, ~ {
               .x %>% filter(between(delta_brain, 0,4.5)) %>% 
                 group_by(xage) %>% 
                 summarise(estimate = weighted.mean(derivative, w, na.rm = T))})
    )
  return(mega)
}


wrapper_hemispheric_differences = function(wd) {
  load(file.path(wd,"mega", "df.mega.rda"))
  beta.estimates = names(mega)[grepl("beta", names(mega))]
  
  dat2 = 
    mega %>% 
    separate(features, c("hemi", "feature2"), remove = F, extra = "merge") %>% 
    filter(!hemi =="SubCortGrayVol") %>% 
    mutate(hemi = if_else(hemi == "Left", "lh", 
                          if_else(hemi == "Right", "rh", hemi))) %>% 
    ungroup() %>%
    select(hemi,
           feature2,
           beta.estimates) %>% 
    pivot_wider(names_from = hemi, 
                values_from = -c(hemi, feature2)) 
  
  
  ## are significantly different the estimates from left and right?
  # yes they are
  (mod.hemi.main.2.0 =t.test(dat2$gamm.w.beta.m2.0.weighted_lh, dat2$gamm.w.beta.m2.0.weighted_rh, paired = T, na.rm = T))
  (mod.hemi.main.5.0 =t.test(dat2$gamm.w.beta.m5.0.weighted_lh, dat2$gamm.w.beta.m5.0.weighted_rh, paired = T, na.rm = T))
  
  x = 
    dat2 %>% 
    select(feature2,
           starts_with("gamm.w.int.beta")) %>% 
    unnest() %>% 
    group_by(xage) %>% 
    nest() %>% 
    mutate(ttest.age20 = map(data, ~ t.test(.x$estimate, .x$estimate1, paired = T, na.rm = T)), 
           ttest.age50 = map(data, ~ t.test(.x$estimate2, .x$estimate3, paired = T, na.rm = T)), 
           t.test.age20.t_statistic = map_dbl(ttest.age20, ~.x$statistic), 
           t.test.age50.t_statistic = map_dbl(ttest.age50, ~.x$statistic), 
           t.test.age20.p.value = map_dbl(ttest.age20, ~.x$p.value), 
           t.test.age50.p.value = map_dbl(ttest.age50, ~.x$p.value), 
           t.test.age20.mean_difference = map_dbl(ttest.age20, ~.x$estimate), 
           t.test.age50.mean_difference = map_dbl(ttest.age50, ~.x$estimate) 
    )
  
  hemi.diff = list()  
  hemi.diff$dat2 = dat2
  hemi.diff$x = x
  hemi.diff$mod.hemi.main.2.0 = mod.hemi.main.2.0
  hemi.diff$mod.hemi.main.5.0 = mod.hemi.main.5.0
  save(hemi.diff,
       file = file.path(wd,"mega","hemispheric_diff", "hemi.rda"))
  
  return(hemi.diff)
}


wrapper_clustering = function(wd) {
  library(M3C)
  load(file.path(wd, "mega", "df.mega.rda"))
  
  sig.features = mega %>% filter(gamm.w.pval.bootstrap.fdr < 0.05)%>% .$features
  
  #load(file.path(wd, "df.change.rda"))
  load(file.path(wd, "merged_reliability.rda"))
  
  M = df.merged %>% 
    unnest(data) %>% 
    ungroup() %>% 
    dplyr::select(rid, features, delta_brain) %>% 
    filter(features %in% sig.features) %>% 
    group_by(features) %>% 
    mutate(delta_brain = scale(delta_brain)%>% .[,1], 
           delta_brain = if_else(between(delta_brain, -4.5,4.5), delta_brain, 0)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = features, values_from = delta_brain)
  
  
  rid = M[,1]
  M = M[,-1]
  
  rownames(M) = rid$rid
  
  # run clustering. different options. suggest M3C
  mod.m3c = M3C(M,
                iters = 100, 
                maxK = 10, 
                pItem = .8,
                repsref = 50,
                repsreal = 50, 
                clusteralg = "spectral",
                seed = 123,
                objective = 'entropy',
                removeplots = F,
                method = 1)
  
  m3c = list()
  m3c$mod = mod.m3c 
  m3c$M = M
  
  save(m3c, 
       file = file.path(wd,"mega","consensus_cluster", "cluster_mod.rda"))
  return(m3c)
}



wrapper_pca = function(wd) {
  load(file.path(wd, "mega", "df.mega.rda"))
  
  sig.features = mega %>% filter(gamm.w.pval.bootstrap.fdr < 0.05)%>% .$features
  
  #load(file.path(wd, "df.change.rda"))
  load(file.path(wd, "merged_reliability.rda"))
  
  M = df.merged %>% 
    ungroup() %>% 
    unnest(data) %>% 
    select(rid, features, delta_brain) %>% 
    filter(features %in% sig.features) %>% 
    group_by(features) %>% 
    mutate(delta_brain = scale(delta_brain)%>% .[,1], 
           delta_brain = if_else(between(delta_brain, -4.5,4.5), delta_brain, 0)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = features, values_from = delta_brain)
  
  rid = M[,1]
  M = M[,-1]
  rownames(M) = rid$rid
  
  #PCA:
  pca_result= prcomp(M, scale = T, center = T )
  
  save(pca_result,
       file = file.path(wd,"mega","PCA", "pca_result.rda"))
  return(pca_result)
}





apoe_plots = function(wd) {
  
  load(file.path(wd,"mega","apoe", "df.mega.apoe.rda"))
  load(file.path(wd,"mega", "df.mega.rda"))
  library(gratia)
  library(patchwork)
  temp_folder = here("data_memory_long/all/mega/apoe/output/temp_models")
  plotdir = here("data_memory_long/all/mega/apoe/plots")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    fileee <- file.path(temp_folder, mylist[i])
    f = gsub("models_","", mylist[i]) %>% gsub(".Rda", "",.)
    if (f %in% c("SubCortGrayVol", "lh_MeanThickness_thickness","rh_MeanThickness_thickness")) {
      next
    }
    load(fileee)
    print(i)
    print(f)
    smoothmain = mega %>% filter(features == f) %>% .$datamainsmooth %>% .[[1]]
    smoothint = mega %>% filter(features == f) %>% .$datainteraction %>% .[[1]]
    
    
    # Predict apoe pos and neg from ApoeO model. Predict also difference in smooths
    mod.gam_main.w.apoeO = models$mod.gam_main.w.apoeO
    grid3 <- expand.grid(delta_brain = seq(-4.5,4.5, length.out = 100),
                         apoe_statusO = c(0,1))
    pred.apoes = predict(mod.gam_main.w.apoeO$gam,grid3, se.fit = T)
    
    predict.apoO.update = 
      data.frame(
        delta_brain = grid3$delta_brain, 
        apoe_statusO = grid3$apoe_statusO %>% as.factor(), 
        .estimate = pred.apoes$fit, 
        .se = pred.apoes$se.fit
      )
    
    x = predictions[[5]] %>% 
      mutate(apoe_statusO = if_else(is.na(apoe_statusO), 0, 1)  %>% as.factor())
    gs1 = ggplot(predict.apoO.update, aes(delta_brain, .estimate, group = apoe_statusO, color = apoe_statusO, fill = apoe_statusO) ) + 
      geom_line() + 
      geom_ribbon(aes(ymin = .estimate - 1.96*.se, ymax = .estimate +1.96*.se ), alpha = .3, lwd = 0) + 
      geom_line(data = x %>% filter(apoe_statusO == 1), color = "brown") + 
      geom_ribbon(data = x %>% filter(apoe_statusO == 1), 
                  mapping = aes(ymin = .estimate - 1.96*.se, ymax = .estimate +1.96*.se ), fill = "brown",  alpha = .3, lwd = 0) +
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") + 
      theme_classic() + 
      theme(legend.position = "bottom")
    
    gs2 = ggplot(smoothmain, aes(delta_brain, .estimate) ) + 
      geom_line() + 
      geom_ribbon(aes(ymin = .estimate - 1.96*.se, ymax = .estimate +1.96*.se ), alpha = .3, lwd = 0) + 
      geom_line(data = predict.apoO.update, aes( group = apoe_statusO, color = apoe_statusO))  +
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") + 
      theme_classic() + 
      theme(legend.position = "bottom")
    
    
    
    #INteraction #2 
    xx =  predictions[[2]] %>% 
      mutate(apoe_statusO = apoe_statusO  %>% as.factor()) %>% 
      filter(xage %in% c(40,50,60,70,80))
    gs3 = ggplot(xx , aes(delta_brain, estimate, group = apoe_statusO, color = apoe_statusO, fill = apoe_statusO) ) + 
      geom_line() + 
      geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate +1.96*se ), alpha = .3, lwd = 0) + 
      facet_grid(cols = vars(xage))
    
    smoothint = smoothint %>%
      filter(xage %in% c(40,50,60,70,80))
    gs4 = ggplot(smoothint , aes(delta_brain, estimate) ) + 
      geom_line() + 
      geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate +1.96*se ), alpha = .3, lwd = 0) + 
      geom_line(data = xx, aes(group = apoe_statusO, color = apoe_statusO)) + 
      facet_grid(cols = vars(xage))
    
    gs5 = gs1 + gs2
    gs6 = gs3/gs4
    
    ggsave(filename = file.path(plotdir, paste0("main.", f, ".png" )), plot = gs5)
    ggsave(filename = file.path(plotdir, paste0("int.", f, ".png" )), plot = gs6, width = 10)
    # ApoeO
    # x = predictions[[5]] %>% 
    #   mutate(apoe_statusO = if_else(is.na(apoe_statusO), 0, 1)  %>% as.factor())
    # ggplot(x, aes(delta_brain, .estimate, group = apoe_statusO, color = apoe_statusO, fill = apoe_statusO) ) + 
    #   geom_line() + 
    #   geom_ribbon(aes(ymin = .estimate - 1.96*.se, ymax = .estimate +1.96*.se ), alpha = .3, lwd = 0)
    
    # ApoeNeg + Pos
    # y = predictions[[3]] %>% 
    #   mutate(apoe_statusO = apoe_statusF %>% as.factor())
    # ggplot(y, aes(delta_brain, .estimate, group = apoe_statusO, color = apoe_statusO, fill = apoe_statusO) ) + 
    #   geom_line() + 
    #   geom_ribbon(aes(ymin = .estimate - 1.96*.se, ymax = .estimate +1.96*.se ), alpha = .3, lwd = 0)
    #ggsave(filename = file.path(wd, "mega", "plots", gsub(".Rda", "",fout)), plot = gs)
  }
}



wrapper_get_predictions_derivative.apoe = function(wd) {
  temp_folder = here("data_memory_long/all/mega/apoe/output/temp_models")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    print(i)
    #i = 4
    fileee <- file.path(temp_folder, mylist[i])
    load(fileee)
    features = gsub(".Rda", "", mylist[i]) %>% gsub("models_", "", .)
    d = density(models$mod.gam_main.w.apoeO$gam$model$delta_brain)
    
    # predict derivative and cumulative distribution in main model
    X = derivatives(models$mod.gam_main.w.apoeO$gam)
    X$w = approx(d$x, d$y, xout = X$delta_brain)$y
    X$features = features
    X = X %>% 
      group_by(features) %>% 
      nest() %>% 
      rename("data.apoeO.int" = "data")
    
    # predict derivative and cumulative distribution in main model
    X1 = derivatives(models$mod.gam_main.w.apoe1$gam)
    X1$w = approx(d$x, d$y, xout = X1$delta_brain)$y
    X1$features = features
    X1 = X1 %>% 
      group_by(features) %>% 
      nest() %>% 
      rename("data.apoe1" = "data")
    
    # predict derivative and cumulative distribution in main model
    X0 = derivatives(models$mod.gam_main.w.apoe0$gam)
    X0$w = approx(d$x, d$y, xout = X0$delta_brain)$y
    X0$features = features
    X0 = X0 %>% 
      group_by(features) %>% 
      nest() %>% 
      rename("data.apoe0" = "data")
    
    
    
    
    # Predict apoe pos and neg from ApoeO model. Predict also difference in smooths
    mod.gam_main.w.apoeO = models$mod.gam_main.w.apoeO
    grid3 <- expand.grid(delta_brain = seq(-4.5,4.5, length.out = 100),
                         apoe_statusO = c(0,1))
    pred.apoes = predict(mod.gam_main.w.apoeO$gam,grid3, se.fit = T)
    
    X2= 
      data.frame(
        delta_brain = grid3$delta_brain, 
        apoe_statusO = grid3$apoe_statusO %>% as.factor(), 
        .estimate = pred.apoes$fit, 
        .se = pred.apoes$se.fit
      ) %>% 
      group_by(apoe_statusO) %>% 
      nest() %>%
      mutate(data = map(data, ~ .x %>%
                          mutate(derivative = c(NA, diff(.$.estimate) / diff(.$delta_brain))))) %>% 
      unnest() %>% ungroup() 
    
    X2$w = approx(d$x, d$y, xout = X2$delta_brain)$y
    X2 = X2 %>% nest() %>% 
      rename("data.apoeO.twogroups" = "data")
    X2$features = features
    
    predictions = reduce(list(X0, X1, X, X2), left_join)
    
    grot[[i]] =predictions
  }
  return(grot)
}




wrapper_draw_smooth_apoe = function(wd) {
  temp_folder = here("data_memory_long/all/mega/apoe/output/temp_models")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[grepl("^models_", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    print(i)
    #i = 4
    fileee <- file.path(temp_folder, mylist[i])
    load(fileee)
    
    predict.apoeO = predictions[[5]]  %>% group_by(features) %>% nest() %>% 
      rename("data.predict.main.apoeO" = "data")
    predict.int.apoeo = predictions[[2]]  %>% group_by(features) %>% nest()%>% 
      rename("data.predict.int.apoeO" = "data")
    predict.int.apoeo.bygroup = predictions[[1]]  %>% group_by(features) %>% nest()%>% 
      rename("data.predict.int.apoeO.bygroup" = "data")
    
    grot[[i]] = reduce(list(predict.apoeO,
                            predict.int.apoeo, 
                            predict.int.apoeo.bygroup), 
                       left_join)
  }
  return(grot)
}



age_plot_apoe_brain = function() {
  load(file.path(wd,"mega","apoe", "df.mega.apoe.brain.rda"))
  library(gratia)
  library(patchwork)
  temp_folder = here("data_memory_long/all/mega/apoe/output/temp_models_mem")
  plotdir = here("data_memory_long/all/mega/apoe/plots_brain")
  mylist = list.files(temp_folder, pattern = ".Rda")
  mylist = mylist[!mylist %in% c("models.Rda", "output.Rda")]
  mylist = mylist[grepl("^models.", mylist)]
  grot = list()
  for (i in 1:length(mylist)) {
    fileee <- file.path(temp_folder, mylist[i])
    f = gsub("models.","", mylist[i]) %>% gsub(".Rda", "",.)
    if (f %in% c("SubCortGrayVol", "lh_MeanThickness_thickness","rh_MeanThickness_thickness")) {
      next
    }
    load(fileee)
    
    
    
    y = predictions[[1]] %>% 
      mutate(apoe_statusO = if_else(is.na(apoe_statusO), 0, 1)  %>% as.factor()) %>% 
      rename(estimate = .estimate, 
             se = .se) %>% 
      filter(apoe_statusO == 1)
    x = predictions[[2]] %>% 
      mutate(apoe_statusO = as.factor(apoe_statusO)) 
    gs1 = ggplot(x, aes(age, estimate, group = apoe_statusO, color = apoe_statusO, fill = apoe_statusO) ) + 
      geom_line() + 
      geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate +1.96*se ), alpha = .3, lwd = 0) + 
      geom_line(data = y, color = "brown") + 
      geom_ribbon(data = y, 
                  mapping = aes(ymin = estimate - 1.96*se, ymax = estimate +1.96*se ), fill = "brown",  alpha = .3, lwd = 0) +
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") + 
      theme_classic() + 
      theme(legend.position = "bottom")
    
    ggsave(filename = file.path(plotdir, paste0("int.", f, ".png" )), plot = gs1)
    
  }
}


retrieve_pboot_cc = function(wd, type ="basic") {
  
  wd = here("data_memory_long/all")
  
  if (!type == "all") {
    if (type == "pca") {
      load(file.path(wd,"mega","consensus_cluster", "pca_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(1:8) %>% as.character())
    } else if (type == "cl1") {
      load(file.path(wd,"mega","consensus_cluster", "cl1_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(2:8) %>% as.character())
    } else if (type == "basic") {
      load(file.path(wd,"mega","consensus_cluster", "basic_consensus_cluster_gam.Rda"))
      cl = paste0("cl", c(1:8) %>% as.character())
      
    }
    
    trues = ls(pattern = "mod.")
    
    pvalues = lapply(trues, function(x) {summary(get(x)$gam)$s.table[1,4]}) %>% 
      simplify2array()
    
    df = 
      data.frame(cl = cl, 
                 pval = pvalues)
    
    df = 
      bootstraps %>%
      filter(!is.na(.))%>% 
      group_by(cl) %>% 
      nest() %>% 
      left_join(df,.)
  } else if (type == "all") {
    load(file.path(wd,"mega","consensus_cluster", "all_consensus_cluster_gam.Rda"))
    df = broom::tidy(mod$gam) %>% 
      rename(cl = term, 
             pval = "p.value")
    bootstraps = bootstraps %>% 
      filter(term == "s(cl.1)" & null == "g.null1" | 
               term == "s(cl.2)" & null == "g.null2" |
               term == "s(cl.3)" & null == "g.null3" |
               term == "s(cl.4)" & null == "g.null4" |
               term == "s(cl.5)" & null == "g.null5" |
               term == "s(cl.6)" & null == "g.null6" |
               term == "s(cl.7)" & null == "g.null7" |
               term == "s(cl.8)" & null == "g.null8") %>% 
      rename(cl = term)
    
    df = 
      bootstraps %>%
      filter(!is.na(.$`p-value`))%>% 
      group_by(cl) %>% 
      nest() %>% 
      left_join(df,.)
    
  }
    
    df = 
      df %>% 
      mutate(pboot = 
               map2_dbl(pval, data, ~sum(.x > .y) / length(.y[[1]]))) %>% 
      select(-data)
  
  
  return(df)
}

wrapper_get_weights_cc_analysis = function(wd, type) {
  library(gratia)
  
   if (type == "pca") {
      load(file.path(wd,"mega","consensus_cluster", "pca_consensus_cluster_gam.Rda"))
      mods = ls(pattern ="^mod\\..*\\.pca$")
      grot = list()
      for (i in mods) {
        model = get(i)
        grot[[i]] = 
          derivatives(model$gam) %>% 
          filter(is.na(pca)) %>% 
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        dt = model$gam$model %>% 
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        d = density(dt$delta_brain)
        grot[[i]]$w = approx(d$x, d$y, xout = grot[[i]]$delta_brain)$y
      }
    } else if (type == "cl1") {
      load(file.path(wd,"mega","consensus_cluster", "cl1_consensus_cluster_gam.Rda"))
      mods = ls(pattern ="^mod\\..*\\.cl1$")
      grot = list()
      for (i in mods) {
        model = get(i)
        grot[[i]] = 
          derivatives(model$gam) %>% 
          filter(is.na(cl.1)) %>% 
          rename(grot ="cl.1") %>% 
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        dt = model$gam$model %>%
          rename("grot" = "cl.1") %>% 
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        d = density(dt$delta_brain)
        grot[[i]]$w = approx(d$x, d$y, xout = grot[[i]]$delta_brain)$y
      }
    } else if (type == "basic") {
      load(file.path(wd,"mega","consensus_cluster", "basic_consensus_cluster_gam.Rda"))
      mods = ls(pattern ="^mod\\..*\\.b$")
      grot = list()
      for (i in mods) {
        model = get(i)
        grot[[i]] = 
          derivatives(model$gam) %>% 
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        dt = model$gam$model %>%
          rename_with(~ "delta_brain", .cols = starts_with("cl."))
        d = density(dt$delta_brain)
        grot[[i]]$w = approx(d$x, d$y, xout = grot[[i]]$delta_brain)$y
      }
    } else if (type == "all") {
      
      load(file.path(wd,"mega","consensus_cluster", "all_consensus_cluster_gam.Rda"))
      model = mod
      X = derivatives(model$gam)
      cls = names(X)[grepl("cl", names(X))]
      grot = list()
      for (i in cls) {
        grot[[i]] = X %>% filter(!is.na(get(i))) %>% 
          rename("delta_brain" = i) %>% 
          dplyr::select(-starts_with("cl."))
        d = density(model$gam$model[[i]])
        grot[[i]]$w = approx(d$x, d$y, xout = grot[[i]]$delta_brain)$y
      }
      dt = data.table::rbindlist(grot, idcol ="model")
      
      load(file.path(wd,"mega","consensus_cluster", "all_consensus_cluster_gam.Rda"))
      model = mod
      X = derivatives(model$gam)
      cls = names(X)[grepl("cl", names(X))]
      grot = list()
      for (i in cls) {
        grot[[i]] = X %>% filter(!is.na(get(i))) %>% 
          rename("delta_brain" = i) %>% 
          dplyr::select(-starts_with("cl."))
        d = density(model$gam$model[[i]])
        grot[[i]]$w = approx(d$x, d$y, xout = grot[[i]]$delta_brain)$y
      }
    }
    
    dt = data.table::rbindlist(grot, idcol = "model")
    dt = dt %>% 
      group_by(model) %>% 
      nest() %>% 
      mutate(gamm.w.beta.m2.0.weighted = 
               map_dbl(data, ~ { 
                 grot = .x %>% filter(between(delta_brain, -2,0))
                 weighted.mean(grot$.derivative, grot$w, na.rm = T)}), 
             gamm.w.beta.m5.0.weighted = 
               map_dbl(data, ~ { 
                 grot = .x %>% filter(between(delta_brain, -4.5,0))
                 weighted.mean(grot$.derivative, grot$w, na.rm = T)}), 
             gamm.w.beta.m5.5.weighted = 
               map_dbl(data, ~ { 
                 grot = .x %>% filter(between(delta_brain, -4.5,4.5))
                 weighted.mean(grot$.derivative, grot$w, na.rm = T)}),
             gamm.w.beta.m0.5.weighted = 
               map_dbl(data, ~ { 
                 grot = .x %>% filter(between(delta_brain, -0,4.5))
                 weighted.mean(grot$.derivative, grot$w, na.rm = T)}))
    
    dt = dt %>% select(-data)
    return(dt)
  }
  
  
#### Simulation functions ####
model_twoModel_discussion = 
  function(n = 1000, x1, sd1, skew1, sd2, r, ysd, tau1=0) {
    library(gratia)
    library(mgcv)
    
    # n number of samples
    # x1 mean aging distribution
    # sd1 sd aging distribution
    # skew1 skew aging dist. 
    # sd noise distribution x2
    # r relationship between aging and mem
    # ysd sd noise memory
    
    X1 =sn::rsn(n, xi=x1, omega=sd1, alpha=skew1, tau=tau1,  dp=NULL)
    X2 =sn::rsn(n, xi=0, omega=sd2, alpha=0, tau=0,  dp=NULL) 
    Ynoise <- rnorm(n, mean = 0, sd = ysd)
    
    Y <- X1 * 0.3 + Ynoise
    Y <- as.numeric(Y)  # Remove attributes from scaling
    
    dataset <- data.frame(X1, X2, Y)
    dataset = dataset %>% mutate(X = X1 + X2)
    
    gam_model <- gam(Y ~ s(X), data = dataset)
    gam_aging <- gam(Y ~ s(X1), data = dataset)
    lm_model <- lm(Y ~ X, data = dataset)
    lm_aging <- lm(Y ~ X1, data = dataset)
    
    grot = list()
    grot$tidy.gam_aging = broom::tidy(gam_aging)
    grot$tidy.gam.model = broom::tidy(gam_model)
    grot$tidy.lm_aging = broom::tidy(lm_aging)
    grot$tidy.lm_model = broom::tidy(lm_model)
    
    grot$smooth.gam.model  = smooth_estimates(gam_model)
    grot$smooth.gam.aging  = smooth_estimates(gam_aging)
    return(grot)
  }

summary_tidy = function(mod) {
  grot = 
    rbind(
      mod %>% 
        summarise_if(is.numeric, mean), 
      mod %>% 
        summarise_if(is.numeric, sd))
  grot$meas = c("mean", "sd")
  return(grot)
}

tidy_stats = function(out) {
  
  tidy.gam_aging = lapply(out, function(x) {x$tidy.gam_aging}) %>% 
    data.table::rbindlist()
  tidy.gam_model = lapply(out, function(x) {x$tidy.gam.model}) %>% 
    data.table::rbindlist()
  tidy.lm_aging = lapply(out, function(x) {x$tidy.lm_aging}) %>% 
    data.table::rbindlist()
  tidy.lm_model = lapply(out, function(x) {x$tidy.lm_model}) %>% 
    data.table::rbindlist()
  df.gam_aging.summary = summary_tidy(tidy.gam_aging)
  df.gam_model.summary = summary_tidy(tidy.gam_model)
  df.lm_aging.summary = summary_tidy(tidy.lm_aging %>% filter(term =="X1"))
  df.lm_model.summary = summary_tidy(tidy.lm_model %>% filter(term =="X"))
  
  df = list()
  df$gam_aging.summary = df.gam_aging.summary
  df$gam_model.summary = df.gam_model.summary
  df$lm_aging.summary = df.lm_aging.summary
  df$lm_model.summary = df.lm_model.summary
  return(df)
}

get_mean_smooths = function(out, smooth = "model") {
  if (smooth == "model") {
    smt = 
      lapply(out, function(x) {x$smooth.gam.model}) 
  } else if (smooth == "aging") {
    smt = 
      lapply(out, function(x) {x$smooth.gam.aging}) 
  }
  smt = 
    smt %>% 
    data.table::rbindlist(idcol = "it") 
  
  if (smooth == "aging") {
    smt = smt %>% 
      rename(X = X1)
  }
  
  
  smt = 
    smt %>% 
    group_by(it) %>% 
    nest() %>% 
    mutate(G = smooth)
  
  range = 
    smt %>% 
    unnest() %>% 
    mutate(minX = min(X), 
           maxX = max(X)) %>% 
    ungroup() %>% 
    summarise(liX = max(minX), 
              upX = min(maxX))
  
  df.grid = 
    data.frame(
      X = seq(range$liX, range$upX, by = .01), 
      G = smooth
    ) %>% 
    group_by(G) %>% 
    nest() %>% 
    rename("predict" = "data")
  
  smt= left_join(smt, df.grid)
  
  smt = 
    smt %>% 
    mutate(predict = map2(data, predict, ~ {
      .y %>% mutate(.estimate = approx(.x$X, .x$.estimate, xout = .y$X)$y)
    }))
  
  df.model = 
    smt %>% 
    dplyr::select(predict) %>% 
    unnest() %>% 
    ungroup() %>% 
    group_by(X) %>% 
    summarise(mean = mean(.estimate), 
              li = quantile(.estimate,probs = .025), 
              ui = quantile(.estimate,probs = .975)) %>%
    mutate(G = smooth)
  
  o = list(smt, 
           df.model)
  return(o)
}

twoModel_distribution_example = function(n = 1000, x1, sd1, skew1, sd2, r, ysd) {
  library(gratia)
  library(mgcv)
  
  
  X1 =sn::rsn(n, xi=x1, omega=sd1, alpha=skew1, tau=0,  dp=NULL)
  X2 =sn::rsn(n, xi=0, omega=sd2, alpha=0, tau=0,  dp=NULL) 
  Ynoise <- rnorm(n, mean = 0, sd = ysd)
  
  Y <- X1 * 0.3 + Ynoise
  Y <- as.numeric(Y)  # Remove attributes from scaling
  
  dataset <- data.frame(X1, X2, Y, Ynoise)
  dataset = dataset %>% mutate(X = X1 + X2)
  
  return(dataset)
}


adjust_omega <- function(alpha, target_sd) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- target_sd * sqrt(pi / (pi - 2 * delta^2))
  return(omega)
}
