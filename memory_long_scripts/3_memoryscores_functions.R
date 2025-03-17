set_links = function(folder) {
  datadir <<- here('data_memory_long', folder, 'other_output')
  outdir <<- here('data_memory_long', folder)
  #if (!dir.exists(outdir)) {dir.create(outdir)}
}

# scale memory score function
scale.mem.var = function(df, mem_var){
  var1 <- enquo(mem_var)
  
  df <- df %>% 
    group_by(subject_id) %>%
    arrange(timepoint) %>%
    mutate(first_mem = first(!!var1)) 
  
  #mean and SD of first timepoint
  first_mem_mean = mean(df$first_mem, na.rm = T)
  first_mem_sd = sd(df$first_mem, na.rm = T)
  
  colnm <- paste(get_expr(var1), "scaled", sep = "_")
  
  df <- df %>% 
    mutate(!! colnm := (!!var1-first_mem_mean)/first_mem_sd) %>% 
    #mutate(!! colnm := !!var1-first_mem_mean/first_mem_sd) %>% 
    select(-first_mem)
  return(df)
}


#uio
memoryscore.uio = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, CVLT_short_delay)
  df = scale.mem.var(df, CVLT_delayed)
  df = scale.mem.var(df, CVLT_total_learning)
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(CVLT_short_delay_scaled, CVLT_delayed_scaled,CVLT_total_learning_scaled)
  memvar_all <- df %>% 
    ungroup() %>% 
    select(CVLT_short_delay_scaled, CVLT_delayed_scaled,CVLT_total_learning_scaled)
  first_tp = imputePCA(first_tp, ncp = 2)$completeObs
  memvar_all = imputePCA(memvar_all, ncp = 2)$completeObs
  
  #PCA on first timepoint
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  #PCA_tp1_2 <- prcomp(first_tp[, c(16,17,19)], center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all )
  
  df$PCA = PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
  
}


#umu
memoryscore.umu = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, VT_free_recall_vtb)
  df$PCA = df$VT_free_recall_vtb_scaled
  save(df,
       file = file.path(outdir, "df.memory.rda"))
}


#ub
memoryscore.ub = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  
  
  df = scale.mem.var(df, RAVLT_total_learning)
  df = scale.mem.var(df, RAVLT_delayed)
  
  # no missing data
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(RAVLT_total_learning_scaled, RAVLT_delayed_scaled)
  memvar_all <- df %>% 
    select(RAVLT_total_learning_scaled, RAVLT_delayed_scaled)
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all[, -1] )
  
  df$PCA = PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
  # library(factoextra)
  # var <- get_pca_var(PCA_tp1)
  # # Contributions of variables to PC1
  # fviz_contrib(PCA_tp1, choice = "var", axes = 1, top = 10)
}


#mpib
memoryscore.mpib = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  #df = scale.mem.var(df, VLMT_learning)
  df = scale.mem.var(df, VLMT_total_learning)
  df = scale.mem.var(df, VLMT_short_delay)
  df = scale.mem.var(df, VLMT_delayed)
  df = scale.mem.var(df, io_hit_minus_fa)
  df = scale.mem.var(df, fp_hit_minus_new_fa)
  df = scale.mem.var(df, fp_hit_minus_rear_fa)
  df = scale.mem.var(df, ol_acc)

  
  #df = scale.mem.var(df, VLMT_recognition_rates)
  #df = scale.mem.var(df, VLMT_true_recognition_rates)
  
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup()%>% 
    select(VLMT_total_learning_scaled,
           VLMT_short_delay_scaled,
           VLMT_delayed_scaled,
           io_hit_minus_fa_scaled, 
           fp_hit_minus_new_fa_scaled, 
           fp_hit_minus_rear_fa_scaled, 
           ol_acc_scaled)
  memvar_all <- df %>% 
    ungroup() %>% 
    select(VLMT_total_learning_scaled,
           VLMT_short_delay_scaled,
           VLMT_delayed_scaled,
           io_hit_minus_fa_scaled, 
           fp_hit_minus_new_fa_scaled, 
           fp_hit_minus_rear_fa_scaled, 
           ol_acc_scaled)
  
  first_tp = imputePCA(first_tp, ncp = 2)$completeObs
  memvar_all = imputePCA(  memvar_all, ncp = 2)$completeObs
  PCA_tp1 <- prcomp(first_tp, center = T, scale = T)
  PCA_all <- predict(PCA_tp1, newdata = memvar_all )
  df$PCA = PCA_all[,1]
  
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
}


#adni
memoryscore.adni = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, ADNI_MEM)
  df$PCA = df$ADNI_MEM_scaled

  save(df,
       file = file.path(outdir, "df.memory.rda"))
}


#aibl
memoryscore.aibl = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  
  df = scale.mem.var(df, LogMem_short_delay)
  df = scale.mem.var(df, LogMem_delayed)
  
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(LogMem_short_delay_scaled,LogMem_delayed_scaled)
  memvar_all <- df %>% 
    ungroup() %>% 
    select(LogMem_short_delay_scaled,LogMem_delayed_scaled)
  
  first_tp = imputePCA(first_tp, ncp = 1)$completeObs
  memvar_all = imputePCA(memvar_all, ncp = 1)$completeObs
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all)
  
  df$PCA =  PCA_all[,1] 
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
  # library(factoextra)
  # var <- get_pca_var(PCA_tp1)
  # # Contributions of variables to PC1
  # fviz_contrib(PCA_tp1, choice = "var", axes = 1, top = 10)
  
}


#ous
memoryscore.ous = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  
  df = scale.mem.var(df, CERAD_immediate)
  df = scale.mem.var(df, CERAD_delayed)
  
  # no missing data in df
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(CERAD_immediate,CERAD_delayed)
  memvar_all <- df %>% 
    ungroup() %>% 
    select(CERAD_immediate,CERAD_delayed)
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all)
  
  df$PCA =  PCA_all[,1] 
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
}


#ukb
memoryscore.ukb = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  
  df = scale.mem.var(df, PAL)
  df$PCA = df$PAL_scaled
  
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
}

#habs
memoryscore.habs = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  
  df = scale.mem.var(df, SRT_delayed)
  df = scale.mem.var(df, SRT_total_recall)
  df = scale.mem.var(df, LogMem_short_delay)
  df = scale.mem.var(df, LogMem_delayed)
  df = scale.mem.var(df, FCsrt_Free)
  
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(SRT_delayed_scaled,SRT_total_recall_scaled,LogMem_short_delay_scaled,LogMem_delayed_scaled, FCsrt_Free_scaled)
  first_tp = imputePCA(first_tp, ncp = 2)$completeObs
  
  memvar_all <- df %>% 
    ungroup() %>% 
    select(SRT_delayed_scaled,SRT_total_recall_scaled,LogMem_short_delay_scaled,LogMem_delayed_scaled, FCsrt_Free_scaled)
  memvar_all = imputePCA(  memvar_all, ncp = 2)$completeObs
  
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  PCA_all <- predict(PCA_tp1, newdata = memvar_all)
  
  df$PCA <- PCA_all[,1] 
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
}


#preventAD
memoryscore.preventad = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, list_learning_total_score)
  df = scale.mem.var(df, list_recall_score)
  df = scale.mem.var(df, story_memory_score)
  df = scale.mem.var(df, story_recall_score)
  df = scale.mem.var(df, figure_recall_total_score)
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup()%>% 
    select(list_learning_total_score_scaled, list_recall_score_scaled, story_memory_score_scaled, story_recall_score_scaled, figure_recall_total_score_scaled)

  memvar_all <- df %>% 
    select(list_learning_total_score_scaled, list_recall_score_scaled, story_memory_score_scaled, story_recall_score_scaled, figure_recall_total_score_scaled)
  
  #memvar_all <- df %>% 
   # select(list_learning_total_score_scaled, list_recall_score_scaled, story_memory_score_scaled, story_recall_score_scaled)
  
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  #PCA_tp1 <- prcomp(first_tp[, c(17,18,19,20)], center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all[, -1] )
  
  df$PCA <- PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
}

#oasis3
memoryscore.oasis3 = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df,LOGIMEM)
  df = scale.mem.var(df,MEMUNITS)
  df = scale.mem.var(df,srtfree)
  df = scale.mem.var(df,srttotal)
  df = scale.mem.var(df,asscmem)
 
  
  # inputing missing data for first time point
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(LOGIMEM_scaled, MEMUNITS_scaled, srtfree_scaled, srttotal_scaled, asscmem_scaled) 
  first_tp = imputePCA(first_tp, ncp = 2)$completeObs
  PCA_tp1 <- prcomp(first_tp, center = T, scale = T)
  
  memvar_all <- df %>% 
    ungroup() %>% 
    select(LOGIMEM_scaled, MEMUNITS_scaled, srtfree_scaled, srttotal_scaled, asscmem_scaled)
  memvar_all = imputePCA(  memvar_all, ncp = 2)$completeObs
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all )
  df$PCA = PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
}

#bbhi
memoryscore.bbhi = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, RAVLT_short_delay)
  df = scale.mem.var(df, RAVLT_delayed)
  df = scale.mem.var(df, RAVLT_total_learning)
  
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(RAVLT_short_delay_scaled,RAVLT_delayed_scaled,RAVLT_total_learning_scaled)
  memvar_all <- df %>% 
    ungroup() %>% 
    select(RAVLT_short_delay_scaled,RAVLT_delayed_scaled,RAVLT_total_learning_scaled)
  first_tp = imputePCA(first_tp, ncp = 2)$completeObs
  memvar_all = imputePCA(memvar_all, ncp = 2)$completeObs
  
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all)
  
  df$PCA = PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
}

#uio
memoryscore.vetsa = function(site){
  load(file.path(datadir, "df.long_prepared.rda"))
  df <- df.memory
  df = scale.mem.var(df, CVLT_short_delay)
  df = scale.mem.var(df, CVLT_delayed)
  df = scale.mem.var(df, CVLT_total_learning)
  
  # no missing data
  first_tp <- df %>% 
    filter(timepoint == 1) %>% 
    ungroup() %>% 
    select(CVLT_short_delay_scaled, CVLT_delayed_scaled,CVLT_total_learning_scaled)
  
  memvar_all <- df %>% 
    select(CVLT_short_delay_scaled, CVLT_delayed_scaled,CVLT_total_learning_scaled)
  
  
  #PCA on first timepoint
  PCA_tp1 <- prcomp(first_tp, center = TRUE, scale = TRUE)
  
  
  PCA_all <- predict(PCA_tp1, newdata = memvar_all[, -1] )
  
  df$PCA = PCA_all[,1]
  
  save(df,
       file = file.path(outdir, "df.memory.rda"))
  
   # library(factoextra)
   # var <- get_pca_var(PCA_tp1)
   # # Contributions of variables to PC1
   # fviz_contrib(PCA_tp1, choice = "var", axes = 1, top = 10)
}

