#set links
set_links = function(folder) {
  datadir <<- here('data_memory_long', folder, 'other_output')
  mridir <<- here('data_normative_long/df_mri', folder)
  }




#rename function cog data
rename_var_cog = function(df, site) { 
  g1 = df %>% 
    select(df.harmonize_cog[[site]])
  colnms <- names(g1) 
  colnms <- colnms[ !colnms == "subses_id"]
  
  names(g1) = df.harmonize_cog$harmonize
  
  g2 <- df %>% 
    select(-colnms)
  
  df_out <- inner_join(g1, g2)
  
  sex.f = sex.rename$female[grep(site,sex.rename$dataset)]
  sex.m = sex.rename$male[grep(site,sex.rename$dataset)]
  sex.fn = sex.rename$female[grep("normative",sex.rename$dataset)] %>% as.numeric()
  sex.mn = sex.rename$male[grep("normative",sex.rename$dataset)] %>% as.numeric()
  df_out = 
    df_out %>% 
    mutate(sex= 
             if_else(sex == sex.f, sex.fn,
                     if_else(sex == sex.m,sex.mn,NaN)))
  return(df_out)
}

#rename function mri data
rename_var_mri = function(df, site) { 
  g1 = df %>% 
    select(df.harmonize_mri[[site]])
  colnms <- names(g1) 
  colnms <- colnms[ !colnms == "subses_id"]
  
  names(g1) = df.harmonize_mri$harmonize
  
  g2 <- df %>% 
    select(-colnms)
  
  df_out <- inner_join(g1, g2)
  
  sex.f = sex.rename$female[grep(site,sex.rename$dataset)]
  sex.m = sex.rename$male[grep(site,sex.rename$dataset)]
  sex.fn = sex.rename$female[grep("normative",sex.rename$dataset)] %>% as.numeric()
  sex.mn = sex.rename$male[grep("normative",sex.rename$dataset)] %>% as.numeric()
  df_out = 
    df_out %>% 
    mutate(sex= 
             if_else(sex == sex.f, sex.fn,
                     if_else(sex == sex.m,sex.mn,NaN)))
  return(df_out)
}


prepare.uio = function(site){
  load(file.path(datadir, "df.rda"))

  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "CVLT") %>% 
    rename(CVLT_total_learning = cvlt_a_total) %>% 
    rename(CVLT_short_delay = cvlt_5min_free ) %>% 
    rename(CVLT_delayed = cvlt_30min_free) %>% 
    rename(CVLT_recognition = cvlt_recog_hits) %>% 
    mutate(CVLT_learning = (cvlt_a5-cvlt_a1)) %>% 
    mutate(CVLT_true_recognition = (CVLT_recognition-cvlt_recog_fp))
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,CVLT_learning,CVLT_total_learning,CVLT_short_delay,CVLT_delayed,CVLT_recognition,CVLT_true_recognition) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup() %>% 
    mutate(rid = paste0("uio-", subject_id)) 
  
  
  
  save(df.memory,
    file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.umu = function(site){
  load(file.path(datadir, "df.rda"))

  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "Recall of senteces and words")
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,VT_free_recall_vtb) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup() %>% 
    mutate(rid = paste0("umu-", subject_id)) 
  
  #link with MRI data
  
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
    file = file.path(datadir, "df.long_prepared.rda"))
}


prepare.ub = function(site){
  load(file.path(datadir, "df.rda"))
 
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "RAVLT") %>% 
    rename(RAVLT_total_learning = Verbal_memory_RAVLT_Total) %>% 
    rename(RAVLT_delayed = Verbal_Memory_RAVLT_delayed) %>% 
    mutate(RAVLT_learning = (Verbal_memory_RAVLT_trial5-Verbal_memory_RAVLT_trial1))
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,RAVLT_total_learning, RAVLT_delayed,RAVLT_learning) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("ub-", subject_id)) 
  
  #link with MRI data
  
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.mpib = function(site){
  load(file.path(datadir, "df.rda"))
  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "VLMT & Scene encoding & face-proffesion & Object locatio") %>% 
    mutate(VLMT_learning = (trial5-trial1)) %>%
    rename(VLMT_short_delay = trial7) %>% 
    rename(VLMT_delayed = recall1) %>% 
    rename(VLMT_recognition_rates = vlmt_hit_r) %>% 
    rename(VLMT_true_recognition_rates = vlmt_hit_minus_fa) %>% 
    rowwise(subject_id) %>% 
    mutate(VLMT_total_learning = sum(c_across(trial1:trial5))) 
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,
           memory_task, 
           VLMT_total_learning, 
           VLMT_learning, 
           VLMT_short_delay, 
           VLMT_delayed, 
           VLMT_recognition_rates, 
           VLMT_true_recognition_rates, 
           io_hit_minus_fa,
           fp_hit_minus_new_fa,
           fp_hit_minus_rear_fa,
           ol_acc) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup() %>% 
    mutate(rid = paste0("mpib-", subject_id)) 
  
  #link with MRI data
  
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.adni = function(site){
  load(file.path(datadir, "df.rda"))
 
  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "ADNI_MEM") %>% 
    filter(!is.na(ADNI_MEM))

  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,ADNI_MEM) %>% 
    mutate(age = as.numeric(age)) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("adni-", subject_id)) 
  
  #link with MRI data
  
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.aibl = function(site){
  load(file.path(datadir, "df.rda"))

  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "LogMem") %>% 
    rename(LogMem_short_delay = LIMMTOTAL) %>% 
    rename(LogMem_delayed = LDELTOTAL)
  
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,LogMem_short_delay,LogMem_delayed) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("aibl-", subject_id)) 
  
  #link with MRI data
  
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.habs = function(site){
  load(file.path(datadir, "df.rda"))

  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "SRT & LogMem & FCsrt_Free") %>% 
    rename(LogMem_short_delay = LogicMem_IL) %>% 
    rename(LogMem_delayed = LogicMem_DR) %>% 
    rename(SRT_delayed = SRT_dr) %>% 
    rename(SRT_total_recall = SRT_tr )
  
  memory.long = memory.long[!duplicated(memory.long), ]  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,LogMem_short_delay, LogMem_delayed, SRT_delayed, SRT_total_recall,FCsrt_Free) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("habs-", subject_id)) 
  
  #link with MRI data
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  

  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}



prepare.ous = function(site){
  load(file.path(datadir, "df.rda"))

  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "CERAD") %>% 
    rename(CERAD_immediate = Cerad_immediate) %>% 
    rename(CERAD_delayed = Cerad_delayed) %>%
    rename(CERAD_PCA = PCA_CERAD)
  
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,CERAD_immediate,
           CERAD_delayed, CERAD_PCA) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("ous-", subject_id)) 
  
  #link with MRI data
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}


prepare.preventad = function(site){
  load(file.path(datadir, "df.rda"))
  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "wordlist & storyrecall")
  
  

  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,
           memory_task,
           list_learning_total_score, 
           list_recall_score, 
           story_memory_score, 
           story_recall_score, 
           figure_recall_total_score) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("preventad-", subject_id)) 
  
  #link with MRI data
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}

prepare.ukb = function(site){
  load(file.path(datadir, "df.rda"))
  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "word pair") %>% 
    rename(Pairs_matching = number_of_incorrect_matches_in_round_f399_round2) %>% 
    rename(PAL = number_of_word_pairs_correctly_associated_f20197) 
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,Pairs_matching,PAL) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup() %>% 
    mutate(rid = paste0("ukb-", subject_id)) 
  
  #link with MRI data
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
  
  # #for PAL memory task
  # df.memory_PAL =
  #   df.memory_PAL %>%
  #   mutate(site = site) %>%
  #   mutate(memory_task = "PAL") %>%
  #   rename(Pairs_matching = number_of_incorrect_matches_in_round_f399_round2) %>%
  #   rename(PAL = number_of_word_pairs_correctly_associated_f20197)
  # 
  # df.memory_PAL = rename_var_cog(df.memory_PAL,site)
  # df.memory_PAL <- df.memory_PAL %>%
  #   select(site:subses_id,memory_task,Pairs_matching,PAL) %>%
  #   group_by(subject_id) %>%
  #   arrange(age) %>%
  #   mutate(timepoint = row_number()) %>%
  #   add_count(subject_id, name = "total_tp") %>%
  #   ungroup()
  
}

prepare.oasis3 = function(site){
  load(file.path(datadir, "df.rda"))

  df.memory = df.memory[!duplicated(df.memory), ]  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "LogMem & SRT & ASSMEM")
    # rename(LogMem_short_delay = logmem) %>% 
    # rename(LogMem_delayed = lmdelay)
  
  df.memory = rename_var_cog(memory.long,site)
  df.memory <- df.memory %>% 
    select(site:subses_id,
           memory_task,
           LOGIMEM, 
           MEMUNITS, 
           srtfree, 
           srttotal, 
           asscmem) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("oasis3-", subject_id)) 
  
  #link with MRI data
  # df.memory <- df.memory %>% 
  #   filter(rid %in% df.mri$rid)
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}

prepare.bbhi = function(site) {
  load(file.path(datadir, "df.rda"))

  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(project_id = "NA") %>% 
    mutate(memory_task = "RAVLT") %>% 
    rename(RAVLT_short_delay = inm) %>% 
    rename(RAVLT_delayed = delayed) %>% 
    rowwise(subject_id) %>% 
    mutate(RAVLT_total_learning = sum(c_across(inm_recall_1_raw:inm_recall_5_raw)))
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,memory_task,RAVLT_short_delay,RAVLT_delayed,RAVLT_total_learning ) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup()%>% 
    mutate(rid = paste0("bbhi-", subject_id)) 
  
   save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}

prepare.vetsa = function(site){
  load(file.path(datadir, "df.rda"))
  
  memory.long = 
    df.memory %>% 
    mutate(site = site) %>% 
    mutate(memory_task = "CVLT") %>% 
    rename(CVLT_total_learning = cvatot) %>% 
    rename(CVLT_short_delay = CVSDFR ) %>% 
    rename(CVLT_delayed = CVLDFR)
  
  
  df.memory = rename_var_cog(memory.long,site)
  
  df.memory <- df.memory %>% 
    select(site:subses_id,Case,memory_task,CVLT_total_learning,CVLT_short_delay,CVLT_delayed) %>% 
    group_by(subject_id) %>%
    arrange(age) %>% 
    mutate(timepoint = row_number()) %>% 
    add_count(subject_id, name = "total_tp") %>% 
    ungroup() %>% 
    mutate(rid = paste0("vetsa-", subject_id)) 
  
  
  save(df.memory,
       file = file.path(datadir, "df.long_prepared.rda"))
}


