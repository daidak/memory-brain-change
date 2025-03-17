set_links = function(folder) {
  datadir <<- here('data_memory_long', folder)
  outdir <<- here('data_memory_long', folder)
  #if (!dir.exists(outdir)) {dir.create(outdir)}
}


extract_output.full = function(site) {
  load(file.path(datadir, "memory_output.rda"))
  load(file.path(datadir, "df.memory.rda"))
  df <- df %>% 
    mutate(agebl = first(age), #Age baseline
           time = (age-agebl), #time (years) between TP
           xtime = time - mean(time),# subtract mean time per subject time
           total_followup = last(time)) #Total follow-up time (years)
  
  
  # check PCA is in the right sign
  if (site == "uio") {
    if (cor(df$PCA, df$CVLT_total_learning, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "umu") {
    if (cor(df$PCA, df$VT_free_recall_vtb, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "ub") {
    if ( cor(df$PCA, df$RAVLT_total_learning_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  }  else if (site == "mpib") {
    if ( cor(df$PCA, df$VLMT_total_learning_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  }  else if (site == "adni") {
    if ( cor(df$PCA, df$ADNI_MEM, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  }  else if (site == "aibl") {
    if (cor(df$PCA, df$LogMem_delayed_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "ous") {
    if (cor(df$PCA, df$CERAD_delayed, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "ukb") {
    if ( cor(df$PCA, df$PAL, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "habs") {
    if ( cor(df$PCA, df$LogMem_short_delay_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "preventad") {
    if ( cor(df$PCA, df$story_recall_score_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "oasis3") {
    if ( cor(df$PCA, df$LOGIMEM_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "bbhi") {
    if ( cor(df$PCA, df$RAVLT_delayed_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  } else if (site == "vetsa") {
    if ( cor(df$PCA, df$CVLT_total_learning_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  }  else if (site == "ucam") {
    if ( cor(df$PCA, df$storyrecall_immediate_scaled, use =  "pairwise.complete.obs") < 0) {
      df$PCA = -df$PCA
    }
  }
  
  df <- df %>% 
    select(subject_id,
           subses_id,
           site, 
           rid, 
           age,
           sex, 
           memory_task, 
           timepoint, 
           total_tp,
           total_followup,
           agebl, 
           time, 
           PCA) %>% 
    rename(total_tp_memory = total_tp) %>% 
    rename(total_followup_memory = total_followup) %>% 
    rename(subject = subject_id) %>% 
    mutate(subject = as.character(subject),
           site = as.character(site))
  
  out <- out %>% 
    mutate(subject = as.character(subject)) %>% 
    rename("memory_slope" = "slope") %>% 
    rename("memory_endtercept" = "endtercept")
  
  dt <- left_join(df, out)
  
  if (cor(dt$PCA, dt$Int.ranef, use = "pairwise.complete.obs") < 0) {
    dt$Int.ranef = -dt$Int.ranef
    dt$memory_slope = -dt$memory_slope
    dt$intercept = -dt$intercept
    dt$memory_endtercept = -dt$memory_endtercept
  }
  return(dt)
}
