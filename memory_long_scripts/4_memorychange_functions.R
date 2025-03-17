set_links = function(folder) {
  datadir <<- here('data_memory_long', folder)
  outdir <<- here('data_memory_long', folder)
  #if (!dir.exists(outdir)) {dir.create(outdir)}
}

compute_var = function(df) {
  df <- df %>% 
    group_by(subject_id) %>% 
    arrange(timepoint) %>% 
    mutate(agebl = first(age), #Age baseline
           time = (age-agebl), #time (years) between TP
           xtime = time - mean(time),# subtract mean time per subject time
           total_followup = last(time)) #Total follow-up time (years)
  
  df <- df %>% 
    mutate(subject = factor(subject_id)) %>% 
    mutate(sex = factor(sex)) %>% 
    mutate(retest.dummy = if_else(timepoint > 1,1,0), #01111 
           retest.dummy2 = if_else(timepoint > 2,1,0), #00111 
           )
  return(df)
}


#gamms
run_gamms = function(df,type){
  
  #for gamm, only random intercept, if type = 1, add only one retest.dummy(for datasets with mostly 2 TPs)
  #and type = 2 add another retest.dummy (for datasets with more timepoints)
  # retest.dummy : 01111 
  # retest.dummy2: 00111
  
  
  
  if (type == 1){
    model <- gamm4(memory ~ s(age) + sex + retest.dummy, random = ~(1|subject), data = df)
  } else if (type == 2){
    model <- gamm4(memory ~ s(age) + sex + retest.dummy + retest.dummy2, random = ~(1|subject), data = df)
  } else {
    stop("invalid input")
  }
  
  save(model,
       file = file.path(outdir, "gamms.rda"))


  
}

#extract residuals 
extract_residuals = function(model,df) {
  res <- model$gam$residuals
  subject_id <- df$subject_id
  resid <- data.frame(subject_id,res) %>% 
    select(res)
  df <- cbind(df,resid)
   
  # mod_res1 =
  #   df  %>%
  #   group_by(subject_id) %>%
  #   mutate(#agebsl = min(age),
  #     #time = age - agebsl,
  #     xage = mean(age),
  #     xtime = age - xage) %>%
  #   select(subject_id, xtime, res) %>%
  #   filter(!is.na(res)) %>%
  #   group_by(subject_id) %>%
  #   nest() %>%
  #   mutate(mod = map(data, ~ lm(res ~ xtime, data = .x)),
  #          tidy = map(mod, ~broom::tidy(.x))) %>%
  #   unnest(tidy) %>%
  #   select(subject_id, term, estimate) %>%
  #   pivot_wider(names_from = "term",
  #               values_from = c("estimate"))
  #  names(mod_res) = c("subject", "intercept", "slope")
   
  mod_res =
    df  %>%
    group_by(subject_id) %>%
    mutate(#agebsl = min(age),
      Mage = max(age),
      xtime = age - Mage) %>%
    select(subject_id, xtime, res) %>%
    filter(!is.na(res)) %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(res ~ xtime, data = .x)),
           tidy = map(mod, ~broom::tidy(.x))) %>%
    unnest(tidy) %>%
    select(subject_id, term, estimate) %>%
    pivot_wider(names_from = "term",
                values_from = c("estimate"))
  names(mod_res) = c("subject", "endtercept", "slope")
  return(mod_res)
  
}

extract_residuals_v2 = function(model,df) {
  res <- model$gam$residuals
  subject_id <- df$subject_id
  resid <- data.frame(subject_id,res) %>% 
    select(res)
  df <- cbind(df,resid)
  
  ss = ranef(model$mer)
  df.ranef = data.frame(subject = rownames(ss$subject), 
                        Int.ranef =  ss$subject$`(Intercept)`)
  mod_res1 =
    df  %>%
    group_by(subject_id) %>%
    mutate(agebsl = min(age),
      time = age - agebsl) %>% 
      #xage = mean(age),
      #xtime = age - xage) %>%
    select(subject_id, time, res) %>%
    filter(!is.na(res)) %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(res ~ time, data = .x)),
           tidy = map(mod, ~broom::tidy(.x))) %>%
    unnest(tidy) %>%
    select(subject_id, term, estimate) %>%
    pivot_wider(names_from = "term",
                values_from = c("estimate"))  
   names(mod_res1) = c("subject", "intercept", "slope")

  mod_res =
    df  %>%
    group_by(subject_id) %>%
    mutate(#agebsl = min(age),
      Mage = max(age),
      xtime = age - Mage) %>%
    select(subject_id, xtime, res) %>%
    filter(!is.na(res)) %>%
    group_by(subject_id) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(res ~ xtime, data = .x)),
           tidy = map(mod, ~broom::tidy(.x))) %>%
    unnest(tidy) %>%
    select(subject_id, term, estimate) %>%
    pivot_wider(names_from = "term",
                values_from = c("estimate"))
  names(mod_res) = c("subject", "endtercept", "slope")
  
  mod_res$intercept = mod_res1$intercept
  
  if (is.numeric(mod_res$subject)) {
  mod_res = inner_join(df.ranef %>% mutate(subject = as.numeric(subject)), mod_res)
  } else {
    mod_res = inner_join(df.ranef %>% mutate(subject = as.character(subject)), mod_res)
  }
  return(mod_res)
  
}



#memory change, uio
memory.uio = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,1))
  cat("computed gamm, saved as gamms.rda")
  
  load(file.path(datadir, "gamms.rda"))
  
  out = extract_residuals_v2(model,df) 

  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
  
    }
  

#memory change, umu
memory.umu = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA)
  try(run_gamms(df,2))
  
  cat("computed gamm, saved as gamms.rda")
  
  load(file.path(datadir, "gamms.rda"))

  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, ub
memory.ub = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,1))
  cat("computed gamm, saved as gamms.rda")
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}


#memory change, mpib
memory.mpib = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA)
  
  try(run_gamms(df,1))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  
  cat("saving memory_output.rda ")
}


#memory change, adni
memory.adni = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, aibl
memory.aibl = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, ous
memory.ous = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA)
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  
  cat("saving memory_output.rda ")
}

#memory change, ukb
memory.ukb = function(site){
  
  set.seed(123)
  #load memory data (just  longitudinalsince its a lot of data)
  load(file.path(datadir, "df.memory.rda"))
 
  
  df = compute_var(df)
  
  df <- df %>%
    rename(memory = PCA)
  
  try(run_gamms(df,1))
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))

}

#memory change, habs
memory.habs = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA)
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, preventAd
memory.preventad = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, oasis3
memory.oasis3 = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  try(run_gamms(df,2))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  cat("saving memory_output.rda ")
}

#memory change, bbhi
memory.bbhi = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  
  try(run_gamms(df,1))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  
  cat("saving memory_output.rda ")
}

#memory change, vetsa
memory.vetsa = function(site){
  load(file.path(datadir, "df.memory.rda"))
  
  df = compute_var(df)
  df <- df %>%
    rename(memory = PCA) 
  
  df <- df %>% 
    mutate(Case = factor(Case))
  
  #remove sex since its only males
  
  model <- gamm4(memory ~ s(age) + retest.dummy + retest.dummy2, random = ~(1|subject) + (1|Case), data = df)
  
  save(model,
     file = file.path(outdir, "gamms.rda"))
  cat("computed gamm, saved as gamms.rda")
  
  
  load(file.path(datadir, "gamms.rda"))
  out = extract_residuals_v2(model,df) %>%      filter(!is.na(slope))
  
  save(out,
       file = file.path(outdir, "memory_output.rda"))
  
  cat("saving memory_output.rda ")
}
