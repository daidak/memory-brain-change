args = commandArgs(TRUE)

mmodels=as.character(args[1])
model_name=as.character(args[2])
phenotypes=as.character(args[3])


compute_Zchange = function(mmodels, model_name, phenotypes) {
  library(tidyverse)
  library(rio)
  library(broom)
  
  df_te = file.path(mmodels, "df_te.csv")
  df = import(df_te)
  idpdir = file.path(mmodels, model_name, phenotypes)
  
  df.Z = import(file.path(idpdir, "Z_predict.txt"))
  
  if(dim(df.Z)[2] == 1) { # training sets
    names(df.Z) = c("Znorm")
    df = cbind(df, df.Z)  
  }else {
    names(df.Z) = c("sub_id", "Znorm")
    df = left_join(df, df.Z)
  }
  mod = 
    df  %>% 
    group_by(rid) %>% 
    mutate(#agebsl = min(age), 
      #time = age - agebsl,
      xage = mean(age), 
      xtime = age - xage) %>% 
    select(rid, xtime, Znorm) %>%
    filter(!is.na(Znorm)) %>% 
    group_by(rid) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ lm(Znorm ~ xtime, data = .x)), 
           tidy = map(mod, ~broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    select(rid, term, estimate, std.error) %>% 
    pivot_wider(names_from = "term", 
                values_from = c("estimate", "std.error"))
  names(mod) = c("rid", "meanZ", "deltaZ", "se.meanZ", "se.deltaZ")
  df.out = 
    left_join(
      df %>% select(sub_id, rid, Znorm),
      mod)
  save("df.out", 
       file = file.path(idpdir, "Z_predict_extended.Rda"))
}

compute_Zchange(mmodels, model_name, phenotypes)




