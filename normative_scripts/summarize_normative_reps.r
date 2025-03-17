args = commandArgs(TRUE)

ds=as.character(args[1])
model_name=as.character(args[2])
nreps=as.numeric(args[3])
phenotypes=as.character(args[4])
wd=as.character(args[5])

summarize_normative_reps = function(ds,model_name, nreps, phenotypes, wd) {
  
  library(tidyverse)
  library(rio)

  phenotypes = 
    import(file.path(
      wd, 
      "scripts/braincharts-master",
      "docs", 
      phenotypes), 
      header = F) %>% 
    .$V1
  
  outdir = file.path(
      wd, 
      'data_normative_long/df_mri/',
      ds, 
      "normative",
      model_name)
  
  try(dir.create(outdir))
  
  te = 
    list.files(file.path(wd, 
                         'data_normative_long/df_mri/',
                          ds, 
                          "normative",  
                          str_pad(as.character(1:nreps),4,pad = "0")), 
                pattern = "df_te.csv", 
                full.names = T)
  te = lapply(te, function(x) {import(x, drop = "V1")})
  
  df.te = reduce(te, rbind) %>% distinct()
  
  write.csv(df.te, 
            file = file.path(wd, 
                             'data_normative_long/df_mri',
                             ds, 
                             "normative",
                             "df_te.csv"),
            quote = F)
  
    for (i in 1:length(phenotypes)) {
      print(paste(ds, phenotypes[i], i , length(phenotypes)))
      
      
      try(dir.create(file.path(outdir, phenotypes[i])))
      
      Z = 
        list.files(file.path(wd,
                             'data_normative_long/df_mri/',
                              ds, 
                              "normative",  
                              str_pad(as.character(1:nreps),4,pad = "0"), 
                              model_name,
                              phenotypes[i]), 
                   pattern = "Z_predict.txt", 
                   full.names = T)
      
      
      Z = lapply(Z, function(x) {import(x)})
      
      grot = lapply(1:nreps, function(x) { data.frame(sub_id = te[[x]]$sub_id, V1 = Z[[x]]$V1)})
      
      
      df.out = data.table::rbindlist(grot) %>% 
        group_by(sub_id) %>% 
        summarise(M = mean(V1, na.rm = T))
      df.out = 
        left_join(df.te %>% select(sub_id),
                df.out)
      
      write.table(df.out, 
                  file = file.path(outdir, phenotypes[i], "Z_predict.txt"),
                  quote = F, 
                  row.names = F,
                  col.names = F)
      
    }
  }
  

summarize_normative_reps(ds,model_name, nreps, phenotypes, wd)