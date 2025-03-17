args = commandArgs(TRUE)

type=as.character(args[1])
wd=as.character(args[2])


run_cluster_mega = function(wd, type) {
  # all = F run all cluster in a single gam
    
  library(here)  
  library(tidyverse)  
  library(mgcv)
  nrep = 100
  #load cluster data
  load(file.path(wd,"mega","consensus_cluster", "cluster_mod.rda"))
  annon <- m3c$mod$realdataresults[[8]]$ordered_annotation
  old.names = gsub("\\.","-",rownames(annon))
  annon$features = old.names  
  
  # load raw data
  load(file.path(wd, "merged_reliability.rda"))
  
  # load reliability data
  options(bitmapType = "cairo")
  rel_parameters = "/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/all/reliability/reliability_parameters.rda"
  load(rel_parameters)
  
  # load cluster data
  load(file.path(wd,"mega","PCA", "pca_result.rda"))
  
  # prepare matrices (select covariates, cluster data and pca data)
  
  df1 = 
    left_join(annon, df.merged) %>% 
    unnest() 
  
  # get time invariant covariates
  df3 = df1 %>% 
    group_by(rid) %>%  
    summarise(
      dataset = first(dataset), 
      sex = first(sex), 
      total_tp_memory = first(total_tp_memory), 
      total_followup_memory = first(total_followup_memory), 
      n=first(n),
      time  = first(time), 
      xage = first(xage), 
      ageblmri = first(ageblmri), 
      agediff= first(agediff), 
      ageF = first(ageF), 
      ageFmri = first(ageFmri), 
      overlap  = first(overlap), 
      apoe4 =first(apoe4), 
      apoe_status = first(apoe_status),
      icv  = first(icv))
  
  # get delta brain delta memory for cluster
  df2 = df1 %>% 
    group_by(consensuscluster, rid) %>% 
    summarise(memory_slope = mean(memory_slope, na.rm = T), 
              delta_brain = mean(delta_brain, na.rm = T)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = consensuscluster, 
                values_from = delta_brain, 
                names_prefix = "cl.")
  
  
  # get pca
  pca = -pca_result$x[,1]
  df.pca = data.frame(rid = names(pca), pca = pca)
  
  # merge provisional dataframes
  df.cc = inner_join(df2, df.pca)
  df.cc2 = left_join(df.cc,df3) 
  
  #C = df.cc[,-c(1:2)] %>% cor()

  # merge reliability parameters
  df.features.cc = left_join(annon, df.features)
  df.features.cc = 
    rbind(
      df.features.cc %>% 
        group_by(consensuscluster) %>% 
        summarise_if(is.numeric, mean), 
      df.features.cc = 
        df.features.cc %>% 
        summarise_if(is.numeric, mean) %>% 
        mutate(consensuscluster = "pca")) %>% 
    mutate(bs = seD^2, 
           lf = "grot")
  
  
  x = df.cc2 %>% mutate(lf = "grot")          
  x = left_join(x, df.features.cc,relationship = "many-to-many")
  x = x %>% 
    mutate(ws = pct.err.mean^2*((time^2*n*(n+1))/ (12*(n-1)))^-1,
           icc.sub = bs /(ws + bs),
           weight.sub2 = icc.sub^2) %>% 
    mutate(weight.sub2 = if_else(weight.sub2 < .09, .09, weight.sub2)) %>% 
    dplyr::select(-c(pct.err.mean, seD, bs, ws, icc.sub)) %>% 
    pivot_wider(names_from = "consensuscluster", 
                values_from = "weight.sub2",
                names_prefix = "w2.")
  
  
   x = 
     x %>% 
     mutate(pca = scale(pca)[,1])
   
   x = x %>% 
     filter(
       between(cl.1, -4.5, 4.5) &
         between(cl.2, -4.5, 4.5) &
         between(cl.3, -4.5, 4.5) &
         between(cl.4, -4.5, 4.5) &
         between(cl.5, -4.5, 4.5) &
         between(cl.6, -4.5, 4.5) &
         between(cl.7, -4.5, 4.5) &
         between(cl.8, -4.5, 4.5) &
         between(pca, -4.5, 4.5))
  ## run main models. 
  
  if (type == "all") {
    mod = 
      gamm(memory_slope ~ s(cl.1, bs = "cr") +  s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.pca)
    
    g.true = mod
    g.null1 = 
      gamm(memory_slope ~  s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.pca)
     g.null2 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null3 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null4 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null5 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null6 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null7 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.8, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
     g.null8 =
       gamm(memory_slope ~  s(cl.1, bs = "cr") + s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr"),
            data = x, random=list(dataset=~1), weights = w2.pca)
    
    nulls = ls(pattern ="g.null")
    grot.d = list()
    grot = list()
    for (i in 1:length(nulls)) {
      gnull  = get(nulls[i])
      grot = replicate(nrep,  tryCatch({gamm_bootstrap_all(gnull,x)}, error = function(e){ NaN}), simplify = F)
      grot2d = 
        lapply(grot, function(x) {if(any(is.nan(x))) return(NULL)
        as.data.frame(x) %>% rownames_to_column(var = "term")})
      grot2d = data.table::rbindlist(grot2d)
      grot2d$null = nulls[i]
      grot.d[[i]] =  grot2d
    }
    
    bootstraps = data.table::rbindlist(grot.d)
    save(mod, 
         bootstraps,
         file = file.path(wd,"mega","consensus_cluster", "all_consensus_cluster_gam.Rda"))
  }
  
  if (type == "pca") {
    mod.1.vs.pca = 
      gamm(memory_slope ~ s(cl.1, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.1)
    
    mod.2.vs.pca = 
      gamm(memory_slope ~ s(cl.2, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.2)
    
    mod.3.vs.pca = 
      gamm(memory_slope ~ s(cl.3, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.3)
    
    mod.4.vs.pca = 
      gamm(memory_slope ~ s(cl.4, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.4)
    
    mod.5.vs.pca = 
      gamm(memory_slope ~ s(cl.5, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.5)

    mod.6.vs.pca = 
      gamm(memory_slope ~ s(cl.6, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.6)
    mod.7.vs.pca = 
      gamm(memory_slope ~ s(cl.7, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.7)
    mod.8.vs.pca = 
      gamm(memory_slope ~ s(cl.8, bs = "cr") +  s(pca, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.8)
    
    g.null1.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.1)
    
    g.null2.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.2)
    
    g.null3.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.3)
    
    g.null4.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.4)
    
    g.null5.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.5)
    
    g.null6.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.6)
    
    g.null7.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.7)
    
    g.null8.vs.pca = 
      gamm(memory_slope ~  s(pca, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.8)
    
    
    grot = ls(pattern = "vs.pca")
    trues = grot[!grepl("null", grot)]
    nulls = grot[grepl("null", grot)]
    cl = gsub("mod.", "cl", trues) %>% gsub(".vs.pca", "", .)
   
    
    
    grot.d = list()
    for (i in 1:length(nulls)) {
      gnull  = get(nulls[i])
      cc = cl[i]
      grot = replicate(nrep,  tryCatch({gamm_bootstrap_pca(gnull,x,cc)}, error = function(e){ NaN}))
      grot = grot %>% as.data.frame()
      grot$cl = cc
      grot.d[[i]] =  grot
    }
    
    bootstraps = data.table::rbindlist(grot.d)
    save(mod.1.vs.pca, 
         mod.2.vs.pca, 
         mod.3.vs.pca, 
         mod.4.vs.pca, 
         mod.5.vs.pca, 
         mod.6.vs.pca, 
         mod.7.vs.pca, 
         mod.8.vs.pca, 
         bootstraps,
         file = file.path(wd,"mega","consensus_cluster", "pca_consensus_cluster_gam.Rda"))
  }
  
  
  if (type == "cl1") {
    
    mod.2.vs.cl1 = 
      gamm(memory_slope ~ s(cl.2, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.2)
    
    mod.3.vs.cl1 = 
      gamm(memory_slope ~ s(cl.3, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.3)
    
    mod.4.vs.cl1 = 
      gamm(memory_slope ~ s(cl.4, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.4)
    
    mod.5.vs.cl1 = 
      gamm(memory_slope ~ s(cl.5, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.5)
    
    mod.6.vs.cl1 = 
      gamm(memory_slope ~ s(cl.6, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.6)
    mod.7.vs.cl1 = 
      gamm(memory_slope ~ s(cl.7, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.7)
    mod.8.vs.cl1 = 
      gamm(memory_slope ~ s(cl.8, bs = "cr") +  s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.8)
    
    
    g.null2.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.2)
    
    g.null3.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.3)
    
    g.null4.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.4)
    
    g.null5.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.5)
    
    g.null6.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.6)
    
    g.null7.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.7)
    
    g.null8.vs.cl1 = 
      gamm(memory_slope ~  s(cl.1, bs = "cr"), data = x, random=list(dataset=~1), weights = w2.8)
    
    
    grot = ls(pattern = "vs.cl1")
    trues = grot[!grepl("null", grot)]
    nulls = grot[grepl("null", grot)]
    cl = gsub("mod.", "cl", trues) %>% gsub(".vs.cl1", "", .)
    
    grot.d = list()
    for (i in 1:length(nulls)) {
      gnull  = get(nulls[i])
      cc = cl[i]
      grot = replicate(nrep,  tryCatch({gamm_bootstrap_cl(gnull,x,cc)}, error = function(e){ NaN}))
      grot = grot %>% as.data.frame()
      grot$cl = cc
      grot.d[[i]] =  grot
    }
    
    bootstraps = data.table::rbindlist(grot.d)
    save(mod.2.vs.cl1, 
         mod.3.vs.cl1, 
         mod.4.vs.cl1, 
         mod.5.vs.cl1, 
         mod.6.vs.cl1, 
         mod.7.vs.cl1, 
         mod.8.vs.cl1, 
         bootstraps,
         file = file.path(wd,"mega","consensus_cluster", "cl1_consensus_cluster_gam.Rda"))
  }
  
  if (type == "basic") {
    mod.1.vs.b = 
      gamm(memory_slope ~ s(cl.1, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.1)
    
    mod.2.vs.b = 
      gamm(memory_slope ~ s(cl.2, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.2)
    
    mod.3.vs.b = 
      gamm(memory_slope ~ s(cl.3, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.3)
    
    mod.4.vs.b = 
      gamm(memory_slope ~ s(cl.4, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.4)
    
    mod.5.vs.b = 
      gamm(memory_slope ~ s(cl.5, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.5)
    
    mod.6.vs.b = 
      gamm(memory_slope ~ s(cl.6, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.6)
    
    mod.7.vs.b = 
      gamm(memory_slope ~ s(cl.7, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.7)
    mod.8.vs.b = 
      gamm(memory_slope ~ s(cl.8, bs = "cr"),
           data = x, random=list(dataset=~1), weights = w2.8)
    
    g.null1.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.1)
    
    g.null2.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.2)
    
    g.null3.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.3)
    
    g.null4.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.4)
    
    g.null5.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.5)
    
    g.null6.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.6)
    
    g.null7.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.7)
    
    g.null8.vs.b = 
      gamm(memory_slope ~  1,
           data = x, random=list(dataset=~1), weights = w2.8)
    
    
    
    grot = ls(pattern = "vs.b")
    trues = grot[!grepl("null", grot)]
    nulls = grot[grepl("null", grot)]
    cl = gsub("mod.", "cl", trues) %>% gsub(".vs.b", "", .)
   
    
    
    grot.d = list()
    for (i in 1:length(nulls)) {
      gnull  = get(nulls[i])
      cc = cl[i]
      grot = replicate(nrep,  tryCatch({gamm_bootstrap_basic(gnull,x,cc)}, error = function(e){ NaN}))
      grot = grot %>% as.data.frame()
      grot$cl = cc
      grot.d[[i]] =  grot
    }
    
    bootstraps = data.table::rbindlist(grot.d)
    save(mod.1.vs.b, 
         mod.2.vs.b, 
         mod.3.vs.b, 
         mod.4.vs.b, 
         mod.5.vs.b, 
         mod.6.vs.b, 
         mod.7.vs.b, 
         mod.8.vs.b, 
         bootstraps,
         file = file.path(wd,"mega","consensus_cluster", "basic_consensus_cluster_gam.Rda"))
  }
}
  
 
  
  gamm_bootstrap_all = function(g.null, df) { 
    y.null = predict(g.null$lme)
    n = length(y.null)
    wild = sample(c(-1 ,1), size=n, replace=TRUE)
    y.resample = y.null + g.null$gam$residuals*wild           
    df$y_resample = y.resample  %>% as.numeric() 
    
    bs = 
      gamm(y_resample ~ s(cl.1, bs = "cr") +  s(cl.2, bs = "cr") + s(cl.3, bs = "cr") + s(cl.4, bs = "cr") + s(cl.5, bs = "cr") + s(cl.6, bs = "cr") + s(cl.7, bs = "cr") + s(cl.8, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.pca)
    stable = summary(bs$gam)$s.table
    return(stable)
  }
  
  
  gamm_bootstrap_pca = function(g.null, df,cc) { 
    y.null = predict(g.null$lme)
    n = length(y.null)
    wild = sample(c(-1 ,1), size=n, replace=TRUE)
    y.resample = y.null + g.null$gam$residuals*wild           
    df$y_resample = y.resample  %>% as.numeric() 
    
    if (cc == "cl1") {  
      bs = 
      gamm(y_resample ~ s(cl.1, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.1)
    } else if (cc == "cl2") {
    bs = 
      gamm(y_resample ~ s(cl.2, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.2)
    } else if (cc == "cl3") {
    bs = 
      gamm(y_resample ~ s(cl.3, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.3)
    } else if (cc == "cl4") {
    bs = 
      gamm(y_resample ~ s(cl.4, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.4)
    } else if (cc == "cl5") {
    bs = 
      gamm(y_resample ~ s(cl.5, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.5)
    } else if (cc == "cl6") {
    bs = 
      gamm(y_resample ~ s(cl.6, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.6)
    } else if (cc == "cl7") {
    bs = 
      gamm(y_resample ~ s(cl.7, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.7)
    } else if (cc == "cl8") {
    bs = 
      gamm(y_resample ~ s(cl.8, bs = "cr") +  s(pca, bs = "cr"),
           data = df, random=list(dataset=~1), weights = w2.8)
    }
    stable = summary(bs$gam)$s.table[1,4]
    return(stable)
  }
  
  
  gamm_bootstrap_cl = function(g.null, df,cc) { 
    y.null = predict(g.null$lme)
    n = length(y.null)
    wild = sample(c(-1 ,1), size=n, replace=TRUE)
    y.resample = y.null + g.null$gam$residuals*wild           
    df$y_resample = y.resample  %>% as.numeric() 
    
    if (cc == "cl2") {
      bs = 
        gamm(y_resample ~ s(cl.2, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.2)
    } else if (cc == "cl3") {
      bs = 
        gamm(y_resample ~ s(cl.3, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.3)
    } else if (cc == "cl4") {
      bs = 
        gamm(y_resample ~ s(cl.4, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.4)
    } else if (cc == "cl5") {
      bs = 
        gamm(y_resample ~ s(cl.5, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.5)
    } else if (cc == "cl6") {
      bs = 
        gamm(y_resample ~ s(cl.6, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.6)
    } else if (cc == "cl7") {
      bs = 
        gamm(y_resample ~ s(cl.7, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.7)
    } else if (cc == "cl8") {
      bs = 
        gamm(y_resample ~ s(cl.8, bs = "cr") +  s(cl.1, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.8)
    }
    stable = summary(bs$gam)$s.table[1,4]
    return(stable)
  }
  
  
  gamm_bootstrap_basic = function(g.null, df,cc) { 
    y.null = predict(g.null$lme)
    n = length(y.null)
    wild = sample(c(-1 ,1), size=n, replace=TRUE)
    y.resample = y.null + g.null$gam$residuals*wild           
    df$y_resample = y.resample  %>% as.numeric() 
    
    if (cc == "cl1") {
      bs = 
        gamm(y_resample ~ s(cl.1, bs = "cr") ,
             data = df, random=list(dataset=~1), weights = w2.1)
    } else if (cc == "cl2") {
      bs = 
        gamm(y_resample ~ s(cl.2, bs = "cr") ,
             data = df, random=list(dataset=~1), weights = w2.2)
    } else if (cc == "cl3") {
      bs = 
        gamm(y_resample ~ s(cl.3, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.3)
    } else if (cc == "cl4") {
      bs = 
        gamm(y_resample ~ s(cl.4, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.4)
    } else if (cc == "cl5") {
      bs = 
        gamm(y_resample ~ s(cl.5, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.5)
    } else if (cc == "cl6") {
      bs = 
        gamm(y_resample ~ s(cl.6, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.6)
    } else if (cc == "cl7") {
      bs = 
        gamm(y_resample ~ s(cl.7, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.7)
    } else if (cc == "cl8") {
      bs = 
        gamm(y_resample ~ s(cl.8, bs = "cr"),
             data = df, random=list(dataset=~1), weights = w2.8)
    }
    stable = summary(bs$gam)$s.table[1,4]
    return(stable)
  }
  
  
  
  
run_cluster_mega(wd, type)
  