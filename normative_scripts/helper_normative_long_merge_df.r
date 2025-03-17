set_links = function(folder) {
  sitedir = here('data_normative_long','df_mri', folder)
  load(file.path(sitedir, "df.rda"))
  return(df.out)
}

rename_variables = function(df.out, site) { 
  grot = df.out %>% 
    ungroup() %>% 
    select(df.harmonize[[site]])
  names(grot) = df.harmonize$normative_modelling
  sex.f = sex.rename$female[which(sex.rename$dataset == site)]
  sex.m = sex.rename$male[which(sex.rename$dataset == site)]
  sex.fn = sex.rename$female[grep("normative",sex.rename$dataset)] %>% as.numeric()
  sex.mn = sex.rename$male[grep("normative",sex.rename$dataset)] %>% as.numeric()
  grot = 
    grot %>% 
    mutate(sex= 
             if_else(sex == sex.f, sex.fn,
                     if_else(sex == sex.m,sex.mn,NaN)))
  return(grot)
}


squeue = function(user, job) {
  df.squeue = system(paste0("squeue --name=", job," -u ",user), intern = T) %>% 
    strsplit(., " +") %>% 
    simplify2array() %>% 
    t() %>% 
    as.data.frame()
  return(df.squeue)
  
}

prepare.uio.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_")) 
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("uio data renamed. ready to be merged")
}





prepare.mpib.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("mpib data renamed. ready to be merged")
}

prepare.umu.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("umu data renamed. ready to be merged")
}


prepare.ub.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("ub data renamed. ready to be merged")
}



prepare.bbhi.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("bbhi data renamed. ready to be merged")
}

prepare.habs.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("habs data renamed. ready to be merged")
}

prepare.ous.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("ous data renamed. ready to be merged")
}

prepare.aibl.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, model, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("aibl data renamed. ready to be merged")
}


prepare.adni.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    select(-site) %>% 
    mutate(site = paste(site, paste(Site, model,sep ="-"), sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           str_pad(df.out$site %>% 
               as.factor() %>% 
               as.numeric() %>% 
               as.character(),
             3, 
             pad = "0")) %>% 
    as.numeric()
  
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("adni data renamed. ready to be merged")
}

prepare.preventad.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "1", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("preventad data renamed. ready to be merged")
}





prepare.oasis3.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("oasis3 data renamed. ready to be merged")
}

prepare.vetsa.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_")) 
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("vetsa data renamed. ready to be merged")
}

prepare.ukb.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(input = paste(eid, tp, sep = "-"),
           site = paste(site, uk_biobank_assessment_centre_f54, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum) %>% 
    filter(!is.na(EstimatedTotalIntraCranialVol))
  return(df.out)
  
  cat("ukb data renamed. ready to be merged")
}

stats.overview = function(df.merge.long, df.merge) {
  mod = list()
  df.merge.long = 
    df.merge.long %>% 
    group_by(rid) %>% 
    mutate(futime = max(age)- min(age),
           mage = mean(age))  
  
  df.merge = 
    df.merge %>% 
    group_by(rid) %>% 
    mutate(mage = mean(age))
  
  df.long.dmg = 
    df.merge.long %>% 
    group_by(rid) %>% 
    summarise(sex = first(sex), 
              mage = first(mage), 
              futime = first(futime),
              n = n_distinct(sub_id), 
              dataset = first(dataset))
  
  df.all.dmg = 
    df.merge %>% 
    group_by(rid) %>% 
    summarise(sex = first(sex), 
              mage = first(mage), 
              dataset = first(dataset))
  
  # n information
  mod$n$long.data = df.merge.long$sub_id %>% length()
  
  mod$n$all.data = df.merge$sub_id %>% length()
  
  mod$n$long.data.dataset = 
    df.merge.long %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$all.data.dataset = 
    df.merge %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$long.data.site = 
    df.merge.long %>% 
    group_by(site) %>% 
    tally()
  
  mod$n$all.data.site = 
    df.merge %>% 
    group_by(site) %>% 
    tally()
  
  mod$n$long.unique.subs = df.long.dmg$rid %>% length()
  mod$n$all.unique.subs = df.all.dmg$rid %>% length()
  
  mod$n$long.unique.subs.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$all.unique.subs.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    tally()
  
  # age information  
  mod$age$all.mean.subs.age.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$long.mean.subs.age.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$all.mean.subs.age = 
    df.all.dmg %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$long.mean.subs.age = 
    df.long.dmg %>%
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$n$sites = 
    df.merge %>% 
    summarise(sites = n_distinct(site))
  
  mod$n$sites.dataset = 
    df.merge %>% 
    group_by(dataset) %>% 
    summarise(sites = n_distinct(sitenum))
  
  # SEX
  mod$sex$all.sex.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$long.sex.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$all.sex = 
    df.all.dmg %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$long.sex = 
    df.long.dmg %>%
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  # follow up time and observations
  mod$follow.up$long.follow.up.span.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(futime), sd = sd(futime))
  
  mod$follow.up$long.follow.up.span = 
    df.long.dmg %>%
    summarise(mean = mean(futime), sd = sd(futime))
  
  mod$follow.up$long.observations.per.participant.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(n), sd = sd(n))
  
  mod$follow.up$long.observations.per.participant = 
    df.long.dmg %>%
    summarise(mean = mean(n), sd = sd(n))
  
  return(mod)
} 

linechart = function(infile, outdir, filename) {
  x = infile %>% 
    select(dataset, rid, age, sub_id ) %>% 
    group_by(rid) %>% 
    mutate(minAge = min(age),
           maxAge = max(age)) %>% 
    ungroup() %>% 
    group_by(dataset) %>% 
    mutate(nrid = n_distinct(rid)) %>% 
    ungroup() %>% 
    arrange(minAge, maxAge)
  
  xx = data.frame(rid = unique(x$rid) )
  xx = left_join(xx,
                 x %>% select(rid, dataset) %>% group_by(rid,dataset) %>% tally())
  xx = xx %>% group_by(dataset) %>% mutate(i = row_number())
  x = left_join(x,xx)
  
  ggplot(x, aes(x = age,
                y = i, 
                group = as.character(i), 
                color = dataset)) + 
    geom_line(size = .1) + 
    geom_point(size =.1) + 
    #scale_color_viridis_c() +
    theme_classic() + 
    theme(legend.position = 'none') +
    facet_wrap(.~dataset, scales = "free_y")
  ggsave( file = file.path(outdir, filename))
} 
  
filter_sites_with_low_n_v2 = function(df.merge.long,df.merge, thr.obs, thr.subs) {
  
  df.merge.subs = df.merge %>% 
    group_by(site) %>% 
    mutate(nsubs = n_distinct(sub_id)) %>% 
    tally() %>% 
    filter(n >= thr.obs)
  
  df.merge.sites = df.merge %>% 
    group_by(site) %>% 
    summarise(nsite = n_distinct(rid)) %>% 
    filter(nsite >= thr.subs)
  
  sites <<- inner_join(df.merge.subs, df.merge.sites)
  
  df.merge <<- 
    df.merge %>% filter(site %in% sites$site)
  
  df.merge.long <<- 
    df.merge.long %>% filter(site %in% sites$site)
  
}

filter_sites_with_low_n = function(df.merge.long,df.merge, thr.obs) {
  i = 0
  c = 0
  # step filters scans with low numbers of subs, and the removes subs with n == 1 obserbation
  while (i == 0) {
    c = c +1
    print(c)
    grot.long = 
      df.merge.long %>% 
      group_by(site) %>% 
      mutate(nscans = n_distinct(sub_id)) %>% 
      filter(nscans >= thr.obs) %>% 
      ungroup() %>% 
      group_by(rid) %>% 
      mutate(nsub = n_distinct(age)) %>% 
      filter(nsub > 1) %>% 
      ungroup()
    
    # removes scanners where n given tp == 1  < 25
    grot.all = 
      df.merge %>% 
      group_by(rid) %>% 
      mutate(minage = min(age)) %>% 
      filter(age == minage) %>% 
      ungroup() %>% 
      group_by(site) %>% 
      mutate(nscans = n_distinct(sub_id)) %>% 
      ungroup() %>% 
      filter(nscans >= thr.obs) 
    
    # filter df based on both condition
    sites = intersect(grot.all$site, grot.long$site)
    
    grot.long = 
      grot.long %>% 
      filter(site %in% sites)
    
    grot.all2 = 
      df.merge %>% 
      filter(site %in% sites)
    
    if  (is_empty(setdiff(unique(df.merge.long$site),sites))) {
      df.merge <<- grot.all2
      df.merge.long <<- grot.long
      i = 1
    } else {
      df.merge = grot.all2
      df.merge.long = grot.long
    }
  }
}


save_sample_from_model = function(df.long, modeldir, train.datasets) {
  
  try(dir.create(modeldir))
  
  df.norm =
    df.long %>% 
    mutate(site = if_else(site == "ukb_11025", "ukb-11025.0",
                          if_else(site == "ukb_11027", "ukb-11027.0", 
                                  if_else(site == "ucam_01", "cam", site)))) %>% 
    filter(dataset %in% train.datasets)
  
  write.csv(df.norm, 
            file = file.path(modeldir, "df_te.csv"),
            quote = F)
  print("data overlapping in normative model saved")
}

save_sample_first_tp_in_calibration = function(df.long, df.all, modeldir, train.datasets) {
  try(dir.create(modeldir))
  
  df.norm =
    df.long %>% 
    filter(!dataset %in% train.datasets)
  
  write.csv(df.norm, 
            file = file.path(modeldir, "df_te.csv"),
            quote = F)
  
  df.calibration = 
    df.all %>% 
    group_by(rid) %>% 
    mutate(minage = min(age)) %>% 
    filter(age == minage) %>% 
    filter(!dataset%in% train.datasets) %>% 
    ungroup() %>% 
    mutate(sub_id = paste0("ca_", sub_id))
  
  write.csv(df.calibration, 
            file = file.path(modeldir,  "df_ca.csv"),
            quote = F)
  print("data overlapping not with normative model saved")
}


save_sample_min_tp_in_calibration = function(df.long, df.all, modeldir, train.datasets,n.calibration){
  try(dir.create(modeldir)) 
  
  df.norm =
    df.long %>% 
    filter(!dataset %in% train.datasets)
  
  df.sites = data.frame(site = df.norm$site %>% unique())
  df.sites$criteria = df.sites$n.calibration = df.sites$n.tp = 0
  
  grot = 
    df.all %>% 
    group_by(rid) %>% 
    mutate(n = n_distinct(age),
           minage = min(age)) %>% 
    filter(!dataset%in% train.datasets, 
           age == minage) %>%
    ungroup() 
  
  # save n timepoints x site 
  for (i in 1:max(grot$n)) {
    grot2 = grot %>% 
      filter(n == i) %>% 
      group_by(site) %>% 
      summarise( "nsub.{i}" := n_distinct(age))
    
    df.sites = left_join(df.sites, grot2)
    
    if (i == 1) {
      grot3 = df.sites$nsub.1
    } else {
      grot3 = rowSums(df.sites[,grepl("nsub", names(df.sites))], na.rm = T)  
    }
    idx = df.sites$criteria == 0 & grot3 >= n.calibration & !is.na(grot3)
    df.sites[idx,"criteria"] = 1 
    df.sites[idx,"n.calibration"] = grot3[idx] 
    df.sites[idx,"n.tp"] = i 
  }
  
  ## selection of participants for calibration sample
  subs.calibration = c()
  for (i in 1:length(df.sites$site)) {
    #print(df.sites$site[i])
    
    if ( df.sites$n.tp[i] == 1) {
      subs.calibration =   
        c(subs.calibration, 
          grot %>% filter(site == df.sites$site[i], n == 1) %>% .$rid)
    } else {
      grot4 = 
        unique(
          grot %>% 
            filter(site == df.sites$site[i] & 
                   n < df.sites$n.tp[i]) %>% 
            .$rid)
      grot5 = 
        unique(
          grot %>% 
            filter(site == df.sites$site[i] & 
                   n == df.sites$n.tp[i]) %>% 
            .$rid)
      subs.calibration =   
        c(subs.calibration,
          grot4, 
          sample(grot5, n.calibration - length(grot4)))
    }
  }
  
  df.calibration =
    grot  %>% 
    filter(rid %in% subs.calibration)
  
  df.norm =
    df.long %>% 
    filter(!dataset %in% train.datasets & 
           !rid %in% subs.calibration) 
  
  write.csv(df.norm, 
            file = file.path(modeldir, "df_te.csv"),
            quote = F)
  
  write.csv(df.calibration, 
            file = file.path(modeldir, "df_ca.csv"),
            quote = F)
  
  write.csv(df.sites, 
            file = file.path(modeldir, "sites_x_timepoints.csv"),
            quote = F)
  
}

save_sample_random_in_calibration = function(df.long, df.all, modeldir, train.datasets,n.calibration) {
  try(dir.create(modeldir))
  grot = 
    df.all %>% 
    group_by(rid) %>% 
    mutate(minage = min(age)) %>% 
    filter(!dataset%in% train.datasets, 
           minage == age) %>% 
    ungroup()
  
  df.calibration =  
    grot %>% 
    group_by(site) %>% 
    slice_sample(n = n.calibration)
  
  subs.calibration = df.calibration$rid  
  
  df.norm =
    df.long %>% 
    filter(!dataset %in% train.datasets, 
           !rid %in% subs.calibration) 
  
  write.csv(df.norm, 
            file = file.path(modeldir, "df_te.csv"),
            quote = F)
  
  write.csv(df.calibration, 
            file = file.path(modeldir, "df_ca.csv"),
            quote = F)  
}


save_sample_random_in_calibration_v2 = function(df.long, df.merge, modeldir, n.calibration, sites.site, ds) {
  try(dir.create(modeldir))
  
  df.grot =  
    df.merge %>%  
    filter(dataset == ds) %>% 
    group_by(site, rid) %>% 
    tally() %>% 
    ungroup() %>% 
    nest_by(site)
  df.grot = left_join(df.grot, sites.site) %>% 
    mutate(ncalibr = if_else(nsite >= n.calibration, n.calibration, nsite -2)) %>% 
    mutate(ca = list(slice_sample(data, n = ncalibr))) %>% 
    select(site, ca) %>%
    unnest(ca) %>% 
    select(-n)
  
  df.ca = left_join(df.grot,
                    df.merge)  
  
  df.long = df.long %>%  
    filter(dataset == ds)
  
  df.te =
    anti_join(df.long, df.grot)
  
  write.csv(df.te, 
            file = file.path(modeldir, "df_te.csv"),
            quote = F)
  
  write.csv(df.ca, 
            file = file.path(modeldir, "df_ca.csv"),
            quote = F)  
}




merge_normative_mri_data = function(analysis, phenotypes, df.harmonize, outdir, model_name) {
  thr.miss=.2
  grot = list()
  er = c()
  for (i in 1:length(phenotypes)) {
    idpdir = file.path(outdir, analysis, model_name, phenotypes[i])
    if (file.exists(file.path(idpdir, "data.rda"))) {
      load(file.path(idpdir, "data.rda")) 
      if (i == 1){
        df.base = df.pheno$df$df.all.rid.filt  %>% 
          select(-c(meanZ, 
                    deltaZ,
                    se.meanZ, 
                    se.deltaZ))
      }
      grot$idx.delta[[i]] = df.pheno$df$idx.rid.deltaZ
      grot$idx.mean[[i]] = df.pheno$df$idx.rid.meanZ
      grot$deltaZ[[i]] = df.pheno$df$df.all.rid.filt$deltaZ
      grot$meanZ[[i]] = df.pheno$df$df.all.rid.filt$meanZ
      grot$se.deltaZ[[i]] = df.pheno$df$df.all.rid.filt$se.deltaZ
      grot$se.meanZ[[i]] = df.pheno$df$df.all.rid.filt$se.meanZ
    } else {
      er = c(er,i)
    }
  }
  
  idx.delta = as.data.frame(do.call(cbind, grot$idx.delta))
  idx.mean = as.data.frame(do.call(cbind, grot$idx.mean))
  df.deltaZ = as.data.frame(do.call(cbind, grot$deltaZ))
  df.meanZ = as.data.frame(do.call(cbind, grot$meanZ))
  df.se.deltaZ = as.data.frame(do.call(cbind, grot$se.deltaZ))
  df.se.meanZ = as.data.frame(do.call(cbind, grot$se.meanZ))
  
  
  idx.delta[is.na(idx.delta)] = F
  idx.mean[is.na(idx.mean)] = F
  # create output dataframe
  df = list()
  df$outliers$delta$idx = idx.delta
  df$outliers$mean$idx = idx.mean
  df$outliers$delta$roi = colSums(idx.delta == F)
  df$outliers$delta$rid = rowSums(idx.delta == F)
  df$outliers$mean$roi = colSums(idx.mean == F)
  df$outliers$mean$rid = rowSums(idx.mean == F)
  df$outliers$delta$rmsubs.miss.data = df$outliers$delta$rid < 187*thr.miss
  df$outliers$mean$rmsubs.miss.data = df$outliers$mean$rid < 187*thr.miss
  df$outliers$thr.miss = thr.miss
  
  df$df$delta$df = df.deltaZ
  df$df$mean$df = df.meanZ
  df$df$delta$df.se = df.se.deltaZ
  df$df$mean$df.se = df.se.meanZ
  df$df$base$df = df.base
  df.ggeg = rename_variables_ggseg(phenotypes, df.harmonize)
  
  df$misc$phenotypes = df.ggeg$phenotypes
  df$misc$phenotypes.atlas = df.ggeg$atlas
  
  return(df)
}





filter_out_outliers = function(dataset, phenos, thr.time, thr.mean, thr.delta, model_name) {
  
  datasetdir = here("data_normative_long/df_mri/", dataset)
  idpdir = file.path(datasetdir,"normative", model_name, phenos)
  
  sub_vars = c("sub_id", 
               "rid", 
               "dataset", 
               "age", 
               "sex", 
               "site", 
               "sitenum")
  
  df =  import(file.path(file.path(datasetdir, 
                                   "normative", 
                                   "df_te.csv"))) %>% 
    select(sub_vars) %>% 
    group_by(rid) %>% 
    mutate(xage = mean(age),
           agebsl = min(age),
           time = max(age)- min(age))
  
  df.rid= 
    df %>% 
    group_by(rid) %>% 
    summarise(n = n_distinct(age),
              sex = first(sex),
              dataset = first(dataset),
              site = first(site),
              sites = n_distinct(site),
              time = first(time),
              xage = first(age))
  
  
  
  df.p = lapply(idpdir, 
                function(x) {data.frame(t(import(file.path(x, "Z_predict.txt"))))}) %>% 
    data.table::rbindlist() %>% 
    t() %>% 
    data.frame()
  names(df.p) = phenos
  
  df.norm = 
    cbind(df, df.p) %>% 
    filter(time > thr.time)
  
  df.filt = 
    df %>% 
    filter(time > thr.time)
  
  df.rid.filt = 
    df.rid %>% 
    filter(time > thr.time)
  
  grot = list()
  for(j in 1:length(idpdir)) {
    load(file.path(idpdir[j],"Z_predict_extended.Rda"))
    grot[[j]] = df.out
  }
  names(grot) = basename(idpdir)
  
  grot = data.table::rbindlist(grot, idcol = "feature") %>% 
    group_by(feature, rid) %>% 
    summarise(meanZ = mean(meanZ), 
              deltaZ = mean(deltaZ)) %>% 
    left_join(df.rid.filt,.)
  
  
  x.mean = find_outliers(grot, "meanZ", thr.mean, phenos)
  x.delta = find_outliers(grot, "deltaZ", thr.delta, phenos)
  
  
  
  # saving data
  df = list()
  df$outliers$delta$idx = x.delta[[2]]
  df$outliers$mean$idx = x.mean[[2]]
  df$outliers$delta$roi = colSums(x.delta[[2]] == F, na.rm = T)
  df$outliers$delta$rid = rowSums(x.delta[[2]] == F, na.rm = T)
  df$outliers$mean$roi = colSums(x.mean[[2]] == F, na.rm = T)
  df$outliers$mean$rid = rowSums(x.mean[[2]] == F, na.rm = T)
  df$outliers$delta$rmsubs.miss.data = df$outliers$delta$rid < dim(x.delta[[1]])[[2]]*thr.miss
  df$outliers$mean$rmsubs.miss.data = df$outliers$mean$rid < dim(x.mean[[1]])[[2]]*thr.miss
  df$outliers$thr.miss = thr.miss
  
  df$df$delta$df = x.delta[[1]]
  df$df$mean$df = x.mean[[1]]
  #df$df$delta$df.se = df.se.deltaZ
  #df$df$mean$df.se = df.se.meanZ
  df$df$base$df = df.rid.filt
  df.ggeg = rename_variables_ggseg(phenos, df.harmonize)
  
  df$misc$phenotypes = df.ggeg$phenotypes
  df$misc$phenotypes.atlas = df.ggeg$atlas
  
    save(df,
         file = file.path(datasetdir,
                          "normative",
                          paste("data",
                                model_name,
                                "rda", 
                                sep = ".")))
  
  return(df)
}



filter_out_outliers_gamm = function(dataset, phenos, thr.time, thr.mean, thr.delta, model_name) {
  
    datasetdir = here("data_normative_long/df_mri/", dataset)  
    idpdir = file.path(datasetdir, "GammHarm", phenos)  
  
  sub_vars = c("sub_id", 
               "rid", 
               "dataset", 
               "age", 
               "sex", 
               "site", 
               "sitenum")
  
  load(file.path(file.path(datasetdir, 
                                   "df.all.filt.Rda" ))) 
  
  
  df = df.merge.long  %>% 
    select(sub_vars) %>% 
    group_by(rid) %>% 
    mutate(xage = mean(age),
           agebsl = min(age),
           time = max(age)- min(age))
  
  df.rid= 
    df %>% 
    group_by(rid) %>% 
    summarise(n = n_distinct(age),
              sex = first(sex),
              dataset = first(dataset),
              site = first(site),
              sites = n_distinct(site),
              time = first(time),
              xage = first(age)) %>%  
    filter(time > thr.time)
  
  df.rid.filt = df.rid
  
  grot = list()
  for(j in 1:length(idpdir)) {
    load(file.path(idpdir[j],"Z_predict_extended.Rda"))
    grot[[j]] = df.out
  }
  names(grot) = basename(idpdir)
  
  grot = data.table::rbindlist(grot, idcol = "feature") %>% 
    group_by(feature, rid) %>% 
    summarise(meanZ = mean(meanZ), 
              deltaZ = mean(deltaZ)) 
  
  grot = 
    grot %>% filter(rid %in% df.rid.filt$rid)
  
  
  x.mean = find_outliers(grot, "meanZ", thr.mean, phenos)
  x.delta = find_outliers(grot, "deltaZ", thr.delta, phenos)
  
  
  
  # saving data
  df = list()
  df$outliers$delta$idx = x.delta[[2]]
  df$outliers$mean$idx = x.mean[[2]]
  df$outliers$delta$roi = colSums(x.delta[[2]] == F, na.rm = T)
  df$outliers$delta$rid = rowSums(x.delta[[2]] == F, na.rm = T)
  df$outliers$mean$roi = colSums(x.mean[[2]] == F, na.rm = T)
  df$outliers$mean$rid = rowSums(x.mean[[2]] == F, na.rm = T)
  df$outliers$delta$rmsubs.miss.data = df$outliers$delta$rid < dim(x.delta[[1]])[[2]]*thr.miss
  df$outliers$mean$rmsubs.miss.data = df$outliers$mean$rid < dim(x.mean[[1]])[[2]]*thr.miss
  df$outliers$thr.miss = thr.miss
  
  df$df$delta$df = x.delta[[1]]
  df$df$mean$df = x.mean[[1]]
  #df$df$delta$df.se = df.se.deltaZ
  #df$df$mean$df.se = df.se.meanZ
  df$df$base$df = df.rid.filt
  df.ggeg = rename_variables_ggseg(phenos, df.harmonize)
  
  df$misc$phenotypes = df.ggeg$phenotypes
  df$misc$phenotypes.atlas = df.ggeg$atlas
  
    save(df,
         file = file.path(datasetdir,
                          "GammHarm",
                          paste("data",
                                model_name,
                                "rda", 
                                sep = ".")))
  return(df)
}


find_outliers = function(grot, var, thr.mads, phenos) {
  x = grot %>% 
    select(rid, feature, var) %>% 
    pivot_wider(rid, 
                names_from = "feature", 
                values_from = var) 
  x = x[,match(phenos,names(x))]  
  
  mm = apply(x, 2,median, na.rm = T)
  mads = apply(x, 2,mad, na.rm = T)
  xx = sweep(x, 2, mm, "-")
  xx = sweep(xx, 2, mads, "/")
  idx = abs(xx) < thr.mads
  out = list(x, 
             idx)
  return(out)
}


rename_variables_ggseg = function(phenotypes, df.harmonize) { 
  library(ggsegDesterieux)
  library(ggseg)
  df.p = data.frame(phenotypes = phenotypes)
  gg.dest = unique(desterieux$data$label)
  gg.aseg = unique(aseg$data$label) 
  
  for (i in 1:length(phenotypes)) {
    library(ggseg)
    library(ggsegDesterieux)
    idx = grep(paste0("\\b",phenotypes[i],"\\b"),df.harmonize$normative_modelling)
    df.p$phenotypes[i] = df.harmonize$ggseg[idx]
    df.p$atlas[i] = ifelse(df.p$phenotypes[i] %in% gg.dest, "dest", 
                           ifelse(df.p$phenotypes[i] %in% gg.aseg, "aseg", "global"))
  }
  
  return(df.p)
}

select_vars_change_matrix = function(x, pheno) {
  names(x) = pheno
  rm.rois = 
    c("rh_MeanThickness_thickness",
      "lh_MeanThickness_thickness",
      "EstimatedTotalIntraCranialVol",
      "SubCortGrayVol", 
      "TotalGrayVol",
      "SupraTentorialVol",
      "SupraTentorialVolNotVent",
      "CSF")
  change.sign = c(names(x)[grepl("entricl",names(x))],
                  names(x)[grepl("Lat-Vent",names(x))])
  
  x = 
    x %>% select(-rm.rois)
  x[,change.sign] = -x[,change.sign]
  return(x)
}

merge_normative_mri_data_v2 = function(datasets, prefixin, prefixout) {
  db = list()
  for (ds in 1:length(datasets)) {
    datasetdir = here("data_normative_long/df_mri/", datasets[ds], "normative")
    filename = file.path(datasetdir, 
                         paste("data",
                               prefixin,
                               "rda", 
                               sep = "."))
    load(filename)
    db[[ds]] = df
  }
  
  df = list()
  
  # outlier data merging
  x = lapply(db, function(x) {x$outliers$delta$idx})
  df$outliers$delta$idx = reduce(x, rbind) %>% data.frame()
  
  
  x = lapply(db, function(x) {x$outliers$delta$roi})
  df$outliers$delta$roi = 
    reduce(x, rbind) %>% 
    data.frame() %>% 
    colSums()
  
  
  x = lapply(db, function(x) {x$outliers$delta$rid})
  df$outliers$delta$rid = 
    x %>% 
    simplify_all() %>% 
    simplify()
  
  x = lapply(db, function(x) {x$outliers$delta$rmsubs.miss.data})
  df$outliers$delta$rmsubs.miss.data = 
    x %>% 
    simplify_all() %>% 
    simplify()
  
  
  x = lapply(db, function(x) {x$outliers$mean$idx})
  df$outliers$mean$idx = reduce(x, rbind) %>% data.frame()
  
  
  x = lapply(db, function(x) {x$outliers$mean$roi})
  df$outliers$mean$roi = 
    reduce(x, rbind) %>% 
    data.frame() %>% 
    colSums()
  
  
  x = lapply(db, function(x) {x$outliers$mean$rid})
  df$outliers$mean$rid = 
    x %>% 
    simplify_all() %>% 
    simplify()
  
  x = lapply(db, function(x) {x$outliers$mean$rmsubs.miss.data})
  df$outliers$mean$rmsubs.miss.data = 
    x %>% 
    simplify_all() %>% 
    simplify()
  
  df$outliers$thr.miss = db[[1]]$outliers$thr.miss
  
  
  # misc
  df$misc$phenotypes = db[[1]]$misc$phenotypes 
  df$misc$phenotypes.atlas = db[[1]]$misc$phenotypes.atlas
  
  # df
  x = lapply(db, function(x) {x$df$base$df})
  df$df$base$df = data.table::rbindlist(x)
  
  x = lapply(db, function(x) {x$df$delta$df})
  df$df$delta$df = data.table::rbindlist(x)
  
  
  x = lapply(db, function(x) {x$df$mean$df})
  df$df$mean$df = data.table::rbindlist(x)
  
  
  m = db[[1]]$df$imputed.delta$m
  for (i in 1:m) {
    x = lapply(db, function(x) {mice::complete(x$df$imputed.mean, i)})
    df$df$imputed.mean[[i]] = x %>% data.table::rbindlist()
    x = lapply(db, function(x) {mice::complete(x$df$imputed.delta, i)})
    df$df$imputed.delta[[i]] = x %>% data.table::rbindlist()
  }
  
  save(df, file = file.path(outdir, 
                            paste("df.merged",
                                  prefixout,
                                  "Rda", 
                                  sep = ".")))
  return(df)
}
