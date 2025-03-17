mutate_ub = function(cohort, is.long = F) {
  if (is.long == T) {
    cohort =
      cohort %>% 
      mutate(Sess = if_else(`Round_ id` =="R01",1, 
                            if_else(`Round_ id` == "R02",2,3))) 
  } else {
    cohort = 
      cohort %>% 
      mutate(Sess = 1)   
  }
  return(cohort)
}


df.ub.change.names = function(df.mri) {
  df.mri$Subject_id[ df.mri$Subject_id == "hc001"] = "hc001_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc002"] = "hc002_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc003"] = "hc003_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc005"] = "hc005_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc010"] = "hc010_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc011"] = "hc011_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc014"] = "hc014_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc015"] = "hc015_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc016"] = "hc016_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc017"] = "hc017_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc018"] = "hc018_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc019"] = "hc019_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc020"] = "hc020_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc021"] = "hc021_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc023"] = "hc023_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc024"] = "hc024_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc027"] = "hc027_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc028"] = "hc028_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc029"] = "hc029_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc031"] = "hc031_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc032"] = "hc032_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc033"] = "hc033_nc059"
  df.mri$Subject_id[ df.mri$Subject_id == "hc035"] = "hc035_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc036"] = "hc036_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc039"] = "hc039_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc040"] = "hc040_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc041"] = "hc041_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc042"] = "hc042_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc043"] = "hc043_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc044"] = "hc044_nc060"
  df.mri$Subject_id[ df.mri$Subject_id == "nc052"] = "nc052_2"  
  df.mri$Subject_id[ df.mri$Subject_id == "001_MGS"] = "001_MSG"
  df.mri$Subject_id[ df.mri$Subject_id == "013_RFS"] = "013_RFJ"
  return(df.mri)
}

load_scan_info_adni = function(fname, field.strength) {
  
  mri.info = read.csv(file.path(sitedir, fname ))
  names(mri.info)[10] = "Image_ID"
  
  tmp =mri.info %>% 
    filter(Type =="Original") %>% 
    group_by(Subject.ID, 
             Visit,
             Imaging.Protocol) %>% 
    tally()
  
  mri.info =
    left_join(mri.info %>% 
                select(-Imaging.Protocol),
              tmp) %>% 
    mutate(field.strength = field.strength)
  return(mri.info)  
}

              "Mfg Model=Allegra"                        
#[6] "Mfg Model=Ingenuity"         "Mfg Model=Ingenia Elition X" "Mfg Model=MAGNETOM Vida"  
df.adni.set.scans = function() {
  triotim = c("Mfg Model=TrioTim",
              "Mfg Model=Trio")
  prisma = c("Mfg Model=Prisma_fit",
             "Mfg Model=Prisma")
  verio = c("Mfg Model=Verio")
  skyra = c("Mfg Model=Skyra",
            "Mfg Model=Skyra|DicomCleaner",
            "Mfg Model=Skyra_fit")
  achieva = c("Mfg Model=Achieva",
              "Mfg Model=Achieva dStream")
  ingenia = c("Mfg Model=Ingenia")
  signa = c("Mfg Model=Signa HDxt",
            "Mfg Model=SIGNA HDx",
            "Mfg Model=GENESIS_SIGNA",
            "Mfg Model=SIGNA EXCITE",
            "Mfg Model=SIGNA HDx",
            "Mfg Model=Signa HDxt",
            "Mfg Model=SIGNA Premier")
  discovery = c("Mfg Model=DISCOVERY MR750",
                "Mfg Model=DISCOVERY MR750w")
  intera = c("Mfg Model=Intera",
             "Mfg Model=Gyroscan Intera",
             "Mfg Model=Gyroscan NT",
             "Mfg Model=Intera",
             "Mfg Model=Intera Achieva")
  gemini = c("Mfg Model=GEMINI",
             " Mfg Model=Ingenuity")
  biographmMR = c("Mfg Model=Biograph_mMR")
  avanto=c("Mfg Model=Avanto")
  sonata=c("Mfg Model=Sonata",
           "Mfg Model=SonataVision",
           "Mfg Model=Espree")
  symphony=c("Mfg Model=Symphony",
             "Mfg Model=SymphonyTim",
             "Mfg Model=NUMARIS/4")
  allegra = c("Mfg Model=Allegra")
  models = list("triotim" = triotim,
                "prisma" = prisma,
                "verio" = verio,
                "skyra" = skyra, 
                "achieva" = achieva, 
                "ingenia" = ingenia, 
                "signa" = signa, 
                "discovery" = discovery, 
                "intera" = intera, 
                "gemini" = gemini, 
                "biographmMR" = biographmMR,
                "avanto" = avanto, 
                "sonata" = sonata, 
                "symphony" = symphony,
                "allegra" = allegra)
  
  return(models)
}
set_links = function(folder) {
  tabulateddata <<- here('data-raw/tabulated')
  sitedir <<- file.path(tabulateddata, folder)
  outdir <<- here('data_normative_long','df_mri', folder)
  if (!dir.exists(outdir)) {dir.create(outdir)}
}

check.uio.fs = function(sitedir, outdir, mrifiles) {
  df.uio = import(file.path(sitedir, "noas_query_2022-11-29_21-26-23_041b29e_b0f5e7b.csv"))
  df.uio = df.uio %>% filter(subject_shareable == 1)
  suppressWarnings({
    df.uio = 
      df.uio %>% 
      group_by(subject_id) %>% 
      mutate(mean_age = mean(visit_age)) %>% 
      ungroup()
    
    ss = df.uio %>% 
      mutate(mmse.crit = if_else(mms_score < 26 & !is.na(mms_score), 1,0)) %>% 
      group_by(subject_id) %>% 
      mutate(low.mmse = max(mmse.crit))
    
    df.uio = ss %>% 
      filter(low.mmse == 0) %>% 
      select(-c(low.mmse, mmse.crit))
    
    aa = ss %>% filter(low.mmse == 1) %>% 
      group_by(subject_id) %>% 
      mutate(time.norm =if_else(mmse.crit == 0, visit_age, NaN), 
             time.norm = max(time.norm, na.rm = T)) %>%
      ungroup() %>% 
      mutate(crit.max.normal = if_else(visit_age <= time.norm & !is.infinite(time.norm), 1,0))  %>% 
      filter(crit.max.normal == 1) %>% 
      select(-c(low.mmse, mmse.crit, time.norm, crit.max.normal))
    
    
    df.uio = rbind(df.uio,aa) %>% 
      filter(mean_age > 20)
  })
  
  grot.scanner = import(file.path(sitedir, "noas_query_2022-11-29_21-28-27_94ea888_b0f5e7b.csv"))
  grot.scanner.1 = grot.scanner %>% group_by(mri_info_folder) %>% summarise_if(is.numeric, mean)
  grot.scanner.2 = grot.scanner %>% group_by(mri_info_folder) %>% summarise_if(is.character, first)
  grot.scanner.3 = grot.scanner %>% group_by(mri_info_folder) %>% summarise(mri_info_date = first(mri_info_date))
  df.scanner = 
    left_join(grot.scanner.1, grot.scanner.2) %>% 
    left_join(., grot.scanner.3)
  
  df.scanner = 
    df.scanner %>% filter(subject_id %in% df.uio$subject_id)
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(subject_id = substr(input, 5,11) %>% 
             as.numeric(),
           visit_number = substr(input, 16,16) %>% 
             as.numeric(),
           scanner = substr(input,17,30)) %>% 
    separate(scanner, into = "scanner", sep = ".lo")
  
  df.out = inner_join(df.uio, df.mri)
  
  # create mri-to-database file
  save(df.out, 
       df.scanner, 
       df.uio, 
       file = file.path(outdir, "df.rda"))
  
  
  cat("uio data checked. rda file save with valid data only")
}


check.mpib.fs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  list.files(file.path(sitedir,"behavioural_data"))
  
  df.mpib = import(file.path(sitedir, "sheet_1.tsv")) 
  df.mpib = df.mpib %>% 
    filter(!is.na(do_mrt)) %>% 
    mutate(round_id = if_else(round_id == "R01", 1, 2))
  
  df.MPIB = read.csv(file.path(sitedir,"Data_BASEII_YlvaEniko211011.csv"))
  
  df.id.mpib = df.MPIB %>%
    mutate(ses01 = if_else(!is.na(MRT12_age), 1,NaN),
           ses02 = if_else(!is.na(MRT16_age), 2,NaN)) %>% 
    #select(ID, ses01, ses02) %>% 
    pivot_longer(cols = c(ses01,ses02), 
                 names_to = "grot1",
                 values_to = "SubjectRound_id",
                 values_drop_na = T) %>% 
    mutate(visit_age = if_else(SubjectRound_id ==1, MRT12_age, MRT16_age)) %>% 
    select(-c(Age, 
              MPIB_2012_Date1, 
              MRT12_age, 
              MRT16_age))  %>% 
    select(ID, SubjectRound_id, MMSE)
  names(df.id.mpib) = c("subject_id", "round_id", "MMSE")
  
  df.mpib = left_join(df.mpib, df.id.mpib) %>% 
    filter(MMSE > 25 | is.na(MMSE)) # two-timepoints so if mmse is low there wil be not long data
  
  
  # check if fs data
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  
  df.mri =
    df.mri %>% 
    mutate(subject_id = substr(input, 5,11),
           round_id = substr(input, 15,16) %>% as.numeric())
  
  # missing subjects
  df.out = left_join(df.mpib, df.mri)
  
  save(df.out, 
       df.mpib, 
       file = file.path(outdir, "df.rda"))
  
  
  cat("mpib data checked. rda file save with valid data only")
  
}
check.umu.fs = function(sitedir, outdir, mrifiles) {
  # load databases and get mri_filepaths
  
  grot.t5 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T5.xlsx"))
  grot.t6 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T6.xlsx"))
  grot.t7 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T7.xlsx"))
  
  db.exclusion = 
    import(file.path(sitedir, "Exclusions_Betula_Imaging_Sample_cleaned.xlsx"))
  
  db.exclusion = 
    c(paste0(db.exclusion  %>% 
               filter(`Exclude T5` == 1) %>% 
               .$`Subject code`, 
             "_T5"),
      paste0(db.exclusion  %>% 
               filter(`Exclude T6` == 1) %>% 
               .$`Subject code`, 
             "6"),
      paste0(db.exclusion  %>% 
               filter(`Exclude T7` == 1) %>% 
               .$`Subject code`, 
             "7"))
  
  
  df.Umea = data.table::rbindlist(list(grot.t5, 
                                  grot.t6, 
                                  grot.t7), 
                             use.names = F) %>% 
    rowwise() %>% 
    mutate(age = mean(c(rounded_age_HT, rounded_age_MT), na.rm = T),
           input = paste(Lifebrain_Subject_id, Lifebrain_StudyRound_id, sep = "_T")) %>% 
    filter(!`Demens::dementiaStatusAtEvaluationDate_binary` == 1 & !is.na(age) & !input %in% db.exclusion)
  
  # check if fs data
  df.mri = readxl::read_xlsx(file.path(sitedir, 
                            "mri_fs7_JamesProc", 
                            "FreeSurf711_LongPipe_VolumeAreaThickness_20210423.xlsx"))
  
  df.mri =
    df.mri %>% 
    rename("Lifebrain_Subject_id" = "subject") %>% 
    mutate(Lifebrain_StudyRound_id = if_else(TimePoint == "T5", 5, 
                                             if_else(TimePoint == "T6", 6,7)),
           input = paste(Lifebrain_Subject_id, Lifebrain_StudyRound_id, sep = "_T")) 
  df.mri = 
    df.mri %>% 
    filter(!is.na(`aseg:Left-Lateral-Ventricle.Volume_mm3`))
  # missing subjects
  df.out = inner_join(df.Umea, df.mri)
  

  
  save(df.out, 
       df.Umea,
       file = file.path(outdir, "df.rda")) 
  cat("umu data checked. rda file save with valid data only")
}

check.ub.fs = function(sitedir, outdir, mrifiles) {
  # load databases and get mri_filepaths
  grot.PDcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                    sheet= "PDcohort")
  grot.WAHAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "WAHAcohort")
  grot.MSAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                     sheet= "MSAcohort")
  grot.CRcohort = readxl::read_xlsx(file.path(sitedir,  "UB_Depression_late_Eniko_May2021.xlsx"), 
                                    sheet= "CRcohort")
  grot.GABAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "GABAcohort")
  grot.iTBScohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "iTBScohort")

  grot.WAHAcohort = grot.WAHAcohort %>% 
    rename(c("calculated_age_MRI" = "calculated_age",
             "Sex" = `Sex (1=woman)`,
             "Other_Praxis_Executive_ROCF_Copy" = "Praxis_Executive_ROCF_Copy")) %>% 
    mutate(Subject_id = as.character(Subject_id))
  grot.PDcohort = grot.PDcohort %>% 
    rename(c("calculated_age_MRI" = "calculated_age",
             "Sex" = `Sex (1=woman)`))

    
  grot.PDcohort  = mutate_ub(grot.PDcohort, T)
  grot.WAHAcohort  = mutate_ub(grot.WAHAcohort, T)
  grot.MSAcohort  = mutate_ub(grot.MSAcohort)
  grot.CRcohort = mutate_ub(grot.CRcohort)
  grot.GABAcohort = mutate_ub(grot.GABAcohort)
  grot.iTBScohort = mutate_ub(grot.iTBScohort)
  
  suppressMessages({
  df.UB = 
  list(grot.PDcohort, 
       grot.WAHAcohort, 
       grot.MSAcohort, 
       grot.CRcohort, 
       grot.GABAcohort, 
       grot.iTBScohort) %>% reduce(full_join)
  })
  df.UB = df.UB %>% 
    mutate(MMSE = if_else(MMSE < 0,NaN, MMSE)) %>% 
  drop_na(any_of(c("Subject_id",
               "Sex", 
               "Sess", 
               "calculated_age_MRI")))
  

  df.UB = 
    df.UB %>% filter(MMSE > 25 | is.na(MMSE)) 
  
  # # check if fs data
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    separate(input, 
             c("Subject_id", "sess"), 
             sep = "ses0", 
             remove = F) %>% 
    mutate(Subject_id = substr(Subject_id, 5, 30),
           Sess = substr(sess, 1,1) %>% as.numeric()) 
  
  # change ID strings for CR cohort
  stringi::stri_sub(df.mri$Subject_id[df.mri$Subject_id %>% startsWith("CTR")],7,6) <- "_"
  stringi::stri_sub(df.mri$Subject_id[grepl("^[[:digit:]]+", df.mri$Subject_id) & !grepl("[[:digit:]]+$", df.mri$Subject_id)],4,3)  <- "_"
  df.mri = df.ub.change.names(df.mri)
  
  
  df.out = inner_join(df.UB, df.mri)
 
  save(df.out, 
       df.UB,
       file = file.path(outdir, "df.rda")) 
  
  
  cat("ub data checked. linkage file created")
}



check.habs.fs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  df.habs = import(file.path(sitedir, "Demographics_HABS_DataRelease_2.0.csv"))
  names(df.habs)[1] = "SubjIDshort"
  df.grot1 =import(file.path(sitedir, "ClinicalMeasures_HABS_DataRelease_2.0.csv"))
  df.grot2 =import(file.path(sitedir, "Cognition_HABS_DataRelease_2.0.csv"))
  df.grot3 =import(file.path(sitedir, "PACC_HABS_DataRelease_2.0.csv"))
  #df.grot4 =import(file.path(sitedir, "ADNI_MRI_FS6_XSec_HABS_DataRelease_2.0.csv"))
  x = left_join(df.habs, df.grot1)
  xx = left_join(df.grot2, df.grot3)
  
  df.habs = left_join(x, xx) %>% 
    rename("SubjID" = "SubjIDshort")
  
  df.habs = 
    df.habs %>% 
    group_by(SubjID) %>% 
    mutate(AgeBsl = min(NP_Age),
           SubjectRound_id = if_else(StudyArc == "HAB_1.0", 0,
                                     if_else(StudyArc == "HAB_2.0", 12,
                                             if_else(StudyArc == "HAB_3.0", 24,
                                                     if_else(StudyArc == "HAB_4.0", 36,
                                                             if_else(StudyArc == "HAB_5.0", 48,60)))))) %>% 
    filter(HABS_DX == "CN")
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(SubjID = paste0("P_", substr(input, 9,14)),
           SubjectRound_id = substr(input, 19,20) %>% as.numeric(), 
           SubjectRound_id = if_else(SubjectRound_id == 18, 24,SubjectRound_id)) 
  # use for merging purposes use MR age for modeling. 
  
 # missing subjects
  df.out = inner_join(df.habs, df.mri)
  
  save(df.out, 
       df.habs,
       file = file.path(outdir, "df.rda")) 
  cat("habs data checked. rda file save with valid data only")
}


# check.ous.fs = function(sitedir, outdir, mrifiles) {
#   
#   x = haven::read_sav(file.path(sitedir, "CogNorm_incl_exl_criteria_Anders.SAV"), encoding="latin1")
#   load(file.path(sitedir, "NBM_MRI_long_new.Rda"))
#   df.ous = NBM_MRI_long_new %>% 
#     select(-c(contains("Vol"), 
#               contains("bil"),
#               eTIV))
#   
#   subs.exclude = 
#     c(x$ID[x$Umar_IL8 == 2], 
#       x$ID[x$MMSE_score_baseline < 26 & 
#              x$Cerad_umiddelbar_baseline_15SD %in% c(1,3,9)]) %>% 
#     unique()
#   
#   df.ous = 
#     df.ous %>% filter(!ID %in% subs.exclude)
#   
#   df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
#   for (df in 1:length(df.mri)) {
#     names(df.mri[[df]])[1] <- "input"
#   }
#   df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
#   
#   df.mri =
#     df.mri %>% 
#     mutate(ID = substr(input, 5,11) %>% 
#              as.character(),
#            Subject_Timepoint = substr(input, 16,16) %>% 
#              as.character(),
#            scanner = substr(input,17,30)) %>% 
#     separate(scanner, into = "scanner", sep = ".l")
#   
#   # missing subjects
#   df.out = inner_join(df.ous, df.mri)
#   
#   save(df.out, 
#        df.ous,
#        file = file.path(outdir, "df.rda")) 
#   cat("ous data checked. rda file save with valid data only")
#   
# }


check.ous.fs = function(sitedir, outdir, mrifiles) {
  
  load(file.path(sitedir, "cognNBM_updated.Rda"))
  df.scan.date = import(file.path(sitedir, "scan.date.csv")) %>% 
    separate(acq_time, c("scan_date","tod"), sep = " ") %>% 
    mutate(Subject_Timepoint = substr(ses, 6,6) %>% as.character(),
           scanner = substring(ses, 7), 
           ID = substring(sub,5)) 
  
  df.out = 
    cognNBM_updated %>% 
    group_by(ID) %>% 
    summarise(Birth_Date = first(Birth_Date), 
              Sex = first(Sex), 
              Edu_Years = first(Edu_Years))
  df.out = left_join(df.out, df.scan.date)
  
  df.out = 
    df.out %>% 
    mutate(MRI_Age = (as.Date(scan_date) - as.Date(Birth_Date)) %>% as.numeric(),
           MRI_Age = MRI_Age/365.25)
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(ID = substr(input, 5,11) %>% 
             as.character(),
           Subject_Timepoint = substr(input, 16,16) %>% 
             as.character(),
           scanner = substr(input,17,30)) %>% 
    separate(scanner, into = "scanner", sep = ".l")
  
  # missing subjects
  df.ous =  cognNBM_updated #%>% 
    # mutate(time.mmse = if_else(MMSE < 26, Years, NaN),
    #        time.norm = if_else(MMSE > 25, Years, NaN)) %>% 
    # group_by(ID) %>% 
    # mutate(max.time.mmse = max(time.mmse, na.rm = T),
    #        max.time.norm = max(time.norm, na.rm = T)) %>% 
    # ungroup() %>% 
    # mutate(crit.mmse = if_else(is.infinite(max.time.mmse) | Years <= max.time.norm, 1,0)) %>% 
    # filter(crit.mmse == 1)
    # 
  df.out = inner_join(df.out, df.mri)
  
  
  
  save(df.out, 
       df.ous,
       df.mri,
       file = file.path(outdir, "df.rda")) 
  cat("ous data checked. rda file save with valid data only")
  
}
check.aibl.fs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  
  df.aibl.demog =import(file.path(sitedir, "aibl_ptdemog_01-Jun-2018.csv")) %>% 
    select(-VISCODE)
  df.aibl.mmse =import(file.path(sitedir, "aibl_mmse_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.cdr =import(file.path(sitedir, "aibl_cdr_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.cogn =import(file.path(sitedir, "aibl_neurobat_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.dx = import(file.path(sitedir, "aibl_pdxconv_01-Jun-2018.csv"))  
  
  #df.aibl.visits = import(file.path(sitedir, "aibl_visits_01-Jun-2018.csv"))
  df.grot.mri =import(file.path(sitedir, "aibl_mri3meta_01-Jun-2018.csv")) %>% 
    mutate(Tscan = 3)
  df.grot.mri3 =import(file.path(sitedir, "aibl_mrimeta_01-Jun-2018.csv")) %>% 
    mutate(Tscan = 1.5)
  df.aibl.mri =
    rbind(df.grot.mri, 
          df.grot.mri3) %>% 
    filter(MMSMPRAGE == 1) %>% 
    select(-EXAMDATE)
  
  df.aibl = purrr::reduce(
    list(df.aibl.demog, 
         df.aibl.mmse, 
         df.aibl.cdr, 
         df.aibl.cogn, 
         df.aibl.dx,
         df.aibl.mri), 
    left_join)
  
  df.aibl = df.aibl %>% 
    mutate(visit = if_else(VISCODE == "bl", 0, 
                           if_else(VISCODE == "m18", 18, 
                                   if_else(VISCODE == "m36", 36, 
                                           if_else(VISCODE == "m54", 54, 
                                                   if_else(VISCODE == "m72", 72, NaN))))))
  
  df.aibl = df.aibl %>%
    group_by(RID) %>% 
    mutate(time.norm =if_else(DXCURREN == 1, visit, NaN), 
           time.norm = max(time.norm, na.rm = T)) %>%
    ungroup() %>% 
    mutate(crit.max.normal = if_else(visit <= time.norm & !is.nan(time.norm) & !is.infinite(time.norm), 1,0))  %>% 
    filter(crit.max.normal == 1)
  
  # df.aibl = df.aibl %>% 
  #   filter(DXCURREN == 1 & MMSMPRAGE == 1)
  
  df.aibl.idaSearch = import(file.path(sitedir, "idaSearch_12_05_2022.csv")) %>% 
    separate(`Imaging Protocol`, 
             c("grot1","field_strength", "grot3", "model", "grot4", "slice_thicknes"), 
             sep="=|;") %>% 
    mutate(Session_ID = if_else(Visit =="Baseline", "bl",
                                if_else(Visit =="18 Month follow-up", "m18",
                                        if_else(Visit =="36 Month follow-up", "m36",
                                                if_else(Visit =="54 Month follow-up", "m54",
                                                        if_else(Visit =="72 Month follow-up", "m72","")))))) %>% 
    rename("Subjects_ID" = "Subject ID") %>% 
    select(-c(contains("grot"),
              Description, 
              Sex,
              Visit)) 
  
  df.aibl.idaSearch = 
    df.aibl.idaSearch %>% 
    group_by(Subjects_ID, Session_ID) %>% 
    summarise(field_strength = first(field_strength), 
              slice_thicknes = first(slice_thicknes), 
              Session_ID = first(Session_ID), 
              model = first(model), 
              Age = mean(Age, na.rm = T))
  
  df.aibl.paths = import(file.path(sitedir, "t1_paths_aibl.tsv")) %>% 
    mutate(Session_ID = if_else(Session_ID == "M00", "bl", Session_ID)) %>% 
    group_by(Subjects_ID, Session_ID) %>% 
    tally()
  
  
  df.aibl.paths = 
    inner_join(df.aibl.paths, df.aibl.idaSearch) %>% 
    rename(c("VISCODE" = "Session_ID", 
             "RID" = "Subjects_ID")) 
  
  df.aibl = left_join(df.aibl, df.aibl.paths)
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    mutate(grot1 = gsub("sub-AIBL", "", input)) %>% 
    separate(grot1, c("RID",
                      "grot2"),
             sep = "sesM") %>% 
    mutate(visit = substr(grot2, 1,2) %>% as.numeric(),
           RID = as.numeric(RID))
  
  
  # missing subjects
  df.out = inner_join(df.aibl, df.mri)
  
  save(df.out, 
       df.aibl,
       file = file.path(outdir, "df.rda")) 
  cat("aibl data checked. rda file save with valid data only")
}

check.adni.fs = function(sitedir, outdir, mrifiles) {
  #ADNI.stripped = read.csv(file.path(sitedir, "ADNI_stripped.csv"))
  df.adni.participants =import(file.path(sitedir, "participants.tsv")) %>% 
    rename("RID" = "adni_rid") %>%
    mutate(RID = as.numeric(RID))
  df.adni.conversioninfo =import(file.path(sitedir, "t1_paths_new.tsv")) 
  
  load(file.path(sitedir, "adni.sessions.Rda"))
  df.adni.sessions[df.adni.sessions == "n/a"] = NA 
  idx = !apply (is.na(df.adni.sessions), 2, all)
  df.adni.sessions = df.adni.sessions[, ..idx ]
  df.adni = inner_join(df.adni.sessions, df.adni.participants)
  
  mri.info15 = load_scan_info_adni("idaSearch_1_25_2022_1.5T.csv", "1.5")
  mri.info3 = load_scan_info_adni("idaSearch_1_25_2022_3T.csv", "3.0")
  
  mri.info = rbind(mri.info15, 
                   mri.info3)
  
  df.adni.mri = 
    inner_join(df.adni.conversioninfo,mri.info) %>% 
    mutate(Subject.ID = as.character(Subject.ID)) %>% 
    separate(Subject.ID, c("Site", "ID"), sep = "_S_", remove = F)
  
  df.adni.mri$model = ""
  models = df.adni.set.scans()
  
  for (x in 1:length(models)) {
    xx = names(models)[[x]]
    df.adni.mri = 
      df.adni.mri %>% 
      mutate(model = if_else(Imaging.Protocol %in% models[[x]], xx, model))
  }
  
  df.adni.mri = 
    df.adni.mri %>% 
    mutate(RID = ID %>% as.numeric(),
           session_id = gsub("m", "ses-M", VISCODE),
           session_id = if_else(session_id =="bl", "ses-M00", session_id))
  
  
  df.adni = left_join(df.adni,df.adni.mri)
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    mutate(grot1 = gsub("sub-ADNI", "", input)) %>% 
    separate(grot1, c("grot4",
                      "grot2"),
             sep = "ses") %>% 
    separate(grot2, 
             c("session_id", "grot3"), 
             sep = ".long") %>%
    separate(grot4, 
             c("grot5", "RID"), 
             sep = "S") %>% 
    mutate(RID = as.numeric(RID),
           session_id = paste("ses", session_id, sep = "-")) %>% 
    select(-starts_with("grot"))
  
  
  df.adni.hc = df.adni %>%
    group_by(RID) %>% 
    mutate(time.norm =if_else(diagnosis == "CN", Age, NaN), 
           time.norm = max(time.norm, na.rm = T)) %>%
    ungroup() %>% 
    mutate(crit.max.normal = if_else(Age <= time.norm & !is.nan(time.norm) & !is.infinite(time.norm), 1,0))  %>% 
    filter(crit.max.normal == 1)
  
  df.out = inner_join(df.adni.hc, df.mri)
  
  
  save(df.out, 
       df.adni,
       df.adni.hc,
       df.mri,
       file = file.path(outdir, "df.rda")) 
  cat("adni data checked. rda file save with valid data only")
}

check.preventad.fs = function(sitedir, outdir, mrifiles) {
  datadir = file.path(sitedir, "2022_11_23", "tabular")
  
  df.preventAD.main = import(file.path(sitedir, "data-2022-12-11T20_48_55.373Z.csv")) %>% 
    rename(c("CONP_ID" = "PSCID",
             "CONP_CandID" = "DCCID",
             "Study_visit_label" = "Visit Label"))
  
  df.preventAD.demographics = import(file.path(datadir, "Demographics_Registered_PREVENTAD.csv"))
  df.preventAD.genetics = import(file.path(datadir, "Genetics_Registered_PREVENTAD.csv"))
  df.preventAD.DKEFS = import(file.path(datadir, "Neuropsych_DKEFS-CWIT_Registered_PREVENTAD.csv"))
  df.preventAD.RAVLT = import(file.path(datadir, "Neuropsych_RAVLT_Registered_PREVENTAD.csv"))
  df.preventAD.TMT = import(file.path(datadir, "Neuropsych_TMT_Registered_PREVENTAD.csv" ))
  df.preventAD.RBANS = import(file.path(datadir, "RBANS_Registered_PREVENTAD.csv" ))
  df.preventAD.AD8 = import(file.path(datadir, "AD8_Registered_PREVENTAD.csv" ))
  df.preventAD.APS = import(file.path(datadir, "APS_Registered_PREVENTAD.csv" ))
  df.preventAD.CDR = import(file.path(datadir, "EL_CDR_MoCA_Registered_PREVENTAD.csv" ))
  df.preventAD.MedHist = import(file.path(datadir, "EL_Medical_history_Registered_PREVENTAD.csv"))                                                   
  df.preventAD.Reg = import(file.path(datadir, "Lab_Registered_PREVENTAD.csv"))
  #df.preventAD.siblings = import(file.path(datadir, "List_of_participants_with_only_1_sibling.txt"))
  #df.preventAD.Switch = import(file.path(datadir, "List_of_participants_switched_back_to_cohort.txt"))
  
  
  df.preventAD.main = 
    inner_join(df.preventAD.demographics, 
               df.preventAD.main)
  df.preventAD.main$age_at_mci = NaN
  df.preventAD.main$age_at_mci[df.preventAD.main$probable_MCI_visit == df.preventAD.main$Study_visit_label ] = df.preventAD.main$`Age At MRI In Months`
  
  
  df.preventAD.main = 
    df.preventAD.main %>% 
    group_by(CONP_ID) %>% 
    mutate(age_at_mci = min(age_at_mci, na.rm = T), 
           age = `Age At MRI In Months`/12,
           mci_visit = if_else(age_at_mci <= `Age At MRI In Months`, 1, 0))
  
  
  df.preventAD = full_join(df.preventAD.Reg, df.preventAD.main)
  
  df.preventAD = 
    reduce(
      list(df.preventAD, 
           df.preventAD.AD8, 
           df.preventAD.APS, 
           df.preventAD.CDR, 
           df.preventAD.DKEFS,
           #df.preventAD.genetics, 
           df.preventAD.MedHist, 
           df.preventAD.RAVLT, 
           df.preventAD.RBANS, 
           df.preventAD.TMT),
      left_join,
      by = c("CONP_ID", "CONP_CandID",  "Visit_label")
    )
  
  df.preventAD = 
    left_join(df.preventAD,
              df.preventAD.genetics)
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    mutate(grot1 = gsub("sub-", "", input)) %>% 
    separate(grot1, c("CONP_CandID",
                      "grot2"),
             sep = "ses") %>% 
    separate(grot2, 
             c("Study_visit_label.x", "grot3"), 
             sep = ".long") %>% 
    mutate(CONP_CandID = as.numeric(CONP_CandID)) %>% 
    select(-starts_with("grot"))
  
  df.preventAD = df.preventAD %>% filter(mci_visit == 0 | is.na(mci_visit))
  df.out = inner_join(df.preventAD, df.mri) 
  
  save(df.out, 
       df.preventAD,
       df.mri,
       file = file.path(outdir, "df.rda")) 
  cat("preventad data checked. rda file save with valid data only")
}




check.ukb.fs = function(sitedir, outdir, mrifiles) {
  # load databases and get mri_filepaths
  #list.files(file.path(sitedir, "45249", "normative_modelling"))
  dfdir=file.path(sitedir, "45249", "normative_modelling")
  
  load(file.path(dfdir, "basic_demographics.Rda"))
  #load(file.path(dfdir, "mri_data_tp23.Rda"))
  load(file.path(dfdir, "demographics_expanded.Rda"))
  
  df.ukb = df.dmg %>% filter(!is.na(age_3)) %>% 
    group_by(eid) %>% 
    select(-(c(age_0, age_1))) %>% 
    mutate(agebsl = min(c(age_2, age_3)))
  
  df.ukb = 
    left_join(df.ukb, 
              df.dmg.extra %>% 
                select(
                  eid,
                  uk_biobank_assessment_centre_f54_2_0,
                  uk_biobank_assessment_centre_f54_3_0
                ))
  grot.idx = grepl("age_",names(df.ukb))
  names(df.ukb)[grot.idx] = paste0(names(df.ukb)[grot.idx], "_0")
  
  df.ukb = 
    df.ukb %>% 
    mutate(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0),
           t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0),
           uk_biobank_assessment_centre_f54_2_0 = as.numeric(uk_biobank_assessment_centre_f54_2_0),
           uk_biobank_assessment_centre_f54_3_0 = as.numeric(uk_biobank_assessment_centre_f54_3_0)) %>% 
    pivot_longer(-c(eid, Sex, agebsl), 
                 names_to = "names", 
                 values_to = "values") %>% 
    mutate(grot1 = stringi::stri_reverse(names)) %>% 
    separate(grot1, 
             c("grot2", "tp", "grot3"),
             sep ="_",
             extra = "merge") %>% 
    mutate(names = stringi::stri_reverse(grot3)) %>% 
    select(-starts_with("grot")) %>% 
    pivot_wider(names_from = names, 
                values_from = values)
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(dfdir, "tsd_processed_fs.7.1.0", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(eid = substr(input, 5,11) %>% as.integer(),
           tp = substr(input, 15,15)) 
  
  df.out = inner_join(df.ukb, df.mri)
  
  save(df.out, 
       df.ukb,
       file = file.path(outdir, "df.rda")) 
  cat("ukb data checked. rda file save with valid data only")
}




check.oasis.fs = function(sitedir, outdir, mrifiles) {
  datadir = file.path(sitedir, "data", "data_files")
  
  #load(file.path(sitedir,"sessions.info.extended.Rda"))
  load(file.path(sitedir,"sessions.info.Rda"))
  
  sessions = sessions %>% 
    mutate(Subject = substr(subs, 5,12), 
           days_to_visit = gsub("ses-d", "", ses) %>% as.numeric(), 
           scanner = paste(ManufacturersModelName, DeviceSerialNumber, sep = "_"))
  
  
  df.oasis.clinical = import(file.path(datadir, "UDSb4/csv/OASIS3_UDSb4_cdr.csv"))
  df.oasis.demo = import(file.path(datadir, "OASIS3_demographics.csv"))
  df.oasis.mri.info = import(file.path(datadir, "MRI_info.csv"))
  df.oasis.psych = import(file.path(datadir, "pychometrics/csv/OASIS3_UDSc1_cognitive_assessments.csv"))
  df.oasis.cognorm = import(file.path(sitedir,"data/cohort_files/CogNorm/csv/OASIS3_unchanged_CDR_cognitively_healthy.csv"))
  
  
  idx = !colSums(is.na(df.oasis.clinical)) == length(df.oasis.clinical$OASISID)
  df.oasis.clinical = df.oasis.clinical[,idx] %>% 
    mutate(session_clinical = OASIS_session_label) %>% 
    select(-OASIS_session_label) 
  
  idx = !colSums(is.na(df.oasis.demo)) == length(df.oasis.demo$OASISID)
  df.oasis.demo = df.oasis.demo[,idx]
  
  idx = !colSums(is.na(df.oasis.psych)) == length(df.oasis.psych$OASISID)
  df.oasis.psych = df.oasis.psych[,idx] %>% 
    mutate(session_psych = OASIS_session_label) %>% 
    select(-OASIS_session_label) 
  
  
  df.oasis.mri.info = 
    df.oasis.mri.info %>% 
    separate(`MR ID`, c("OASISID", "days_to_visit"), sep = "_MR_d", remove = F) %>% 
    filter(!is.na(days_to_visit)) %>% 
    mutate(days_to_visit = as.numeric(days_to_visit),
           `age at visit` = Age)
  
  df.oasis.cognorm= 
    df.oasis.cognorm %>%  
    mutate(OASISID = OASIS3_id, 
           cognorm = 1) %>% 
    select(OASISID, cognorm)
  
  
  
  df.oasis3 =
    reduce(list(
      df.oasis.psych, 
      df.oasis.mri.info,
      df.oasis.clinical,
      df.oasis.demo, 
      df.oasis.cognorm),
      full_join
    )
  
  df.oasis3 = 
    df.oasis3 %>% 
    mutate(time = days_to_visit/366.25, 
           AGE = AgeatEntry + time) %>% 
    filter(!time < 0)
  
  
  df.oasis3 = 
    df.oasis3 %>% 
    mutate(time.dx1.normal =if_else(dx1 %in% c("Cognitively normal","No dementia"), time, NaN),
           time.dx1.exclude =if_else(!dx1 %in% c("Cognitively normal","No dementia") & !is.na(dx1), time, NaN)) %>% 
    group_by(OASISID) %>% 
    mutate(time.dx1.normal.max = max(time.dx1.normal, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(time.dx1.exclude.filt = if_else(time.dx1.exclude > time.dx1.normal.max, time.dx1.exclude, NaN)) %>% 
    group_by(OASISID) %>% 
    mutate(time.dx1.exclude.min = min(time.dx1.exclude.filt, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(time.dx1.criteria = (time.dx1.exclude.min + time.dx1.normal.max)/2,
           crit.dx1.normal.max = if_else(time <= time.dx1.criteria, 1,0 ))
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  
  
  df.mri =
    df.mri %>% 
    mutate(Subject = substr(input, 5,12),
           days_to_visit = substr(input, 17,20) %>% 
             as.numeric())
  
  df.mri = inner_join(df.mri, sessions)
  
  df.out = inner_join(df.mri, df.oasis3)
  
  # df.out = 
  #   df.out %>% filter(cognorm == 1)
  # 
  df.out = 
    df.out %>% filter(crit.dx1.normal.max == 1)
  
  
  save(df.out, 
       df.oasis3,
       file = file.path(outdir, "df.rda")) 
  cat("oasis3 data checked. rda file save with valid data only")
}

check.bbhi.fs = function(sitedir, outdir, mrifiles) {
  
  db.bbhi = import(file.path(sitedir, "BBHI_NPS-2tps_sexAgeYoE.csv"))
  db.bbhi.apoe = import(file.path(sitedir, "BBHI_APOE.xlsx"))
  names(db.bbhi.apoe)[1] = "id_user"  
  db.bbhi = left_join(db.bbhi, db.bbhi.apoe)
  
  df.bbhi = db.bbhi %>% 
    select(-c(contains("comment"), V1)) %>%
    rename("w2_tmtb" = "w2_tmt_b_raw",
           "w1_tmtb" = "w1_tmt_b_raw") %>% 
    pivot_longer(-c(id_user, birth_date, sex, YoE, apoe), 
                 names_to = c("variable", ".value"),
                 names_sep = "_") %>% 
    drop_na(any_of("nps")) %>% 
    mutate(sess = gsub("w","", variable) %>% as.numeric())
  
  
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    separate(input, c("grot1",
                      "id_user"),
             sep = ".long.",
             remove = F) %>% 
    separate(grot1, 
             c("grot2", "sess"), 
             sep = "_") %>% 
    select(-starts_with("grot")) %>% 
    mutate(sess = as.numeric(sess),
           id_user = as.numeric(id_user))
  
  
  df.bbhi = 
    df.bbhi %>% 
    mutate(age = (as.Date(nps)-as.Date(birth_date, "%d/%m/%Y")) %>% as.numeric, 
           age = age/365.25) %>% 
    filter(age > 40)
  
  # missing subjects
  df.out = inner_join(df.bbhi, df.mri) 
  
  #df.notinmri = anti_join(db.bbhi, df.mri) 
  #df.notinbehav = anti_join(df.mri, db.bbhi) 
  
  save(df.out, 
       df.bbhi,
       df.mri, 
       file = file.path(outdir, "df.rda")) 
  cat("bbhi data checked. rda file save with valid data only")
}



check.vetsa.fs = function(sitedir, outdir, mrifiles) {
  
  db = import(file.path(sitedir, "Cognitive",  "admin_2023_07_07.sas7bdat"))
  # db.cogn1 = import(file.path(sitedir, "Cognitive",  "cid_v1_2023_07_07.csv"))
  # db.cogn2 = import(file.path(sitedir, "Cognitive",  "cid_v2_2023_07_07.csv"))
  # db.cogn3 = import(file.path(sitedir, "Cognitive",  "cid_v3_2023_07_07.csv"))
  db = db %>% select(CID, 
                     Case, 
                     Twin, 
                     AGE_V1, 
                     AGE_V2, 
                     AGE_V3, 
                     DECEASED_STATUS, 
                     TEDALL, 
                     ETHNALL, 
                     RACEALL, 
                     DADOCCAVG, 
                     apoe2019, 
                     year_born, 
                     LT_OCC, 
                     HANDSW, 
                     OLDTOT, 
                     ZHAND)
  
  db = db %>% 
    pivot_longer(c("AGE_V1", 
                   "AGE_V2", 
                   "AGE_V3"), 
                 names_to = "wave", 
                 values_to = "age") %>% 
    filter(!is.na(age)) %>% 
    mutate(sex = "male", 
           wave = if_else(wave == "AGE_V1", 1, 
                          if_else(wave == "AGE_V2", 2,3)))
  
  
  db.mri_v1 = import(file.path(sitedir, "Cognitive", "MRI_readme", "vetsa1_MRI_admin_cid_2023-06-12.csv"))
  db.mri_v2 = import(file.path(sitedir, "Cognitive", "MRI_readme", "vetsa2_MRI_admin_cid_2023-06-12.csv"))
  db.mri_v3 = import(file.path(sitedir, "Cognitive", "MRI_readme", "vetsa3_MRI_admin_cid_2023-06-12.csv"))
  df.inclusion = import(file.path(sitedir, "Cognitive", "MRI_readme", "VETSA123_healthMRIexclusions_cid_2023-06-28.csv"))
  
  vars_db_mri = 
    c("CID",
      "Case",
      "Twin",
      "scanner", 
      "acq_T1",
      "QC_flag")
  
  db.mri_v1 = db.mri_v1 %>%
    rename("scanner" = "scanner_v1", 
           "acq_T1" = "acq_T1_v1",
           "QC_flag" = "QC_flag_v1") %>% 
    select(vars_db_mri) %>% 
    mutate(wave = 1)
  
  db.mri_v2 = db.mri_v2 %>%
    rename("scanner" = "scanner_v2", 
           "acq_T1" = "acq_SagT1_v2",
           "QC_flag" = "QC_flag_v2") %>% 
    select(vars_db_mri) %>% 
    mutate(wave = 2)
  
  db.mri_v3 = db.mri_v3 %>%
    rename("scanner" = "scanner_v3", 
           "acq_T1" = "acq_SagT1_v3",
           "QC_flag" = "QC_flag_v3") %>% 
    select(vars_db_mri)  %>% 
    mutate(wave = 3)
  
  db.mri = 
    list(db.mri_v1, 
         db.mri_v2, 
         db.mri_v3) %>% 
    data.table::rbindlist() %>% 
    filter(!QC_flag == 2)
  
  df.inclusion = df.inclusion %>% 
    filter(!MRI_addtl_exclusion_v3 == 1, 
           !SCHZYN_V1 == 1, 
           !SCHZYN_V2 == 1,
           !SCHZYN_V3 == 1,
           !AIDSYN_V1 == 1,
           !AIDSYN_V2 == 1,
           !AIDSYN_V3 == 1,
           !MSYN_V1 == 1,
           !MSYN_V2 == 1,
           !MSYN_V3 == 1,)
  
  db.vetsa = inner_join(db.mri, df.inclusion)
  db.vetsa = inner_join(db.vetsa, db)
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    separate(input, c("grot1",
                      "CID"),
             sep = ".long.",
             remove = F) %>% 
    separate(grot1, 
             c("grot2", "wave"), 
             sep = "_v") %>% 
    select(-starts_with("grot")) %>% 
    mutate(wave = as.numeric(wave),
           CID = as.numeric(CID))
  
  
  df.out = inner_join(db.vetsa, df.mri)
  save(df.out, 
       db,
       df.mri, 
       file = file.path(outdir, "df.rda")) 
  cat("vetsa data checked. rda file save with valid data only")
  
}
