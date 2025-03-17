library(dplyr)
library(eeptools)

## set_links function ##
set_links = function(folder) {
  tabulateddata <<- here('data-raw/tabulated')
  sitedir <<- file.path(tabulateddata, folder)
  outdir <<- here('data_memory_long', folder)
  if (!dir.exists(outdir)) {dir.create(outdir)}
}

## function to filter longitudinal memory and mri data 
filter.long = function(df1, df2){
  #find longitudinal data in df1 (n>1)
  subs.df1 = df1 %>% 
    group_by(subject_id) %>% 
    tally() %>% 
    filter(n>1)
  #find longitudinal data in df2 (n>1)
  subs.df2 = df2 %>% 
    group_by(subject_id) %>% 
    tally() %>% 
    filter(n>1)
  #filter in all longitudinal data in df1
  df.long <- df1 %>% 
    filter(subject_id %in% subs.df1$subject_id)
  #filter data that also have longitudinal data in df2
  df.long <- df.long %>% 
    filter(subject_id %in% subs.df2$subject_id)
  return(df.long)
}

##mutate_ub function
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

## adni scan info
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



extract.longdata.uio = function(sitedir, outdir, mrifiles) { 
  #load tabular data
  df.uio = import(file.path(sitedir, "noas_query_2022-11-29_21-26-23_041b29e_b0f5e7b.csv")) %>% 
    mutate(subses_id = paste(subject_id, visit_number, sep = '_'))
  # exclude non-shareable, mean age under 20y, max mmse score under 25, projects with memory training
  df.uio <- df.uio %>% filter(subject_shareable == 1)
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
    
   ## subjects that have a low mms score in some timepoints, filter in all timepoints up to the timepoint with low score
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
    
    df.uio <- df.uio %>% 
      filter(!project_id %in% c("Loci", "NCP")) #filter out projects with memory training

      
  })
  
  
  mem_vars <- 
    c("cvlt_a5",
      "cvlt_a_total", 
      "cvlt_5min_free", 
      "cvlt_30min_free")
  
  # Remove rows where all specified columns are NA
  df.memory <-  df.uio[!apply( df.uio[mem_vars], 1, function(row) all(is.na(row))), ]
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  cat("uio data extracted. saved as df.rda in /other_output/")
}


## for umu ##
extract.longdata.umu = function(sitedir, outdir, mrifiles) {
  # load tab data
  grot.t1 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T1.xlsx"))
  grot.t2 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T2.xlsx"))
  grot.t3 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T3.xlsx"))
  grot.t4 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T4.xlsx"))
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
  df.Umea = data.table::rbindlist(list(grot.t1,
                                       grot.t2,
                                       grot.t3,
                                       grot.t4,
                                       grot.t5, 
                                       grot.t6, 
                                       grot.t7), 
                                  use.names = F) %>% 
    rowwise() %>% 
    mutate(age = mean(c(rounded_age_HT, rounded_age_MT), na.rm = T),
           input = paste(Lifebrain_Subject_id, Lifebrain_StudyRound_id, sep = "_T")) %>% 
    filter(!`Demens::dementiaStatusAtEvaluationDate_binary` == 1 & !is.na(age) & !input %in% db.exclusion) 
  

  df.memory =
    df.Umea %>% 
    rename("subject_id" = "Lifebrain_Subject_id")%>% 
    select(Lifebrain_Site_id:Lifebran_Sex, "SPT_VTFreeRecall::vtb", age) %>%
    mutate(subses_id = paste(subject_id, Lifebrain_StudyRound_id, sep = '_')) %>% 
    rename(VT_free_recall_vtb = "SPT_VTFreeRecall::vtb") %>% 
    filter(!is.na(VT_free_recall_vtb))
  
  
  
  # # #get a sublist with subject with more than 2 TPs
  # sub.mem <- df.memory %>% 
  #   group_by(subject_id) %>% 
  #   tally() %>% 
  #   filter(n>1)
  # #filter in all longitudinal data in df1
  # only.memory.long <- df.memory %>% 
  #   filter(subject_id %in% sub.mem$subject_id)
  # 
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  cat("umu data extracted. saved as df.rda in /other_output/")
  
  
}


## for ub ##
extract.longdata.ub = function(sitedir, outdir, mrifiles) {
  # load databases and get mri_filepaths
  grot.PDcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                    sheet= "PDcohort")
  grot.WAHAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "WAHAcohort")
  
  grot.WAHAcohort = grot.WAHAcohort %>% 
    #rename("calculated_age_MRI" = "calculated_age") %>% 
    rename("Sex" = `Sex (1=woman)`) %>% 
    rename("Other_Praxis_Executive_ROCF_Copy" = "Praxis_Executive_ROCF_Copy") %>% 
    mutate(Subject_id = as.character(Subject_id)) %>% 
    mutate(project_id = "WAHA")
  
  grot.PDcohort = grot.PDcohort %>%
    #rename("calculated_age_MRI" = "calculated_age") %>%
    rename("Sex" = `Sex (1=woman)`) %>% 
    mutate(Subject_id = as.character(Subject_id)) %>% 
    mutate(project_id = "PD")
  
  grot.PDcohort  = mutate_ub(grot.PDcohort, T)
  grot.WAHAcohort  = mutate_ub(grot.WAHAcohort, T)
  
  
  df.UB <- full_join(grot.PDcohort,grot.WAHAcohort)
  
  df.UB <- df.UB %>% 
    filter(MMSE > 25)  %>% 
    drop_na(any_of(c("Subject_id",
                     "Sex", 
                     "Sess", 
                     "calculated_age")))
  
  #df.UB %>% filter(MMSE > 25 | is.na(MMSE))
  
  df.UB <- df.UB %>% 
    rename("subject_id" = "Subject_id")%>% 
    rename("Round_id" = "Round_ id") %>% 
    mutate(Sess = substr(Round_id,3,3) %>% as.numeric()) %>% 
    mutate(subses_id = paste(subject_id, Sess, sep = '_'))
  
  ##memory data ##
  df.memory <- df.UB %>% 
    select(subject_id:calculated_age,project_id,subses_id, contains(c("Verbal_memory", "Verbal_Memory"))) %>% 
    filter(!if_all(Verbal_memory_RAVLT_trial1:Verbal_Memory_RAVLT_delayed, is.na)) %>% 
    filter(!Verbal_memory_RAVLT_trial1 == -9999)
  

  
  # # #get a sublist with subject with more than 2 TPs
  # sub.mem <- df.memory %>% 
  #   group_by(subject_id) %>% 
  #   tally() %>% 
  #   filter(n>1)
  # #filter in all longitudinal data in df1
  # only.memory.long <- df.memory %>% 
  #   filter(subject_id %in% sub.mem$subject_id)
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  cat("ub data extracted. saved as df.rda in /other_output/")
}

## for mpib ## 
extract.longdata.mpib = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING datals
  df.mem1 <- read_tsv(file.path(sitedir,"verbal_learning.tsv")) %>% 
    select(!"test_id")
  df.mem2 <- read_tsv(file.path(sitedir,"vl_recall.tsv"))%>% 
    select(!"test_id")
  df.mem.scene.encoding = rio::import(file.path(sitedir, "behavioural_data", "indoor_outdoor.tsv"))%>% 
    select(!c(test_id, do_t))
  df.mem.face_profession = rio::import(file.path(sitedir, "behavioural_data", "face_profession.tsv"))%>% 
    select(!c(test_id, do_t))
  df.mem.object_location = rio::import(file.path(sitedir, "behavioural_data", "object_location.tsv"))%>% 
    select(!c(test_id, do_t, age, tod))
  
  "object_location.tsv"
  df.mem <- full_join(df.mem1, df.mem2) %>% 
    full_join(., df.mem.scene.encoding) %>% 
    full_join(., df.mem.face_profession) %>% 
    full_join(., df.mem.object_location) %>% 
    rename(round_id = study_round_id)
  
  df.MPIB = read.csv(file.path(sitedir,"Data_BASEII_YlvaEniko211011.csv")) 
  df.id.mpib = df.MPIB %>%
    mutate(ses01 = if_else(!is.na(MRT12_age), 1,NaN),
           ses02 = if_else(!is.na(MRT16_age), 2,NaN)) %>% 
    pivot_longer(cols = c(ses01,ses02), 
                 names_to = "grot1",
                 values_to = "round_id",
                 values_drop_na = T) %>% 
    mutate(visit_age = if_else(round_id ==1, MRT12_age, MRT16_age)) %>% 
    select(-c(Age, 
              MPIB_2012_Date1, 
              MRT12_age, 
              MRT16_age)) %>% 
    filter(MMSE > 25) %>% 
    rename(subject_id = ID)
  # select(ID,round_id,sex ,visit_age) %>%
  # rename(SubjectRound_id = round_id)
  
  
  
  df.mpib<- left_join(df.mem, df.id.mpib)
  
  demo <- read_tsv(file.path(sitedir,"sheet_1.tsv")) %>% 
    #select(!age_group, edu) %>% 
    rename(sex_1 = sex) %>% 
    mutate(sex_1 = if_else(sex_1 == "male", 1, 2)) %>% 
    mutate(round_id = if_else(round_id == "R01", 1,
                              if_else(round_id =="R02", 2, 3))) %>% 
    #rename(SubjectRound_id = round_id) %>% 
    distinct()
  
  
  df.mpib<- left_join(df.mpib, demo) %>% 
    mutate(sex_new = coalesce(sex,sex_1)) %>% 
    select(-sex, -sex_1) %>% 
    rename(sex = sex_new)
  
  df.mpib <- df.mpib %>% 
    group_by(subject_id) %>% 
    fill(sex, .direction ="downup") %>% 
    ungroup()
  
  mem_vars = c(
    "io_hit_minus_fa",
    "fp_hit_minus_new_fa",
    "fp_hit_minus_rear_fa",
    "ol_acc", 
    "recall1",
    "trial1",
    "trial2",
    "trial3",
    "trial4",
    "trial5",
    "trial7"
  )
  
  #filter out those with no memory data
  df.memory <-  df.mpib[!apply( df.mpib[mem_vars], 1, function(row) all(is.na(row))), ]
 
 
  df.memory <-df.memory %>% 
    filter(!is.na(sex)) %>% 
    mutate(subses_id = paste(subject_id, round_id, sep = '_'))

  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  
  cat("mpib data extracted. saved as df.rda in /other_output/")
  
}

## for adni ##
extract.longdata.adni = function(sitedir, outdir, mrifiles) {
  #ADNI.stripped = read.csv(file.path(sitedir, "ADNI_stripped.csv"))
  df.adni.participants =import(file.path(sitedir, "participants.tsv")) %>% 
    rename("RID" = "adni_rid") %>%
    mutate(RID = as.numeric(RID))
  
  load(file.path(sitedir, "adni.sessions.Rda"))
  df.adni.sessions[df.adni.sessions == "n/a"] = NA 
  idx = !apply (is.na(df.adni.sessions), 2, all)
  df.adni.sessions = df.adni.sessions[, ..idx ]
  df.adni = inner_join(df.adni.sessions, df.adni.participants)
  
  df.psych <- import(file.path(sitedir, "studydata", "UWNPSYCHSUM_03_09_21.csv")) %>% 
    mutate(RID = RID %>% as.numeric(),
           session_id = gsub("m", "ses-M", VISCODE2),
           session_id = if_else(session_id =="bl", "ses-M00", session_id))
  
  #df.adni <- left_join(df.adni, df.psych)
  df.adni <- full_join(df.adni, df.psych)
  
  adni.merge <- import(file.path(sitedir, "studydata", "ADNIMERGE.csv")) %>% 
    mutate(RID = RID %>% as.numeric(),
           session_id = gsub("m", "ses-M", VISCODE),
           session_id = if_else(session_id =="bl", "ses-M00", session_id)) %>% 
    select(-MMSE,-FAQ,-VISCODE)
  
  #df.adni <- left_join(df.adni, adni.merge)
  df.adni <- full_join(df.adni, adni.merge) %>% 
    mutate(age = age %>% as.numeric())
  
  
  df.adni <- df.adni %>% 
    rename(subject_id = RID) %>% 
    mutate(subses_id = paste(subject_id, session_id, sep = '_'))
  
  
  
  df.adni.hc = df.adni %>%
    group_by(subject_id) %>% 
    mutate(time.norm =if_else(diagnosis == "CN", age, NaN), 
           time.norm = max(time.norm, na.rm = T)) %>%
    ungroup() %>% 
    mutate(crit.max.normal = if_else(age <= time.norm & !is.nan(time.norm) & !is.infinite(time.norm), 1,0))  %>% 
    filter(crit.max.normal == 1)
  
  
  
  df.memory <- df.adni.hc %>% 
    filter(!is.na(ADNI_MEM))
  
  
  
  load("/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/adni/df.rda")
  df.out <- df.out %>% 
    rename(subject_id = RID)
  
  mricheck <- df.out %>%
    group_by(subject_id) %>%
    mutate(last_mri = max(Scan_Date)) %>%
    ungroup() %>%
    select(subject_id, session_id, Scan_Date, last_mri)
  
  
  merge <- left_join(df.memory, mricheck) %>%
    group_by(subject_id) %>%
    fill(last_mri, .direction = "down") %>%
    ungroup()
  
  exc.data <- merge %>%
    filter(!is.na(last_mri)) %>%
    mutate(differ = difftime(examdate, last_mri, units = "days"))
  
  exc.data <- exc.data %>%
    filter(differ>183)
  
  
  #exclude memory data after 6 months after last MRI
  df.memory <- df.memory %>%
    filter(!subses_id %in% exc.data$subses_id)
 
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  cat("adni data extracted. saved as df.rda in /other_output/")
}



## for aibl
extract.longdata.aibl = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  df.aibl.demog =import(file.path(sitedir, "aibl_ptdemog_01-Jun-2018.csv")) %>% 
    select(-VISCODE)
  df.aibl.mmse =import(file.path(sitedir, "aibl_mmse_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.cdr =import(file.path(sitedir, "aibl_cdr_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.dx = import(file.path(sitedir, "aibl_pdxconv_01-Jun-2018.csv")) 
  df.aibl.cogn =import(file.path(sitedir, "aibl_neurobat_01-Jun-2018.csv")) # %>% 
  #select(-EXAMDATE)
  
  
  df.aibl = purrr::reduce(
    list(df.aibl.demog, 
         df.aibl.mmse, 
         df.aibl.cdr, 
         df.aibl.cogn, 
         df.aibl.dx), 
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
  
  
  df.aibl <- df.aibl %>% 
    mutate(PTDOB = substr(PTDOB, 2,5) %>%  as.numeric(),
           EXAMDATE = substr(EXAMDATE, 7, 10) %>%  as.numeric(),
           age_Cog = EXAMDATE-PTDOB) 

 
   
  
  #rename subject_id and create subses_id (subject_id and ses_id)
  df.aibl <- df.aibl %>% 
    rename(subject_id = RID) %>% 
    mutate(subses_id = paste(subject_id, visit, sep = '_'))
  
  #filter out memory NAs
  df.memory <- df.aibl %>% 
    mutate(LIMMTOTAL = if_else(LIMMTOTAL<0, NaN, LIMMTOTAL), 
           LDELTOTAL = if_else(LDELTOTAL<0, NaN, LDELTOTAL)) %>% 
    filter(!if_all(LIMMTOTAL:LDELTOTAL, is.na)) %>%
    filter(!is.na(age_Cog))
  

  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  cat("aibl data extracted. saved as df.rda in /other_output/")
}

## for habs ## 
extract.longdata.habs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  df.habs = import(file.path(sitedir, "Demographics_HABS_DataRelease_2.0.csv")) %>% 
    rename(SubjIDshort = SubjID)
  df.grot1 =import(file.path(sitedir, "ClinicalMeasures_HABS_DataRelease_2.0.csv"))
  df.grot2 =import(file.path(sitedir, "Cognition_HABS_DataRelease_2.0.csv"))
  df.grot3 =import(file.path(sitedir, "PACC_HABS_DataRelease_2.0.csv"))
  #df.grot4 =import(file.path(sitedir, "ADNI_MRI_FS6_XSec_HABS_DataRelease_2.0.csv"))
  df.HABS <-
    list(df.habs, 
         df.grot1,
         df.grot2,
         df.grot3) %>% reduce(left_join, by=c("SubjIDshort", "StudyArc"))
  
  df.habs = 
    df.HABS %>% 
    rename(SubjID = SubjIDshort) %>% 
    group_by(SubjID) %>% 
    mutate(AgeBsl = min(NP_Age),
           SubjectRound_id = if_else(StudyArc == "HAB_1.0", 0,
                                     if_else(StudyArc == "HAB_2.0", 12,
                                             if_else(StudyArc == "HAB_3.0", 24,
                                                     if_else(StudyArc == "HAB_4.0", 36,
                                                             if_else(StudyArc == "HAB_5.0", 48, 60)))))) %>%
    #if_else(StudyArc == "HAB_5.0", 48, 60)))))) %>% 
    filter(HABS_DX == "CN")
  
  

  
  
  #rename subject_id and create subses_id (subject_id and ses_id)
  df.habs <- df.habs %>% 
    rename(subject_id = SubjID) %>% 
    mutate(subses_id = paste(subject_id, SubjectRound_id, sep = '_'))
  

  
  df.memory <- df.habs %>% 
    select(subject_id:HABS_DX, AgeBsl:SubjectRound_id, FCsrt_FNC:FCsrt_Free, LogicMem_IL:SRT_mc, 
           FCsrt_FNC:FCsrt_Free, LogicMem_IL:SRT_tr, LogicMem_DR_Z:FCsrt96_Z, subses_id) #%>% 
  #filter(!if_all(FCsrt_FNC:FCsrt96_Z, is.na))
  
  mem_vars <- 
    c("SRT_dr",
      "SRT_tr", 
      "LogicMem_DR", 
      "LogicMem_IL", 
      "FCsrt_Free")
  
  # Remove rows where all specified columns are NA
  df.memory <-  df.memory[!apply( df.memory[mem_vars], 1, function(row) all(is.na(row))), ]
  
  
  
  load("/ess/p274/cluster/projects/p039_image_brain/data_normative_long/df_mri/habs/df.rda")
  
  df.out <- df.out %>% 
    rename(subject_id = SubjID)
  
  mricheck <- df.out %>%
    group_by(subject_id) %>%
    mutate(last_mri = max(NP_Age)) %>%
    ungroup() %>%
    select(subject_id, StudyArc, NP_SessionDate, last_mri, NP_Age)
  
  
  merge <- left_join(df.memory, mricheck) %>%
    group_by(subject_id) %>%
    fill(last_mri, .direction = "down") %>%
    ungroup()
  
  exz <- merge %>% 
    group_by(subject_id) %>% 
    mutate(ex = if_else(NP_Age > last_mri, 1,0)) %>% 
    filter(ex == 1)
  

  
  #exclude memory data after 6 months after last MRI
  df.memory1 <- df.memory %>%
    filter(!subses_id %in% exz$subses_id)
  
  #include subjects tht have memory data but not mri)
  y <- anti_join(df.memory,df.out, by="subject_id")

  df.memory <- rbind(y,df.memory1)
    
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  cat("habs data extracted. saved as df.rda in /other_output/")
}

## for ous ##
extract.longdata.ous = function(sitedir, outdir, mrifiles) {
  
  load(file.path(sitedir, "cognNBM_updated.Rda")) 
  cognNBM_updated <- cognNBM_updated %>% 
    mutate(project_id = "CogNorm") %>% 
    rename(subject_id = ID) %>% 
    mutate(subses_id = paste(subject_id, timepoint, sep = '_'))
  
  
  #memory data
  df.memory <- cognNBM_updated %>% 
    filter(!(is.na(PCA_CERAD))) 
  
  
  ss <- df.memory %>% 
    mutate(mmse.crit = if_else(MMSE < 26 & !is.na(MMSE), 1,0)) %>% 
    group_by(subject_id) %>% 
    mutate(low.mmse = max(mmse.crit))
  
  df.grot <- ss %>% 
    filter(low.mmse == 0) %>% 
    select(-c(low.mmse, mmse.crit))
  
  aa = ss %>% filter(low.mmse == 1) %>% 
    group_by(subject_id) %>% 
    mutate(time.norm =if_else(mmse.crit == 0, Years, NaN), 
           time.norm = max(time.norm, na.rm = T)) %>%
    ungroup() %>% 
    mutate(crit.max.normal = if_else(Years <= time.norm & !is.infinite(time.norm), 1,0))  %>% 
    filter(crit.max.normal == 1) %>% 
    select(-c(low.mmse, mmse.crit, time.norm, crit.max.normal))
  
  
  df.memory = rbind(df.grot,aa) 

  
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  cat("ous data extracted. saved as df.rda in /other_output/")
  
}

## for preventAD ##
extract.longdata.preventad = function(sitedir, outdir, mrifiles) {
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
  
 
  df.preventAD <- df.preventAD %>% 
    rename(subject_id = CONP_CandID) %>% 
    mutate(subses_id = paste(subject_id, Study_visit_label.x, sep = '_')) %>% 
    mutate(project_id = "preventAD") %>% 
    select(-c(Study_visit_label.x.x,Study_visit_label.x.x.x,Study_visit_label.x.x.x.x,Study_visit_label.y,Study_visit_label.y.y,Study_visit_label.y.y.y,Study_visit_label.y.y.y.y))
  
  
  df.preventAD = 
    df.preventAD %>% 
    mutate(age = 
             if_else(!is.na(Candidate_Age.x), Candidate_Age.x, 
                     if_else(!is.na(Candidate_Age.y), Candidate_Age.y,
                             if_else(!is.na(Candidate_Age.x.x), Candidate_Age.x.x,
                                     if_else(!is.na(Candidate_Age.y.y), Candidate_Age.y.y,
                                             if_else(!is.na(Candidate_Age.x.x.x), Candidate_Age.x.x.x,
                                                     if_else(!is.na(Candidate_Age.y.y.y), Candidate_Age.y.y.y,
                                                             if_else(!is.na(Candidate_Age.x.x.x.x), Candidate_Age.x.x.x.x,
                                                                     if_else(!is.na(Candidate_Age.y.y.y.y), Candidate_Age.y.y.y.y,`Age At MRI In Months`)))))))),
           age = age/12)

  
  
  
  df.preventAD = 
    df.preventAD %>% 
    group_by(CONP_ID) %>% 
    mutate(age_at_mci = min(age_at_mci, na.rm = T)/12, 
           mci_visit = if_else(age_at_mci <= age, 1, 0)) %>% 
    filter(mci_visit == 0)
  
  
  mem_vars = c("list_learning_total_score", 
               "list_recall_score", 
               "story_memory_score", 
               "story_recall_score",
               "figure_recall_total_score")
  
  # Remove rows where all specified columns are NA
  df.memory <-  df.preventAD[!apply( df.preventAD[mem_vars], 1, function(row) all(is.na(row))), ]

  
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.long.rda"))
  
  
  
  cat("preventad data extracted. saved as df.rda in /other_output/")
}

## for ukb ##
extract.longdata.ukb = function(sitedir, outdir) {
  #load databases and get mri_filepaths
  #list.files(file.path(sitedir, "45249", "ukb_new"))
  dfdir=file.path(sitedir, "45249", "ukb_new")
  
  load(file.path(dfdir, "basic_demographics.Rda"))
  load(file.path(dfdir, "mri_data_tp23.Rda"))
  load(file.path(dfdir, "demographics_expanded.Rda"))
  load(file.path(dfdir, "cognition.Rda"))
  
  # df.ukb2 = df.dmg %>% filter(!is.na(age_3)) %>% 
  #   group_by(eid) %>% 
  #   select(-(c(age_0, age_1))) %>% 
  #   mutate(agebsl = min(c(age_2, age_3)))
  
  df.ukb = 
    left_join(df.dmg, 
              df.dmg.extra %>% 
                select(
                  eid,
                  uk_biobank_assessment_centre_f54_2_0,
                  uk_biobank_assessment_centre_f54_3_0
                ))
  
  df.ukb =
    left_join(df.ukb,
              df.cognition %>%
                select(
                  eid,
                  number_of_incorrect_matches_in_round_f399_0_1,
                  number_of_incorrect_matches_in_round_f399_0_2,
                  number_of_incorrect_matches_in_round_f399_0_3,
                  number_of_incorrect_matches_in_round_f399_1_1,
                  number_of_incorrect_matches_in_round_f399_1_2,
                  number_of_incorrect_matches_in_round_f399_1_3,
                  number_of_incorrect_matches_in_round_f399_2_1,
                  number_of_incorrect_matches_in_round_f399_2_2,
                  number_of_incorrect_matches_in_round_f399_2_3,
                  number_of_word_pairs_correctly_associated_f20197_2_0,
                  number_of_word_pairs_correctly_associated_f20197_3_0
                  
                ))
  
  
  df.ukb <- df.ukb %>% 
    rename(number_of_incorrect_matches_in_round_f399_round1_0_1 = number_of_incorrect_matches_in_round_f399_0_1) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round2_0_2 = number_of_incorrect_matches_in_round_f399_0_2) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round3_0_3 = number_of_incorrect_matches_in_round_f399_0_3) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round1_1_1 = number_of_incorrect_matches_in_round_f399_1_1) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round2_1_2 = number_of_incorrect_matches_in_round_f399_1_2) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round3_1_3 = number_of_incorrect_matches_in_round_f399_1_3) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round1_2_1 = number_of_incorrect_matches_in_round_f399_2_1) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round2_2_2 = number_of_incorrect_matches_in_round_f399_2_2) %>% 
    rename(number_of_incorrect_matches_in_round_f399_round3_2_3 = number_of_incorrect_matches_in_round_f399_2_3) 
  
  
  
  
  ##exclude those missing memory data on one or both tp
  #df.ukb <- df.ukb %>% 
  # filter(!if_any( prospective_memory_result_f20018_2_0:prospective_memory_result_f20018_3_0, is.na))
  
  
  
  grot.idx = grepl("age_",names(df.ukb))
  names(df.ukb)[grot.idx] = paste0(names(df.ukb)[grot.idx], "_0")
  
  df.ukb <- df.ukb %>%
    mutate(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0),
           t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0),
           uk_biobank_assessment_centre_f54_2_0 = as.numeric(uk_biobank_assessment_centre_f54_2_0),
           uk_biobank_assessment_centre_f54_3_0 = as.numeric(uk_biobank_assessment_centre_f54_3_0),)%>%
    pivot_longer(-c(eid, Sex),
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
  
  df.ukb <- df.ukb %>% 
    filter(!is.na(age))
  
  # df.ukb_long <- df.ukb %>% 
  #   group_by(eid) %>% 
  #   tally() %>% 
  #   filter(n>1)
  
  
  
  df.ukb <- df.ukb %>% 
    rename(subject_id = eid) %>%
    mutate(subses_id = paste(subject_id, tp, sep = '_'))%>% 
    mutate(project_id = "ukb")
  
  ## memory data ##
  #pairs matching task, round 2 
  df.memory.pm <- df.ukb %>% 
    filter(!is.na(number_of_incorrect_matches_in_round_f399_round2))
  
  #Paired associated learning (PAL) - word pairs
  df.memory <- df.ukb %>% 
    filter(!is.na(number_of_word_pairs_correctly_associated_f20197))
  
  
  
 
  
  # sub.mem <- df.memory %>% 
  #   group_by(subject_id) %>% 
  #   tally() %>% 
  #   filter(n>1)
  # #filter in all longitudinal data in df1
  # only.memory.long <- df.memory %>% 
  #   filter(subject_id %in% sub.mem$subject_id)
  
  #for PAL task
  # sub.mem2 <- df.memory_PAL %>% 
  #   group_by(subject_id) %>% 
  #   tally() %>% 
  #   filter(n>1)
  # #filter in all longitudinal data in df1
  # only.memory.long_PAL <- df.memory_PAL %>% 
  #   filter(subject_id %in% sub.mem2$subject_id)
  
  
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  cat("ukb data extracted. saved as df.rda in /other_output/")
}


## for oasis ## 
extract.longdata.oasis3 = function(sitedir, outdir, mrifiles) {
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
  
  df.oasis.cognorm= 
    df.oasis.cognorm %>%  
    mutate(OASISID = OASIS3_id, 
           cognorm = 1) %>% 
    select(OASISID, cognorm)
  
  
  
  df.oasis3 =
    reduce(list(
      df.oasis.psych, 
      #df.oasis.mri.info,
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
  
  
  df.oasis3 <- df.oasis3 %>% 
    rename(subject_id = OASISID) %>% 
    mutate(subses_id = paste(subject_id,days_to_visit, sep = '_')) %>% 
    mutate(project_id = "oasis3") %>%
    filter(crit.dx1.normal.max == 1)
  
  
  
  
  df.memory <- df.oasis3 %>%
    select(subject_id,
           project_id,
           days_to_visit,
           subses_id, 
           AGE, 
           GENDER, 
           MMSE,
           LOGIMEM,
           MEMUNITS,
           craftvrs,
           crafturs,
           craftdvr,
           craftdre,
           srtfree, 
           srttotal,
           asscmem) 
  
  # craft story has much less valud values than the remaing test. so it won'tbe integrated
  # other variables assessing memmory not included due to low n of valid values
  mem_vars <- 
    c("LOGIMEM",
      "MEMUNITS",
      #"craftvrs",
      #"crafturs",
      #"craftdvr",
      #"craftdre",
      "srtfree",
      "srttotal",
      "asscmem")
  
  # Remove rows where all specified columns are NA
  df.memory <-  df.memory[!apply( df.memory[mem_vars], 1, function(row) all(is.na(row))), ]
  
  
  # ## find longitudinal data in memory df (n>1) with longitudinal mri data
  # memory.long <- filter.long(df.memory,df.mri)
  # ## find longitudinal data in mri df (n>1) with longitudinal memory data
  # mri.long <- filter.long(df.mri, df.memory)
  # 
  
  # sub.mem <- df.memory %>% 
  #   group_by(subject_id) %>% 
  #   tally() %>% 
  #   filter(n>1)
  # #filter in all longitudinal data in df1
  # only.memory.long <- df.memory %>% 
  #   filter(subject_id %in% sub.mem$subject_id)
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  cat("oasis3 data extracted. saved as df.rda in /other_output/")
}



## for bbhi ## 
extract.longdata.bbhi = function(sitedir, outdir, mrifiles) {
  bbhi.data1 = read.csv(file.path(sitedir, "BBHI_NPS-2tps_sexAgeYoE.csv"))
  bbhi.data2 = read.csv(file.path(sitedir, "RAVLT_periode_aprenentatge_BBHI.csv"))
  
  
  df.bbhi1 = bbhi.data1 %>% 
    select(-c(contains("comment"))) %>%
    rename("w2_tmtb" = "w2_tmt_b_raw",
           "w1_tmtb" = "w1_tmt_b_raw") %>% 
    rename("w1_date" = "w1_nps_date") %>% 
    rename("w2_date" = "w2_nps_date") %>% 
    pivot_longer(-c(id_user, birth_date, sex, YoE,X), 
                 names_to = c("variable", ".value"),
                 names_sep = "_") %>% 
    #drop_na(any_of("nps")) %>% 
    mutate(sess = gsub("w","", variable) %>% as.numeric()) %>% 
    select(-X)
  
  #other variables
  df.bbhi2 <- bbhi.data2 %>% 
    pivot_longer(-c(id_user,X), 
                 names_to = "names", values_to = "values")
  df.bbhi2$sess <- substring(df.bbhi2$names, 1,2)
  df.bbhi2$v2 <- substring(df.bbhi2$names,4)

  df.bbhi2 <-df.bbhi2 %>%
    select(-names) %>% 
    pivot_wider(names_from = v2,
                values_from = values)
  
  
  df.bbhi2 <- df.bbhi2 %>% 
    mutate(sess = gsub("w","", sess) %>% as.numeric()) %>% 
    select(-X)
  #merge the dfs
  df.bbhi <- left_join(df.bbhi1,df.bbhi2)
  
  
  
  #reformat date and calculate age
  df.bbhi$date <- format(as.Date(df.bbhi$date, format = "%Y-%m-%d"), "%d/%m/%Y")
  df.bbhi$date <- as.Date(df.bbhi$date, format = "%d/%m/%Y" )
  
  df.bbhi$birth_date <- as.Date(df.bbhi$birth_date,"%d/%m/%Y")
  
  #filter out NA age
  df.bbhi <- df.bbhi %>% 
    filter(!is.na(date)) %>% 
    filter(!is.na(birth_date))
  
  #calculate age
  df.bbhi$age <- age_calc(df.bbhi$birth_date, df.bbhi$date, units = "years")
  
  #filter out mmse score >25
  # df.bbhi <- df.bbhi %>% 
  #   filter(mmse>25)
  
  mem_vars = c("inm",
               "delayed",
               "inm_recall_1_raw",
               "inm_recall_2_raw",
               "inm_recall_3_raw",
               "inm_recall_4_raw",
               "inm_recall_5_raw")
  
  df.memory <-  df.bbhi[!apply( df.bbhi[mem_vars], 1, function(row) all(is.na(row))), ]
  
  df.memory <- df.memory %>% 
    rename("subject_id" = "id_user") %>% 
    mutate(subses_id = paste(subject_id, sess, sep = '_')) 
   
  
  
  # create output
  save(df.memory,
       #only.memory.long,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  cat("bbhi data extracted. saved as df.rda in /other_output/")
}






extract.longdata.vetsa = function(sitedir, outdir, mrifiles) {
  db = import(file.path(sitedir, "Cognitive",  "admin_2023_07_07.sas7bdat"))
  db.cogn1 = import(file.path(sitedir, "Cognitive",  "cid_v1_2023_07_07.csv")) %>% 
    mutate(wave = 1) 
  db.cogn2 = import(file.path(sitedir, "Cognitive",  "cid_v2_2023_07_07.csv")) %>% 
    mutate(wave = 2)%>% 
    rename_with(~str_replace(., "_v2", ""), everything())
  db.cogn3 = import(file.path(sitedir, "Cognitive",  "cid_v3_2023_07_07.csv")) %>% 
    mutate(wave = 3)%>% 
    rename_with(~str_replace(., "_V3", ""), everything())
  
  df.cogn <- db.cogn1 %>% 
    full_join(db.cogn2) %>% 
    full_join(db.cogn3) %>% 
    select(CID, wave, everything())
  
  
  
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
  
  
  
  db <- left_join(db, df.cogn)
  df.inclusion = import(file.path(sitedir, "Cognitive", "MRI_readme", "VETSA123_healthMRIexclusions_cid_2023-06-28.csv"))
  df.inclusion1 = df.inclusion %>% 
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
  
  exc <- anti_join(df.inclusion, df.inclusion1)
  
  db = db %>% 
    filter(!CID %in% exc$CID)

  mem_vars = c("cvatot",
               "CVSDFR", 
               "CVLDFR")
  
  
  df.memory <-  db[!apply( db[mem_vars], 1, function(row) all(is.na(row))), ]
  df.memory <- df.memory %>% 
    rename(subject_id =  CID) %>% 
    mutate(subses_id = paste(subject_id, wave, sep = '_')) %>% 
    mutate(project_id = "vetsa")

  
  
  # create output
  save(df.memory,
       file = file.path(outdir,'other_output', "df.rda"))
  
  
  cat("vetsa data extracted. saved as df.rda in /other_output/")
}
