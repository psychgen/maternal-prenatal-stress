# Read in and prepare data for mat stress gxe analysis

library(foreign)
library(tidyverse)



## PRIMARY ANALYSIS VARIABLES

# Stress at work - these need reintegrating
work_stress15g <- c("AA1172","AA1177","AA1180", "AA1183")
work_stress30g <- c("CC923", "CC928")

# Select response vars - parent
# Relationship problems between parents
relprobs15g <- c("AA1533","AA1536","AA1538")
relprobs15g_r <- c("AA1532","AA1534","AA1535","AA1537","AA1539","AA1540","AA1541")
relprobs30g <- c("CC1193","CC1196","CC1198")
relprobs30g_r <- c("CC1192","CC1194","CC1195","CC1197","CC1199","CC1200","CC1201")
relprobs6 <- c("DD785","DD788","DD790")
relprobs6_r <- c("DD784","DD786","DD787","DD789","DD791","DD792","DD793")
relprobs18 <- c("EE611","EE614","EE616")
relprobs18_r <- c("EE610","EE612","EE613","EE615","EE617","EE618","EE619")
relprobs36 <- c("GG509")
relprobs36_r <- c("GG510","GG511","GG512","GG513" )
relprobs60 <- c("LL371")
relprobs60_r <- c("LL372","LL373","LL374","LL375")


relprobs = c(relprobs15g, relprobs30g,relprobs6, relprobs18, relprobs36, relprobs60,
             relprobs15g_r, relprobs30g_r,relprobs6_r, relprobs18_r, relprobs36_r, relprobs60_r)

# Life events
# each item first coded yes/no, then coded for severity, 1-3
# should be coded so that no=zero, yes=1 (not using severity scores)
LEyn30g<- c("CC1233","CC1235","CC1237","CC1239","CC1241","CC1243","CC1245","CC1247","CC1249") # yes/no responses
LEsev30g<- c("CC1234","CC1236","CC1238","CC1240","CC1242","CC1244","CC1246","CC1248","CC1250") # severity

LEyn6<- c("DD805","DD807","DD809","DD811","DD813","DD815","DD817","DD819","DD821","DD823","DD825") # yes/no responses
LEsev6<- c("DD806","DD808","DD810","DD812","DD814","DD816","DD818","DD820","DD821","DD824","DD826") # severity

LEyn18<- c("EE649","EE651","EE653","EE655","EE657","EE659","EE661","EE663","EE665","EE667","EE669") # yes/no responses
LEsev18<- c("EE650","EE652","EE654","EE656","EE658","EE660","EE662","EE664","EE666","EE668","EE670") # severity

LEyn36<- c("GG522","GG524","GG526","GG528","GG530","GG532","GG534","GG536","GG538","GG540") # yes/no responses
LEsev36<- c("GG523","GG525","GG527","GG529","GG531","GG533","GG535","GG537","GG539","GG541") # severity

LEn60<- c("LL388","LL391","LL394","LL397","LL400","LL403","LL406","LL409","LL412","LL415","LL418") # NO responses
LEy60<- c("LL389","LL392","LL395","LL398","LL401","LL404","LL407","LL410","LL413","LL416","LL419") # YES, during the last year


# Child outcomes

int <- c("EE908", "EE438", "EE439", "EE909", "EE963",
         "GG321", "GG335", "GG317", "GG328", "GG336","GG318", "GG323", "GG334", "GG337",
         "LL309", "LL305", "LL504", "LL317", "LL315","LL505", "LL321", "LL322", "LL310", "LL320", "LL323")

ext <- c("EE435", "EE903", "EE961", "EE446", "EE447", "EE962", "EE442","EE448",
         "GG314", "GG315", "GG330","GG332", "GG316", "GG319", "GG320", "GG324", "GG326", "GG329", "GG331",
         "LL303", "LL302", "LL324", "LL319", "LL304", "LL307", "LL308", "LL311","LL313", "LL316", "LL318")


q1 <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q1_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, relprobs15g,relprobs15g_r,work_stress15g)
q3 <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q3_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306,LEyn30g, relprobs30g,relprobs30g_r,work_stress30g)
q4 <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q4_6months_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, LEyn6,LEsev6, relprobs6,  relprobs6_r)
q5 <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5_18months_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, relprobs18,relprobs18_r, LEyn18, LEsev18,
         int[grepl("EE", int)], ext[grepl("EE", ext)])
q6 <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q6_3yrs_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, relprobs36,relprobs36_r, LEyn36, LEsev36,
         int[grepl("GG", int)], ext[grepl("GG", ext)])
q5yr <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_Q5yrs_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR,relprobs60_r,relprobs60, int[grepl("LL", int)], ext[grepl("LL", ext)],LEn60,LEy60)

# recode LE60 to a single binary yes/no variable

q5yr[LEn60]<-q5yr[LEn60]-1 
q5yr$LL388[1:20]
q5yr$LL389[1:20]

q5yr$LL389yn=ifelse(!is.na(q5yr$LL389),q5yr$LL389, q5yr$LL388)
q5yr$LL388[1:20]
q5yr$LL389[1:20]
q5yr$LL389yn[1:20]
q5yr$LL392yn=ifelse(!is.na(q5yr$LL392),q5yr$LL392, q5yr$LL391)
q5yr$LL395yn=ifelse(!is.na(q5yr$LL395),q5yr$LL395, q5yr$LL394)
q5yr$LL398yn=ifelse(!is.na(q5yr$LL398),q5yr$LL398, q5yr$LL397)
q5yr$LL401yn=ifelse(!is.na(q5yr$LL401),q5yr$LL401, q5yr$LL400)
q5yr$LL404yn=ifelse(!is.na(q5yr$LL404),q5yr$LL404, q5yr$LL403)
q5yr$LL407yn=ifelse(!is.na(q5yr$LL407),q5yr$LL407, q5yr$LL406)
q5yr$LL410yn=ifelse(!is.na(q5yr$LL410),q5yr$LL410, q5yr$LL409)
q5yr$LL413yn=ifelse(!is.na(q5yr$LL413),q5yr$LL413, q5yr$LL412)
q5yr$LL416yn=ifelse(!is.na(q5yr$LL416),q5yr$LL416, q5yr$LL415)
q5yr$LL419yn=ifelse(!is.na(q5yr$LL419),q5yr$LL419, q5yr$LL418)

LEyn60<- c("LL389yn","LL392yn","LL395yn","LL398yn","LL401yn","LL404yn","LL407yn","LL410yn","LL413yn","LL416yn","LL419yn") # YES/no, during the last year

LE = c(LEyn30g, LEyn6, LEyn18, LEyn36, LEyn60)


#

mbrn <- read.spss("//ess01/P471/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_MBRN_541_v12.sav", to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, KJONN, VEKT, SVLEN_DG, DODKAT_G, PARITET_5)



alldata <- mbrn %>% 
  left_join(q1) %>% 
  left_join(q3 ) %>%
  left_join(q4 ) %>%
  left_join(q5 ) %>%
  left_join(q6  ) %>%
  left_join(q5yr ) 

rm(q1, q3, q4, q5, q6, q5yr )

save(alldata, file="./scratch_data/alldata.RData")

###

load("./scratch_data/alldata.RData")

wsdata <- alldata %>%
  select( PREG_ID_2306,BARN_NR, work_stress15g, work_stress30g) %>%
  gather(item,val,-PREG_ID_2306, -BARN_NR) %>%
  mutate(true_val = case_when(val %in% c("Yes, every day more than half of the working day","Agree")~ 3,
                              val %in% c("Yes, every day less than half of the working day","Agree mostly")~ 2,
                              val %in% c("Yes, periodically but not daily","Disagree mostly")~ 1,
                              val %in% c("Seldom or never","Disagree completely")~ 0),
         age=str_sub(item, end=2),
         age=ifelse(age=="CC","AA",age)) %>%
  group_by(PREG_ID_2306, BARN_NR, age) %>%
  dplyr::summarize(items_present = sum(!is.na(val)),
            items_possible = n(),
            score = round(mean(true_val, na.rm=T)*n(),0)) %>%
  ungroup() %>%
  mutate(score = ifelse(items_present<(0.5*items_possible), NA, score),
         std_score = scale(score)) %>%
  select(-items_present, -items_possible) %>%
  gather(var, val, -PREG_ID_2306, -BARN_NR, -age) %>%
  mutate(age = case_when(age == "AA" ~ "wstressPRE"))%>%
  unite(temp, age, var) %>%
  spread(key = temp, value = val) 


rpdata <- alldata %>%
  select( PREG_ID_2306,BARN_NR, relprobs) %>%
  gather(item,val,-PREG_ID_2306, -BARN_NR) %>%
  mutate(true_val = case_when(val %in% c("Strongly agree") ~ 5,
                              val %in% c("Agree")~ 4,
                              val %in% c("Agree somewhat")~ 3,
                              val %in% c("Disagree somewhat")~ 2,
                              val %in% c("Disagree")~ 1,
                              val %in% c("Strongly disagree")~ 0),
         true_val = ifelse(item %in% c(relprobs15g_r, relprobs30g_r,relprobs6_r, relprobs18_r, relprobs36_r,
                                       relprobs60_r), (-1*true_val)+5, true_val ),
         age=str_sub(item, end=2)) %>%
  group_by(PREG_ID_2306, BARN_NR, age) %>%
  dplyr::summarize(items_present = sum(!is.na(val)),
            items_possible = n(),
            score = round(mean(true_val, na.rm=T)*n(),0)) %>%
  ungroup() %>%
  mutate(score = ifelse(items_present<(0.5*items_possible), NA, score)) %>% 
  group_by(age) %>%
  mutate(std_score = scale(score)) %>%
  ungroup() %>% 
  select(-items_present, -items_possible) %>%
  gather(var, val, -PREG_ID_2306, -BARN_NR, -age) %>%
  mutate(age = case_when(age == "AA" ~ "rprobsAA",
                         age == "CC" ~ "rprobsCC",
                         age == "DD" ~ "rprobs6mo",
                         age == "EE" ~ "rprobs18m",
                         age == "GG" ~ "rprobs3yr",
                         age == "LL" ~ "rprobs5yr"))%>%
  unite(temp, age, var) %>%
  spread(key = temp, value = val) %>% 
  mutate(rprobsPRE_score = ifelse(!is.na(rprobsCC_score), (rprobsAA_score+rprobsCC_score)/2, rprobsAA_score ),
         rprobsPRE_std_score = scale(rprobsPRE_score)) %>% 
  select(-matches("AA|CC"))

ledata <- alldata %>%
  select( PREG_ID_2306,BARN_NR, LE) %>%
  gather(item,val,-PREG_ID_2306, -BARN_NR) %>%
  mutate(age=str_sub(item, end=2),
         age=ifelse(age=="CC","AA",age),
         true_val = case_when(val %in% c("Yes") ~ 1,
                              val %in% c("No")~ 0,
                              val %in% c("1") ~ 1,
                              val %in% c("0")~ 0)) %>%
  filter(true_val %in% c(0,1)) %>%
  group_by(PREG_ID_2306, BARN_NR, age) %>%
  dplyr::summarize(items_present = sum(!is.na(val)),
            items_possible = n(),
            score = round(mean(true_val, na.rm=T)*n(),0)) %>%
  ungroup() %>%
  mutate(score = ifelse(items_present<(0.5*items_possible), NA, score),
         std_score = score) %>% #Not necessary to standardise this as a count variable
  select(-items_present, -items_possible) %>%
  gather(var, val, -PREG_ID_2306, -BARN_NR, -age) %>%
  mutate(age = case_when(age == "AA" ~ "levntsPRE",
                         age == "CC" ~ "levntsPRE",
                         age == "DD" ~ "levnts6mo",
                         age == "EE" ~ "levnts18m",
                         age == "GG" ~ "levnts3yr",
                         age == "LL" ~ "levnts5yr"))%>%
  unite(temp, age, var) %>%
  spread(key = temp, value = val) 
  
outdata <- alldata %>%
  select(PREG_ID_2306,BARN_NR, int, ext) %>%
  rename_at(vars(int), funs(paste0("int",.))) %>%
  rename_at(vars(ext), funs(paste0("ext",.))) %>%
  gather(item, val, -PREG_ID_2306,-BARN_NR) %>%
  mutate(age=str_sub(item, end=5),
         true_val = case_when(val %in% c("Not true","Rarely/never") ~ 0,
                              val %in% c("Sometimes","Somewhat or sometimes true")~ 1,
                              val %in% c("Often/ typical","Very true or often true") ~ 2)) %>%
  filter(true_val %in% c(0,1,2)) %>%
  group_by(PREG_ID_2306, BARN_NR, age) %>%
  dplyr::summarize(items_present = sum(!is.na(val)),
            items_possible = n(),
            score = round(mean(true_val, na.rm=T)*n(),0)) %>%
  ungroup() %>%
  mutate(score = ifelse(items_present<(0.5*items_possible), NA, score)) %>% 
  group_by(age) %>% 
  mutate(std_score = scale(score)) %>%
  ungroup() %>% 
  select(-items_present, -items_possible) %>%
  gather(var, val, -PREG_ID_2306, -BARN_NR, -age) %>%
  mutate(age = case_when(age == "intEE" ~ "cbclin18m",
                         age == "intGG" ~ "cbclin3yr",
                         age == "intLL" ~ "cbclin5yr",
                         age == "extEE" ~ "cbclex18m",
                         age == "extGG" ~ "cbclex3yr",
                         age == "extLL" ~ "cbclex5yr"))%>%
  unite(temp, age, var) %>%
  spread(key = temp, value = val) 
  

## JOIN, restricting to live births only

all <- mbrn %>% 
  left_join(wsdata) %>%
  left_join(rpdata) %>% 
  left_join(ledata) %>%
  left_join(outdata) %>% 
  filter(!DODKAT_G=="Stillborn (died before/during delivery or unknown time of death) or Abortion (requiring approval (§2.3c))")

# curate a phenotools dataset to get m_id and f_id for each child

ids <- phenotools::curate_dataset("asq_lang_c_18m", out_format = "merged_df")
ids <- ids %>% 
  select("PREG_ID_2306"=preg_id,BARN_NR,m_id,f_id) %>% 
  mutate(PREG_ID_2306=as.integer(PREG_ID_2306))

all <- all %>% 
  left_join(ids)

save(all, file= "./scratch_data/all_raw.RData")

############################################################

## ADD POLYGENIC SCORES

#library(genotools)

#Fetch

geno <- fetch_pgs(c("neurot2018","ptsd2019","adhd2022","height2"),
                  maf = "0.01",
                  clump = "250_1_0.1")  
#Process (regress on genotyping batch, imputation batch, and 20PCs)
geno_procd <- process_pgs(geno)

#Run PGS-PCA separately for kids/Mums/dads
geno_pcs_kids <-  geno_procd %>% 
  select(IID, preg_id, BARN_NR, matches("_res")) %>% 
  drop_na(preg_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd2022","height2")) %>% 
  select('PREG_ID_2306'= preg_id, BARN_NR, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".child"), matches("pgs.pc")) %>% 
  mutate(BARN_NR=as.numeric(BARN_NR))

geno_pcs_dads <-  geno_procd %>% 
  select(IID, f_id, matches("_res")) %>% 
  drop_na(f_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd2022","height2")) %>% 
  select( f_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".father"), matches("pgs.pc"))

geno_pcs_mums <-  geno_procd %>% 
  select(IID, m_id, matches("_res")) %>% 
  drop_na(m_id) %>% 
  pgs_pca(indid = "IID",
          pgs_var_stem = c("neurot2018","ptsd2019","adhd2022","height2")) %>% 
  select( m_id, matches("pgs.pc")) %>% 
  rename_with(~paste0(.x,".mother"), matches("pgs.pc")) 

geno_ids <- geno_procd %>% 
               select(IID,f_id,m_id,'PREG_ID_2306'=preg_id,BARN_NR)
dat_final <- all  %>% 
  left_join(geno_pcs_kids%>% 
              mutate(PREG_ID_2306=as.integer(PREG_ID_2306))) %>% 
  left_join(geno_pcs_dads) %>% 
  left_join(geno_pcs_mums) %>% 
  select(PREG_ID_2306,BARN_NR,m_id, KJONN, matches("std"),VEKT,SVLEN_DG, PARITET_5,matches("pgs")) %>% 
  mutate(indid= paste0(PREG_ID_2306,BARN_NR),
         VEKT= as.numeric(scale(VEKT)),
         SVLEN_DG= as.numeric(scale(SVLEN_DG)) ) %>%
  rename("bwt"=VEKT,"ges"=SVLEN_DG,"parity"=PARITET_5) %>% 
  mutate(sex= case_when(KJONN=="Male"~0,
                        KJONN=="Female"~1)) %>% 
  select(indid, m_id,
         matches("stress"),
         matches("probs"),
         matches("evnts"),
         bwt,ges,parity,sex,
         matches("cbclin"),
         matches("cbclex"),
         matches(".mother"),
         matches(".father"),
         matches(".child")) %>% 
  filter(!is.na(sex))

lkp<-tibble("old"=names(dat_final),
            "new"= c("indid","m_id",
                     "wsPRE","rp18m","rp3yr","rp5yr","rp6mo","rpPRE","le18m",
                     "le3yr","le5yr","le6mo","lePRE","bwt","ges","par","sex","emo1","emo2","emo3","beh1",
                     "beh2","beh3","neuM","ptsM","adhM","heiM","neuF","ptsF","adhF","heiF",
                     "neuC","ptsC","adhC","heiC"))

dat_final <- dat_final %>% 
  `colnames<-`(lkp$new) %>% 
  mutate(rpPRE=as.numeric(rpPRE))

save(dat_final, file = './scratch_data/prepped_00_vars.RData')

rawdat_final <- all  %>% 
  left_join(geno_pcs_kids%>% 
              mutate(PREG_ID_2306=as.integer(PREG_ID_2306))) %>% 
  left_join(geno_pcs_dads) %>% 
  left_join(geno_pcs_mums) %>% 
  select(PREG_ID_2306,BARN_NR,m_id, KJONN, matches("score"),VEKT,SVLEN_DG, PARITET_5,matches("pgs")) %>% 
  mutate(indid= paste0(PREG_ID_2306,BARN_NR),
         bwt_std= as.numeric(scale(VEKT)),
         ges_std= as.numeric(scale(SVLEN_DG)) ) %>%
  rename("bwt"=VEKT,"ges"=SVLEN_DG,"parity"=PARITET_5) %>% 
  mutate(sex= case_when(KJONN=="Male"~0,
                        KJONN=="Female"~1)) %>% 
  select(indid, m_id,
         matches("stress"),
         matches("probs"),
         matches("evnts"),
         bwt,ges,parity,sex,
         matches("cbclin"),
         matches("cbclex"),
         matches(".mother"),
         matches(".father"),
         matches(".child"))%>% 
  filter(!is.na(sex))

save(rawdat_final, file= "./scratch_data/all_raw.RData")


MplusAutomation::prepareMplusData(dat_final, filename = paste0("./data/all_for_mplus.dat"),
                 inpfile = TRUE)

sibids = dat_final %>% 
  filter(duplicated(m_id))

sibsmplusdat_final = dat_final %>% 
  filter(m_id %in% sibids$m_id)

#Load sib weighted data from 00.1 and join weights

load("./data/smoothed_weights.RData")

sibsmplusdat_final <- sibsmplusdat_final %>% 
  left_join(smoothed_weights) %>% 
  drop_na(ipsw)

message("Preparing data... \n")
MplusAutomation::prepareMplusData(sibsmplusdat_final, filename = paste0("./data/sibs_for_mplus.dat"),
                 inpfile = TRUE)


