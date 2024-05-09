#03_mr.R

#library(genotools)
library(tidyverse)

#Check SNPs/proxies in moba

geno<-readr::read_delim("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bim", col_names = c("CHR","rsid","a", "BP","a1","a2"))

snps<- c("rs9989237", "rs2736898", "rs11620763", "rs7146221")

geno %>% filter(rsid%in%snps)

#all but rs11620763 are present; find a proxy (the following code only runs outside of TSD:)

# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")

# LDLINK token must be set as an environment variable

# get_proxies("rs11620763", pop="CEU", results_dir="./", r2_threshold= 0.8)

# rs7141205 (r2 = 1) and rs11620777 (r2 = .967) are proxies - see if they are in MoBa:

snps<- c( "rs7141205", "rs11620777")

geno %>% filter(rsid%in%snps)

# Use rs7141205

#Use plink to code raw genetic data for canonical SNPs (see cortSNPs.sh)
#Read in raw genetic data for canonical SNPs

snps <- readr::read_delim("./data/canon.raw")


#Get effect allele frequencies

snps <- snps %>% 
  select(-contains("_HET"))
colnames(snps)

## Get allele frequencies for all SNPs
cols = c(7:10) #change to the columns of snps that contain SNPs

AA = lapply(cols, function(x) length(which(snps[,x] == 2)))
AA = data.frame(matrix(unlist(AA)))
AA$matrix.unlist.AA..= as.numeric(AA$matrix.unlist.AA..)
AB = lapply(cols, function(x) length(which(snps[,x] == 1)))
AB = data.frame(matrix(unlist(AB)))
AB$matrix.unlist.AB..= as.numeric(AB$matrix.unlist.AB..)
BB = lapply(cols, function(x) length(which(snps[,x] == 0)))
BB = data.frame(matrix(unlist(BB)))
BB$matrix.unlist.BB..= as.numeric(BB$matrix.unlist.BB..)
total = lapply(1:4,function(x) c(2*(AA[x,]+AB[x,]+BB[x,]))) #change to the length of AA/AB/BB
total = data.frame(matrix(unlist(total)))
total$matrix.unlist.total..= as.numeric(total$matrix.unlist.total..)

FRQ = cbind(AA,AB,BB,total)
colnames(FRQ) = c("AA","AB","BB","total")

FRQ["FRQ_A1"]  = (((2*FRQ$AA)+FRQ$AB)/FRQ$total)
FRQ = FRQ[,c("FRQ_A1")]

SNP <- colnames(snps[c(7:10)]) #should be the same as cols above
FRQ1 <- data.frame(SNP, FRQ)
FRQ1

rm(AA,AB,BB,total)

# Fetch and process a PGS using genotools to get an id list and covariates

geno <- fetch_pgs(c("neurot2018"),
                  maf = "0.01",
                  clump = "250_1_0.1") %>% 
  process_pgs() 

#get unrelated trios for analyses (using unrelate from the development version of genotools)

unrelated <- unrelate()

#Join snps and geno, restricting to IDs in unrelated

snpdat<- snps %>% 
  left_join(geno %>% 
              select(-matches("neurot"))) %>% 
  filter(IID %in% unrelated$IID) %>% 
  group_by(FID) %>% 
  fill(c(preg_id,BARN_NR,f_id,m_id),.direction="downup") %>% 
  ungroup()

#Read in and join to analytic dataset ahead of performing snp models

load('./scratch_data/prepped_00_vars.RData')

#Reshape snpdat to wide with mother / father / child SNP cols

snpdatwide = snpdat  %>% 
  pivot_longer(cols=matches("rs"), names_to="snp") %>% 
  mutate("snp_role"= paste(snp, Role, sep="_") )%>% 
  select(-snp) %>% 
  pivot_wider(names_from = snp_role, values_from= value ) %>%
  group_by(preg_id) %>% 
  fill(matches("rs"),.direction="downup") %>% 
  ungroup() %>% 
  filter(Role=="Mother") # The PC values that are retained are from the mothers


snpmodelsdat <- snpdatwide %>% 
  mutate(indid=paste0(preg_id,BARN_NR)) %>% 
  select(-m_id)%>% 
  left_join(dat_final) 

save(snpmodelsdat, file="./scratch_data/snpmodelsdat.RData")

for (out in c("bwt","ges","emo1","emo2","emo3","beh1","beh2","beh3")){
  
  snpmodelsdat$yvar = unlist(snpmodelsdat[,out])
  
  triogenodat = snpmodelsdat %>% 
    filter(!is.na(`rs7146221_A(/G)_Mother`),
           !is.na(`rs7146221_A(/G)_Father`),
           !is.na(`rs7146221_A(/G)_Child`))
  
  
  covars <- paste0(" + ", paste0("PC",seq(1,20),collapse=" + ")," + genotyping_batch ")
  
  trio1 <- lm(paste0('yvar ~ `rs7146221_A(/G)_Mother` + `rs7146221_A(/G)_Father`+ `rs7146221_A(/G)_Child`',covars) , data=triogenodat) %>%  summary()
  mum1 <- lm(paste0('yvar ~ `rs7146221_A(/G)_Mother` ',covars) , data=triogenodat) %>%  summary()
  
  trio2 <- lm(paste0('yvar ~ `rs9989237_T(/C)_Mother` + `rs9989237_T(/C)_Father`+ `rs9989237_T(/C)_Child`',covars), data=triogenodat) %>%  summary()
  mum2 <- lm(paste0('yvar ~ `rs9989237_T(/C)_Mother`',covars), data=triogenodat) %>%  summary()
  
  trio3 <- lm(paste0('yvar ~ `rs2736898_T(/C)_Mother` + `rs2736898_T(/C)_Father`+ `rs2736898_T(/C)_Child`',covars), data=triogenodat) %>%  summary()
  mum3 <- lm(paste0('yvar ~ `rs2736898_T(/C)_Mother` ',covars), data=triogenodat) %>%  summary()
  
  trio4 <- lm(paste0('yvar ~ `rs7141205_G(/A)_Mother` + `rs7141205_G(/A)_Father`+ `rs7141205_G(/A)_Child`',covars), data=triogenodat) %>%  summary()
  mum4 <- lm(paste0('yvar ~ `rs7141205_G(/A)_Mother` ',covars), data=triogenodat) %>%  summary()
  
  modlist <- list(trio1,trio2,trio3,trio4,mum1,mum2,mum3,mum4)
  
  tmpres <- map(modlist, function(x){
    broom::tidy(x) %>% 
      filter(str_detect(term,"rs")) %>% 
      cbind(broom::glance(x) %>% 
              select('N'=nobs)) %>% 
      mutate(across(c(term), ~ str_remove_all(.x,"`"))) %>% 
      separate(term, into=c("SNP","Alleles","Role"), sep="_") %>%
      left_join(FRQ1 %>% 
                  separate(SNP, into=c("SNP","Alleles"), sep="_")) %>% 
      select(SNP,Alleles,Role, "beta"=estimate, "SE"=std.error,"pval"=p.value, N, FRQ) 
  }) %>% 
    purrr::reduce(bind_rows) %>% 
    mutate(effect_allele=str_sub(Alleles,1,1),
           other_allele=str_sub(Alleles,4,4)) %>% 
    select(-Alleles)
  
  write.csv(tmpres, paste0("./output/sstats/",out,"_all.csv"), quote=F, row.names=F)
  
  #maternal effects only
  
  matres = tmpres %>% 
    filter(Role == "Mother") %>% 
    mutate(Model = rep(c("Trio","Single"),each=4))
  
  write.csv(matres %>% filter(Model=="Trio"), paste0("./output/sstats/",out,"_trio.csv"), quote=F, row.names=F)
  write.csv(matres %>% filter(Model=="Single"), paste0("./output/sstats/",out,"_single.csv"), quote=F, row.names=F)
}



