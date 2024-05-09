#03_2smr.R
##Adapted from script by Robyn Wootton 

# In this script, which should be run offline (i.e., outside of the TSD environment), we get
# Wald ratios for the single SNP associations and perform IVW meta-analyses for the MR component
# accounting for the correlations between SNPs

library(TwoSampleMR)
library(MendelianRandomization)
library(tidyverse)

##############################################################################################
#Single SNP analysis for canonical SNPs
##############################################################################################
#Read in the exposure effects for each of the canonical SNPs

canon <-readr::read_delim("./gwama_maf0.01_n10k.txt",delim="\t")

canon_hits_only <- canon %>% 
  filter(SNP %in% c("rs9989237", "rs2736898",  "rs7146221","rs7141205"))

write_tsv(canon_hits_only, file="./sstats/cort_canon.tsv")

cort_canon <- read_exposure_data("./sstats/cort_canon.tsv", sep="\t", 
                                 snp_col = "SNP",
                                 beta_col = "Effect",
                                 se_col = "StdErr",
                                 eaf_col = "Freq1",
                                 effect_allele_col = "Allele1",
                                 other_allele_col = "Allele2")
head(cort_canon)



trio_data_list <- list.files("./sstats")[ grep(pattern="trio.csv" ,list.files("./sstats"))] 
single_data_list <- list.files("./sstats")[ grep(pattern="single.csv" ,list.files("./sstats"))] 

##Canonical - cort
res_ssnp_trio <- c()
res_ssnp_single <- c()

for(i in 1:length(single_data_list)){
  tname <- trio_data_list[i]
  sname <- single_data_list[i]
  tout_dat = read_outcome_data(snps=cort_canon$SNP, 
                              filename = paste0("./sstats/",trio_data_list[i]),
                              sep =",",
                              snp_col="SNP",
                              beta_col="beta",
                              se_col="SE",
                              effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele",
                              eaf_col = "FRQ",
                              pval_col = "pval",
                              samplesize_col = "N")
  sout_dat = read_outcome_data(snps=cort_canon$SNP, 
                               filename = paste0("./sstats/",single_data_list[i]),
                               sep =",",
                               snp_col="SNP",
                               beta_col="beta",
                               se_col="SE",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               eaf_col = "FRQ",
                               pval_col = "pval",
                               samplesize_col = "N")
  tout_dat["outcome"] = tname
  tharm_dat = harmonise_data(cort_canon, tout_dat)
  #tharm_dat = subset(tharm_dat, mr_keep == "TRUE")
  tharm_dat$exposure <- "cort Canonical SNPs"
  tmr_ssnp = mr_singlesnp(tharm_dat, all_method = c("mr_ivw"))
  
  ##replace IVW estimate with one from MendelianRandomization package that accounts for inter-SNP correlations
  
  tIVWcordat <- dat_to_MRInput(tharm_dat, get_correlations = T)
  tIVWcor <-MendelianRandomization::mr_ivw(tIVWcordat[[1]], correl = TRUE)
  tmr_ssnp <- tmr_ssnp %>% 
    mutate(b=case_when( SNP == "All - Inverse variance weighted"~ tIVWcor$Estimate,
                        TRUE ~ b),
           se=case_when( SNP == "All - Inverse variance weighted"~ tIVWcor$StdError,
                         TRUE ~ se),
           p=case_when( SNP == "All - Inverse variance weighted"~ tIVWcor$Pvalue,
                        TRUE ~ p))

  
  tp2 <- mr_forest_plot(tmr_ssnp)
  filename <- paste0("./plots/cort_singleSNP", tname, ".png")
  ggsave(tp2[[1]], file=filename, width=7, height=5)
  
  sout_dat["outcome"] = sname
  sharm_dat = harmonise_data(cort_canon, sout_dat)
  #sharm_dat = subset(sharm_dat, mr_keep == "TRUE")
  sharm_dat$exposure <- "cort Canonical SNPs"
  smr_ssnp = mr_singlesnp(sharm_dat, all_method = c("mr_ivw"))
  
  ##replace IVW estimate with one from MendelianRandomization package that accounts for inter-SNP correlations
  
  sIVWcordat <- dat_to_MRInput(sharm_dat, get_correlations = T)
  sIVWcor <-MendelianRandomization::mr_ivw(sIVWcordat[[1]], correl = TRUE)
  smr_ssnp <- smr_ssnp %>% 
    mutate(b=case_when( SNP == "All - Inverse variance weighted"~ sIVWcor$Estimate,
                        TRUE ~ b),
           se=case_when( SNP == "All - Inverse variance weighted"~ sIVWcor$StdError,
                         TRUE ~ se),
           p=case_when( SNP == "All - Inverse variance weighted"~ sIVWcor$Pvalue,
                        TRUE ~ p))
  
  sp2 <- mr_forest_plot(smr_ssnp)
  filename <- paste0("./plots/cort_singleSNP", sname, ".png")
  ggsave(sp2[[1]], file=filename, width=7, height=5)
  
  
  #save all of the effects
  
  res_ssnp_trio <- rbind(res_ssnp_trio, tmr_ssnp)
  res_ssnp_single <- rbind(res_ssnp_single, smr_ssnp)
}

#Save
write.csv(res_ssnp_trio, "./output/mum_trio_2024.csv", row.names=FALSE, quote=FALSE)
write.csv(res_ssnp_single, "./output/mum_single_2024.csv", row.names=FALSE, quote=FALSE)


#Plot trio
res_trio<- readr::read_csv("./output/mum_trio_2024.csv")

plot_res <- res_trio %>% 
  mutate(outcome = factor(str_remove_all(outcome,"_trio.csv"),
                          levels=rev(c("bwt","ges","beh1","beh2","beh3","emo1","emo2","emo3")),
                          labels=rev(c("Birthweight","Gestational\nage","Behavioral\n(1.5yrs)",
                                       "Behavioral\n(3yrs)","Behavioral\n(5yrs)","Emotional\n(1.5yrs)",
                                       "Emotional\n(3yrs)","Emotional\n(5yrs)"))),
         est_type = ifelse(str_detect(SNP,"All"),"IVW","SNP"),
         SNP = str_replace_all(SNP,"All - Inverse variance weighted","IVW"))

ggplot(plot_res, aes(x=b,y=SNP, colour=est_type, shape= est_type))+
  geom_vline(xintercept=0,linetype=2,colour="grey80",size=1.2)+
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0)+
  geom_point(size=2)+
  scale_colour_manual("Estimate",values=c("red","black"))+
  scale_shape_manual("Estimate",values=c(23,16))+
  theme(
    text=element_text(size=12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=14, margin= margin(t=20,r=0,b=0,l=0)),
    panel.background = element_rect(fill="grey98", colour="white"),
    panel.border= element_rect(colour="black", fill=NA),
    panel.grid.major = element_line(colour="white"),
    strip.text.y = element_text(angle=0),
    strip.background = element_blank(),
    strip.text = element_text(face="bold"),
    legend.position = "bottom",
    legend.direction="horizontal"
  )+
  facet_grid(outcome~.)  +
  scale_x_continuous("Two-Sample MR\n effect estimate")

ggsave("./plots/mr_trio.tiff",device="tiff",width=9, height= 17, unit="cm", dpi=420,bg="white")

# PLot single

res_single<- readr::read_csv("./output/mum_single_2024.csv")

plot_res <- res_single %>% 
  mutate(outcome = factor(str_remove_all(outcome,"_single.csv"),
                          levels=rev(c("bwt","ges","beh1","beh2","beh3","emo1","emo2","emo3")),
                          labels=rev(c("Birthweight","Gestational\nage","Behavioral\n(1.5yrs)",
                                       "Behavioral\n(3yrs)","Behavioral\n(5yrs)","Emotional\n(1.5yrs)",
                                       "Emotional\n(3yrs)","Emotional\n(5yrs)"))),
         est_type = ifelse(str_detect(SNP,"All"),"IVW","SNP"),
         SNP = str_replace_all(SNP,"All - Inverse variance weighted","IVW"))

ggplot(plot_res, aes(x=b,y=SNP, colour=est_type, shape= est_type))+
  geom_vline(xintercept=0,linetype=2,colour="grey80",size=1.2)+
  geom_errorbarh(aes(xmin=b-1.96*se, xmax=b+1.96*se), height=0)+
  geom_point(size=2)+
  scale_colour_manual("Estimate",values=c("red","black"))+
  scale_shape_manual("Estimate",values=c(23,16))+
  theme(
    text=element_text(size=12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=14, margin= margin(t=20,r=0,b=0,l=0)),
    panel.background = element_rect(fill="grey98", colour="white"),
    panel.border= element_rect(colour="black", fill=NA),
    panel.grid.major = element_line(colour="white"),
    strip.text.y = element_text(angle=0),
    strip.background = element_blank(),
    strip.text = element_text(face="bold"),
    legend.position = "bottom",
    legend.direction="horizontal"
  )+
  facet_grid(outcome~.)  +
  scale_x_continuous("Two-Sample MR\n effect estimate")

ggsave("./plots/mr_single.tiff",device="tiff",width=9, height= 17, unit="cm", dpi=420,bg="white")


