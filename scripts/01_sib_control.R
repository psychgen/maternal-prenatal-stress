#01_sib_control.R

# Load packages

library(tidyverse)
library(MplusAutomation)
library(tictoc)

# Run models

filepath1 <- "./scripts/mplus/sib_control"

tic()
runModels(filepath1, recursive =F,replaceOutfile="modifiedDate", Mplus_command = "C:/Program Files/Mplus/Mplus" )
toc()

# Move the results files to the output file and clean up the scripts folder
file.copy(from=paste0(filepath1,"/",list.files(filepath1)[str_detect(list.files(filepath1),".inp",negate=T)]),
          to="./output/mplus/sib_control",
          overwrite = TRUE, recursive = F,
          copy.mode = TRUE)

junk <- dir(path=filepath1, pattern=".out|.dat") 
file.remove(paste0(filepath1,"/",junk))

# Read outputs
# read in output
mplusOutput <- readModels("./output/mplus/sib_control", recursive=FALSE,
                          what = c("input", "warn_err", "data_summary", "sampstat", "covariance_coverage", "summaries", "parameters", "class_counts", "indirect", "mod_indices", "residuals",  "tech3", "tech4", "tech7", "tech8", "tech9", "tech10", "tech12", "fac_score_stats", "lcCondMeans", "gh5", "output"))


# Check whether emo/beh models constrained across time are acceptable 

#Calculate test statistics

beh1l_2LLdiff <- -2 * (as.numeric(mplusOutput$beh_1l_cons.out$summaries$LL)-as.numeric(mplusOutput$beh_1l.out$summaries$LL))
emo1l_2LLdiff <- -2 * (as.numeric(mplusOutput$emo_1l_cons.out$summaries$LL)-as.numeric(mplusOutput$emo_1l.out$summaries$LL))
beh2l_2LLdiff <- -2 * (as.numeric(mplusOutput$beh_2l_cons.out$summaries$LL)-as.numeric(mplusOutput$beh_2l.out$summaries$LL))
emo2l_2LLdiff <- -2 * (as.numeric(mplusOutput$emo_2l_cons.out$summaries$LL)-as.numeric(mplusOutput$emo_2l.out$summaries$LL))

#Calcuated df difference (same for both)
dfdiff1l <- mplusOutput$emo_1l.out$summaries$Parameters - mplusOutput$emo_1l_cons.out$summaries$Parameters
dfdiff2l <- mplusOutput$emo_2l.out$summaries$Parameters - mplusOutput$emo_2l_cons.out$summaries$Parameters
#Run test

if(pchisq(emo1l_2LLdiff, df = dfdiff1l, lower.tail = FALSE)<0.05/4){
  emo_1l_mod=mplusOutput$emo_1l.out
}else{
  emo_1l_mod=mplusOutput$emo_1l_cons.out
}
if(pchisq(emo2l_2LLdiff, df = dfdiff2l, lower.tail = FALSE)<0.05/4){
emo_2l_mod=mplusOutput$emo_2l.out
}else{
emo_2l_mod=mplusOutput$emo_2l_cons.out
}

if(pchisq(beh1l_2LLdiff, df = dfdiff1l, lower.tail = FALSE)<0.05/4){
  beh_1l_mod=mplusOutput$beh_1l.out
}else{
  beh_1l_mod=mplusOutput$beh_1l_cons.out
}
if(pchisq(beh2l_2LLdiff, df = dfdiff2l, lower.tail = FALSE)<0.05/4){
  beh_2l_mod=mplusOutput$beh_2l.out
}else{
  beh_2l_mod=mplusOutput$beh_2l_cons.out
}

# Summarise results

bwt_1l_mod=mplusOutput$bwt_1l.out
bwt_2l_mod=mplusOutput$bwt_2l.out
ges_1l_mod=mplusOutput$ges_1l.out
ges_2l_mod=mplusOutput$ges_2l.out
bwt_full_mod=mplusOutput$bwt_full.out
ges_full_mod=mplusOutput$ges_full.out
emo_full_mod=mplusOutput$emo_full.out
beh_full_mod=mplusOutput$beh_full.out

mods <- list(bwt_1l_mod, bwt_2l_mod,
             ges_1l_mod, ges_2l_mod,
             emo_1l_mod,emo_2l_mod,
             beh_1l_mod, beh_2l_mod,
             bwt_full_mod, ges_full_mod,
             emo_full_mod, beh_full_mod)

results<- purrr::map(mods, function(x){
   cis <- x$parameters$ci.stdyx.standardized %>% 
    as_tibble() %>% 
    filter(param %in% c("WSPRE","RPPRE","LEPRE")) %>% 
    mutate(model=x$input$title[1]) %>% 
     filter(!str_detect(paramHeader, "WITH|Variances"),
            BetweenWithin=="Within") %>% 
     select(model,paramHeader,param,est,"lci"=`low2.5`, "uci"=`up2.5`)
   ps <- x$parameters$stdyx.standardized %>% 
     as_tibble() %>% 
     filter(param %in% c("WSPRE","RPPRE","LEPRE")) %>% 
     mutate(model=x$input$title[1])%>% 
     filter(!str_detect(paramHeader, "WITH|Variances"),
            BetweenWithin=="Within") %>% 
     select(model,paramHeader,param,est,se,pval,model)
   left_join(ps,cis)
})

# Save out

save(results, file="./output/obs_sib_mods_res.RData")


# Plot

plotres <- results %>% 
  purrr::reduce(bind_rows) %>% 
  separate(model, into = c("outcome","model")) %>% 
  mutate(exposure= factor(tolower(param),levels=c("lepre","rppre","wspre"),
                          labels=c("Life events","Relationship problems","Work stress")),
         outcome_specific=factor(tolower(str_remove_all(paramHeader,".ON")), 
                                 levels=rev(c("bwt","ges","beh1","beh2","beh3","emo1","emo2","emo3")),
                                 labels=rev(c("","","1.5yrs","3yrs","5yrs","1.5yrs","3yrs","5yrs"))),
         outcome= factor(outcome,levels=rev(c("bwt","ges","beh","emo")),
                         labels=rev(c("Birthweight","Gestational age","Behavioral","Emotional"))),
         Wave= factor(ifelse(outcome%in%c("Birthweight","Gestational age"),"Birth",outcome_specific),
                      levels=c("Birth","1","2","3"),
                      labels=c("Birth", "1.5yrs","3yrs","5yrs")))
  

pd=0.4
ggplot(plotres %>% 
         filter(model %in% c("1l","2l")), aes(x=est, y=outcome_specific, fill= Wave,colour=Wave, shape =model))+
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+
  geom_errorbarh(aes(xmin=lci,xmax=uci),position=position_dodge(pd), size=1.2, height=0, alpha=0.4) +
  geom_point(position = position_dodge(pd),colour= "grey40", size=3.4, stroke=1.2)+  
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_shape_manual(values=c( 21, 24), name = "Familial confounders...", breaks=c("1l","2l"), labels=c("Undjusted", "Adjusted")   )+
  theme(text = element_text(size=12),
        axis.title.y = element_text(size=14,margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=14,margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "grey98", colour = "white"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0)) +
  facet_grid(outcome~exposure, scales= "free_y",  space="free", switch ="y")+
  scale_y_discrete("Outcome")+
  scale_x_continuous("Standardised estimate of effect in multilevel SEM")+
  guides(fill=guide_legend(override.aes=list(shape=21)))


ggsave("./output/plots/sib_control.tiff", device="tiff", width=30,height=14,units="cm",dpi=320,bg="white")


#Supp plot including obs estimates

ggplot(plotres %>% 
         filter(model %in% c("1l","full")), aes(x=est, y=outcome_specific, fill= Wave,colour=Wave, shape =model))+
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+
  geom_errorbarh(aes(xmin=lci,xmax=uci),position=position_dodge(pd), size=1.2, height=0, alpha=0.4) +
  geom_point(position = position_dodge(pd),colour= "grey40", size=3.4, stroke=1.2)+  
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  scale_shape_manual(values=c( 21, 22), name = "Familial confounders...", breaks=c("1l", "full"), labels=c("Undjusted (sibling sample)" ,"Unadjusted (full sample)" )   )+
  theme(text = element_text(size=12),
        axis.title.y = element_text(size=14,margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=14,margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "grey98", colour = "white"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0)) +
  facet_grid(outcome~exposure, scales= "free_y",  space="free", switch ="y")+
  scale_y_discrete("Outcome")+
  scale_x_continuous("Standardised estimate of effect in multilevel SEM")+
  guides(fill=guide_legend(override.aes=list(shape=21)))


ggsave("./output/plots/SUPP_obs_sib_control.tiff", device="tiff", width=30,height=14,units="cm",dpi=320,bg="white")



