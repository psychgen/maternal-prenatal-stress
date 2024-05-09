#02_pgs_gxe.R

#Load packages

library(tidyverse)
library(lavaan)

#Load data

load('./scratch_data/prepped_00_vars.RData')

dat = dat_final %>% 
  mutate(par=as.integer(par))

#Source functions

source("./scripts/02.1_pgs_gxe_funs.R")


# Create a data frame with inputs and to hold results

linear_mods <- expand.grid(c("bwt","ges","emo","beh"),
                           c("neu","pts","adh","hei")) %>% 
  `colnames<-`(c("out","mod")) %>% 
  as_tibble()

# Run the models across each combination of outcome and moderator, saving results in nested dataframe

linear_mods_res <- linear_mods %>% 
  mutate(results = pmap(., run_test_models, data=dat),
         ests = map(results, function(x) x[["ests"]]),)

head(linear_mods_res)

# Save the result

save(linear_mods_res, file="./output/gxe_linear_mods_res.RData")

# Plot the results

load( file="./output/gxe_linear_mods_res.RData")

allests <- do.call(rbind, linear_mods_res$ests) %>% 
  filter(str_detect(label, "b7|b8|b9")) %>% 
  mutate(exposure= factor(str_sub(tolower(rhs),end=5),levels=c("lepre","rppre","wspre"),
                          labels=c("Life evts.","R'ship probs.","Work stress")),
         outcome_specific=factor(tolower(lhs), 
                                 levels=rev(c("bwt","ges","beh1","beh2","beh3","emo1","emo2","emo3")),
                                 labels=rev(c("","","1.5yrs","3yrs","5yrs","1.5yrs","3yrs","5yrs"))),
         outcome= factor(str_sub(lhs, end=3),
                         levels=rev(c("bwt","ges","beh","emo")),
                         labels=rev(c("Birthweight","Gest.age","Behavioural","Emotional"))),
         Wave= factor(ifelse(outcome%in%c("Birthweight","Gest.age"),"Birth",outcome_specific),
                      levels=c("Birth","1","2","3"),
                      labels=c("Birth", "1.5yrs","3yrs","5yrs")),
         `PGS moderator`=factor(mod, levels=c("neu","pts","adh","hei"),
                          labels=c("Neurot.","PTSD","ADHD","Height")),
         fdrp= p.adjust(pvalue, method="fdr")) %>% 
  group_by(outcome,`PGS moderator`) %>% 
  summarise(est= mean(est.std),
            se= mean(se))
pd=0.8
ggplot(allests,aes(x=outcome,y=est, colour=`PGS moderator`))+
  geom_hline(yintercept=0, linetype=2, colour="grey80", size=1.2)+
  geom_errorbar(aes(ymin=est-1.96*se,ymax=est+1.96*se),position=position_dodge(pd), size=1.2, width=0, alpha=0.4) +
  geom_point(position = position_dodge(pd), size=3.4, stroke=1.2)+  
  scale_colour_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
 theme(text = element_text(size=12),
        axis.title.y = element_text(size=14,margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=14,margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.text.y  = element_blank(),
        #axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "grey98", colour = "white"),
        panel.border = element_rect(colour="black", fill=NA),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="white"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank()) +
 scale_x_discrete("Outcome")+
  coord_cartesian(ylim=c(-0.05,0.05))+
  scale_y_continuous("Average interaction effect between\n PGS moderator and exposures ")


ggsave("./output/plots/gxe.tiff", device="tiff", width=14,height=12,units="cm",dpi=320,bg="white")


