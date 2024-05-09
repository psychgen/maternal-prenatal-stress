#04_negative_ctrl.R

library(tidyverse)
library(lavaan)

#Load data

load('./scratch_data/prepped_00_vars.RData')

dat = dat_final %>% 
  mutate(par=as.integer(par))
#Source functions

source("./scripts/04.1_neg_ctrl_funs.R")


# Create a data frame with inputs and to hold results

linear_mods <- expand.grid(c("bwt","ges","emo1","emo2","emo3","beh1","beh2","beh3"),
                           c("PRE","18m","3yr","5yr")) %>% 
  `colnames<-`(c("out","expwave")) %>% 
  as_tibble()

# Run the models across each combination of outcome and wave

linear_mods_res <- linear_mods %>% 
  mutate(results = pmap(., run_test_models, data=dat))

head(linear_mods_res)

# Save the result

save(linear_mods_res, file="./output/neg_ctrl_linear_mods_res.RData")
load(file="./output/neg_ctrl_linear_mods_res.RData")
# Plot the results

allests <- do.call(rbind, linear_mods_res$results) %>% 
  mutate(model = ifelse(out=="ges"&model=="invalid","validnc",model )) %>% 
  filter(model!="invalid",
         !str_detect(out,"3")) %>% 
  mutate(exposure=factor(str_sub(exp,end=2),
                         levels=c("rp","le"),
                         labels=c("Relationship problems", "Life events")),
         outcome=factor(str_sub(out,end=3),
                        levels=rev(c("bwt","ges","emo","beh")),
                        labels=rev(c("Birthweight","Gest.age","Behavioral","Emotional"))),
         `Analysis type`=factor(model, levels=c("expout","validnc"),
                                labels=c("Exposure-outcome","Negative control"))) %>%
  group_by(exposure, outcome, `Analysis type`) %>% 
  summarise(est = mean(est.std), se=mean(se))
           
                          
                        

pd=0.4
ggplot(allests,aes(x=outcome,y=est, colour=`Analysis type`,shape=`Analysis type`))+
  geom_hline(yintercept=0, linetype=2, colour="grey80", size=1.2)+
  geom_errorbar(aes(ymin=est-1.96*se,ymax=est+1.96*se),position=position_dodge(pd), size=1.2, width=0, alpha=0.4) +
  geom_point(position = position_dodge(pd), size=3.4, stroke=1.2)+  
  scale_colour_brewer(palette="Dark2")+
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
  scale_y_continuous("Average effect in full sample exposure-outcome\n and negative control analyses ")+
  coord_flip()+
  facet_grid(.~exposure, scales="free", space="free")



ggsave("./output/plots/neg_ctrl.tiff", device="tiff", width=16,height=11,units="cm",dpi=320,bg="white")

