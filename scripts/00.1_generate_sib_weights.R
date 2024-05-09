# 00.1_generate_sib_weights.R

# Within-level estimates in MLMs rely substantially on families with more than one
# participating child; for estimates pertaining to fathers' variables, the fathers
# also need to have participated multiple times - this is likely a selected sample
# In this script, we explore this issue and seek adjustments

library(tidyverse)
library(ggplot2)
library(mice)
library(weights)
library(tictoc)

# source a script with modified functions needed
source("./scripts/99_utils.R")

############ 1 - EXPLORE BIASES ASSOCITED WITH SIB STATUS ###############

# read in processed data from 00_data_preparation.R
load(file = './scratch_data/prepped_00_vars.RData')

dat <- dat_final %>% 
  mutate(par=as.integer(par))


# separate out siblings and singletons
sib_ids= dat %>% 
  filter(duplicated(m_id)) %>% 
  .$m_id
sibs_dat = dat %>% 
  filter(m_id %in% sib_ids)
sing_dat = dat %>% 
  filter(!m_id %in% sibs_dat$m_id)

# make function to compare using t-tests
compare = function(sibs, sings, vars){
  all=data.frame()
  for(var in vars){
    temp=t.test(sibs[,var], sings[,var]) %>% broom::tidy() %>% 
      mutate(pred=var) %>% 
      select(pred, "mean_sibs"=estimate1,"mean_sings"=estimate2,"diff"=estimate, everything())
    all=rbind(all,temp)
  }
  return(all)
}

# run function
comparison = compare(sibs_dat,sing_dat,names(dat %>% 
                                               select(wsPRE:heiC) %>% 
                                               select(-sex,-par)))
comparison %>%  print(n=Inf)

ggplot(comparison %>% 
         filter(pred%in%names(dat %>% select(wsPRE:heiC)))
         , aes(x=diff, y=pred))+
  geom_vline(aes(xintercept=0),colour="grey60", linetype=2,size=1.1)+
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0, size=1.1)+
  geom_point(shape=21, size=4, stroke=1.1, fill="blue")+
  theme_minimal()+
  scale_x_continuous("Mean difference between siblings and singletons (SDs)")+
  scale_y_discrete("Variable")

ggsave("./output/plots/sib_single_comp.tiff", device="tiff", width=15,height=14,units="cm",dpi=320,bg="white")

  

###############################################################################  
############ 2 - MI and IPW to set up adjustment for sib biases ###############
###############################################################################

### MI - run on cluster
dat_for_mi <- dat %>% 
  mutate(indid=as.character(indid),
         m_id=as.character(m_id))

#save(dat_for_mi,file="./data/dat_for_MI.RData")

# example code: check runtime, RAM
# tic()
# MI <- parlmice(dat_for_mi,n.core=parallel::detectCores(), n.imp.core=2, maxit = 2)#for linux, cl.type= "FORK")
# toc()

load("./data/MI_60imp_15its.RData")

# obj name is MI

### IPW
all = complete(MI, "all")
all_weights = dat_for_mi %>% 
  select(indid)

for (i in 1:length(all)){
  
    comp_dat <-  all[[i]]  %>%
    mutate(sibs=ifelse(m_id%in%sib_ids,1,0))
  

  # use to generate ipw
  # stabilized
  ippws <- myipwpoint(exposure=sibs,
                      family="binomial",
                      link="logit",
                      numerator = ~ 1,
                      denominator = ~ 1 +wsPRE + rpPRE + rp18m + rp3yr + rp5yr + rp6mo + le18m + lePRE + le3yr + le5yr + le6mo + bwt + ges + par + emo1 + emo2 + emo3 + beh1 + beh2 + beh3 + neuM + ptsM + adhM + heiM + neuF + ptsF + adhF + heiF + neuC + ptsC + adhC + heiC,
                      data= comp_dat)
  
  range(ippws$ipw.weights)
 
  # check accuracy of predictions
  test <- comp_dat %>% 
    mutate(pred = round(predict(ippws$den.mod, newdata=., type="response"))) 
  
  accuracy <- table(test$pred, test$sibs)
  print(sum(diag(accuracy))/sum(accuracy))
  
  all_weights =cbind(all_weights, ippws$ipw.weights)
}

# rename and save raw weights
colnames(all_weights)<-c("indid", paste0("ipw",seq(1,60)))

save(all_weights, file="./data/all_weights.RData")

load(file="./data/all_weights.RData")


##########################################################
### smooth weights by averaging across all imputations ###
##########################################################

smoothed_weights <- all_weights %>% 
  rowwise() %>% 
  group_by(indid) %>% 
  summarise(ipsw= mean(c_across(ipw1:ipw60)))

save(smoothed_weights, file="./data/smoothed_weights.RData")

load(file="./data/smoothed_weights.RData")


# check performance with weighted t-tests
wdat <- dat %>% 
  left_join(smoothed_weights)

# show relationship between weights and variables
demo = wdat %>% 
  mutate(sib= ifelse(m_id %in% sib_ids,"Sibling","Singleton")) %>% 
  group_by(sib) %>% 
  slice_sample(n=2000) %>% 
  ungroup() %>% 
  select(wsPRE:sib) %>% 
  select(-sex,-par) %>% 
  pivot_longer(wsPRE:heiC, names_to ="Variable")

ggplot(demo , aes(x=ipsw,y=value,size=ipsw, fill=sib)) +
  geom_point(shape=21, alpha=0.6)+
  geom_smooth(aes(color=sib),method="lm", size=1)+
  facet_wrap(vars(Variable))+
  theme_minimal()

ggsave("./output/plots/weights_by_vars.tiff", device="tiff", width=20,height=14,units="cm",dpi=320,bg="white")



sibs_wdat = wdat %>% 
  filter(m_id %in% sib_ids)
sing_wdat = wdat %>% 
  filter(!m_id %in% sib_ids)

# make function to compare using weighted t-tests
wtcompare = function(sibs, sings, vars, wts){
  all=data.frame()
  for(var in vars){
    temp=weights::wtd.t.test(x=sibs[,var], y= sings[,var], weight=sibs[,wts]) 
    temp = c(temp$additional,temp$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("wtd.",.))
    base= weights::wtd.t.test(x=sibs[,var], y= sings[,var]) 
    base = c(base$additional,base$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("unw.",.))
    comp= cbind(base,temp) %>% 
      mutate(pred=var) %>% 
      select(pred,everything())
    all=rbind(all,comp)
  }
  return(as_tibble(all))
}

# run function
wtcomparison = wtcompare(sibs_wdat,sing_wdat,names(dat %>% 
                                                     select(wsPRE:heiC) %>% 
                                                     select(-sex,-par)), "ipsw")
                         
wtcomparison %>% print(n=Inf)
                                                 
wtcomplong = wtcomparison %>% 
  filter(pred%in%names(dat %>% select(wsPRE:heiC))) %>% 
  select(pred, unw.Difference,`unw.Std. Err`,wtd.Difference,`wtd.Std. Err`) %>% 
  pivot_longer(cols=-pred) %>% 
  mutate(name = str_replace_all(name,"Std\\. ","Std")) %>% 
  separate(name, into=c("model","param"), sep="\\." ) %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate(conf.low=Difference-1.96*StdErr,
         conf.high=Difference+1.96*StdErr,
         pred=as.factor(pred),
         model=as.factor(model))


# create plot
ggplot(wtcomplong, aes(x=Difference, y=pred, fill=factor(model)))+
  geom_vline(aes(xintercept=0),colour="black", linetype=2,size=1)+
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0, size=1.1, position=position_dodge(0.4))+
  geom_point(shape=21, size=4, stroke=1.1, position=position_dodge(0.4))+
  theme_minimal(base_size=16)+
  scale_x_continuous("Mean difference between siblings and singletons (SDs)")+
  scale_fill_manual(values=c("blue","red"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))


ggsave("./output/plots/sib_single_comp_weighted.tiff", device="tiff", width=18,height=23,units="cm",dpi=320,bg="white")


# save weighted dataset
save(wdat, file= "./data/sib_weighted_data.RData")

###############################################################################




