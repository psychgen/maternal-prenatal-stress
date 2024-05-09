#01a_linear_model_funs.R

# The purpose of this script is to make functions for running models in 
# 01_run_linear_models.R

library(tictoc)

run_test_models <- function(outname, #the outcome name (root only if emo/beh)
                            modname, #the root of the moderating PGS, before .pgs.pc and the role-identifying suffix
                            data){ #the dataset
  
  
  message(paste0("Running models for ",outname,", with moderator ",modname))
  
  
  #Code for simple models
  if(outname %in% c("bwt","ges")){
    
    #Specify
    
    model <-
      paste0(outname,' ~ b1*wsPRE + b2*rpPRE +b3*lePRE + b4*',modname,'M +',
             'b5*sex + b6*par + b7*wsPRE:',modname,'M + b8*rpPRE:',modname,
             'M +b9*lePRE:',modname,'M + b10*wsPRE:sex + b11*rpPRE:sex +', 
             'b12*lePRE:sex + b13*wsPRE:par + b14*rpPRE:par +', 
             'b15*lePRE:par + b16*',modname,'M:sex + b17*',modname,'M:sex +',
             'b18*',modname,'M:par')
    
    #Fit the model with the tier 1 covariates
    
    fit <-lavaan::sem(model,
                      data=data,
                      estimator="MLR",
                      missing="fiml.x",
                      se="robust",
                      cluster= "m_id")
    lrt=NULL
  }else{
    #Code for multi-wave outcome models
    model <-
      paste0('#Regression 18m
             ',outname,'1 ~ b11*wsPRE + b21*rpPRE +b31*lePRE + b41*',modname,'M +',
             'b51*sex + b61*par + b71*wsPRE:',modname,'M + b81*rpPRE:',modname,
             'M +b91*lePRE:',modname,'M + b101*wsPRE:sex + b111*rpPRE:sex +', 
             'b121*lePRE:sex + b131*wsPRE:par + b141*rpPRE:par +', 
             'b151*lePRE:par + b161*',modname,'M:sex + b171*',modname,'M:sex +',
             'b181*',modname,'M:par
             #Regression 3yr
             ',outname,'2 ~ b12*wsPRE + b22*rpPRE +b32*lePRE + b42*',modname,'M +',
             'b52*sex + b62*par + b72*wsPRE:',modname,'M + b82*rpPRE:',modname,
             'M +b92*lePRE:',modname,'M + b102*wsPRE:sex + b112*rpPRE:sex +', 
             'b122*lePRE:sex + b132*wsPRE:par + b142*rpPRE:par +', 
             'b152*lePRE:par + b162*',modname,'M:sex + b172*',modname,'M:sex +',
             'b182*',modname,'M:par
             #Regression 5yr
             ',outname,'3 ~ b13*wsPRE + b23*rpPRE +b33*lePRE + b43*',modname,'M +',
             'b53*sex + b63*par + b73*wsPRE:',modname,'M + b83*rpPRE:',modname,
             'M +b93*lePRE:',modname,'M + b103*wsPRE:sex + b113*rpPRE:sex +', 
             'b123*lePRE:sex + b133*wsPRE:par + b143*rpPRE:par +', 
             'b153*lePRE:par + b163*',modname,'M:sex + b173*',modname,'M:sex +',
             'b183*',modname,'M:par
             # Residual variances
           ',outname,'1 ~~ outres1*',outname,'1
           ',outname,'2 ~~ outres2*',outname,'2
           ',outname,'3 ~~ outres3*',outname,'3')  
    
    
    #Fit the model and a constrained (across wave) version
    
    fit <-lavaan::sem(model,
                      data=data,
                      estimator="MLR",
                      missing="fiml.x",
                      se="robust",
                      cluster= "m_id") 
    fit_cons <-lavaan::sem(paste0(model,'
                                              # Constraints
                                              b11==b12
                                              b12==b13
                                              b21==b22
                                              b22==b23
                                              b31==b32
                                              b32==b33
                                              b71==b72
                                              b72==b73
                                              b81==b82
                                              b82==b83
                                              b91==b92
                                              b92==b93'),
                           data=data,
                           estimator="MLR",
                           missing="fiml.x",
                           se="robust",
                           cluster= "m_id") 
    
    #LRT to select a best fitting model
    lrt=lavaan::lavTestLRT(fit,fit_cons)
    if(lrt$`Pr(>Chisq)`[2]<0.05/2){
      fit=fit
    }else{
      fit=fit_cons
    }
    
  }
  #Now we have a fit object, regardless of the outcome type
  
  message("All model-fitting for this outcome-moderator combination is complete.")
  #Create a function to extract standardized estimates in the desired format
  get_ests <- function(x){
    
    model_est <- standardizedsolution(x, type="std.nox") %>% 
      filter (op == "~" ) %>% 
      as_tibble() %>% 
      filter(rhs%in%c("wsPRE","rpPRE","lePRE")|str_detect(rhs,paste0(":",modname) )) %>% 
      mutate(mod=modname)
    return(model_est)
    
    
  }
  
  #Return results in a nested list
  
  return(list("fit"=fit,"LRT"=lrt,"ests"=get_ests(fit) ))
  
  
}
