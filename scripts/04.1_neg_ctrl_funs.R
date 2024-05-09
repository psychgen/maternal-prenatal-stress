#04.1_neg_ctrl_funs.R

# The purpose of this script is to make functions for running models in 
# 04_negative_control.R

library(tictoc)

run_test_models <- function(outname, #the outcome name (root only if emo/beh)
                            expwave, #the wave of the exposure
                            data){ #the dataset
  
  
  message(paste0("Running models for ",outname," with ", expwave))
  
  modtype <- if(expwave=="PRE"){
    "expout"
  }else if(expwave%in%c("18m","3yr","5yr") & outname %in% c("bwt","ges")){
    "validnc"
  }else if(expwave%in%c("3yr","5yr") &  str_detect(outname,"1") ){
    "validnc"
  }else if(expwave%in%c("5yr") &  str_detect(outname,"2")){
    "validnc"
  }else { "invalid"}

    #Specify
    
    model <-
      paste0(outname,' ~ b1*rp',expwave,' +b2*le',expwave,' + b10*sex + b13*par ')
    
    #Fit the model 
    
    fit <-lavaan::sem(model,
                      data=data,
                      estimator="MLR",
                      missing="fiml.x",
                      se="robust",
                      cluster= "m_id")

  
  message("All model-fitting for this outcome-exposure combination is complete.")
  #Create a function to extract standardized estimates in the desired format
  get_ests <- function(x){
    
    model_est <- standardizedsolution(x, type="std.nox") %>% 
      filter (op == "~" ) %>% 
      as_tibble() %>% 
      filter(label%in%c("b1","b2") ) %>% 
      rename('out'=lhs,'exp'=rhs) %>% 
      mutate(model= modtype)
    return(model_est)
    
    
  }
 
  #Return results in a df
  
  return(get_ests(fit))
         
  
}
