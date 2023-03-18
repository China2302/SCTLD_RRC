# # Function_geno_cast returns a cast data frame with Sample_ID as columns plus the columns 
#describing if the samples have control pair scrapes (C_PE_C_EE) or diseased pair scrapes(D_PE_D_EE)

function_sample_period<- function(Exp) {
  
  geno<- unique(Exp$ID)
  
  temp_data<- data.frame()
  
  for(i in 1:length(geno)) {
    
    geno_select<- dcast(Exp, ID + Resistance ~ Sample_Period , value.var = 'uniqsample', fun.aggregate = sum) %>% 
      
      rename(SP1= '1', SP2="2", SP3="3") %>% 
      
      filter(ID == geno[i]) %>% 
      
      mutate(SP1_SP2_temp = case_when((SP1 == 1 & SP2 == 1) ~ 1)) %>% 
      
      mutate(SP1_SP2 = case_when((SP1_SP2_temp == 1 & SP3 == 0) ~ 1)) %>% 
      
      mutate(SP1_SP3_temp = case_when( (SP1 == 1 & SP3 == 1) ~ 1)) %>% 
      
      mutate(SP1_SP3 = case_when((SP1_SP3_temp == 1 & SP2 == 0) ~ 1)) %>% 
      
      mutate(SP2_SP3_temp = case_when((SP2 == 1 & SP3 == 1) ~ 1)) %>%
      
      mutate(SP2_SP3 = case_when((SP2_SP3_temp == 1 & SP1 == 0) ~ 1)) %>%
      
      mutate(SP1_SP2_SP3= case_when((SP1_SP2_temp == 1 & SP3 == 1) ~ 1)) %>% 

      dplyr::select(ID, SP1_SP2, SP1_SP3, SP1_SP2_SP3, SP2_SP3 )
    
    geno_select<- as.list(geno_select) 
    
    geno_all<- data.frame(ID = geno[i],
                          SP1_SP2 =  geno_select$SP1_SP2,
                          SP1_SP3 = geno_select$SP1_SP3,
                          SP1_SP2_SP3 = geno_select$SP1_SP2_SP3,
                          SP2_SP3 = geno_select$SP2_SP3, 
                          stringsAsFactors = F)
    
    temp_data<- rbind(temp_data, geno_all) 
    
    Pairs_sample_period<- dcast(Exp,Resistance + ID ~ Sample_ID) %>% 
      
      left_join(temp_data, by="ID") 
  }
  return(Pairs_sample_period)
}

#write_rds(function_sample_period,'function_sample_period.rds')
