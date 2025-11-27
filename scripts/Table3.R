######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# load in the test results for B cell, Immuno, genetic, cluster
br <- fread( "../result/B_cell_HPO_tests_ratio_pvalue.csv") # B cell
ir <- fread( "../result/Immuno_HPO_tests_ratio_pvalue.csv") # immuno
gr <- fread( "../result/genetic_HPO_tests_ratio_pvalue.csv") # genetic


######################## analyse ##########################################
# make one table per test; subset to the biologically important
bat <- br %>% 
  dplyr::filter(HPO_code %in% c("HP:0001744",
                                "HP:0001973",
                                "keyHPOautoimmunity",
                                "keyHPOsplenomegaly")) %>% 
  dplyr::select(-c(biological_measure))


gat <- gr %>% 
  dplyr::select(-c(biological_measure)) %>% 
  dplyr::filter(HPO_code %in% c("HP:0001890", # Autoimmune hemolytic anemia
                                "HP:0005401",#Recurrent candida infections
                                "HP:0001904", # Autoimmune neutropenia /Neutropenia in presence of anti-neutropil antibodies
                                "HP:0001973"))  # Autoimmune thrombocytopenia 


iat <- ir %>% 
  dplyr::select(-c(biological_measure)) %>% 
  dplyr::filter(HPO_code %in% c("HP:0033582", # Pulmonary interstitial lymphocyte infiltration
                                "HP:0410301")) # Partial absence of specific antibody response to unconjugated pneumococcus vaccine 

############ save all
ig_cohort_size <- 500
ig_cohort_size - iat$patients_with_hpo
fwrite(bat, "../result/Table3/Table3A.csv")
fwrite(gat, "../result/Table3/Table3C.csv")
fwrite(iat, "../result/Table3/Table3B.csv")

