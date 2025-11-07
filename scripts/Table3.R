######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# load in the test results for B cell, Immuno, genetic, cluster
br <- fread( "../result/B_cell_HPO_tests.csv") # B cell
ir <- fread( "../result/Immuno_HPO_tests.csv") # immuno
gr <- fread( "../result/genetic_HPO_tests.csv") # genetic
cr <- fread( "../result/bio_cluster_tests.csv") # cluster

# load in the ratio
brr <- fread( "../result/B_cell_HPO_tests_summary_contingency_with_ratio.csv")
irr <- fread( "../result/Immuno_HPO_tests_summary_contingency_with_ratio.csv")
grr <- fread( "../result/genetic_HPO_tests_summary_contingency_with_ratio.csv")
crr <- fread( "../result/cluster_HPO_tests_summary_contingency_with_ratio.csv")

######################## analyse ##########################################
# make one table per test
bat <- brr %>% 
  mutate(adjP = br$BH_adjust_pvalue[match(brr$HPO_code, br$phenotype_code)],
         phenotype_name = br$phenotype_name[match(brr$HPO_code,
                                                  br$phenotype_code)]) %>% 
  dplyr::select(-c(patients_with_hpo, total_patients_in_group))

fwrite(bat, "../result/Table3/Table3A.csv")

gat <- grr %>% 
  mutate(adjP = gr$BH_adjust_pvalue[match(grr$HPO_code, gr$phenotype_code)],
         phenotype_name = gr$phenotype_name[match(grr$HPO_code,
                                                  gr$phenotype_code)]) %>% 
  dplyr::select(-c(patients_with_hpo, total_patients_in_group))

fwrite(gat, "../result/Table3/Table3B.csv")

# immuno
irr$biological_measure <- ir$biological_measure[match(irr$HPO_code, 
                                                      ir$phenotype_code)]

iat <- irr %>% 
  mutate(adjP = ir$BH_adjust_pvalue[match(irr$HPO_code, ir$phenotype_code)],
         phenotype_name = ir$phenotype_name[match(irr$HPO_code,
                                                  ir$phenotype_code)]) %>% 
  group_by(HPO_code,	phenotype) %>% 
  dplyr::select(-c(patients_with_hpo, total_patients_in_group))

fwrite(iat, "../result/Table3/Table3C.csv")


crr$biological_category <- crr$biological_measure %>% 
  gsub(., pattern = "^no", replacement = "") %>% 
  gsub(., pattern = "_.*$", replacement = "")


catr <- crr %>% 
  mutate(adjP = cr$BH_adjust_pvalue[match(crr$biological_category, 
                                          cr$biological_measure)]) %>% 
  dplyr::select(-c(patient_count, total_patient_in_biological_measure))

fwrite(catr, "../result/Table3/Table3D.csv")
