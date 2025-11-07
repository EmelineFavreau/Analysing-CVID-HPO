######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")


############# Suppl Table Complete HPO list ####################################

df <- patient_hpo_bio_mat[grep(rownames(patient_hpo_bio_mat), pattern = "^HP"),]

HPO_code = row.names(df)

df1 <- data.table(HPO_code = HPO_code,
                 HPO_freq = rowSums(df),
                 HPO_name = hpo$name[match(HPO_code, names(hpo$name))])
fwrite(df1, "../result/SI/Suppl_Table_Complete_HPO_list.csv")
