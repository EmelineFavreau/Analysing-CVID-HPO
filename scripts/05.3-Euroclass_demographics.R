######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")


####### EUROClass groups Demographics ########################################
# Following the EUROCLASS decision tree from the paper

# If subsets of B cells ( CD21low, smB and transitional ) are present, 
# then the patient has more than 1% of total B cells
# if there is a row of 0, 
# this means patients have less than 1% and  are not included
i <- c("smb_normal", "smb_minus",
       "cd21_low", "cd21_plus", # cd21_low >10, cd21_plus < 10
       "tr_norm", "tr_high" )
t <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i,]
t <- t[,colSums(t) != 0]

# it is imperative to have a measure for smb; 
# if not the data are excluded from this analysis
tt <- t[,t[rownames(t) %in% c("smb_minus",  "smb_normal" ), ] %>%
          colSums(.) == 1]

# create a final dataset
euroclass_STUDYID_df <- as.data.frame(rbind(
  cbind("smBplus_Euroclass",
        colnames(tt)[tt["smb_normal",] == 1]),
  
  cbind("smBplus21norm_Euroclass",
        colnames(tt)[tt["smb_normal",] == 1 & tt["cd21_plus",] == 1]),
  
  cbind("smBplus21lo_Euroclass", 
        colnames(tt)[tt["smb_normal",] == 1 & tt["cd21_low",] == 1]),
  
  cbind("smBminus_Euroclass",
        colnames(tt)[tt["smb_minus",] == 1]),
  
  cbind("smBminusTrhi_Euroclass",
        colnames(tt)[tt["smb_minus",] == 1 & tt["tr_high",] == 1]),
  
  cbind("smBminusTrnorm_Euroclass",
        colnames(tt)[tt["smb_minus",] == 1 & tt["tr_norm",] == 1]),
  
  cbind("smBminusCD21norm_Euroclass",
        colnames(tt)[tt["smb_minus",] == 1 & tt["cd21_plus",] == 1]),
  
  cbind("smBminusCD21lo_Euroclass",
        colnames(tt)[tt["smb_minus",] == 1 & tt["cd21_low",] == 1])))


# update colnames
colnames(euroclass_STUDYID_df) <- c("euroclass_group", "STUDY_ID")

# save the table
fwrite(euroclass_STUDYID_df, "../result/euroclass_STUDYID.csv")
