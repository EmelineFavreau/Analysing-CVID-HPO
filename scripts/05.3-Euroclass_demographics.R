######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

phbm <- readRDS("../result/patient_hpo_bio_mat.RDS")


####### EUROClass groups Demographics ########################################
# patient subsets with B cells  present, 
i <- c("smb_normal", "smb_minus",
       "cd21_low", "cd21_plus",
       "tr_norm", "tr_high" )
t <- phbm[row.names(phbm) %in% i,]
t <- t[,colSums(t) != 0]

# EUROClass filter
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
