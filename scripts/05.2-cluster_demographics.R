######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
#nbr <- readRDS("../result/tidy_data")
patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")
# load patient clusters 
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
####### Cluster Demographics ############################################

mat_long <- matrix(data = 0, ncol = 3)
for (cl in c("infection_cluster", "complex_cluster" )){
  for (bio in c("smb_normal", "smb_minus",
                "cd21_low", "cd21_plus",
                "tr_norm","tr_high",
                "llh", "lln", "lll", "lnl",
                "canonicalTNFRSF13B",
                "rareTNFRSF13B", 
                "NFKB1",
                "anyPathogenic")){
    i <- c(cl, bio)
    id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
    # count the number of patients
    s <- sum(colSums(id) == 2) 
    mat_long <- rbind(mat_long, c(bio, cl,s))
  } 
}

# tidy the table
mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
colnames(mat_long) <- c("name", "cluster", "patient_number")

# add three custom queries
# (1) N SmB+ Tr norm
for (cl in c("infection_cluster","complex_cluster" )){
  bio = c("smb_normal", "tr_norm")
  i <- c(cl, bio)
  id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
  # count the number of patients
  s <- sum(colSums(id) == 3) 
  mat_long <- rbind(mat_long, c("smb_normal & tr_norm", cl,s))
}

# (2) N SmB Tr high
for (cl in c("infection_cluster","complex_cluster" )){
  bio=c("smb_normal","tr_high")
  i <- c(cl, bio)
  id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
  # count the number of patients
  s <- sum(colSums(id) == 3) 
  mat_long <- rbind(mat_long, c("smb_normal & tr_high", cl,s))
}


# save the table
fwrite(mat_long, "../result/cluster_bio_demographics.csv")

cluster_demographics <- as.data.frame(mat_long)
cluster_demographics$patient_number <- as.numeric(cluster_demographics$patient_number)
# The y axis is a "proportion of patients within the cluster" 
cluster_demographics$cluster_proportion[cluster_demographics$cluster ==
                                          "infection_cluster"] <- 
  cluster_demographics$patient_number[
    cluster_demographics$cluster == "infection_cluster"]/nrow(infection_cluster)


cluster_demographics$cluster_proportion[cluster_demographics$cluster == 
                                          "complex_cluster"] <- 
  cluster_demographics$patient_number[
    cluster_demographics$cluster == "complex_cluster"]/nrow(complex_cluster)

cluster_demographics$cluster_proportion <- 
  unlist(cluster_demographics$cluster_proportion)*100

# readable name
cluster_demographics$long_name[cluster_demographics$name == "smb_normal"] <- "SmB+"
cluster_demographics$long_name[cluster_demographics$name == "smb_minus"] <- "SmB-"
cluster_demographics$long_name[cluster_demographics$name == "cd21_low"] <- "CD21low high" # cd21_low >= 10
cluster_demographics$long_name[cluster_demographics$name == "cd21_plus"] <- "CD21low normal" # cd21_plus < 10
cluster_demographics$long_name[cluster_demographics$name == "tr_norm"] <- "Transitional B normal"
cluster_demographics$long_name[cluster_demographics$name == "tr_high"] <- "Transitional B high"
cluster_demographics$long_name[cluster_demographics$name == "llh"] <- "lowIgG-lowIgA-highIgM"
cluster_demographics$long_name[cluster_demographics$name == "lln"] <- "lowIgG-lowIgA-normalIgM"                     
cluster_demographics$long_name[cluster_demographics$name == "lll"] <- "lowIgG-lowIgA-lowIgM"
cluster_demographics$long_name[cluster_demographics$name == "lnl"] <- "lowIgG-normalIgA-lowIgM"
cluster_demographics$long_name[cluster_demographics$name == "canonicalTNFRSF13B"] <- "canonical TNFRSF13B"
cluster_demographics$long_name[cluster_demographics$name == "rareTNFRSF13B"] <- "rare TNFRSF13B"
cluster_demographics$long_name[cluster_demographics$name == "NFKB1"] <- "NFKB1" 
cluster_demographics$long_name[cluster_demographics$name == "anyPathogenic"] <- "any Pathogenic variant"
cluster_demographics$long_name[cluster_demographics$name == "smb_normal & tr_high"] <- "SmB+ & Transitional B high"
cluster_demographics$long_name[cluster_demographics$name == "smb_normal & tr_norm"] <- "SmB+ & Transitional B normal"

cluster_demographics$cluster <- gsub(pattern = "_cluster",
                                     x = cluster_demographics$cluster,
                                     replacement = "")

# split between bio and gene
cluster_demographics$category <- "phenotype"
i <- c("anyPathogenic","NFKB1","rareTNFRSF13B","canonicalTNFRSF13B" )
cluster_demographics$category[cluster_demographics$name %in% i] <- "genotype"
fwrite(cluster_demographics, "../result/cluster_demographics.csv")

