######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load presence/absence matrix for all patients
phbm <- readRDS("../result/patient_hpo_bio_mat.RDS")

# load patient clusters 
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")

####### Cluster Demographics ############################################

# count the number of patients
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
    id <- phbm[row.names(phbm) %in% i, ]
    s <- sum(colSums(id) == 2) 
    mat_long <- rbind(mat_long, c(bio, cl,s))
  } 
}

# tidy the table
mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
colnames(mat_long) <- c("name", "cluster", "patient_number")

# add count for SmB+ Tr norm
for (cl in c("infection_cluster","complex_cluster" )){
  bio = c("smb_normal", "tr_norm")
  i <- c(cl, bio)
  id <- phbm[row.names(phbm) %in% i, ]
  s <- sum(colSums(id) == 3) 
  mat_long <- rbind(mat_long, c("smb_normal & tr_norm", cl,s))
}

# add count for SmB Tr high
for (cl in c("infection_cluster","complex_cluster" )){
  bio = c("smb_normal","tr_high")
  i <- c(cl, bio)
  id <- phbm[row.names(phbm) %in% i, ]
  s <- sum(colSums(id) == 3) 
  mat_long <- rbind(mat_long, c("smb_normal & tr_high", cl,s))
}

# save the table
fwrite(mat_long, "../result/cluster_bio_demographics.csv")

# add number and proportions
cdem <- as.data.frame(mat_long)
cdem$patient_number <- as.numeric(cdem$patient_number)
cdem$cluster_proportion[cdem$cluster == "infection_cluster"] <- 
  cdem$patient_number[cdem$cluster == "infection_cluster"]/nrow(infection_cluster)
cdem$cluster_proportion[cdem$cluster == "complex_cluster"] <- 
  cdem$patient_number[cdem$cluster == "complex_cluster"]/nrow(complex_cluster)

cdem$cluster_proportion <- unlist(cdem$cluster_proportion)*100

# named vector of traits
trait_vec <- c("smb_normal", "smb_minus",
               "cd21_low", "cd21_plus",
               "tr_norm", "tr_high",
               "llh", "lln", "lll", "lnl",
               "canonicalTNFRSF13B", "rareTNFRSF13B", "NFKB1",
               "any Pathogenic variant", 
               "smb_normal & tr_high", "smb_normal & tr_norm")
names(trait_vec) <- c("SmB+", "SmB-",
                      "CD21low high", "CD21low normal",
                      "Transitional B normal", "Transitional B high",
                      "lowIgG-lowIgA-highIgM", "lowIgG-lowIgA-normalIgM",
                      "lowIgG-lowIgA-lowIgM", "lowIgG-normalIgA-lowIgM",
                      "canonical TNFRSF13B", "rare TNFRSF13B", "NFKB1",
                      "anyPathogenic",
                      "SmB+ & Transitional B high", "SmB+ & Transitional B normal")
# readable name
cdem$long_name <- names(trait_vec)[match(cdem$name, trait_vec)]


cdem$cluster <- gsub(pattern = "_cluster",
                                     x = cdem$cluster,
                                     replacement = "")

# split between bio and gene
cdem$category <- "phenotype"
i <- c("anyPathogenic","NFKB1","rareTNFRSF13B","canonicalTNFRSF13B" )
cdem$category[cdem$name %in% i] <- "genotype"

fwrite(cdem, "../result/cluster_demographics.csv")
