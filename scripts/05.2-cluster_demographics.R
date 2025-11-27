######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load presence/absence matrix for all patients
phbm <- readRDS("../result/patient_hpo_bio_mat.RDS")

# load patient clusters 
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")

# clinical lab label
LC <- fread("../input/HPO_freq_name_labORclinical_LC.csv")
####### Cluster Demographics ############################################
genetic_groups <- c("canonicalTNFRSF13B",       
                    "rareTNFRSF13B",
                    "NFKB1",
                    "anyOtherPathogenic",
                    "anyPathogenic")
# count the number of patients
mat_long <- matrix(data = 0, ncol = 3)
for (cl in c("infection_cluster", "complex_cluster")){
  for (bio in c("smb_normal", "smb_minus",
                "cd21_low", "cd21_plus",
                "tr_norm","tr_high",
                "llh", "lln", "lll", "lnl",
                genetic_groups)){
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
               genetic_groups, 
               "smb_normal & tr_high", "smb_normal & tr_norm")
names(trait_vec) <- c("SmB+", "SmB-",
                      "CD21low high", "CD21low normal",
                      "Transitional B normal", "Transitional B high",
                      "lowIgG-lowIgA-highIgM", "lowIgG-lowIgA-normalIgM",
                      "lowIgG-lowIgA-lowIgM", "lowIgG-normalIgA-lowIgM",
                      "canonical TNFRSF13B", "rare TNFRSF13B", "NFKB1",
                      "any other Pathogenic", "any Pathogenic",
                      "SmB+ & Transitional B high", "SmB+ & Transitional B normal")
# readable name
cdem$long_name <- names(trait_vec)[match(cdem$name, trait_vec)]


cdem$cluster <- gsub(pattern = "_cluster",
                                     x = cdem$cluster,
                                     replacement = "")

# split between bio and gene
cdem$category <- "phenotype"
cdem$category[cdem$name %in% genetic_groups] <- "genotype"


# HPO frequency in infection/complex cohort
HPO_code <- rownames(phbm)[grep("^HP", rownames(phbm))]
HPO_freq <- rowSums(phbm)[grep("^HP", rownames(phbm))]
HPO_name <- hpo$name[match(HPO_code, names(hpo$name))]
Hfnl <- tibble(HPO_code = HPO_code,
                                      HPO_freq = HPO_freq,
                                      HPO_name = HPO_name)

# filter each term in 10 or more patients
Hfnl <- Hfnl %>%dplyr::filter(HPO_freq >= 10)
Hfnl$Category <- LC$Category[match(Hfnl$HPO_code, LC$HPO_code)]
  
############## save all
# main
 
#Splenomegaly  in the complex group,
h = "Splenomegaly"
names(h) <- names(hpo$name)[match(h, hpo$name)]
i <- c("complex_cluster", names(h))
id <- phbm[row.names(phbm) %in% i, ]
s <- sum(colSums(id) == 2) 
(s/nrow(complex_cluster))*100 

#Autoimmune cytopenias  in the complex group, respectively
i <- c("complex_cluster", "HP:0001973")
id <- phbm[row.names(phbm) %in% i, ]
s <- sum(colSums(id) == 2) 
(s/nrow(complex_cluster))*100 # 44.77


fwrite(mat_long, "../result/cluster_bio_demographics.csv")
fwrite(cdem, "../result/cluster_demographics.csv")
fwrite(Hfnl, "../result/HPO_freq_name_labORclinical.csv")