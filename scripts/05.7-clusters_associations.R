######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")
unique_codes <- readRDS("../result/unique_codes.RDS")

hpo_vec <- rownames(patient_hpo_bio_mat)[grep("^HP",
                                              rownames(patient_hpo_bio_mat))]

key_hpo_vec <- rownames(patient_hpo_bio_mat)[grep("^key",
                                              rownames(patient_hpo_bio_mat))]
####### Hypothesis D: Association between patient cluster and biology ##########
# Hypothesis D underlying all those tests
# e.g. clustering patients in groups are associated with a given phenotype or genotype

# result p values list
pvalue_cluster_result_list <- list()
contingency_cluster_result_list <- list()
summary_contingency_cluster_result_list <- list()

##### PART 1: B cell
# setting parameters
biological_names <- c("smb", "cd21", "tr")

biological_cat <- c("smb_minus", "smb_normal",
                    "cd21_low", "cd21_plus",
                    "tr_norm", "tr_high")

names(biological_cat) <- c(rep(biological_names[1:3], each = 2))

# looping through biological categories
for(k in 1:length(biological_names)){
  
  # select measure
  focus_name <- biological_names[k]
  focus_cat <- biological_cat[names(biological_cat) == focus_name]
  
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    # loop through the B cell
    for(b in focus_cat){
      # loop through the cluster membership
      for(g in c("infection_cluster", "complex_cluster")){
        # collect the categories
        i <- c(b, g, c) 
        # subset the data to those categories
        id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
        # count the number of patients in this centre and this cluster
        s <- sum(colSums(id) == 3) 
        mat_long <- rbind(mat_long, c(i,s))
      }
    }
  }
  
  # remove rows with 0
  mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
  
  # remove rows if the whole centre has no data for it
  for(c in names(centres_vec)){
    if(sum(as.integer(mat_long[,4][mat_long[,3] == c])) == 0){
      mat_long <- mat_long[mat_long[,3] !=c, ]
    }
  }
  
  # name the columns
  colnames(mat_long) <- c(focus_name, "phenotype", "centre", "count")
  
  # make it into a dataframe
  mat_long <- as.data.frame(mat_long)
  mat_long$count <- as.integer(mat_long$count)
  
  # make a three-dimensional table
  contingency_table <- xtabs(
    as.formula(paste("count ~ phenotype +", focus_name, "+ centre")),
    data = mat_long
  )
  contingency_table1 <- xtabs(
    as.formula(paste("count ~ phenotype +", focus_name)),
    data = mat_long
  )
  
  # test for association, save p value
  cluster_result <- mantelhaen.test(contingency_table)
  pvec <- cluster_result$p.value
  thisname <- paste(focus_name, "_cluster", sep = "")
  names(pvec) <- thisname
  contingency_cluster_result_list[[thisname]] <- contingency_table
  summary_contingency_cluster_result_list[[thisname]] <- contingency_table1
  
  # some results are NaN, because the contingency tables are so imbalanced
  pvalue_cluster_result_list[[k]] <- pvec[!is.nan(pvec)]
  names(pvalue_cluster_result_list)[k] <- focus_name
}


# PART 2: Ig and genetics and key HPO groups
# setting parameters
biological_names <- c("lowIgGlowIgAhighIgM",
                      "lowIgGlowIgAnormalIgM", 
                      "lowIgGlowIgAlowIgM",
                      "lowIgGnormalIgAlowIgM",
                      "canonicalTNFRSF13B",
                      "rareTNFRSF13B",
                      "NFKB1", 
                      "anyPathogenic",
                      key_hpo_vec)

biological_cat <- c("llh",
                    "lln",
                    "lll",
                    "lnl",
                    "canonicalTNFRSF13B",
                    "rareTNFRSF13B",
                    "NFKB1",
                    "anyPathogenic",
                    key_hpo_vec)

names(biological_cat) <- biological_names


# looping through biological categories
for(k in 1:length(biological_names)){
  # select measure
  focus_name <- biological_names[k]
  focus_cat <- biological_cat[names(biological_cat) == focus_name]
  
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    # loop through the cluster membership
    for(g in c("infection_cluster", "complex_cluster")){
      # collect the categories
      i <- c(focus_cat, g, c) 
      # subset the data to those categories
      id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
      # count the number of patients of that centre in the cluster and with the bio
      s <- sum(colSums(id) == 3) 
      mat_long <- rbind(mat_long, c(i,s))
      # count the number of patients of that centre in the cluster and without the bio
      s1 <- sum(id[focus_cat, ] == 0 &
                  id[c, ] == 1 &
                  id[g, ] == 1)
      mat_long <- rbind(mat_long,
                        c(paste("no", focus_cat, sep =""),
                          i[2:3],
                          s1))
      
    }
  }
  
  # remove rows with 0
  mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
  
  # remove rows if the whole centre has no data for it
  for(c in names(centres_vec)){
    if(sum(as.integer(mat_long[,4][mat_long[,3] == c])) == 0){
      mat_long <- mat_long[mat_long[,3] !=c, ]
    }
  }
  
  # name the columns
  colnames(mat_long) <- c(focus_name, "phenotype", "centre", "count")
  
  # make it into a dataframe
  mat_long <- as.data.frame(mat_long)
  mat_long$count <- as.integer(mat_long$count)
  
  # insert edge case about no variable within a stratum
  # remove rows if the whole Ig group has no data for it
  for(c in names(centres_vec)){
    if(sum(mat_long[,4][mat_long[,3] == c][1:2] == 0) == 2 |
       sum(mat_long[,4][mat_long[,3] == c][3:4] == 0) == 2){
      #print("the patients in this centre present no variation for this hpo") 
      mat_long <- mat_long[mat_long[,3] !=c, ]}
  }
  
  # make a three-dimensional table
  contingency_table <- xtabs(
    as.formula(paste("count ~ phenotype +", focus_name, "+ centre")),
    data = mat_long
  )
  contingency_table1 <- xtabs(
    as.formula(paste("count ~ phenotype +", focus_name)),
    data = mat_long
  )
  
  
  # test for association, save p value
  cluster_result <- mantelhaen.test(contingency_table)
  pvec <- cluster_result$p.value
  thisname <- paste(focus_name, "_cluster", sep = "")
  names(pvec) <- thisname
  contingency_cluster_result_list[[thisname]] <- contingency_table
  summary_contingency_cluster_result_list[[thisname]] <- contingency_table1
  
  # some results are NaN, because the contingency tables are so imbalanced
  # B cell occupies 1 to 3, so add
  result_row <- 3
  pvalue_cluster_result_list[[k+result_row]] <- pvec[!is.nan(pvec)]
  names(pvalue_cluster_result_list)[k+result_row] <- focus_name
  
}


# adjustment for multiple comparison Benjamini and Hochberg
BH_cluster_pvalue_result <- p.adjust(unlist(pvalue_cluster_result_list),
                                     method = "BH")
biological_measure <- gsub("\\..*$",
                           names(BH_cluster_pvalue_result[BH_cluster_pvalue_result < 0.05]),
                           replacement = "")
BH_adjust_pvalue <- unname(BH_cluster_pvalue_result[BH_cluster_pvalue_result < 0.05])

BH_cluster_CHM_table <- data.frame(
  biological_measure,BH_adjust_pvalue)

# human friendly p value. 
BH_cluster_CHM_table$BH_adjust_pvalue <- 
  formatC(BH_cluster_CHM_table$BH_adjust_pvalue,
          format = "g")


# save table
fwrite(BH_cluster_CHM_table , "../result/bio_cluster_tests.csv")

# qualify result
c1 <- matrix(NA, ncol = 3)
df <- BH_cluster_CHM_table
rl <- summary_contingency_cluster_result_list
for(i in 1:nrow(df)){
  thisname <- paste(df[i, 1],
                    "cluster", sep = "_")
  c2 <- as.data.table(rl[[thisname]])
  colnames(c2) <- NULL
  c1 <- rbind(c1, c2)
}
c1 <- c1[complete.cases(c1), ]
colnames(c1) <- c("cluster",
                  "biological_measure",
                  "patient_count")

# find the num of patients in each biological group
bg <- unique(c1$biological_measure)
jj <- c()
for (i in 1:length(bg)){
  j <- sum(patient_hpo_bio_mat[rownames(patient_hpo_bio_mat) == bg[i], ])
  jj <- c(jj, j)
}
names(jj) <- bg

c1$total_patient_in_biological_measure <- 
  jj[match(c1$biological_measure, names(jj) )]

c1$percentage_in_biological_measure <- 
  (c1$patient_count/c1$total_patient_in_biological_measure)*100



fwrite(c1,
       "../result/cluster_HPO_tests_summary_contingency_with_ratio.csv")


