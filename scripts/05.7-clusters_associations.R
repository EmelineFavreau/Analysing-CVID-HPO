######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# presence/absence matrix
phbm <- readRDS("../result/patient_hpo_bio_mat.RDS")

# HPO codes
unique_codes <- readRDS("../result/unique_codes.RDS")

# present HPO codes
hpo_vec <- rownames(phbm)[grep("^HP", rownames(phbm))]

# present Key HPO groups
key_hpo_vec <- rownames(phbm)[grep("^key", rownames(phbm))]

######################## association tests #####################################
# association between infection/conplex and biomarkers

 
# result p values list
p_list <- list()
c_list <- list()
s_list <- list()

##### PART 1: B cell
# setting parameters
bn <- c("smb", "cd21", "tr")

bc <- c("smb_minus", "smb_normal",
                    "cd21_low", "cd21_plus",
                    "tr_norm", "tr_high")

names(bc) <- c(rep(bn[1:3], each = 2))

# looping through biological categories
for(k in 1:length(bn)){
  
  # select measure
  fn <- bn[k]
  fc <- bc[names(bc) == fn]
  
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    # loop through the B cell
    for(b in fc){
      # loop through the cluster membership
      for(g in c("infection_cluster", "complex_cluster")){
        # collect the categories
        i <- c(b, g, c) 
        # subset the data to those categories
        id <- phbm[row.names(phbm) %in% i, ]
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
  colnames(mat_long) <- c(fn, "phenotype", "centre", "count")
  
  # make it into a dataframe
  mat_long <- as.data.frame(mat_long)
  mat_long$count <- as.integer(mat_long$count)
  
  # make a three-dimensional table
  contingency_table <- xtabs(
    as.formula(paste("count ~ phenotype +", fn, "+ centre")),
    data = mat_long
  )
  contingency_table1 <- xtabs(
    as.formula(paste("count ~ phenotype +", fn)),
    data = mat_long
  )
  
  # test for association, save p value
  cluster_result <- mantelhaen.test(contingency_table)
  pvec <- cluster_result$p.value
  thisname <- paste(fn, "_cluster", sep = "")
  names(pvec) <- thisname
  c_list[[thisname]] <- contingency_table
  s_list[[thisname]] <- contingency_table1
  
  # some results are NaN, because the contingency tables are so imbalanced
  p_list[[k]] <- pvec[!is.nan(pvec)]
  names(p_list)[k] <- fn
}


# PART 2: Ig and genetics and key HPO groups
# setting parameters
genetic <- c("canonicalTNFRSF13B",       
              "rareTNFRSF13B",
              "NFKB1",
              "anyOtherPathogenic",
              "anyPathogenic")
bn <- c("lowIgGlowIgAhighIgM",
        "lowIgGlowIgAnormalIgM", 
        "lowIgGlowIgAlowIgM",
        "lowIgGnormalIgAlowIgM",
        genetic, key_hpo_vec)

bc <- c("llh",
        "lln",
        "lll",
        "lnl",
        genetic, key_hpo_vec)

names(bc) <- bn


# looping through biological categories
for(k in 1:length(bn)){
  # select measure
  fn <- bn[k]
  fc <- bc[names(bc) == fn]
  
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    # loop through the cluster membership
    for(g in c("infection_cluster", "complex_cluster")){
      # collect the categories
      i <- c(fc, g, c) 
      # subset the data to those categories
      id <- phbm[row.names(phbm) %in% i, ]
      # count the number of patients of that centre in the cluster and with the bio
      s <- sum(colSums(id) == 3) 
      mat_long <- rbind(mat_long, c(i,s))
      # count the number of patients of that centre in the cluster and without the bio
      s1 <- sum(id[fc, ] == 0 &
                  id[c, ] == 1 &
                  id[g, ] == 1)
      mat_long <- rbind(mat_long,
                        c(paste("no", fc, sep =""),
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
  colnames(mat_long) <- c(fn, "phenotype", "centre", "count")
  
  # make it into a dataframe
  mat_long <- as.data.frame(mat_long)
  mat_long$count <- as.integer(mat_long$count)
  
  # insert edge case about no variable within a stratum
  # remove rows if the whole Ig group has no data for it
  for(c in names(centres_vec)){
    if(sum(mat_long[,4][mat_long[,3] == c][1:2] == 0) == 2 |
       sum(mat_long[,4][mat_long[,3] == c][3:4] == 0) == 2){
      mat_long <- mat_long[mat_long[,3] !=c, ]}
  }
  
  # make a three-dimensional table
  contingency_table <- xtabs(
    as.formula(paste("count ~ phenotype +", fn, "+ centre")),
    data = mat_long
  )
  contingency_table1 <- xtabs(
    as.formula(paste("count ~ phenotype +", fn)),
    data = mat_long
  )
  
  
  # test for association, save p value
  cluster_result <- mantelhaen.test(contingency_table)
  pvec <- cluster_result$p.value
  thisname <- paste(fn, "_cluster", sep = "")
  names(pvec) <- thisname
  c_list[[thisname]] <- contingency_table
  s_list[[thisname]] <- contingency_table1
  
  # some results are NaN, because the contingency tables are so imbalanced
  # B cell occupies 1 to 3, so add
  result_row <- 3
  p_list[[k+result_row]] <- pvec[!is.nan(pvec)]
  names(p_list)[k+result_row] <- fn
  
}


# adjustment for multiple comparison Benjamini and Hochberg
BH_cluster_pvalue_result <- p.adjust(unlist(p_list),
                                     method = "BH")
biological_measure <- gsub("\\..*$",
             names(BH_cluster_pvalue_result[BH_cluster_pvalue_result < 0.05]),
             replacement = "")
BH_adjust_pvalue <- 
  unname(BH_cluster_pvalue_result[BH_cluster_pvalue_result < 0.05])

BH_cluster_CHM_table <- data.frame(
  biological_measure,BH_adjust_pvalue)

# human friendly p value. 
BH_cluster_CHM_table$BH_adjust_pvalue <- 
  formatC(BH_cluster_CHM_table$BH_adjust_pvalue,
          format = "g")

# qualify result
c1 <- matrix(NA, ncol = 3)
df <- BH_cluster_CHM_table
rl <- s_list
for(i in 1:nrow(df)){
  thisname <- paste(df[i, 1],
                    "cluster", sep = "_")
  c2 <- as.data.table(rl[[thisname]])
  colnames(c2) <- NULL
  c1 <- rbind(c1, c2)
}
c1 <- c1[complete.cases(c1), ]
colnames(c1) <- c("cluster",
                  "biological_category",
                  "patient_count")

# find the num of patients in each biological group
bg <- unique(c1$biological_category)
jj <- c()
for (i in 1:length(bg)){
  j <- sum(phbm[rownames(phbm) == bg[i], ])
  jj <- c(jj, j)
}
names(jj) <- bg

c1$total_patient_in_biological_measure <- 
  jj[match(c1$biological_category, names(jj) )]

c1$percentage_in_biological_measure <- 
  (c1$patient_count/c1$total_patient_in_biological_measure)*100

c1$biological_measure <- ifelse(grepl("_", c1$biological_category), 
            gsub(pattern = "_.*", c1$biological_category, replacement = ""),
            gsub(pattern = "^no", c1$biological_category, replacement = ""))

# merge info
merged_tibble <- c1 %>%
  full_join(BH_cluster_CHM_table, by = "biological_measure")

############## save all
fwrite(merged_tibble, "../result/cluster_biomarker_tests_ratio_pvalue.csv")

