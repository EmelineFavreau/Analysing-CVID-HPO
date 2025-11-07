######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")
unique_codes <- readRDS("../result/unique_codes.RDS")
####### Hypothesis A: association between B cell and phenotype #################
# a measure of B cell and HPO are associated, controlling for within-centre variability
# Cochran-Mantel-Haenszel test for SmB/CD21/Tr level and each phenotype (HPO, GLILD)
hpo_vec <- rownames(patient_hpo_bio_mat)[grep("^HP",
                                              rownames(patient_hpo_bio_mat))]

key_hpo_vec <- rownames(patient_hpo_bio_mat)[grep("^key",
                                              rownames(patient_hpo_bio_mat))]


# result p values list
pvalue_result_list <- list()
contingency_Bcell_result_list <- list()
summary_contingency_Bcell_result_list <- list()

# setting parameters
biological_names <- c("smb", "cd21", "tr")

biological_cat <- c("smb_minus", "smb_normal",
                    "cd21_low", "cd21_plus",
                    "tr_norm", "tr_high")

names(biological_cat) <- rep(biological_names, each = 2)

# looping through biological categories
for(k in 1:length(biological_names)){
  #print(k)
  # for GLILD:
  focus_name <- biological_names[k]
  focus_cat <- biological_cat[names(biological_cat) == focus_name]
  #print(paste("starting for ", focus_name, " and ", focus_cat, sep =""))
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    #print(paste("starting for ", c, sep=""))
    # loop through the B cell
    for(b in focus_cat){
      # loop through the glild
      for(g in c("reported_glild", "reported_no_glild")){
        # collect the categories
        i <- c(b, g, c) 
        # subset the data to those categories
        id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
        # count the number of patients
        s <- sum(colSums(id) == 3) 
        mat_long <- rbind(mat_long, c(i,s))
      }
    }
  }
  
  # remove rows with 0
  mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
  
  # remove rows if the whole centre has no data for it
  for(c in names(centres_vec)){
    if(sum(as.integer(mat_long[,4][mat_long[,3] == c])) <= 1){
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
  glild_result <- mantelhaen.test(contingency_table)
  pvec <- glild_result$p.value
  thisname <- paste(focus_name, "_glild", sep = "")
  names(pvec) <- thisname
  contingency_Bcell_result_list[[thisname]] <- contingency_table
  summary_contingency_Bcell_result_list[[thisname]] <- contingency_table1
  
  # for each HPO:
  for(h in 1:length(hpo_vec)){
    # get the HPO code
    g <- hpo_vec[h]
    
    # initiate matrix
    mat_long <- matrix(data = 0,
                       ncol = 4)
    
    # loop through the centres
    for(c in names(centres_vec)){
      # loop through the B cell
      for(b in focus_cat){
        # collect the categories
        i <- c(b, g, c) 
        # subset the data to those categories
        id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
        # count the number of patients
        s <- sum(colSums(id) == 3) 
        if(is.integer(s)){
          mat_long <- rbind(mat_long, c(b, paste(g, 1, sep = "_"), c, s))
        } else {
          mat_long <- rbind(mat_long, rep(0, 4))
        }
        # count the number of patients without HPO
        t <- sum(id[g, ] == 0 &
                   id[b, ] == 1 &
                   id[c, ] == 1)
        if(is.integer(t)){
          mat_long <- rbind(mat_long, c(b,paste(g, 0, sep = "_"),c,t))
        } else {
          mat_long <- rbind(mat_long, rep(0, 4))
        }
      }
    }
    
    # remove rows with 0
    mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
    
    # remove rows if the whole centre has no data for it
    for(c in names(centres_vec)){
      if(sum(as.integer(mat_long[,4][mat_long[,3] == c])) <= 1){
        mat_long <- mat_long[mat_long[,3] !=c, ]
      }
    }
    
    # name the columns
    colnames(mat_long) <- c(focus_name, "phenotype", "centre", "count")
    
    # make it into a dataframe
    mat_long <- as.data.frame(mat_long)
    mat_long$count <- as.integer(mat_long$count)
    
    # make a three-dimensional table
    hpo_contingency_table <- xtabs(
      as.formula(paste("count ~ phenotype +", focus_name, "+ centre")),
      data = mat_long
    )
    
    # make a summary table
    hpo_contingency_table1 <- xtabs(
      as.formula(paste("count ~ phenotype +", focus_name)),
      data = mat_long
    )
    
    # test for association
    hpo_result <- mantelhaen.test(hpo_contingency_table)
    cmh_p <- hpo_result$p.value
    thisname <- paste(focus_name, g, sep = "_")
    names(cmh_p) <- thisname
    # combine p vec of all SmB tests (glild and HPOs)
    pvec <- c(pvec, cmh_p)
    # save the contingency table
    contingency_Bcell_result_list[[thisname]] <- hpo_contingency_table
    # save the summary contingency table
    summary_contingency_Bcell_result_list[[thisname]] <- hpo_contingency_table1
  } # end of the loop for each HPO term
  
  
  
  # for each key hpo group:
  for(z in 1:length(key_hpo_vec)){
    #print(z)
    # select key group
    g <- key_hpo_vec[z]
    
    # initiate matrix
    mat_long <- matrix(data = 0,
                       ncol = 4)
    
    # loop through the centres
    for(c in names(centres_vec)){
      # loop through the B cell
      for(b in focus_cat){
        # collect the categories
        i <- c(b, g, c) 
        # subset the data to those categories
        id <- patient_hpo_bio_mat[row.names(patient_hpo_bio_mat) %in% i, ]
        # count the number of patients with key HPO group
        s <- sum(colSums(id) == 3) 
        if(is.integer(s)){
          mat_long <- rbind(mat_long, c(b, paste(g, 1, sep = "_"), c, s))
        } else {
          mat_long <- rbind(mat_long, rep(0, 4))
        }
        # count the number of patients without key HPO group
        t <- sum(id[g, ] == 0 &
                   id[b, ] == 1 &
                   id[c, ] == 1)
        if(is.integer(t)){
          mat_long <- rbind(mat_long, c(b,paste(g, 0, sep = "_"),c,t))
        } else {
          mat_long <- rbind(mat_long, rep(0, 4))
        }
      }
    }
    
    
    # remove rows with 0
    mat_long <- mat_long[rowSums(mat_long != 0) > 0, ]
    
    # sample size in each stratum must be > 1
    
    # remove rows if the whole centre has 0 or 1 
    for(c in names(centres_vec)){
      if(sum(as.integer(mat_long[,4][mat_long[,3] == c])) <= 1){
        mat_long <- mat_long[mat_long[,3] !=c, ]
      }
    }
    
    # name the columns
    colnames(mat_long) <- c(focus_name, "phenotype", "centre", "count")
    
    # make it into a dataframe
    mat_long <- as.data.frame(mat_long)
    mat_long$count <- as.integer(mat_long$count)
    
    # make a three-dimensional table
    hpo_contingency_table <- xtabs(
      as.formula(paste("count ~ phenotype +", focus_name, "+ centre")),
      data = mat_long
    )
    
    # make a summary table
    hpo_contingency_table1 <- xtabs(
      as.formula(paste("count ~ phenotype +", focus_name)),
      data = mat_long
    )
    
    # test for association
    hpo_result <- mantelhaen.test(hpo_contingency_table)
    cmh_p <- hpo_result$p.value
    thisname <- paste(focus_name, g, sep = "_")
    names(cmh_p) <- thisname
    
    # combine p vec of all SmB tests (glild and HPOs)
    pvec <- c(pvec, cmh_p)
    # save the contingency table
    contingency_Bcell_result_list[[thisname]] <- hpo_contingency_table
    # save the summary contingency table
    summary_contingency_Bcell_result_list[[thisname]] <- hpo_contingency_table1
    
    # end of loop for each key hpo group
    #print("end of loop for each key hpo group")
  }
  
  #print(k)
  # some results are NaN, because the contingency tables are so imbalanced
  pvalue_result_list[[k]] <- pvec[!is.nan(pvec)]
  names(pvalue_result_list)[k] <- focus_name
}

# adjustment for multiple comparison Benjamini and Hochberg
BH_pvalue_result <- p.adjust(unlist(pvalue_result_list), method = "BH")

biological_measure <- gsub("\\..*$",
                           names(BH_pvalue_result[BH_pvalue_result < 0.05]),
                           replacement = "")
phenotype_code <- gsub("^.*_",
                       names(BH_pvalue_result[BH_pvalue_result < 0.05]),
                       replacement = "")

phenotype_name <- names(unique_codes)[match(phenotype_code, unique_codes)]

BH_adjust_pvalue <- unname(BH_pvalue_result[BH_pvalue_result < 0.05])

BH_Bcell_CHM_table <- data.frame(
  biological_measure,
  phenotype_code,
  phenotype_name,
  BH_adjust_pvalue)

# human friendly p value. 
BH_Bcell_CHM_table$BH_adjust_pvalue <- 
  formatC(BH_Bcell_CHM_table$BH_adjust_pvalue,
          format = "g")

# save table
fwrite(BH_Bcell_CHM_table, "../result/B_cell_HPO_tests.csv")

# qualify result
c1 <- matrix(NA, ncol = 3)
for(i in 1:nrow(BH_Bcell_CHM_table)){
  thisname <- paste(BH_Bcell_CHM_table[i, 1],
                    BH_Bcell_CHM_table[i, 2], sep = "_")
  c2 <- as.data.table(summary_contingency_Bcell_result_list[[thisname]])
  colnames(c2) <- NULL
  c1 <- rbind(c1, c2)
}
c1 <- c1[complete.cases(c1), ]
colnames(c1) <- c("HPO_presence1_absence0",
                  "phenotype",
                  "patient_count")


# HPO_code | phenotype | patients_with_hpo | total_patients_in_group | percentage_with_hpo
result <- c1 %>%
  # new column for the HPO_code by removing _0 or _1
  mutate(HPO_code = str_remove(HPO_presence1_absence0, "_[01]$")) %>%
  
  # group by the HPO_code and phenotype
  group_by(HPO_code, phenotype) %>%
  
  # summarise the data for each group: sum of patients with HPO, 
  # sum of patients in the group
  summarise(
 patients_with_hpo = sum(patient_count[str_ends(HPO_presence1_absence0, "1")]),
    total_patients_in_group = sum(patient_count),
    .groups = 'drop') %>%
  
  # get the percentage based on the group totals
  mutate(percentage_with_hpo = 
           (patients_with_hpo / total_patients_in_group) * 100) 

fwrite(result,
       "../result/B_cell_HPO_tests_summary_contingency_with_ratio.csv")



