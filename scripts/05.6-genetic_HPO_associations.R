######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# presence/absence matrix
phbm <- readRDS("../result/patient_hpo_bio_mat.RDS")

# HPO codes
unique_codes <- readRDS("../result/unique_codes.RDS")

# present HPO codes
hpo_vec <- rownames(phbm)[grep("^HP", rownames(phbm))]

# present key HPO codes
key_hpo_vec <- rownames(phbm)[grep("^key", rownames(phbm))]


######################## association tests #####################################
# association between genetic variants and phenotypes

# setting parameters
bn <- c("canonicalTNFRSF13B",
                      "rareTNFRSF13B",
                      "NFKB1",
                      "anyPathogenic")

# result p values list
p_list <- list()
c_list <- list()
s_list <- list()

# looping through biological categories
for(k in 1:length(bn)){
  ####### for GLILD ##########
  fn <- bn[k]
  
  # initiate matrix
  mat_long <- matrix(data = 0,
                     ncol = 4)
  
  # loop through the centres
  for(c in names(centres_vec)){
    # loop through glild
    for(g in c("reported_glild", "reported_no_glild")){
      # collect the categories
      i <- c(fn, g, c) 
      # subset the data to those categories
      id <- phbm[row.names(phbm) %in% i, ]
      # count the number of patients with the variant
      s <- sum(colSums(id) == 3) 
      mat_long <- rbind(mat_long, c(i,s))
      # count the number of patients without the variant
      s1 <- sum(id[g, ] == 1 &
                  id[c, ] == 1 &
                  id[fn, ] == 0)
      mat_long <- rbind(mat_long,
                        c(paste("without", i[1], sep =""),
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
  
  # make a three-dimensional table 
  ctw <- xtabs(
    as.formula(paste("count ~ phenotype +", fn, "+ centre")),
    data = mat_long
  )
  ctwo <- xtabs(
    as.formula(paste("count ~ phenotype +", fn)),
    data = mat_long
  )
  
  
  # test for association, save p value
  glild_result <- chisq.test(ctwo)
  pvec <- glild_result$p.value
  thisname <- paste(fn, "_glild", sep = "")
  names(pvec) <- thisname
  c_list[[thisname]] <- ctw
  s_list[[thisname]] <- ctwo
  
  
  ####### for HPO ##########
  for(h in 1:length(hpo_vec)){
    # get the HPO code
    g <- hpo_vec[h]
    
    # initiate matrix
    mat_long <- matrix(data = 0,
                       ncol = 4)
    
    # loop through the centres
    for(c in names(centres_vec)){
      # collect the categories
      i <- c(fn, g, c) 
      
      # subset the data to those categories
      id <- phbm[row.names(phbm) %in% i, ]
      
      # count the number of patients with the variant and the HPO
      s <- sum(colSums(id) == 3) 
      if(is.integer(s)){
        mat_long <- rbind(mat_long, c(fn,
                                      paste(g, 1, sep = "_"),
                                      c,
                                      s))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients without the variant but the HPO
      s1 <- sum(id[g, ] == 1 &
                  id[c, ] == 1 &
                  id[fn, ] == 0)
      if(is.integer(s1)){
        mat_long <- rbind(mat_long,
                          c(paste("without", i[1], sep =""),
                            paste(g, 1, sep = "_"),
                            c,
                            s1))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients with the variant but not the HPO
      s2 <- sum(id[g, ] == 0 &
                  id[c, ] == 1 &
                  id[fn, ] == 1)
      if(is.integer(s2)){
        mat_long <- rbind(mat_long, c(fn,
                                      paste(g, 0, sep = "_"),
                                      c,
                                      s2))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients without the variant and without the HPO
      s3 <- sum(id[g, ] == 0 &
                  id[c, ] == 1 &
                  id[fn, ] == 0)
      if(is.integer(s1)){
        mat_long <- rbind(mat_long,
                          c(paste("without", i[1], sep =""),
                            paste(g, 0, sep = "_"),
                            c,
                            s3))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
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
    hpo_ctw <- xtabs(
      as.formula(paste("count ~ phenotype +", fn, "+ centre")),
      data = mat_long
    )
    
    # make a summary table
    hpo_ctwo <- xtabs(
      as.formula(paste("count ~ phenotype +", fn)),
      data = mat_long
    )
    
    # test for association
    hpo_result <- chisq.test(hpo_ctwo)
    cmh_p <- hpo_result$p.value
    thisname <- paste(fn, g, sep = "_")
    names(cmh_p) <- thisname
    # combine p vec of all genetic tests (glild and HPOs)
    pvec <- c(pvec, cmh_p)
    # save the contingency table
    c_list[[thisname]] <- hpo_ctw
    # save the summary contingency table
    s_list[[thisname]] <- hpo_ctwo
  } 
  
  
  ####### for KEY HPO group ##########
  for(z in 1:length(key_hpo_vec)){
    # select key group
    g <- key_hpo_vec[z]
    
    # initiate matrix
    mat_long <- matrix(data = 0,
                       ncol = 4)
    
    # loop through the centres
    for(c in names(centres_vec)){
      # collect the categories
      i <- c(fn, g, c) 
      
      # subset the data to those categories
      id <- phbm[row.names(phbm) %in% i, ]
      
      # count the number of patients with the variant and the HPO
      s <- sum(colSums(id) == 3) 
      if(is.integer(s)){
        mat_long <- rbind(mat_long, c(fn,
                                      paste(g, 1, sep = "_"),
                                      c,
                                      s))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients without the variant but the HPO
      s1 <- sum(id[g, ] == 1 &
                  id[c, ] == 1 &
                  id[fn, ] == 0)
      if(is.integer(s1)){
        mat_long <- rbind(mat_long,
                          c(paste("without", i[1], sep =""),
                            paste(g, 1, sep = "_"),
                            c,
                            s1))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients with the variant but not the HPO
      s2 <- sum(id[g, ] == 0 &
                  id[c, ] == 1 &
                  id[fn, ] == 1)
      if(is.integer(s2)){
        mat_long <- rbind(mat_long, c(fn,
                                      paste(g, 0, sep = "_"),
                                      c,
                                      s2))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
      }
      
      # count the number of patients without the variant and without the HPO
      s3 <- sum(id[g, ] == 0 &
                  id[c, ] == 1 &
                  id[fn, ] == 0)
      if(is.integer(s1)){
        mat_long <- rbind(mat_long,
                          c(paste("without", i[1], sep =""),
                            paste(g, 0, sep = "_"),
                            c,
                            s3))
      } else {
        mat_long <- rbind(mat_long, rep(0, 4))
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
    hpo_ctw <- xtabs(
      as.formula(paste("count ~ phenotype +", fn, "+ centre")),
      data = mat_long
    )
    
    # make a summary table
    hpo_ctwo <- xtabs(
      as.formula(paste("count ~ phenotype +", fn)),
      data = mat_long
    )
    
    # test for association
    hpo_result <- chisq.test(hpo_ctwo)
    cmh_p <- hpo_result$p.value
    thisname <- paste(fn, g, sep = "_")
    names(cmh_p) <- thisname
    # combine p vec of all genetic tests (glild and HPOs)
    pvec <- c(pvec, cmh_p)
    # save the contingency table
    c_list[[thisname]] <- hpo_ctw
    # save the summary contingency table
    s_list[[thisname]] <- hpo_ctwo
  } 
  
  # some results are NaN, because the contingency tables are so imbalanced
  p_list[[k]] <- pvec[!is.nan(pvec)]
  names(p_list)[k] <- fn
}

# adjustment for multiple comparison Benjamini and Hochberg
BH_genetic_pvalue_result <- p.adjust(unlist(p_list),
                                     method = "BH")
biological_measure <- gsub("\\..*$",
             names(BH_genetic_pvalue_result[BH_genetic_pvalue_result < 0.05]),
             replacement = "")
phenotype_code <- gsub("^.*_",
           names(BH_genetic_pvalue_result[BH_genetic_pvalue_result < 0.05]),
           replacement = "")
phenotype_name <- names(unique_codes)[match(phenotype_code, unique_codes)]
BH_adjust_pvalue <- 
  unname(BH_genetic_pvalue_result[BH_genetic_pvalue_result < 0.05])

BH_genetic_CHM_table <- data.frame(
  biological_measure,
  phenotype_code,
  phenotype_name,
  BH_adjust_pvalue)

# human friendly p value. 
BH_genetic_CHM_table$BH_adjust_pvalue <- 
  formatC(BH_genetic_CHM_table$BH_adjust_pvalue,
          format = "g")

# qualify result
c1 <- matrix(NA, ncol = 3)
df <- BH_genetic_CHM_table
rl <- s_list
for(i in 1:nrow(df)){
  thisname <- paste(df[i, 1],
                    df[i, 2], sep = "_")
  c2 <- as.data.table(rl[[thisname]])
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
  summarise(patients_with_hpo = 
              sum(patient_count[str_ends(HPO_presence1_absence0, "1")]),
            total_patients_in_group = sum(patient_count),
            .groups = 'drop') %>%
  
  # get the percentage based on the group totals
  mutate(percentage_with_hpo = 
           (patients_with_hpo / total_patients_in_group) * 100
  ) 


# merge info
merged_tibble <- result %>%
  full_join(BH_genetic_CHM_table, by = c("HPO_code" = "phenotype_code"))

############## save all
fwrite(merged_tibble, "../result/genetic_HPO_tests_ratio_pvalue.csv")
