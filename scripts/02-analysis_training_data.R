######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# training data: same patient, HPO codes by 10 different clinicians
downloadname <- "../input/pheno-capture-data-all-2022-10-5.tsv"

# make a dataframe from training data
training <- read_tsv(downloadname)

# pre-training data: many patients with HPO codes by untrained clinicians
#2021-11-16: only one centre
t1 <- "../input/pheno-capture-data-2021-11-16.tsv"
d1 <- read_tsv(t1)

# pre-training data: many patients with HPO codes by untrained clinicians
# (except for d1's centre with trained clinicians)
t2 <- "../input/pheno-capture-data-all-2022-9-26.tsv"
d2 <- read_tsv(t2)

# post-training data: many patients with HPO codes by trained clinicians
d3 <- readRDS("../result/tidy_data")

######################## HPO Similarity HPO between training sets ##############
# remove incomplete data
training <- training[!is.na(training$HPO_CODE),]

# keep records of the same patient
training <- training[grep(pattern = ".*atient1.*", x = training$STUDY_ID),]
training <- training[grep(pattern = ".*1E|.*1n|.*1o", x = training$STUDY_ID,
                          invert = T),]

# vector of anonymised names of same patient
patient1_vec <- training$STUDY_ID[grep(pattern = ".*atient1.*",
                                       x = training$STUDY_ID)]

# anonymise the clinicians
training$ClinicianAnon <- 1:10

# list of HPO code per training
patient1_hpolist <- training$HPO_CODE

# vector of HPO codes
all_hpo_vec <- unlist(strsplit(patient1_hpolist, split = "; "))

# vector of unique HPO codes
unique_hpo_vec <- unique(all_hpo_vec) 


# name of clinician in the rows training$ClinicianAnon
# name of HPO in the columns unique_hpo_vec
# 1 for presence, 0 for absence
HPO_pres_abs_mat <- matrix(NA, nrow = nrow(training),
                           ncol = length(unique_hpo_vec))
colnames(HPO_pres_abs_mat) <- unique_hpo_vec
rownames(HPO_pres_abs_mat) <- training$ClinicianAnon

# for each clinician c
for(c in 1:nrow(training)){
  # for each HPO h
  for(h in 1:length(unique_hpo_vec)){
    # add a 1 if present
    HPO_pres_abs_mat[c, h] <- 
      ifelse(unique_hpo_vec[h] %in%
               unlist(strsplit(training$HPO_CODE[c],
                               split = "; ")),
             1, 0)
  }
}

# obtain HPO code from the name
unique_hpo_names_vec <- lapply(unique_hpo_vec, 
           function(x) unname(hpo$name[which(hpo$id %in% x)]))

# 2 columns: HPO names, HPO codes
HPO_name_code <- cbind(unique_hpo_vec, unique_hpo_names_vec)


# frequency of each HPO
HPO_df <- data.frame(HPO_code = unique_hpo_vec,
                     HPO_name = unlist(unique_hpo_names_vec),
                     num_times_coded = colSums(HPO_pres_abs_mat))


# obtain clinician scores
clinician_scores <- list()
for(i in 1:10){
  clinician_scores[[i]] <- names(HPO_pres_abs_mat[i,][HPO_pres_abs_mat[i,] == 1])
}

# obtain minimal sets for each clinician score
minimal_term_sets <- list()
for(i in 1:10){
  minimal_term_sets[[i]] <- minimal_set(hpo, clinician_scores[[i]])
}


# calculate a similarity matrix, containing pairwise term-set similarities:
information_content <- descendants_IC(hpo)
sim_mat <- get_sim_grid(ontology = hpo, term_sets = minimal_term_sets)

# keep informative data
sim_mat[lower.tri(sim_mat, diag = T)] <- NA

# row names have anonymised clinicians names
rownames(sim_mat) <- 1:10

fwrite(sim_mat, "../result/Fig1/sim_mat.csv")


######################## HPO Similarity HPO pre/post training ##################

##### task 1: combine two dataframes

# keep study_id and HPO code 
d1f <- d1[, c("STUDY_ID", "HPO_CODE")]
d2f <- d2[, c("STUDY_ID", "HPO_CODE")]

# remove RFH from d2f, the focus data are in d1f
d2f <- d2f[!d2f$STUDY_ID %in% d1f$STUDY_ID, ]

# combine all info
before_training <- rbind(d1f, d2f)

# remove NA
before_training <- before_training[!is.na(before_training$HPO_CODE), ]


#### task 2: subset after training for the same study_id
# keep "STUDY_ID" and "hpo"  
d3f <- d3[d3$STUDY_ID %in% before_training$STUDY_ID, c("STUDY_ID", "hpo" )]

# remove NA
after_training <- d3f[!is.na(d3f$hpo), ] 

# for each study_id, calculate number of HPO before and after training
num_hpo_comp <- data.frame(STUDY_ID = after_training$STUDY_ID,
                           num_hpo_before = NA,
                           num_hpo_after = NA)

# calculate number of HPO before and after
for(i in 1:nrow(after_training)){
  p <- after_training$STUDY_ID[i]
  num_hpo_after <- length(unlist(after_training$hpo[i]))
  num_hpo_before <-length(unlist(strsplit(
    before_training$HPO_CODE[before_training$STUDY_ID == p], split = "; ")))
  num_hpo_comp[i, 2:3] <- c(num_hpo_before, num_hpo_after)
}

#### task 3: plot distribution
# transform out data into long
num_hpo_compL <- pivot_longer(data = num_hpo_comp, 
                              cols = c("num_hpo_before" , "num_hpo_after" ),
                              names_to = "Timepoint",
                              values_to = "Records",
                              names_prefix ="num_hpo_")

# assessing if distributions are drawn from same continuous 
after_record_vec <- num_hpo_compL %>% dplyr::filter(Timepoint == "after") %>% 
  dplyr::select(Records) %>% unlist() 
before_record_vec <- num_hpo_compL %>% dplyr::filter(Timepoint == "before") %>% 
  dplyr::select(Records) %>% unlist() 
ks.test(before_record_vec, after_record_vec) # p-value < 2.2e-16

fwrite(num_hpo_compL, "../result/Fig1/num_hpo_compL.csv")





#### task 4: Similarity analysis
# for a given patient, how similar are the hpo before and after?

# step 4.1: subset to common STUDY_ID 
common_STUDY_ID_vec <- before_training$STUDY_ID[before_training$STUDY_ID %in%
                                                  after_training$STUDY_ID]

# step 4.2: get terms set 
term_sets <- list()
for(i in 1:length(common_STUDY_ID_vec)){
  p <- common_STUDY_ID_vec[i]
  before_set <- strsplit(
    before_training$HPO_CODE[before_training$STUDY_ID == p], split = "; ")
  after_set <-  after_training$hpo[after_training$STUDY_ID == p] # a list
  term_sets <- c(term_sets,
                 before_set,
                 after_set)
}

# name the elements of the list
n <- c()

# names of list
for(i in 1:length(common_STUDY_ID_vec)){
  ni <- c(paste(common_STUDY_ID_vec[i], "_before", sep =""),
          paste(common_STUDY_ID_vec[i], "_after", sep =""))
  n <- c(n, ni)
}
names(term_sets) <- n

# call minimal_set on each term set to guarantee the 
# similarity expressions are faithfully evaluated
term_sets <- lapply(term_sets,
                    function(x) minimal_set(hpo, x))

# step 4.3: Calculate similarity matrix
sim_mat <- get_sim_grid(ontology = hpo,
                        term_sets = term_sets)

# extract the similarity value before and after for a given STUDY_ID
# every other row, starting from 2 (even)
important_rows <- seq(from = 2, to = nrow(sim_mat), by = 2)

# every other column, starting from 1 (odds)
important_cols <- seq(from = 1, to = nrow(sim_mat), by = 2)

# for a given STUDY_ID, get the similarity of HPO before and after training
similarity_value_vec <- c()
for(i in 1:length(common_STUDY_ID_vec)){
  before_after_sim <- sim_mat[important_rows[i], important_cols[i]]
  similarity_value_vec <- c(similarity_value_vec, before_after_sim)
}

# give STUDY_ID as name of vector of similarity
names(similarity_value_vec) <- common_STUDY_ID_vec

# step 4.4: calculate difference of number of HPO before and after training 
num_hpo_comp$differential <- num_hpo_comp$num_hpo_after - 
  num_hpo_comp$num_hpo_before

# combine results
index <- match(num_hpo_comp$STUDY_ID, common_STUDY_ID_vec)
num_hpo_comp$before_after_similarity_value <- similarity_value_vec[index]

nrow(num_hpo_comp)# for main text
sum(num_hpo_comp$differential < 0)# for main text
summary(num_hpo_comp$differential[num_hpo_comp$differential >0])# for main text
fwrite(num_hpo_comp, "../result/Fig1/num_hpo_comp.csv")

