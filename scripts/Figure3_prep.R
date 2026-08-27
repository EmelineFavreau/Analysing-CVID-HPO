######################## source common libraries ###############################
source("common.R")
library(ggsignif)
library(cluster)
library(ggpubr)
library(ggrepel)
library(vegan)
library(ggtext)

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"


######################## analyse #########################

# obtain minimal sets for each patient
minimal_term_sets <- list()
for(i in 1:nrow(nbr)){
  minimal_term_sets[[i]] <- minimal_set(hpo, nbr$hpo[[i]])
}

names(minimal_term_sets) <- nbr$STUDY_ID

######################## Step 1: calculate similarities ########################
# calculate a similarity matrix, containing pairwise term-set similarities:
calc_sim_mat <- function(patient){
  # subset to the group of patients
  minimal_term_sets <- minimal_term_sets[names(minimal_term_sets) %in% patient]
  
  # calculate IC from large corpus
  information_content <- descendants_IC(hpo)
  
  # pairwise Lin similarity 
  sim_mat <- get_sim_grid(ontology = hpo, term_sets = minimal_term_sets)
  
  return(sim_mat)
}

# similarity for subgroups
infectionOnly_sim <- calc_sim_mat(infection_cluster$STUDY_ID)
complexOnly_sim <- calc_sim_mat(complex_cluster$STUDY_ID)
infectionVsComplex_sim <- calc_sim_mat(c(infection_cluster$STUDY_ID,
                                         complex_cluster$STUDY_ID))
all_sim <- calc_sim_mat(c(infection_cluster$STUDY_ID,
                                         complex_cluster$STUDY_ID))

# subset infectionVsComplex_sim to antagonist pairs
infectionVsComplex_simS <- 
  infectionVsComplex_sim[rownames(infectionVsComplex_sim) %in%
                           complex_cluster$STUDY_ID,
          colnames(infectionVsComplex_sim) %in% infection_cluster$STUDY_ID]

# save
saveRDS(infectionOnly_sim, "infectionOnly_sim.RDS")
saveRDS(complexOnly_sim, "complexOnly_sim.RDS")
saveRDS(infectionVsComplex_simS, "infectionVsComplex_simS.RDS")
saveRDS(infectionVsComplex_sim, "infectionVsComplex_sim.RDS")
saveRDS(all_sim, "all_sim.RDS")

######################## Step 2: prep violin plot ##############################
# plot a violin with x is the groups, y is Lin's similarity
inf <- as.vector(infectionOnly_sim[upper.tri(infectionOnly_sim)])
comp <- as.vector(complexOnly_sim[upper.tri(complexOnly_sim)])
cross <- as.vector(infectionVsComplex_simS)
all <- as.vector(all_sim[upper.tri(all_sim)])

# save 
saveRDS(inf, "inf.RDS")
saveRDS(comp, "comp.RDS")
saveRDS(cross, "cross.RDS")
saveRDS(all, "all.RDS")

# make a df
similarity_df <- tibble(
  Similarity = c(inf, comp, cross),
  Group = c(rep("Within-Infection", times = length(inf)),
            rep("Within-Complex", times = length(comp)),
            rep("Between-Groups", times = length(cross))))

# define the fixed group order
group_order <- c("Within-Infection", "Within-Complex", "Between-Groups")

# create groups and pairs
df <- similarity_df
groups <- unique(df$Group)
pairs <- combn(groups, 2, simplify = FALSE)

results <- lapply(pairs, function(pair) {
  x <- df$Similarity[df$Group == pair[1]]
  y <- df$Similarity[df$Group == pair[2]]
  test <- ks.test(x, y)
  data.frame(
    group1    = pair[1],
    group2    = pair[2],
    statistic = test$statistic,
    p_value   = test$p.value,
    pooled    = FALSE
  )
})
ks_results <- bind_rows(results)

# also between and within
pooled_test <- ks.test(cross, c(inf, comp))

pooled_row <- data.frame(
  group1    = "Within-Infection + Within-Complex (pooled)",
  group2    = "Between-Groups",
  statistic = pooled_test$statistic,
  p_value   = pooled_test$p.value,
  pooled    = TRUE
)

# combine results
ks_results <- bind_rows(ks_results, pooled_row)

saveRDS(ks_results, "ks_results.RDS")

####### permutation of labels in KS test #######################################
# make a dataset STUDYID, group(inf or comp)
input_df <- tibble(STUDY_ID = c(infection_cluster$STUDY_ID, complex_cluster$STUDY_ID),
                   group = c(rep("infection", times = nrow(infection_cluster)),
                             rep("complex", times = nrow(complex_cluster))))
random.seed(42)

# number of permutation
num_perm = 1000

# make output dataset
perm_KS <- matrix(NA_real_,
                  nrow = num_perm,
                  ncol = 4,
                  dimnames = list(
                    NULL,
                    c("inf_vs_comp",
                      "inf_vs_cross",
                      "comp_vs_cross",
                      "pooled")
                  ))

# shuffle the labels
for (perm in seq_len(num_perm)){
  input_dfp <- tibble(STUDY_ID = input_df$STUDY_ID,
                      group = sample(input_df$group))
  
  # similarity for subgroups
  inf_STUDY_IDp <- input_dfp$STUDY_ID[input_dfp$group == "infection"]
  comp_STUDY_IDp <- input_dfp$STUDY_ID[input_dfp$group == "complex"]
  infectionOnly_sim <- calc_sim_mat(inf_STUDY_IDp)
  complexOnly_sim <- calc_sim_mat(comp_STUDY_IDp)
  infectionVsComplex_sim <- calc_sim_mat(c(inf_STUDY_IDp,
                                           comp_STUDY_IDp))
  
  # subset infectionVsComplex_sim to antagonist pairs
  infectionVsComplex_simS <- 
    infectionVsComplex_sim[rownames(infectionVsComplex_sim) %in%
                             comp_STUDY_IDp,
                           colnames(infectionVsComplex_sim) %in% inf_STUDY_IDp]
  
  # Lin's similarity
  inf <- as.vector(infectionOnly_sim[upper.tri(infectionOnly_sim)])
  comp <- as.vector(complexOnly_sim[upper.tri(complexOnly_sim)])
  cross <- as.vector(infectionVsComplex_simS)
  
  # input for ks test
  similarity_dfp <- tibble(
    Similarity = c(inf, comp, cross),
    Group = c(rep("Within-Infection", times = length(inf)),
              rep("Within-Complex", times = length(comp)),
              rep("Between-Groups", times = length(cross))))
  
  # define the fixed group order
  group_order <- c("Within-Infection", "Within-Complex", "Between-Groups")
  
  # create groups and pairs
  groups <- unique(similarity_dfp$Group)
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # KS test
  resultsp <- lapply(pairs, function(pair) {
    x <- similarity_dfp$Similarity[similarity_dfp$Group == pair[1]]
    y <- similarity_dfp$Similarity[similarity_dfp$Group == pair[2]]
    test <- ks.test(x, y)
    data.frame(
      group1    = pair[1],
      group2    = pair[2],
      statistic = test$statistic,
      p_value   = test$p.value,
      pooled    = FALSE
    )
  })
  
  # combine rows
  ks_resultsp <- bind_rows(resultsp)
  
  # also between and within
  pooled_test <- ks.test(cross, c(inf, comp))
  
  pooled_row <- data.frame(
    group1    = "Within-Infection + Within-Complex (pooled)",
    group2    = "Between-Groups",
    statistic = pooled_test$statistic,
    p_value   = pooled_test$p.value,
    pooled    = TRUE
  )
  
  # combine results
  ks_resultsp <- bind_rows(ks_resultsp, pooled_row)
  
  # save D statistics for each test, to create 4 vectors
  perm_KS[perm,] <- ks_resultsp$statistic
}
#### end of loop

# set p value col
ks_results$perm_p_val <- 0

# for each pairwise test
for (i in 1:4){
  ks_results$perm_p_val[i] <- mean(perm_KS[,i] >= ks_results$statistic[i])
}

# save
saveRDS(ks_results, "ks_results_perm.RDS")
saveRDS(perm_KS, "perm_KS.RDS")

####### prepare figure component #######################################
# update P value for plotting
ks_results <- readRDS("ks_results_perm.RDS")

# tidy df
ks_results <- ks_results %>%
  dplyr::mutate(label = case_when(
    perm_p_val == 0 ~ "p < 0.001",
    perm_p_val > 0 ~ paste0("p = ", perm_p_val)
    ))

# set the coordinates for plot
y_max   <- max(df$Similarity, na.rm = TRUE)
y_range <- diff(range(df$Similarity, na.rm = TRUE))
step    <- y_range * 0.12

# to plot Ks result in situ
get_span <- function(g1, g2, order) {
  abs(match(g1, order) - match(g2, order))
}

# KS values
pairwise_ks <- ks_results %>%
  dplyr::filter(!pooled) %>%
  dplyr::mutate(span = mapply(get_span, group1, group2,
                       MoreArgs = list(order = group_order))) %>%
  dplyr::arrange(span)

# create annotations
annotations <- mapply(c,
                      pairwise_ks$group1,
                      pairwise_ks$group2,
                      SIMPLIFY = FALSE)

# position on y axis
y_positions <- seq(y_max + step, by = step, length.out = nrow(pairwise_ks))

# make fig label
pooled_label <- paste0(
  "Between-Groups vs. pooled (Within-Infection + Within-Complex): p = ",
  ks_results$perm_p_val[ks_results$pooled]
)

# vector of colours
infCompBetween <- c(infComp, "Grey")
names(infCompBetween) <- c("Within-Complex", "Within-Infection", "Between-Groups")

# save all components for plotting
saveRDS(similarity_df, "similarity_df.RDS")
saveRDS(annotations, "annotations.RDS")
saveRDS(pairwise_ks, "pairwise_ks.RDS")
saveRDS(y_positions, "y_positions.RDS")
saveRDS(pooled_label, "pooled_label.RDS")
saveRDS(infCompBetween, "infCompBetween.RDS")
saveRDS(group_order, "group_order.RDS")
saveRDS(step, "step.RDS")


######################## Step 3: prep qq plot ##################################
pooled <- c(inf, comp)
n      <- min(length(cross), length(pooled))
q_cross  <- quantile(cross,  probs = seq(0, 1, length.out = n))
q_pooled <- quantile(pooled, probs = seq(0, 1, length.out = n))

qq_df <- data.frame(y = q_pooled, x = q_cross)

# save
saveRDS(qq_df, "qq_df.RDS")

######################## Step 4: MDS ########################################
# take the similarity matrix and plot another type of dimensionality reduction
infectionVsComplex_sim <- readRDS("infectionVsComplex_sim.RDS")
# HPO similarity matrix: infection vs complex patients
ic <- infectionVsComplex_sim

# Perform MDS analysis on converted dist class object
mds_ic <- cmdscale(as.dist(1-ic))

# Add cluster information to the MDS results
mds_df <- as.data.frame(mds_ic)
mds_df$groups <- ifelse(rownames(mds_df) %in% infection_cluster$STUDY_ID, 
                         "Infection", "Complex")
mds_df$centres <- nbr$CENTRE[match(rownames(mds_df), nbr$STUDY_ID) ]

# save
saveRDS(mds_df, "mds_df.RDS")

# colours
InfComp <- infComp
names(InfComp) <- c("Complex", "Infection")

# save
saveRDS(InfComp, "InfComp.RDS")



######################## Step 7: permutation ###################################
# take the current mds, permute the groups (ie MDS coordinates do not move) 
# and calculate the variance between each group. 
# repeat 999 times to test the null hypothesis of no variance between groups

set.seed(42)

# from: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
permanova_result <- adonis2(mds_df[, c("V1", "V2")] ~ groups, 
        data = mds_df, 
        permutations = 999,
        method = "euclidean")

