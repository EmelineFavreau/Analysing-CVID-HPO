######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# Import training data
# name the file
downloadname <- "../input/pheno-capture-data-2026-7-24.tsv"

# make a dataframe 
nbr0 <- read_tsv(downloadname)

######################## analysis ##############################################

# aim: compare semantic similarity between clinicians' HPO assessment
# 5 vignette: each represents a different patient with symptoms
# 5 clinicians: each produced a set of HPO terms

# tidy STUDY_ID
cases <- list(c1 = c("1LC", "1SG", "1JP", "1SE", "1JB"),
              c2 = c("2LC", "2SG", "2JP", "2SE", "2JB"),
              c3 = c("3LC", "3SG", "3JP", "3SE", "3JB"),
              c4 = c("4LC", "4SG", "JP4", "4SE", "4JB"),
              c5 = c("5LC", "5SG", "5JP", "5SE", "5JB")
)

# subset to those cases and their phenotypes
nbr0 <- nbr0 |> dplyr::filter(STUDY_ID %in% unlist(cases)) |>
  dplyr::select(c("STUDY_ID", "PHENOTYPES"))

# tidy phenotypes
patient_hpo_terms <- list()
# for each patient with B cell data
for(i in 1:nrow(nbr0)) {
  
  # assign all HPO related to this patient, to the list's element
  patient_hpo_terms[[i]] <- as.vector(unname(
    sapply(
      strsplit(nbr0$PHENOTYPES[i], split = "; "),
      substring, first = 6, last = 15
    )
  ))
}

# assign a new column, with the HPO terms
nbr0$hpo <- patient_hpo_terms

# analysis: pairwise semantic similarity between clinician, within vignettes
dat <- tibble(vignette = as.character(),
              similarity = as.numeric())

# loop for each case
for(j in 1:length(cases)){

  # subset to the case
  nbr <- nbr0 |> dplyr::filter(STUDY_ID %in% unlist(cases[j])) 
  
  # number of clinicians
  cn <- nrow(nbr)
  
  # list of HPO code per training
  hpolist <- nbr$hpo
  
  # vector of HPO codes
  all_hpo_vec <- unlist(hpolist)
  
  # vector of unique HPO codes
  unique_hpo_vec <- unique(all_hpo_vec) 

  # name of clinician in the rows
  # name of HPO in the columns unique_hpo_vec
  # 1 for presence, 0 for absence
  HPO_pres_abs_mat <- matrix(NA, nrow = cn,
                             ncol = length(unique_hpo_vec))
  colnames(HPO_pres_abs_mat) <- unique_hpo_vec
  rownames(HPO_pres_abs_mat) <- nbr$STUDY_ID
  
  # for each clinician c
  for(c in 1:cn){
    # for each HPO h
    for(h in 1:length(unique_hpo_vec)){
      # add a 1 if phenotype present
      HPO_pres_abs_mat[c, h] <- 
        ifelse(unique_hpo_vec[h] %in%
                 unlist(nbr$hpo[c]),
               1, 0)
    }
  }

  # obtain HPO code from the name
  unique_hpo_names_vec <- lapply(unique_hpo_vec, 
                      function(x) unname(hpo$name[which(hpo$id %in% x)]))
  
  # for each clinician, get all phenotypes they annotated
  clinician_scores <- list()
  for(i in 1:cn){
    tf = HPO_pres_abs_mat[i,] == 1
    clinician_scores[[i]] <- colnames(HPO_pres_abs_mat)[tf]
  }
  
  # obtain minimal sets for each clinician score
  minimal_term_sets <- list()
  for(i in 1:cn){
    minimal_term_sets[[i]] <- minimal_set(hpo, clinician_scores[[i]])
  }

  # calculate a similarity matrix, containing pairwise term-set similarities:
  information_content <- descendants_IC(hpo)
  sim_mat <- get_sim_grid(ontology = hpo, term_sets = minimal_term_sets)
  
  # keep informative data
  sim_mat[lower.tri(sim_mat, diag = T)] <- NA
  
  # row names have anonymised clinicians names
  rownames(sim_mat) <- nbr$STUDY_ID
  
  # a vector of similarities
  sim_vec <- as.vector(sim_mat)
  
  # two columns: the name of the case, the similarity value
  dat0 <- tibble(vignette = as.character(j),
                similarity = sim_vec[!is.na(sim_vec)])
  
  # concatenate similarity values for all cases
  dat <- dat %>% add_row(dat0)
  
  
}

# similarity values grouped by cases, 
# obtain summary stats
stat_sum <- dat %>%
  group_by(vignette) %>%
  reframe(as.list(summary(similarity)) %>% as_tibble())

# plot: dot plot five vignettes on the x axis, 
# the Lin similarity on the y axis and 
# each dot representing a pair of clinicians
ggplot(dat, aes(x = vignette, y = similarity)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = vignette)) +
  ylim(0,1) + theme(legend.position = "none")

############# save figure, stats ###############################################
ggsave("../result/SupplFigE4/FigE4.tiff",
       width = 6,
       height = 4,
       units = "in",
       dpi = 1000)

ggsave("../result/SupplFigE4/FigE4.jpg",
       width = 6,
       height = 4,
       units = "in",
       dpi = 1000)

write_csv(stat_sum, "../result/SupplFigE4/FigE4_summary_stats.csv")
