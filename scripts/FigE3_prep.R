######################## source common libraries ###############################
library(plyr)
source("common.R")
library(ggsignif)
library(cluster)
library(ggpubr)
library(ggrepel)
library(vegan)
library(ggtext)
# downloaded from https://sourceforge.net/projects/hposim/
#path_to_file = "HPO.db_1.9.tar.gz"
#install.packages(path_to_file, repos = NULL, type = "source")
library(HPO.db)
library(igraph)
# downloaded from https://cran.r-project.org/src/contrib/Archive/HPOSim/
#path_to_file = "HPOSim_1.2.tar.gz"
#install.packages(path_to_file, repos = NULL, type = "source")
library(HPOSim)
library(furrr)
library(future)
library(patchwork)

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"

# get current annotations once (uncomment if needed)
#download.file(
#  "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa",
#  destfile = "phenotype.hpoa"
#)
# read disease annotations
annot <- read.delim("../input/phenotype.hpoa", comment.char = "#")

############ calculate Deng similarities ######################################


####### Step 1: update the environment #########################################

# compare how the Lin values compare to those of an alternative semantic 
# similarity metric that is less sensitive to the choice of corpus, 
# such as simGraphIC (see 10.1371/journal.pone.0115692).

## this requires updating the environment to recent HPO download

# from recent HPO download
active_terms <- hpo$id

# build ancestor list that HPOSim expects
Ancestors <- hpo$ancestors
names(Ancestors) <- active_terms

# count number of diseases for each term
term_counts <- table(annot$hpo_id)

# vector of 0, named with HP codes
all_counts <- setNames(numeric(length(active_terms)), active_terms)

for(term in active_terms){
  # find descendants of given HPO term
  desc <- names(which(sapply(hpo$ancestors, function(x) term %in% x)))
  # disease frequency for all descendants of given HPO term
  all_counts[term] <- sum(term_counts[names(term_counts) %in% desc],
                          na.rm = TRUE)
}

# convert to Information Content (IC)
# -log(freq / max_freq)
# avoid -Inf by adding a small value to log(0)
total <- max(all_counts)
IC_values <- -log(all_counts/total + 1e-10)  

# Build IC as a 3-column data frame matching what HPOSim::GIC expects
IC <- data.frame(
  term    = names(IC_values),
  count   = as.numeric(all_counts),
  IC      = as.numeric(IC_values),
  stringsAsFactors = FALSE
)

# create a new environment
HPOSimEnv <- new.env(parent = emptyenv())

# create term IC and parents
termIC <- setNames(hpo$name[active_terms], active_terms)
Parents <- hpo$parents[active_terms]

# assign all
assign("Ancestors", Ancestors, envir = HPOSimEnv)
assign("termIC",    IC,        envir = HPOSimEnv)
assign("termNames", setNames(hpo$name[active_terms], active_terms), envir = HPOSimEnv)
assign("Parents",   hpo$parents[active_terms],                      envir = HPOSimEnv)

saveRDS(Ancestors, "../result/Ancestors.RDS")
saveRDS(IC, "../result/IC.RDS")
saveRDS(termIC, "../result/termIC.RDS")
saveRDS(Parents, "../result/Parents.RDS")


my_initialize <- function() {
  # Prevent re-initialisation if already loaded
  if (exists("Ancestors", envir = HPOSimEnv, inherits = FALSE)) return(invisible(NULL))
  
  # Load
  Ancestors  <- readRDS("../result/Ancestors.RDS")
  IC         <- readRDS("../result/IC.RDS")
  termIC     <- readRDS("../result/termIC.RDS")
  Parents    <- readRDS("../result/Parents.RDS")
  
  assign("Ancestors", Ancestors, envir = HPOSimEnv)
  assign("termIC",    IC,        envir = HPOSimEnv)
  assign("termNames", termNames, envir = HPOSimEnv)
  assign("Parents", Parents, envir = HPOSimEnv)
  
  invisible(NULL)
}

# Inject into HPOSim's namespace
assignInNamespace(
  x     = ".initialize",
  value = my_initialize,
  ns    = "HPOSim"
)

####### Step 2: calculate GIC similarities #####################################

# vector of IC terms
build_ic_vec <- function(IC) {
  setNames(IC[, 3], IC[, 1])
}

# function to calculate GIC similarity between two terms
getGIC2_fast <- function(term1, term2, ic_vec, ancestor) {
  if (term1 == term2) return(1)
  an1 <- ancestor[[term1]]
  an2 <- ancestor[[term2]]
  ancommon <- intersect(an1, an2)
  if (length(ancommon) == 0) return(0)
  a <- sum(ic_vec[ancommon], na.rm = TRUE)
  b <- sum(ic_vec[base::union(an1, an2)], na.rm = TRUE)
  if (b == 0) return(0)
  a / b
}

# function to apply method
getTermListSim_GIC_funSimMax <- function(anno1, anno2, ic_vec, ancestor) {
  if (length(anno1) * length(anno2) == 0) {
    warning(paste("No HPO information for", anno1, anno2,
                  ". Similarity set to NaN."))
    return(NaN)
  }
  # kernel matrix of pairwise GIC scores
  ker <- outer(
    anno1, anno2,
    FUN = Vectorize(function(t1, t2)
      getGIC2_fast(t1, t2, ic_vec, ancestor))
  )
  # funSimMax: max of (mean row-maxima, mean col-maxima)
  max(mean(apply(ker, 1, max)),
      mean(apply(ker, 2, max)))
}

# function to make sim mat
calc_GICsim_mat_fast <- function(patient, minimal_term_sets, IC) {
  ic_vec   <- build_ic_vec(IC)
  .initialize() # this should use the updated environment
  ancestor <- get("Ancestors", envir = HPOSimEnv)
  minimal_term_sets <- minimal_term_sets[patient]
  n       <- length(patient)
  sim_mat <- matrix(NA_real_, nrow = n, ncol = n,
                    dimnames = list(patient, patient))
  diag(sim_mat) <- 1
  pairs       <- which(upper.tri(sim_mat), arr.ind = TRUE)
  row_indices <- split(seq_len(nrow(pairs)), pairs[, 1])
  results <- furrr::future_map(
    row_indices,
    function(idx) {
      i     <- pairs[idx[1], 1]
      js    <- pairs[idx, 2]
      anno1 <- minimal_term_sets[[patient[i]]]
      if (length(anno1) == 0) return(rep(NA_real_, length(js)))
      
      vapply(js, function(j) {
        anno2 <- minimal_term_sets[[patient[j]]]
        if (length(anno2) == 0) return(NA_real_)
        tryCatch(
          getTermListSim_GIC_funSimMax(anno1, anno2, ic_vec, ancestor),
          error = function(e) {
            warning(sprintf("Failed for pair (%s, %s): %s",
                            patient[i], patient[j], e$message))
            NA_real_
          }
        )
      }, FUN.VALUE = numeric(1))
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  for (grp in seq_along(row_indices)) {
    i  <- pairs[row_indices[[grp]][1], 1]
    js <- pairs[row_indices[[grp]], 2]
    for (k in seq_along(js)) sim_mat[i, js[k]] <- results[[grp]][k]
  }
  sim_mat
}

# remove redundant terms in matrix
minimal_term_sets <- setNames(
  lapply(seq_len(nrow(nbr)), function(i) minimal_set(hpo, nbr$hpo[[i]])),
  nbr$STUDY_ID
)

# Deng calculation is intensive, so make use of parallel computing
# spins up a pool of background R processes
plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))

infectionOnly_GICsim <- calc_GICsim_mat_fast(
  patient           = infection_cluster$STUDY_ID,
  minimal_term_sets = minimal_term_sets,
  IC                = IC
)

# close the multisession, R will run normally
plan(sequential)

# spins up a pool of background R processes
plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))

complexOnly_GICsim <- calc_GICsim_mat_fast(
  patient           = complex_cluster$STUDY_ID,
  minimal_term_sets = minimal_term_sets,
  IC                = IC
)

# close the multisession, R will run normally
plan(sequential)

# spins up a pool of background R processes
plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))

cross_GICsim <- calc_GICsim_mat_fast(
  patient           = c(infection_cluster$STUDY_ID,complex_cluster$STUDY_ID),
  minimal_term_sets = minimal_term_sets,
  IC                = IC
)
all_GICsim <- cross_GICsim
# close the multisession, R will run normally
plan(sequential)

# subset cross_GICsim to antagonist pairs
cross_GICsim <- cross_GICsim[rownames(cross_GICsim) %in% 
                               complex_cluster$STUDY_ID,
                             colnames(cross_GICsim) %in% infection_cluster$STUDY_ID]

# save object
saveRDS(infectionOnly_GICsim, "../result/infectionOnly_GICsim.RDS")
saveRDS(complexOnly_GICsim, "../result/complexOnly_GICsim.RDS")
saveRDS(cross_GICsim, "../result/cross_GICsim.RDS")
saveRDS(all_GICsim, "../result/all_GICsim.RDS")