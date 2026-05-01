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
#install.packages('igraph')
library(igraph)
# downloaded from https://cran.r-project.org/src/contrib/Archive/HPOSim/
#path_to_file = "HPOSim_1.2.tar.gz"
#install.packages(path_to_file, repos = NULL, type = "source")
library(HPOSim)
#install.packages("furrr")
#install.packages("future")
library(furrr)
library(future)

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"

# get current annotations
#download.file(
#  "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa",
#  destfile = "phenotype.hpoa"
#)
# read disease annotations
annot <- read.delim("phenotype.hpoa", comment.char = "#")

############################# Major Comment 2 ##################################


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

#saveRDS(Ancestors, "Ancestors.RDS")
#saveRDS(IC, "IC.RDS")
#saveRDS(termIC, "termIC.RDS")
#saveRDS(Parents, "Parents.RDS")


my_initialize <- function() {
  # Prevent re-initialisation if already loaded
  if (exists("Ancestors", envir = HPOSimEnv, inherits = FALSE)) return(invisible(NULL))
  
  # Load
  Ancestors  <- readRDS("Ancestors.RDS")
  IC         <- readRDS("IC.RDS")
  termIC     <- readRDS("termIC.RDS")
  Parents    <- readRDS("Parents.RDS")
  
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

# close the multisession, R will run normally
plan(sequential)

# subset cross_GICsim to antagonist pairs
cross_GICsim <- cross_GICsim[rownames(cross_GICsim) %in% 
                                          complex_cluster$STUDY_ID,
                      colnames(cross_GICsim) %in% infection_cluster$STUDY_ID]

####### Step 3: compare Lin to GIC #############################################

# save object
#saveRDS(infectionOnly_GICsim, "infectionOnly_GICsim.RDS")
#saveRDS(complexOnly_GICsim, "complexOnly_GICsim.RDS")
#saveRDS(cross_GICsim, "cross_GICsim.RDS")
# load Lin's similarities
inf <- readRDS("inf.RDS")
comp <- readRDS("comp.RDS")
cross <- readRDS("cross.RDS")
# Load GIC's similarities
infectionOnly_GICsim <- readRDS("infectionOnly_GICsim.RDS")
complexOnly_GICsim <- readRDS("complexOnly_GICsim.RDS")
cross_GICsim <- readRDS("cross_GICsim.RDS")

# keep only the upper triangle
infGIC <- as.vector(infectionOnly_GICsim[upper.tri(infectionOnly_GICsim)])
compGIC <- as.vector(complexOnly_GICsim[upper.tri(complexOnly_GICsim)])
crossGIC <- as.vector(cross_GICsim[upper.tri(cross_GICsim)])

# # exploring the pairs with similarity == 1
# result <- which(cross_GICsim == 1, arr.ind = TRUE)
# df <- data.frame(
#   row_name = rownames(cross_GICsim)[result[, "row"]],
#   col_name = colnames(cross_GICsim)[result[, "col"]]
# )
# 
# df1 <- data.frame(
#   row_name = rownames(complexOnly_GICsim)[result[, "row"]],
#   col_name = colnames(complexOnly_GICsim)[result[, "col"]]
# )
# 
# df <- df[df$row_name != df$col_name, ]
# unique(df$row_name)
# unique(df$col_name)


# compare within-group similarities:
ks.test(infGIC, inf) # different: D = 0.30851, p-value < 2.2e-16
ks.test(compGIC, comp) # different: D = 0.45818, p-value < 2.2e-16

jpeg("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig1.jpg")
qqplot(infGIC, inf) # GIC includes a lot of 1, why?
abline(a = 0, b = 1)
title("Fig1: Infection-only qqplot (KS test D = 0.30851, p-value < 2.2e-16)")
dev.off()

jpeg("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig2.jpg")
qqplot(compGIC,comp)
abline(a = 0, b = 1)
title("Fig2: Complex-only qqplot (KS test D = 0.45818, p-value < 2.2e-16)")
dev.off()

# compare between-group similarities:
jpeg("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig3.jpg")
ks.test(crossGIC, cross)  # different D = 0.36789, p-value < 2.2e-16
qqplot(crossGIC, cross)
abline(a = 0, b = 1)
title("Fig3: Complex-only qqplot (KS test D = 0.36789, p-value < 2.2e-16)")
dev.off()

####### Step 4: plots ##########################################################
crossGIC <- crossGIC[!is.na(crossGIC)]
dat <- data.frame(SimMethod = c(rep("Lin's", times = length(cross)),
                                rep("HPOSim::GraphIC",
                                    times = length(crossGIC))),
                 Similarity = c(cross, crossGIC))

cdat <- ddply(dat, "SimMethod", dplyr::summarise,
              similarity.mean=mean(Similarity))

compareSimFig <- ggplot(dat, aes(x = Similarity)) +
  geom_histogram(binwidth =.005, colour="grey", fill="white") + 
  facet_grid(SimMethod ~ .)+
  geom_vline(data=cdat,
             aes(xintercept=similarity.mean),
                      linetype="dashed", size=0.7, colour="black") +
  xlab("Infection vs Complex Pairwise Similarity") +
  ggtitle("Fig6: Lin vs simGraphIC between groups")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig6.jpg")
# 
# Infection term sets are on average less similar to complex term sets when 
# all ancestors are taken in account (GraphIC). 

infGIC <- infGIC[!is.na(infGIC)]
datInf = data.frame(SimMethod = c(rep("Lin's", times = length(inf)),
                               rep("HPOSim::GraphIC", times = length(infGIC))),
                 Similarity = c(inf, infGIC))


cdatInf <- ddply(datInf, "SimMethod", dplyr::summarise, similarity.mean=mean(Similarity))

compareInfSimFig <- ggplot(datInf, aes(x = Similarity)) +
  geom_histogram(binwidth =.005, colour="grey", fill="white") + 
  facet_grid(SimMethod ~ .)+
  geom_vline(data=cdatInf,
             aes(xintercept=similarity.mean),
             linetype="dashed", size=0.7, colour="black") +
  xlab("Infection only Pairwise Similarity") +
  ggtitle("Fig4: Lin vs simGraphIC within Infection")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig4.jpg")

compGIC <- compGIC[!is.na(compGIC)]
datComp = data.frame(SimMethod = c(rep("Lin's", times = length(comp)),
                                  rep("HPOSim::GraphIC", times = length(compGIC))),
                    Similarity = c(comp, compGIC))


cdatComp<- ddply(datComp, "SimMethod", dplyr::summarise, similarity.mean=mean(Similarity))

compareCompSimFig <- ggplot(datComp, aes(x = Similarity)) +
  geom_histogram(binwidth =.005, colour="grey", fill="white") + 
  facet_grid(SimMethod ~ .)+
  geom_vline(data=cdatComp,
             aes(xintercept=similarity.mean),
             linetype="dashed", size=0.7, colour="black") +
  xlab("Complex only Pairwise Similarity") +
  ggtitle("Fig5: Lin vs simGraphIC within Complex")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig5.jpg")

compareSimFig
compareInfSimFig
compareCompSimFig 
# similarity is higher for Lin's: could be explained by more diverse HPO terms.
# so more unique ancestor terms (demoninator of GIC ratio), thus reducing GIC value

