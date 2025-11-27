######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

cohort_size <- nrow(nbr)

#################### count patients in each group ##############################

#################### CD21 ######################################################

# size and proportion of records with CD21 low 
cd21_cohort_size <- nbr %>%
  dplyr::filter(!is.na(cd21_low_b_cells_pct)) %>% nrow()

cd21_cohort_prop <- round((cd21_cohort_size/cohort_size)*100, 2)

# subset, size, proportion of records with CD21 low "expansion of low CD21"
cd21low_cohort_list <- nbr %>%
  dplyr::filter(cd21_low_b_cells_pct >= 10) %>% 
  dplyr::select(STUDY_ID)

cd21low_cohort_size <- nrow(cd21low_cohort_list)

cd21low_cohort_prop <- round((cd21low_cohort_size/cd21_cohort_size)*100, 2)


# subset, size, proportion of records with CD21 low "normal CD21".
cd21normal_cohort_list <- nbr %>%
  dplyr::filter(cd21_low_b_cells_pct < 10) %>%
  dplyr::select(STUDY_ID)

cd21normal_cohort_size <- nrow(cd21normal_cohort_list)

cd21normal_cohort_prop <- 
  round((cd21normal_cohort_size/cd21_cohort_size)*100, 2)

#################### SmB ######################################################
# size and proportion of records with SmB. 
smb_cohort_size <- nbr %>%
  dplyr::filter(!is.na(class_switched_b_cells_pct)) %>%
  nrow()

smb_cohort_prop <- round((smb_cohort_size/cohort_size)*100, 2)

# subset, size, proportion of records with "SmB+"
smBplus_cohort_list <- nbr %>%
  dplyr::filter(class_switched_b_cells_pct > 2) %>% dplyr::select(STUDY_ID)

smBplus_cohort_size <- nrow(smBplus_cohort_list)

smBplus_cohort_prop <- round((smBplus_cohort_size/smb_cohort_size)*100, 2)

# subset, size, proportion of records with "SmB-"
smBminus_cohort_list <- nbr %>%
  dplyr::filter(class_switched_b_cells_pct <= 2) %>% dplyr::select(STUDY_ID)

smBminus_cohort_size <- nrow(smBminus_cohort_list)

smBminus_cohort_prop <- round((smBminus_cohort_size/smb_cohort_size)*100, 2)

######################### Tr ###################################################
# size and proportion of records with Tr. 
tr_cohort_size <- nbr %>%
  dplyr::filter(!is.na(transitional_b_cells_pct)) %>% nrow()

tr_cohort_prop <- round((tr_cohort_size/cohort_size)*100, 2)

# subset, size, proportion of records with "TrNorm"
TrNorm_cohort_list <- nbr %>%
  dplyr::filter(transitional_b_cells_pct < 9) %>% dplyr::select(STUDY_ID)

TrNorm_cohort_size <- nrow(TrNorm_cohort_list)

TrNorm_cohort_prop <- round((TrNorm_cohort_size/tr_cohort_size)*100, 2)

# subset, size, proportion of records with "TrHigh"
TrHigh_cohort_list <- nbr %>%
  dplyr::filter(transitional_b_cells_pct >= 9) %>% dplyr::select(STUDY_ID)

TrHigh_cohort_size <- nrow(TrHigh_cohort_list)

TrHigh_cohort_prop <- round((TrHigh_cohort_size/tr_cohort_size)*100, 2)


######################### Ig ###################################################

# Number of patients in each subcategory of Ig
ig_vec <- c("lowIgGlowIgAlowIgM","lowIgGlowIgAnormalIgM",
            "lowIgGlowIgAhighIgM","lowIgGnormalIgAlowIgM")

ig_cohort_size <- sum(nbr$Immuno_group %in% ig_vec )
ig_cohort_prop <- round((ig_cohort_size/cohort_size)*100, 2)


# subset, size, proportion of records with low levels of IgG, IgA, IgM.
lll_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[1]]
lll_cohort_prop <- round((lll_cohort_size/ig_cohort_size)*100, 2)


# subset, size, proportion of records with low levels of IgG and IgA,
#normal levels of IgM.
lln_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[2]]
lln_cohort_prop <- round((lln_cohort_size/ig_cohort_size)*100, 2)

# subset, size, proportion of records with  low levels of IgG and IgA,
#high levels of IgM.
llh_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[3]]
llh_cohort_prop <- round((llh_cohort_size/ig_cohort_size)*100, 2)

# subset, size, proportion of records with low levels of IgG and IgM,
#normal levels of IgA.
lnl_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[4]]
lnl_cohort_prop <- round((lnl_cohort_size/ig_cohort_size)*100, 2)

ig_size <- c(lll_cohort_size, lln_cohort_size,
             llh_cohort_size, lnl_cohort_size)

######################### Any pathogenic #######################################

# assign anyOtherPathogenic to records with a diagnosis gene
t <- nbr %>% 
  dplyr::filter(!is.na(Genetic_variant) | !is.na(`Diagnosis gene`)) %>% 
  dplyr::select(c(STUDY_ID, Genetic_variant, `Diagnosis gene`))

t <- t %>% dplyr::filter(!`Diagnosis gene` == "no genotype")

nbr$Genetic_variant[nbr$STUDY_ID %in% t$STUDY_ID[is.na(t$Genetic_variant)]] <-
  "anyOtherPathogenic"

# Number of patients in each subcategory of genetic variants 
# (IUIS genes with frequencies of variants)
# number of patients with WGS 
sequenced_cohort_size <- sum(nbr$sequenced == 1)
sequenced_cohort_prop <-  (sequenced_cohort_size/cohort_size)*100
# number of patients with confirmed genetic diagnosis
diagnosed_cohort_size <- sum(table(nbr$Genetic_variant))
diagnosed_cohort_prop <- (diagnosed_cohort_size/cohort_size)*100

# subset, size, proportion of records with any diagnostic variant
# (no NFKB1, no TNFRSF13B)
any_diagnosis_cohort_size <- unname(table(nbr$Genetic_variant)[
  names(table(nbr$Genetic_variant)) == "anyOtherPathogenic"])

any_diagnosis_cohort_prop <- 
  round((any_diagnosis_cohort_size/diagnosed_cohort_size)*100, 2)

any_diagnosis_sequenced_cohort_prop <- 
  round((any_diagnosis_cohort_size/sequenced_cohort_size)*100, 2)


######################### NFKB1 ################################################
# subset, size, proportion of records with NFKB1
NFKB1_cohort_size <- sum(table(nbr$Genetic_variant)[grep("NFKB1",
                  names(table(nbr$Genetic_variant)))])
NFKB1_cohort_prop <- round((NFKB1_cohort_size/diagnosed_cohort_size)*100, 2)

NFKB1_sequenced_cohort_prop <- 
  round((NFKB1_cohort_size/sequenced_cohort_size)*100, 2)


####################canonical form of TNFRSF13B ################################
# subset, size, proportion of records with the canonical form of TNFRSF13B
canonical_TNFRSF13B_cohort_size <-
  sum(table(nbr$Genetic_variant)[grep("canonicalTNFRSF13B",
                     names(table(nbr$Genetic_variant)))])
canonical_TNFRSF13B_cohort_prop <- 
  round((canonical_TNFRSF13B_cohort_size/diagnosed_cohort_size)*100, 2)

canonical_TNFRSF13B_sequenced_cohort_prop <- 
  round((canonical_TNFRSF13B_cohort_size/sequenced_cohort_size)*100, 2)


#################### rare form of TNFRSF13B #################################### 
# subset, size, proportion of records with the rare form of TNFRSF13B
rare_TNFRSF13B_cohort_size <- sum(table(nbr$Genetic_variant)[grep(
  "rareTNFRSF13B", names(table(nbr$Genetic_variant)))])
rare_TNFRSF13B_cohort_prop <- 
  round((rare_TNFRSF13B_cohort_size/diagnosed_cohort_size)*100, 2)

rare_TNFRSF13B_sequenced_cohort_prop <- 
  round((rare_TNFRSF13B_cohort_size/sequenced_cohort_size)*100, 2)


######################### Summarise ############################################
# putting all descriptive data together
biological_measure <- c("CD21", "CD21",
               "SmB", "SmB", "Tr","Tr",
               rep("Ig", 4),
               rep("Variants", 4))
biological_category <- c(rep("B cell", 6),
     rep("immunoglobulin", 4),
     rep("Genotype", 4))

biological_level <- c("expansion of low CD21", "CD21 normal",
              "SmB minus", "SmB plus", "Transitional B normal",
              "Transitional B high",
              ig_vec,
              "NFKB1",
              "canonical_TNFRSF13B",
              "rare_TNFRSF13B", 
              "any_other_diagnosis_variant")

number_of_patients <- c(cd21low_cohort_size,
                        cd21normal_cohort_size,
           smBminus_cohort_size,
           smBplus_cohort_size, 
           TrNorm_cohort_size,
           TrHigh_cohort_size,
           ig_size,
           NFKB1_cohort_size,
           canonical_TNFRSF13B_cohort_size,
           rare_TNFRSF13B_cohort_size,
           any_diagnosis_cohort_size)

prop_of_patients <- c(cd21low_cohort_prop,
                       cd21normal_cohort_prop,
                       smBminus_cohort_prop,
                       smBplus_cohort_prop, 
                       TrNorm_cohort_prop,
                       TrHigh_cohort_prop,
                       lll_cohort_prop,
                       lln_cohort_prop,
                       llh_cohort_prop,
                       lnl_cohort_prop,
                       NFKB1_cohort_prop,
                       canonical_TNFRSF13B_cohort_prop,
                       rare_TNFRSF13B_cohort_prop,
                       any_diagnosis_cohort_prop)

# plots will be bar plot, % will be in captions.
data <- data.frame(biological_measure,
                   biological_category,
                   biological_level,
                   number_of_patients,
                   prop_of_patients)

####################### save all
# for main
# Switched memory B cell measurements
c(smb_cohort_size, smb_cohort_prop)

# CD21low B cells data
c(cd21_cohort_size, cd21_cohort_prop)

# Transitional B cells
c(tr_cohort_size, tr_cohort_prop)

# details
data

#Genetic data 
c(sequenced_cohort_size, sequenced_cohort_prop)

#relevant genetic finding 
c(diagnosed_cohort_size, diagnosed_cohort_prop)

# NFKB1 variants 
c(NFKB1_cohort_size, NFKB1_cohort_prop, NFKB1_sequenced_cohort_prop)

#Variants in TNFRSF13B canonical form
c(canonical_TNFRSF13B_cohort_size,
  canonical_TNFRSF13B_sequenced_cohort_prop)

#the rare form
c(rare_TNFRSF13B_cohort_size,
  rare_TNFRSF13B_sequenced_cohort_prop)


# variants in genes other than NFKB1 or TNFRSF13B.  
c(any_diagnosis_cohort_size,
  any_diagnosis_cohort_prop )

fwrite(cd21low_cohort_list, "../result/CD21low_high_patients.csv")
fwrite(cd21normal_cohort_list, "../result/CD21low_normal_patients.csv")
fwrite(smBplus_cohort_list,  "../result/SmB_plus_patients.csv")
fwrite(smBminus_cohort_list,  "../result/SmB_minus_patients.csv")
fwrite(TrNorm_cohort_list,  "../result/TransitionalB_Normal_patients.csv")
fwrite(TrHigh_cohort_list,  "../result/TransitionalB_High_patients.csv")
fwrite(data, "../result/Fig2/num_patients_per_biological_category.csv")
