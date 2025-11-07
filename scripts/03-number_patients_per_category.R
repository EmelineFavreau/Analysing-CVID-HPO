######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# Number of patients in each subcategory of B cell:

# There are 528 patients in our cohort.
cohort_size <- nrow(nbr)
cohort_size


#################### count patients in each group ##############################
# Out of those, 294 (55.68%) have a record for CD21 low. 
cd21_cohort_size <- nbr %>%
  dplyr::filter(!is.na(cd21_low_b_cells_pct)) %>% nrow()

cd21_cohort_prop <- round((cd21_cohort_size/cohort_size)*100, 2)
cd21_cohort_size # 294
cd21_cohort_prop # 55.68

# Out of all records with CD21 low: 
# 149 (50.68%) patients are in the category "expansion of low CD21"
cd21low_cohort_list <- nbr %>%
  dplyr::filter(cd21_low_b_cells_pct >= 10) %>% 
  dplyr::select(STUDY_ID)

fwrite(cd21low_cohort_list, "../result/CD21low_high_patients.csv")

cd21low_cohort_size <- nrow(cd21low_cohort_list)

cd21low_cohort_prop <- round((cd21low_cohort_size/cd21_cohort_size)*100, 2)
cd21low_cohort_size # 149
cd21low_cohort_prop # 50.68

# Out of all records with CD21 low: 
# 145 (49.32%) patients are in the category "normal CD21".
cd21normal_cohort_list <- nbr %>%
  dplyr::filter(cd21_low_b_cells_pct < 10) %>% dplyr::select(STUDY_ID)

fwrite(cd21normal_cohort_list, "../result/CD21low_normal_patients.csv")


cd21normal_cohort_size <- nrow(cd21normal_cohort_list)

cd21normal_cohort_prop <- 
  round((cd21normal_cohort_size/cd21_cohort_size)*100, 2)

cd21normal_cohort_size # 145
cd21normal_cohort_prop # 49.32
#########################


#########################
# Out of those, 364 (68.94%) have a record for SmB. 
smb_cohort_size <- nbr %>%
  dplyr::filter(!is.na(class_switched_b_cells_pct)) %>% nrow()

smb_cohort_prop <- round((smb_cohort_size/cohort_size)*100, 2)
smb_cohort_size # 364
smb_cohort_prop # 68.94


# Out of all records with SmB: 
# 184 (50.55%) patients are in the category "SmB+"
smBplus_cohort_list <- nbr %>%
  dplyr::filter(class_switched_b_cells_pct > 2) %>% dplyr::select(STUDY_ID)

fwrite(smBplus_cohort_list,  "../result/SmB_plus_patients.csv")

smBplus_cohort_size <- nrow(smBplus_cohort_list)

smBplus_cohort_prop <- round((smBplus_cohort_size/smb_cohort_size)*100, 2)
smBplus_cohort_size # 184
smBplus_cohort_prop # 50.55

# Out of all records with SmB: 
# 180 (49.45%) patients are in the category "SmB-"
smBminus_cohort_list <- nbr %>%
  dplyr::filter(class_switched_b_cells_pct <= 2) %>% dplyr::select(STUDY_ID)

smBminus_cohort_size <- nrow(smBminus_cohort_list)

fwrite(smBminus_cohort_list,  "../result/SmB_minus_patients.csv")

smBminus_cohort_prop <- round((smBminus_cohort_size/smb_cohort_size)*100, 2)

smBminus_cohort_size  # 180
smBminus_cohort_prop # 49.45
#########################


#########################
# Out of those, 281 (53.02%) have a record for Tr. 
# TrNorm (<9%) and 20 TrHigh (>=9%).
tr_cohort_size <- nbr %>%
  dplyr::filter(!is.na(transitional_b_cells_pct)) %>% nrow()

tr_cohort_prop <- round((tr_cohort_size/cohort_size)*100, 2)
tr_cohort_size # 290
tr_cohort_prop # 54.92

# Out of all records with Tr: 
# 264 (91.03%) patients are in the category "TrNorm"
TrNorm_cohort_list <- nbr %>%
  dplyr::filter(transitional_b_cells_pct < 9) %>% dplyr::select(STUDY_ID)

fwrite(TrNorm_cohort_list,  "../result/TransitionalB_Normal_patients.csv")

TrNorm_cohort_size <- nrow(TrNorm_cohort_list)

TrNorm_cohort_prop <- round((TrNorm_cohort_size/tr_cohort_size)*100, 2)
TrNorm_cohort_size # 264
TrNorm_cohort_prop # 91.03


# Out of all records with Tr: 
# 26 (8.97%) patients are in the category "TrHigh"
TrHigh_cohort_list <- nbr %>%
  dplyr::filter(transitional_b_cells_pct >= 9) %>% dplyr::select(STUDY_ID)

TrHigh_cohort_size <- nrow(TrHigh_cohort_list)

fwrite(TrHigh_cohort_list,  "../result/TransitionalB_High_patients.csv")

TrHigh_cohort_prop <- round((TrHigh_cohort_size/tr_cohort_size)*100, 2)
TrHigh_cohort_size # 26
TrHigh_cohort_prop # 8.97
#########################








#########################
# Number of patients in each subcategory of Ig
# sort(table(nbr$Immuno_group))
ig_vec <- c("lowIgGlowIgAlowIgM","lowIgGlowIgAnormalIgM",
            "lowIgGlowIgAhighIgM","lowIgGnormalIgAlowIgM")

ig_cohort_size <- sum(nbr$Immuno_group %in% ig_vec )
ig_cohort_prop <- round((ig_cohort_size/cohort_size)*100, 2)
ig_cohort_size # 496
ig_cohort_prop # 93.94

#There are all patients with complete sets of immunoglobulin levels.
#Out of those 496 patients, 
# there are 408 (82.26%) patients with low levels of IgG, IgA, IgM.
lll_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[1]]
lll_cohort_prop <- round((lll_cohort_size/ig_cohort_size)*100, 2)

lll_cohort_size # 408
lll_cohort_prop # 82.26

#Out of all patients with complete sets of immunoglobulin levels.
# there are 67 (13.51%) patients with low levels of IgG and IgA,
#normal levels of IgM.
lln_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[2]]
lln_cohort_prop <- round((lln_cohort_size/ig_cohort_size)*100, 2)

lln_cohort_size # 67
lln_cohort_prop  # 13.51


#Out of all patients with complete sets of immunoglobulin levels.
# there are 8 (1.61 %) patients with low levels of IgG and IgA,
#high levels of IgM.
llh_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[3]]
llh_cohort_prop <- round((llh_cohort_size/ig_cohort_size)*100, 2)

llh_cohort_size # 8
llh_cohort_prop # 1.61

#Out of all patients with complete sets of immunoglobulin levels.
# there are 13 (2.62 %) patients with low levels of IgG and IgM,
#normal levels of IgA.
lnl_cohort_size <- table(nbr$Immuno_group)[names(table(nbr$Immuno_group)) ==
                                             ig_vec[4]]
lnl_cohort_prop <- round((lnl_cohort_size/ig_cohort_size)*100, 2)
lnl_cohort_size # 13
lnl_cohort_prop # 2.62 
ig_size <- c(lll_cohort_size, lln_cohort_size, llh_cohort_size, lnl_cohort_size)

#########################


#########################
# unsure about this, I will assign anypathogenic but need to check with Matt
t <- nbr %>% 
  dplyr::filter(!is.na(Genetic_variant) | !is.na(`Diagnosis gene`)) %>% 
  dplyr::select(c(STUDY_ID, Genetic_variant, `Diagnosis gene`))
nbr$Genetic_variant[nbr$STUDY_ID %in% t$STUDY_ID[is.na(t$Genetic_variant)]] <-
  "anyPathogenic"

# Number of patients in each subcategory of genetic variants 
# (IUIS genes with frequencies of variants)
# in the cohort, some have been sequenced, some not. 
sequenced_cohort_size <- sum(nbr$sequenced == 1)
sequenced_cohort_size  # 434 sequenced out of 528
sequenced_cohort_size /nrow(nbr) *100 # 82.2%

# out of the sequenced, 95 has received a confirmed genetic diagnosis
diagnosed_cohort_size <- sum(table(nbr$Genetic_variant))
diagnosed_cohort_size
 
# There are 95 patients with confirmed genetic variants 
diagnosed_cohort_size /sequenced_cohort_size *100 # 21.89%

# Out of those 95 patients, 
# there are 41 patients (43.16%) with 
# any diagnostic variant (no NFKB1, no TNFRSF13B)
any_diagnosis_cohort_size <- unname(table(nbr$Genetic_variant)[
  names(table(nbr$Genetic_variant)) == "anyPathogenic"])
any_diagnosis_diagnosed_cohort_prop <- 
  round((any_diagnosis_cohort_size/diagnosed_cohort_size)*100, 2)
any_diagnosis_cohort_size  # 41
any_diagnosis_diagnosed_cohort_prop # 43.16
any_diagnosis_sequenced_cohort_prop <- 
  round((any_diagnosis_cohort_size/sequenced_cohort_size)*100, 2)
any_diagnosis_sequenced_cohort_prop # 9.45

# Out of those 95 patients, 
# there are 19 (20%) patients with NFKB1
NFKB1_cohort_size <- sum(table(nbr$Genetic_variant)[grep("NFKB1",
                  names(table(nbr$Genetic_variant)))])
NFKB1_cohort_prop <- round((NFKB1_cohort_size/diagnosed_cohort_size)*100, 2)
NFKB1_cohort_size # 19
NFKB1_cohort_prop # 20
NFKB1_sequenced_cohort_prop <- 
  round((NFKB1_cohort_size/sequenced_cohort_size)*100, 2)
NFKB1_sequenced_cohort_prop # 4.38

# Out of those 95 patients, 
# there are 27 (28.42%) patients with the canonical form of TNFRSF13B
canonical_TNFRSF13B_cohort_size <-
  sum(table(nbr$Genetic_variant)[grep("canonicalTNFRSF13B",
                     names(table(nbr$Genetic_variant)))])
canonical_TNFRSF13B_cohort_prop <- 
  round((canonical_TNFRSF13B_cohort_size/diagnosed_cohort_size)*100, 2)
canonical_TNFRSF13B_cohort_size # 27
canonical_TNFRSF13B_cohort_prop # 28.42
canonical_TNFRSF13B_sequenced_cohort_prop <- 
  round((canonical_TNFRSF13B_cohort_size/sequenced_cohort_size)*100, 2)
canonical_TNFRSF13B_sequenced_cohort_prop # 6.22

# Out of those 95 patients, 
# there are 9 (9.47%) patients with the rare form of TNFRSF13B
rare_TNFRSF13B_cohort_size <- sum(table(nbr$Genetic_variant)[grep("rareTNFRSF13B",
                      names(table(nbr$Genetic_variant)))])
rare_TNFRSF13B_cohort_prop <- 
  round((rare_TNFRSF13B_cohort_size/diagnosed_cohort_size)*100, 2)
rare_TNFRSF13B_cohort_size # 9
rare_TNFRSF13B_cohort_prop # 9.47
rare_TNFRSF13B_sequenced_cohort_prop <- 
  round((rare_TNFRSF13B_cohort_size/sequenced_cohort_size)*100, 2)
rare_TNFRSF13B_sequenced_cohort_prop # 2.07

# Out of those 95 patients, 
# there are 1 (1.06%) patients with both the canonical and
# the rare form of TNFRSF13B
can_rare_TNFRSF13B_cohort_size <- 
  sum(table(nbr$Genetic_variant)[grep("NFKB1canonicalTNFRSF13B",
                                    names(table(nbr$Genetic_variant)))])
can_rare_TNFRSF13B_cohort_prop <- 
  round((can_rare_TNFRSF13B_cohort_size/diagnosed_cohort_size)*100, 2)
can_rare_TNFRSF13B_cohort_size # 1
can_rare_TNFRSF13B_cohort_prop # 1.05
can_rare_TNFRSF13B_sequenced_cohort_prop <- 
  round((can_rare_TNFRSF13B_cohort_size/sequenced_cohort_size)*100, 2)
can_rare_TNFRSF13B_sequenced_cohort_prop # 0.23
#########################


#########################
# putting all descriptive data together
biological_measure <- c("CD21", "CD21",
               "SmB", "SmB", "Tr","Tr",
               "Ig", "Ig", "Ig", "Ig",
               "Variants","Variants","Variants", "Variants","Variants")
biological_category <- c(rep("B cell", 6),
     rep("immunoglobulin", 4),
     rep("Genotype", 5))

biological_level <- c("expansion of low CD21", "CD21 normal",
              "SmB minus", "SmB plus", "Transitional B normal",
              "Transitional B high",
              ig_vec,
              "NFKB1",
              "canonical_TNFRSF13B",
              "rare_TNFRSF13B", 
              "NFKB1canonicalTNFRSF13B",
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
           can_rare_TNFRSF13B_cohort_size,
           any_diagnosis_cohort_size)

# plots will be bar plot, % will be in captions.
data <- data.frame(biological_measure,
                   biological_category,
                   biological_level,
                   number_of_patients)

fwrite(data, "../result/Fig2/num_patients_per_biological_category.csv")
