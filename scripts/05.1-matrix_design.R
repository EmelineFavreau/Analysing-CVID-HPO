######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"

# load patients grouped by key hpo groups 
khip <- fread("../result/key_hpo_infection_patients.csv")
khbp <- fread("../result/key_hpo_bronchiectasis_patients.csv")
khgp <- fread("../result/key_hpo_granuloma_patients.csv")
khsp <- fread("../result/key_hpo_splenomegaly_patients.csv")
khnp <- fread("../result/key_hpo_NRH_of_liver_patients.csv")
khlp <- fread("../result/key_hpo_lymphoproliferation_patients.csv")
khap <- fread("../result/key_hpo_autoimmunity_patients.csv")

############## Matrix design ###########################
# matrix design: columns are patients, rows are HPO/glild/pheno/geno/centre

# vector of hpo
hpo_codes <- unlist(nbr$hpo)
unique_codes <- unique(hpo_codes)
names(unique_codes) <- hpo$name[match(unique_codes, names(hpo$name))]
saveRDS(unique_codes, "../result/unique_codes.RDS")

# create a matching patient id vector
patient <- rep(nbr$STUDY_ID, sapply(nbr$hpo, length))

# make empty matrix - default value is 0
phm <- matrix(data = 0,
                          nrow = length(unique_codes),
                          ncol = nrow(nbr),
                          dimnames = list(unique_codes, nbr$STUDY_ID))

# add 1 where the term matches a patient
phm[cbind(hpo_codes, patient)] = 1

# Get the patient order from the matrix columns
patient_order <- colnames(phm)

# Create all phenotype/group rows
phenotype_df <- nbr %>%
  dplyr::transmute(
    STUDY_ID,
    smb_minus = dplyr::if_else(class_switched_b_cells_pct <= 2,
                               1L, 0L, missing = 0L),
    smb_normal = dplyr::if_else(class_switched_b_cells_pct > 2,
                                1L, 0L, missing = 0L),
    cd21_low = dplyr::if_else(cd21_low_b_cells_pct >= 10,
                                1L, 0L, missing = 0L),
    cd21_plus = dplyr::if_else(cd21_low_b_cells_pct < 10,
                                1L, 0L, missing = 0L),
    tr_norm = dplyr::if_else(transitional_b_cells_pct < 9,
                                1L, 0L, missing = 0L),
    tr_high = dplyr::if_else(transitional_b_cells_pct >= 9,
                                1L, 0L, missing = 0L),
    reported_glild = dplyr::if_else(glild == "Yes",
                                1L, 0L, missing = 0L),
    reported_no_glild = dplyr::if_else(glild == "No",
                                1L, 0L, missing = 0L),
    lll = dplyr::if_else(Immuno_group == "lowIgGlowIgAlowIgM",
                                1L, 0L, missing = 0L),
    lln = dplyr::if_else(Immuno_group == "lowIgGlowIgAnormalIgM",
                                1L, 0L, missing = 0L),
    lnl = dplyr::if_else(Immuno_group == "lowIgGnormalIgAlowIgM",
                                1L, 0L, missing = 0L),
    llh = dplyr::if_else(Immuno_group == "lowIgGlowIgAhighIgM",
                                1L, 0L, missing = 0L),
    canonicalTNFRSF13B = dplyr::if_else(Genetic_variant %in% 
                                c("anyPathogeniccanonicalTNFRSF13B",
                                  "anyPathogenicNFKB1canonicalTNFRSF13B",
                                  "canonicalTNFRSF13B"),
                                1L, 0L, missing = 0L),
    rareTNFRSF13B = dplyr::if_else(Genetic_variant ==
                            "anyPathogenicrareTNFRSF13B",
                            1L, 0L, missing = 0L),
    NFKB1 = dplyr::if_else(Genetic_variant %in% 
                             c("NFKB1",
                               "anyPathogenicNFKB1canonicalTNFRSF13B",
                               "anyPathogenicNFKB1"),
                           1L, 0L, missing = 0L),
    anyPathogenic = dplyr::if_else(!is.na(Genetic_variant),
                                   1L, 0L),
  infection_cluster = as.integer(STUDY_ID %in% infection_cluster$STUDY_ID),
  complex_cluster   = as.integer(STUDY_ID %in% complex_cluster$STUDY_ID),
  keyHPOinfection   = as.integer(STUDY_ID %in% khip$infection_patients),
  keyHPObronchiectasis = as.integer(STUDY_ID %in% khbp$bronchiectasis_patients),
  keyHPOgranuloma = as.integer(STUDY_ID %in% khgp$granuloma_patients),
  keyHPOsplenomegaly = as.integer(STUDY_ID %in% khsp$splenomegaly_patients),
  keyHPONRHofliver = as.integer(STUDY_ID %in% khnp$NRH_of_liver_patients),
  keyHPOlymphoproliferation = as.integer(STUDY_ID %in% khlp$lymphoproliferation_patients),
  keyHPOautoimmunity = as.integer(STUDY_ID %in% khap$autoimmunity_patients)) %>%
  dplyr::arrange(match(STUDY_ID, patient_order)) %>%
  dplyr::select(-STUDY_ID)

# transpose
phenotype_mat <- t(as.matrix(phenotype_df))

# centre mapping
centres_vec <- c(CUH      = "Cambridge University Hospitals",
                 Glasgow  = "Queen Elizabeth University Hospital Glasgow",
                 Midlands = "University Hospital of North Midlands",
                 BH       = "Barts Health",
                 RF       = "Royal Free",
                 FP       = "Frimley Park",
                 Papworth = "Papworth",
                 SR       = "Salford Royal",
                 HE       = "Heart of England",
                 LTH      = "Leeds Teaching Hospitals",
                 ICH      = "Imperial College Healthcare")

# get STUDY_ID index
all_patient_ids <- colnames(phm)

# create a list of presence/absence for centre
centre_list <- lapply(names(centres_vec), function(sc) {
  full_centre_name <- centres_vec[sc]
  patients <- nbr %>%
    dplyr::filter(CENTRE == full_centre_name) %>%
    dplyr::pull(STUDY_ID)
  as.integer(all_patient_ids %in% patients)
  }
)

# make a matrix
centre_mat <- do.call(rbind, centre_list)

# name rows
rownames(centre_mat) <- names(centres_vec)

# combine
patient_hpo_bio_mat <- rbind(phm,
                              phenotype_mat,
                              centre_mat)

# save matrix
saveRDS(patient_hpo_bio_mat, "../result/patient_hpo_bio_mat.RDS")
