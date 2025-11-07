######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")


######################## analysis ##########################################
# transform the spc_tbl_ into a data.table 
x <- as.data.table(nbr)

# focus on column: the HPO code from x$hpo

# new column: the HPO term
# find the clinical description for each code
x$hpo_long <- sapply(x$hpo,
    function(y) hpo$name[match(unlist(y), names(hpo$name))])

## make this a vector
hpo_long <- unlist(x$hpo_long)
unique_long <- unique(hpo_long)
hpo_codes <- unlist(x$hpo)
unique_codes <- unique(hpo_codes)

## create a matching patient id vector
patient <- rep(x$STUDY_ID, 
               sapply(x$hpo, length))

## make empty matrix - default value is 0
patient_hpo_mat <- matrix(data = 0,
                          nrow = length(unique_codes),
                          ncol = nrow(x),
                          dimnames=list(unique_long, x$STUDY_ID))

## add 1 where the term matches a patient
patient_hpo_mat[cbind(hpo_long, patient)] <- 1

fwrite(patient_hpo_mat, "../result/patient_hpo_mat.csv", row.names = T)

# shorten name
z <- patient_hpo_mat

## * look at key terms for complex/noncomplex CVID
# find in unique_long all terms related to Luiza's HPO list
# Note: ‘>=’means the term in the brackets + its descendants

# Infection
# Unusual infection [HP:0032101>=]. Comment: this should include the terms you already have on your text file.
# Respiratory tract infection [HP:0011947>=]
# Abscess [HP:0025615>=]. Comment: it can be a sterile abscess but it’s mainly infectious.
keys=list(infection=as.vector(unname(hpo$name[names(hpo$name) %in% 
                                                c(get_descendants(ontology = hpo, "HP:0032101"),
                                                  get_descendants(ontology = hpo, "HP:0011947"),
                                                  get_descendants(ontology = hpo, "HP:0025615"))])),
          
          # Bronchiectasis [HP:0002110>=]. Comment: as a result of recurrent infections
          bronchiectasis= as.vector(unname(hpo$name[names(hpo$name) %in% 
                                                      get_descendants(ontology = hpo, "HP:0002110")])),
          
          # Granulomas
          # Granuloma [HP:0032252>=]
          # Granulomatosis [HP:0002955>=]
          # Pulmonary granulomatosis [HP:0030250>=]
          granuloma=as.vector(unname(hpo$name[names(hpo$name) %in% 
                                                c(get_descendants(ontology = hpo, "HP:0032252"),
                                                  get_descendants(ontology = hpo, "HP:0002955"),
                                                  get_descendants(ontology = hpo, "HP:0030250"))])),
          
          # splenomegaly and its two descendants
          # Hepatosplenomegaly [HP:0001433] as another sign of lymphoproliferation.
          # "Fluctuating splenomegaly"
          splenomegaly=as.vector(unname(hpo$name[names(hpo$name) %in%
                                 get_descendants(ontology = hpo, "HP:0001744")])),
          
          # Liver disease
          # Nodular Regenerative hyperplasia of liver [HP:0011954]
          # comment: not sure whether to include other terms here yet
          NRH_of_liver=as.vector(unname(hpo$name[names(hpo$name) == "HP:0011954"])),
          
          
          
          # Lymphoproliferation
          # It’s worth including Intestinal Lymphoid nodular hyperplasia [HP:0011956] 
          # in addition to Lymphadenopathy [HP:0002716] 
          # Hematological neoplasm [HP:0004377>=] (separate investigation: clinical curiosity, run again with a split here) )
          lymphoproliferation=as.vector(unname(hpo$name[names(hpo$name) %in% 
                            c("HP:0011956", # Intestinal Lymphoid nodular hyperplasia
                              "HP:0002716", # Lymphadenopathy 
                              #"HP:0001433",
                              get_descendants(ontology = hpo, "HP:0004377"))])), # Hematological neoplasm
          
          # Autoimmunity
          # Autoimmunity [HP:0002960>=] Comment: should include all the terms you have on your text file.
          # Type I diabetes [HP:0100651]
          # Rheumatoid arthritis [HP:0001370]
          # vitiligo [HP:0001045]
          # Autoimmune hypoparathyroidism [HP:0011771
          #Billiary cirrhosis [HP:0002613]
          autoimmunity=as.vector(unname(hpo$name[names(hpo$name) %in% 
          c("HP:0100651", #Type I diabetes
            "HP:0001370", #Rheumatoid arthritis
            "HP:0001045", # vitiligo
            "HP:0011771", # autoimmune hypoparathyroidism
            "HP:0002613", # Billiary cirrhosis
             get_descendants(ontology = hpo, "HP:0002960"))])))

# save keys for figures
#saveRDS(keys, "../result/keys")

# list of 7 key phenotypes, each with a vector of clinical names for associated HPOs
keys=lapply(keys, intersect, rownames(z))

# matrix of presence/absence of HPOs (row) per patient (column)
newz=lapply(seq_along(keys), function(i) {
  apply(z[keys[[i]],,drop=FALSE], 2, max)}) %>% do.call("rbind",.)

# name the rows with the key phenotypes
rownames(newz)=names(keys)
#rowSums(newz)

                                                                                                                                                                     
## check frequency 
#rowMeans(newz)  %>% sort()


# some patients have no focus HPO (n = 1)
#no_HPO_patients <- colSums(newz)[colSums(newz) == 0]

# t <- nbr %>% filter(STUDY_ID %in% names(no_HPO_patients)) %>% 
# select(STUDY_ID, CENTRE, hpo, 
#        class_switched_b_cells_pct, cd21_low_b_cells_pct, 
#        transitional_b_cells_pct, glild)
# 
# t$hpo_long <- sapply(t$hpo,
#                      function(y) hpo$name[match(unlist(y), names(hpo$name))])
# 
# fwrite(t, "20250910-patients_no_focus_HPO.csv")
# restrict to patients with 1+ key term  (n=528 patients)
# so we can construct distances between them
newz = newz[, colSums(newz)>0]




# list of patients for each HPO group:
infection_patients <- colnames(newz)[newz[row.names(newz) == "infection", ] == 1]
bronchiectasis_patients <- colnames(newz)[newz[row.names(newz) == "bronchiectasis", ] == 1]
granuloma_patients <- colnames(newz)[newz[row.names(newz) == "granuloma", ] == 1]
splenomegaly_patients <- colnames(newz)[newz[row.names(newz) == "splenomegaly", ] == 1]
NRH_of_liver_patients <- colnames(newz)[newz[row.names(newz) == "NRH_of_liver", ] == 1]
lymphoproliferation_patients <- colnames(newz)[newz[row.names(newz) == "lymphoproliferation", ] == 1]
autoimmunity_patients <- colnames(newz)[newz[row.names(newz) == "autoimmunity", ] == 1]
# length(infection_patients) #481
# length(bronchiectasis_patients) #218
# length(granuloma_patients)#53
# length(splenomegaly_patients)#146
# length(NRH_of_liver_patients)#37
# length(lymphoproliferation_patients)#82
# length(autoimmunity_patients)#152

# save the list for each HPO group
fwrite(as.data.frame(infection_patients), "../result/key_hpo_infection_patients.csv")
fwrite(as.data.frame(bronchiectasis_patients), "../result/key_hpo_bronchiectasis_patients.csv")
fwrite(as.data.frame(granuloma_patients), "../result/key_hpo_granuloma_patients.csv")
fwrite(as.data.frame(splenomegaly_patients), "../result/key_hpo_splenomegaly_patients.csv")
fwrite(as.data.frame(NRH_of_liver_patients), "../result/key_hpo_NRH_of_liver_patients.csv")
fwrite(as.data.frame(lymphoproliferation_patients), "../result/key_hpo_lymphoproliferation_patients.csv")
fwrite(as.data.frame(autoimmunity_patients), "../result/key_hpo_autoimmunity_patients.csv")


# infectionBronchiectasis patients have infection and or bronchiectasis
infectionBronchiectasis <- unique(c(infection_patients,
                                    bronchiectasis_patients))

# patients with other HPO
other_HPO <- unique(c(granuloma_patients,
                      splenomegaly_patients,
         NRH_of_liver_patients,
         lymphoproliferation_patients,
         autoimmunity_patients))

# patients with infection and glild
infection_with_glild <- nbr %>%
  dplyr::filter(STUDY_ID %in% infectionBronchiectasis & glild == "Yes") %>% 
  dplyr::select(STUDY_ID) %>% unlist() %>% unname()

# patients with infection and or bronchiectasis, but no other key HPO
infectionBronchiectasis <- infectionBronchiectasis[!infectionBronchiectasis 
                                                   %in% other_HPO]
# 230 infection and or bronchiectasis patients (no GLILD)
infectionBronchiectasis <- infectionBronchiectasis[!infectionBronchiectasis 
                                                   %in% infection_with_glild]

fwrite(as.data.frame(infectionBronchiectasis),
       "../result/InfectionBronchiectasisPatients.csv")



# 298 complex patients are all remaining patients
complexPatients <- colnames(newz)[!colnames(newz) %in% infectionBronchiectasis]

fwrite(as.data.frame(complexPatients),
       "../result/complexPatients.csv")
