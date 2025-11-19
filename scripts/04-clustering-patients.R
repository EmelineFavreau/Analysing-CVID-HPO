######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# load nbr
nbr <- readRDS("../result/tidy_data")

######################## analysis ##########################################
# transform the spc_tbl_ into a data.table 
x <- as.data.table(nbr)

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

# shorten name
z <- patient_hpo_mat

# Infection
# Unusual infection [HP:0032101>=].
# Respiratory tract infection [HP:0011947>=]
# Abscess [HP:0025615>=].
keys=list(infection = as.vector(unname(hpo$name[names(hpo$name) %in% 
                        c(get_descendants(ontology = hpo, "HP:0032101"),
                          get_descendants(ontology = hpo, "HP:0011947"),
                          get_descendants(ontology = hpo, "HP:0025615"))])),
          
          # Bronchiectasis [HP:0002110>=].
          bronchiectasis = as.vector(unname(hpo$name[names(hpo$name) %in% 
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
          # Hepatosplenomegaly [HP:0001433] another sign of lymphoproliferation.
          # "Fluctuating splenomegaly"
          splenomegaly=as.vector(unname(hpo$name[names(hpo$name) %in%
                             get_descendants(ontology = hpo, "HP:0001744")])),
          
          # Liver disease
          # Nodular Regenerative hyperplasia of liver [HP:0011954]
          NRH_of_liver=as.vector(unname(hpo$name[names(hpo$name) ==
                                                   "HP:0011954"])),
          
          
          
          # Lymphoproliferation
          # Intestinal Lymphoid nodular hyperplasia [HP:0011956] 
          # Lymphadenopathy [HP:0002716] 
          # Hematological neoplasm [HP:0004377>=]
          lymphoproliferation=as.vector(unname(hpo$name[names(hpo$name) %in% 
                            c("HP:0011956", 
                              get_descendants(ontology = hpo, "HP:0002716"),
                              get_descendants(ontology = hpo, "HP:0004377"))])), 
          
          # Autoimmunity
          # Autoimmunity [HP:0002960>=] 
          # Type I diabetes [HP:0100651]
          # Rheumatoid arthritis [HP:0001370]
          # vitiligo [HP:0001045]
          # Autoimmune hypoparathyroidism [HP:0011771
          #Billiary cirrhosis [HP:0002613]
          autoimmunity=as.vector(unname(hpo$name[names(hpo$name) %in% 
          c("HP:0100651",
            "HP:0001370",
            "HP:0001045",
            "HP:0011771",
            "HP:0002613",
             get_descendants(ontology = hpo, "HP:0002960"))])))

# list of 7 key phenotypes, 
# each with a vector of clinical names for associated HPOs
keys=lapply(keys, intersect, rownames(z))

# matrix of presence/absence of HPOs (row) per patient (column)
newz=lapply(seq_along(keys), function(i) {
  apply(z[keys[[i]],,drop=FALSE], 2, max)}) %>% do.call("rbind",.)

# name the rows with the key phenotypes
rownames(newz)=names(keys)

# restrict to patients with 1+ key term 
# to construct distances between them
newz = newz[, colSums(newz)>0]

# list of patients for each HPO group:
key_vec <- c("infection",
             "bronchiectasis",
             "granuloma",
             "splenomegaly",
             "NRH_of_liver",
             "lymphoproliferation",
             "autoimmunity")

# save list of patients for each key
for(k in key_vec){
  kk <- paste(k, "_patients", sep="")
  kkk <- colnames(newz)[newz[row.names(newz) == k, ] == 1]
  assign(kk, kkk)    
  fwrite(as.data.frame(kkk),
         paste("../result/key_hpo_",k,"_patients.csv", sep =""))
}

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
# infection and or bronchiectasis patients (no GLILD)
infectionBronchiectasis <- infectionBronchiectasis[!infectionBronchiectasis 
                                                   %in% infection_with_glild]


# complex patients are all remaining patients
complexPatients <- colnames(newz)[!colnames(newz) %in% infectionBronchiectasis]

####################### save all
fwrite(patient_hpo_mat, "../result/patient_hpo_mat.csv", row.names = T)

fwrite(as.data.frame(infectionBronchiectasis),
       "../result/InfectionBronchiectasisPatients.csv")
fwrite(as.data.frame(complexPatients),
       "../result/complexPatients.csv")
