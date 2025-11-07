# source common.R, which contains common function, data, libraries
source("common.R")


### EPIDEMIOLOGY TABLE

# load nbr
nbr <- readRDS("../result/tidy_data")

# number of participants 528
number_records <- nrow(nbr) 
#number_records
# save list of CVID patients
#write.csv(nbr$STUDY_ID, "../result/CVID_patients.csv")

# total number of patients 
# 528 across 11 centers are CVID and have complete HPO phenotype annotation

# name the focus phenotypes for the table
hpo_names <- list(
  # Splenomegaly
  "Splenomegaly",
  # Lymphadenopathy
  "Lymphadenopathy",
  # Granulomatous disease
  c("Granulomatosis", "Granuloma","Pulmonary granulomatosis"),
  # Autoimmunity
  c("Autoimmunity", "Vitiligo", "Hashimoto thyroiditis",
    "Atrophic gastritis", "Malabsorption of Vitamin B12",
    "Type I diabetes mellitus", "Autoimmune hypoparathyroidism",
    "Billiary cirrhosis"),
  # .... of which are Autoimmune cytopaenias
  c("Autoimmune hemolytic anemia",
    "Autoimmune thrombocytopenia",
    "Neutropenia in presence of anti-neutropil antibodies"),
  # NRH liver
  "Nodular regenerative hyperplasia of liver",
  # Gastrointestinal inflammation
  "Gastrointestinal inflammation",
  # Hematological neoplasm
  "Hematological neoplasm"
)

# obtain HPO code from the name (list of 8 elements)
hpo_codes <- lapply(hpo_names, 
                    function(x) unname(hpo$id[which(hpo$name %in% x)]))

# a vector with one element (number of terms) 
n_terms <- length(hpo_codes)

# empty matrix for results (the number of HPO per category)
counts <- rep(NA, n_terms)

# looping through each phenotype
for(i in 1:n_terms) {
  
  # Get set of terms containing all descendants of terms in a given set
  term_descendants <- get_descendants(hpo, hpo_codes[[i]])
  
  # count patients with at least one descendant term for the phenotype
  counts[i] <- sum(sapply(lapply(nbr$hpo,
                                 function(x) term_descendants %in% x),
                          any))
  
}

# number of records without sex
no_sex <- nbr$STUDY_ID[is.na(nbr$SEX)]
#fwrite(as.data.frame(no_sex), "../result/STUDY_ID_no_sex.csv")
no_sex_centre <- nbr[is.na(nbr$SEX), c("STUDY_ID","CENTRE")]
#fwrite(no_sex_centre, "../result/STUDY_ID_Centre_no_sex.csv")
num_no_sex_recorded <- length(no_sex)

# number of records without age
no_yob <- nbr$STUDY_ID[is.na(nbr$year_of_birth)]
#fwrite(as.data.frame(no_yob), "../result/STUDY_ID_no_yob.csv")
no_yob_centre <- nbr[is.na(nbr$year_of_birth), c("STUDY_ID","CENTRE")]
#fwrite(no_yob_centre, "../result/STUDY_ID_Centre_no_yob.csv")
num_no_yob_record <- sum(is.na(nbr$year_of_birth))

# make the table
epidemiologic_data_table <- data.frame(
  Number_of_Records = c(nbr %>% dplyr::filter(!is.na(SEX)) %>% nrow(),
                        nbr %>% dplyr::filter(!is.na(year_of_birth)) %>% nrow(),
                        counts),
  Clinical_data = c(paste(table(nbr$SEX)[1], "females, ",
                        table(nbr$SEX)[2], "males, and ",
                        num_no_sex_recorded,
                        " with no sex recorded"),
                    paste(ceiling(median(nbr$year_of_birth, na.rm = TRUE)),
                          " (",
                          ceiling(sd(nbr$year_of_birth, na.rm = TRUE)),
                          "), excluding ",
                          num_no_yob_record,
                          " patients with no age recorded",
                          sep = ""),
                paste(round((counts/number_records) * 100, 1), "%", sep = "")),
                  row.names = c("Sex",
                                "Year of birth (SD)",
                                "Splenomegaly, %",
                                "Lymphadenopathy, %",
                                "Granulomatous disease, %",
                                "Autoimmunity, %",
                                ".... of which are Autoimmune cytopaenias, %",
                                "NRH liver, %",
                                "Gastrointestinal inflammation, %",
                                "Hematological neoplasm, %")
                  
)

# save table
write.csv(epidemiologic_data_table, "../result/epidemiologic_data_table.csv")



