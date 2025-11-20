# source common.R, which contains common function, data, libraries
source("common.R")


### EPIDEMIOLOGY TABLE

# load nbr
# col: STUDY_ID,FAMILY_ID,CREATED,CENTRE,SEX,AFFECTED,AETIOLOGY_CONFIRMED,
# AETIOLOGY_REPORT,VARIANTS,GENES,PHENOTYPES,	
# abnormal_t_cell_proliferation_to_pha,age_of_onset_of_symptoms,	
# b_cells_uL,baseline_or_pre_ig_replacement_igg_gL,cd21_low_b_cells_pct,
# cd3_cells_uL,cd4_cells_uL,cd8_cells_uL,class_switched_b_cells_pct,
# consanguinity,diagnosis,Diagnosis gene,family_history,Gene role confirmed,	
# glild	iga_gL,ige_kUL,igm_gL,marginal_zone_b_cells_pct,memory_b_cells_pct,
# naive_b_cells_pct,nk_cells_uL,normal_immunology_test_results,
# on_ig_replacement,on_ig_replacement_igg_gL,	
# on_immunosuppression_at_the_time_of_laboratory_testing,
# other_immunology_test_results,other_information,paediatric_case,
# plasmablasts_pct,transitional_b_cells_pct,year_of_birth,
# BRIDGE,LAB,PACK,PLATEKEY,RNA_Z,SP,WGS,
# Each row describes a participant with PCT entryhpo,
# sequenced,Immuno_group,Genetic_variant
nbr <- readRDS("../result/tidy_data")

# number of participants
number_records <- nrow(nbr) 

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

# number of records without sex or year of birth
no_sex <- nbr$STUDY_ID[is.na(nbr$SEX)]
no_sex_centre <- nbr[is.na(nbr$SEX), c("STUDY_ID","CENTRE")]
num_no_sex_recorded <- length(no_sex)
no_yob <- nbr$STUDY_ID[is.na(nbr$year_of_birth)]
no_yob_centre <- nbr[is.na(nbr$year_of_birth), c("STUDY_ID","CENTRE")]
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
write.csv(epidemiologic_data_table, "../result/Table1/epidemiologic_data_table.csv")
