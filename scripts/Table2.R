# source common.R, which contains common function, data, libraries
source("common.R")


############# prep ############################################################
# load nbr
nbr <- readRDS("../result/tidy_data")

#a table with sex + gene + hpo for the patients with a genetic diagnosis?
tt <- nbr %>% 
  dplyr::filter(!is.na(Genetic_variant)) %>% 
  dplyr::select(SEX, Genetic_variant, hpo)

tt$Genetic_variant[tt$Genetic_variant == "anyPathogenic"] <- "any pathogenic variant"
tt$Genetic_variant[tt$Genetic_variant == "anyPathogenicNFKB1"] <- "NFKB1"
tt$Genetic_variant[tt$Genetic_variant == "anyPathogeniccanonicalTNFRSF13B"] <- "canonical TNFRSF13B"
tt$Genetic_variant[tt$Genetic_variant == "anyPathogenicrareTNFRSF13B"] <- "rare TNFRSF13B"
tt$Genetic_variant[tt$Genetic_variant == "anyPathogenicNFKB1canonicalTNFRSF13B"] <- "NFKB1 & canonical TNFRSF13B"
tt$Genetic_variant[tt$Genetic_variant == "canonicalTNFRSF13B" ] <- "canonical TNFRSF13B"

fwrite(tt, "../result/Table2/Table2.csv")
