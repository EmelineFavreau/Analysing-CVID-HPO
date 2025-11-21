######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS")


############# Suppl Table Complete HPO list ####################################

df <- patient_hpo_bio_mat[grep(rownames(patient_hpo_bio_mat), pattern = "^HP"),]

HPO_code = row.names(df)

df1 <- data.table(HPO_code = HPO_code,
                 HPO_freq = rowSums(df),
                 HPO_name = hpo$name[match(HPO_code, names(hpo$name))])

cohort_size <- ncol(patient_hpo_bio_mat)

# for main
 
output_proportion <- function(vec){
  pp <- (df1$HPO_freq[df1$HPO_name %in%
                  vec]/cohort_size)*100
  return(pp)
}
#decreased circulating IgG, IgA and IgM levels,
output_proportion(c("Decreased circulating IgG level",
                "Decreased circulating IgA level",
                "Decreased circulating total IgM"))

#decreased proportion of class-switched memory B cells 
output_proportion(c("Decreased proportion of class-switched memory B cells"))

#recurrent bacterial infections
output_proportion(c("Recurrent bacterial infections"))

#decreased DLCO
output_proportion(c("Decreased DLCO"))

#ground-glass opacification
output_proportion(c("Ground-glass opacification"))

#nodular regenerative hyperplasia of liver
output_proportion(c("Nodular regenerative hyperplasia of liver"))

#portal hypertension
output_proportion(c("Portal hypertension"))

#villous atrophy 
output_proportion(c("Villous atrophy"))

fwrite(df1, "../result/SI/Suppl_Table_Complete_HPO_list.csv")
