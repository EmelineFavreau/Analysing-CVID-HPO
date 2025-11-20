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
#Expected immunological features such as decreased circulating IgG, IgA and IgM levels, 

(df1$HPO_freq[df1$HPO_name %in% c("Decreased circulating IgG level",
                                 "Decreased circulating IgA level",
                                 "Decreased circulating total IgM")]/cohort_size)*100

#and a decreased proportion of class-switched memory B cells 
(df1$HPO_freq[df1$HPO_name %in% 
                c("Decreased proportion of class-switched memory B cells")]/cohort_size)*100

#were among the most frequently annotated terms, 
#along with recurrent bacterial infections (42.2%). 
(df1$HPO_freq[df1$HPO_name %in% 
    c("Recurrent bacterial infections")]/cohort_size)*100
#In addition to these diagnostics and defining features, 
#granular phenotypes were also captured, including decreased DLCO (9.65%),
(df1$HPO_freq[df1$HPO_name %in% 
                c("Decreased DLCO")]/cohort_size)*100
#ground-glass opacification (9.1%), "ground-glass opacification "
(df1$HPO_freq[df1$HPO_name %in% 
                c("Ground-glass opacification")]/cohort_size)*100
#nodular regenerative hyperplasia of liver (6.8%) 
(df1$HPO_freq[df1$HPO_name %in% 
                c("Nodular regenerative hyperplasia of liver")]/cohort_size)*100
#portal hypertension (6.4%), 
(df1$HPO_freq[df1$HPO_name %in% 
                c("Portal hypertension")]/cohort_size)*100
#and villous atrophy (3%). 
(df1$HPO_freq[df1$HPO_name %in% 
                c("Villous atrophy")]/cohort_size)*100
fwrite(df1, "../result/SI/Suppl_Table_Complete_HPO_list.csv")
