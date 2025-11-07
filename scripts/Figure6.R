######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")


thisfocus <- "no_Bcell_filtering" # the data includes incomplete data 
#(some patients have no B cell) 

# each row describes the relationship between two HPO terms
# e.g. 1st row: relationship between term 1 (HP:0004315) and term 885 (HP:0000001)
result_with_id <- fread(paste("../result/reduc_dim_input_with_id_",
                              thisfocus, ".csv", sep = ""))


# import the embeddings (capturing the context of the HPO, 
# how often it is associated with other HPOs in the onthology pathway)
HPO_id_embeddings <- 
  fread(paste("../result/HPO_id_", thisfocus,"_embeddings.emb", sep =""), skip = 1)
# dim(HPO_id_embeddings)
#2853  129

# import presence absence of HPO
# no obsolete terms (aka terms with no edge)
# HPO terms that are present at least once
# all patients 
phmnoNA <- 
  readRDS(paste("../result/patient_hpo_bio_mat_for_reduc_dim_", thisfocus,".rds", sep = ""))
# 884 528


# load patient clusters 
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")

############# analyse ###########################################################
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"
complex_patients <- complex_cluster$STUDY_ID
infection_patients <- infection_cluster$STUDY_ID

# vector of term ID in the order of the patient hpo matrix
term_order_in_phm <- result_with_id$term_id[match(row.names(phmnoNA),
                                                  result_with_id$term)]

# reorder the rows of HPO_id_embeddings using the order of the phm
hieo <- HPO_id_embeddings[match(term_order_in_phm,
                                HPO_id_embeddings$V1), -c("V1")]
phm <- phmnoNA
# sanity check
#dim(hieo)
# 884 128
#dim(phm)
# 884 528


# set counts of:
N_HPO_TERMS <- nrow(phm) # 884
D_DIMENSION <- ncol(hieo) # 128
N_PATIENTS <- ncol(phm) # 528

# terms x dimensions
HPO_FEATURES_F <- hieo


denom = matrix(colSums(phm),
               nrow = nrow(phm),
               ncol = ncol(phm),
               byrow = TRUE)


# patients x terms
PATIENT_PHENOTYPE_P <- t(phm/denom)

HPO_FEATURES_F <- as.matrix(HPO_FEATURES_F)

# Patient Embeddings = P %*% F Matrix Multiplication
# The matrix multiplication P (528 x 884) %*% F (800x128) yields P_embeddings (884 x 129).
# This transforms the sparse patient profile into a dense vector reflecting the 
# structural similarity encoded by the HPO embeddings.
patient_embeddings_P_embeddings <- PATIENT_PHENOTYPE_P %*% HPO_FEATURES_F


# 528 x 128  embeddings

## pca
p=prcomp(patient_embeddings_P_embeddings)$x  %>% as.data.frame()
p$cluster <- ifelse(rownames(p) %in% complex_patients,
                    "complex",
             ifelse(rownames(p) %in% infection_patients,
                    "infection", NA))
p$CENTRE <- nbr$CENTRE[match(rownames(p), nbr$STUDY_ID)]

left=ggplot(p,
            aes(x=PC1, y=PC2, col=cluster)) +
  geom_point(size = 0.71) + 
  scale_colour_manual(values = infComp) 
right=ggplot(p, aes(x=PC1, y=PC2, col=CENTRE, shape=CENTRE)) +
    geom_point() +
  scale_fill_manual(values = centrescol) +
    scale_colour_manual(values = centrescol) +
    scale_shape_manual(values = c(1:11)) 


patchwork2 <-  left / right 
patchwork2 + plot_annotation(tag_levels = 'a')
ggsave("../result/Fig6/Fig6.jpeg",
       width = 18,
       height = 18,
       units = "cm")

# sanity check
t.test(p$PC1 ~ p$cluster)
t.test(p$PC2 ~ p$cluster)
summary(lm(PC1 ~ CENTRE + cluster, data=p)) # cluster still significant even after adjusting for centre
summary(lm(PC2 ~ CENTRE + cluster, data=p)) # cluster still significant even after adjusting for centre
