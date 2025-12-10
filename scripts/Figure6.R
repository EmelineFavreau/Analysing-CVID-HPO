######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

tf <- "no_Bcell_filtering"

# each row describes the relationship between two HPO terms
# e.g. 1st row: relationship between term 1 (HP:0004315) & term 885 (HP:0000001)
result_with_id <- fread(paste("../result/reduc_dim_input_with_id_",
                              tf, ".csv", sep = ""))


# import the embeddings (capturing the context of the HPO, 
# how often it is associated with other HPOs in the onthology pathway)
HPO_id_embeddings <- 
  fread(paste("../result/HPO_id_",
              tf,"_embeddings.emb", sep =""), skip = 1)

# import presence absence of HPO
# no obsolete terms (aka terms with no edge)
# HPO terms that are present at least once
# all patients 
phmnoNA <- 
  readRDS(paste("../result/patient_hpo_bio_mat_for_reduc_dim_",
                tf,".rds", sep = ""))

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

# set counts of:
nt <- nrow(phm) 
dd <- ncol(hieo)
np <- ncol(phm)

# number of terms per patient
denom <- matrix(colSums(phm),
               nrow = nrow(phm),
               ncol = ncol(phm),
               byrow = TRUE)


# patients x terms
ppp <- t(phm/denom)

hieo <- as.matrix(hieo)

# transforms the sparse patient profile into a dense vector reflecting the 
# structural similarity encoded by the HPO embeddings.
pemb <- ppp %*% hieo

## PCA
p <- prcomp(pemb)$x %>% as.data.frame()
p$cluster <- ifelse(rownames(p) %in% complex_patients,
                    "complex",
             ifelse(rownames(p) %in% infection_patients,
                    "infection", NA))
p$CENTRE <- nbr$CENTRE[match(rownames(p), nbr$STUDY_ID)]

# plots
left <- ggplot(p,
            aes(x=PC1, y=PC2, col=cluster)) +
  geom_point(size = 0.7) + 
  scale_colour_manual(values = infComp) +
  theme(text            = element_text(size = 8, family = "Times"),
        axis.text       = element_text(size = 8, family = "Times"),
        legend.key.size = unit(0.2, "cm"),
        legend.title    = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.2),
        legend.background = element_rect(fill="white",
                                         size=0.1,
                                         linetype="solid", 
                                         colour ="black")) 

right <- ggplot(p, aes(x = PC1,
                       y = PC2,
                       col = CENTRE,
                       shape = CENTRE)) +
    geom_point() +
  scale_fill_manual(values = centrescol) +
    scale_colour_manual(values = centrescol) +
    scale_shape_manual(values = c(1:11)) +
  theme(text            = element_text(size = 8, family = "Times"),
        axis.text       = element_text(size = 8, family = "Times"),
        legend.title    = element_blank(),
        legend.key.size = unit(0.2, "cm")) 



############# Layout ###########################################################
patchwork2 <-  left / right 
patchwork2 + plot_annotation(tag_levels = 'a')

ggsave("../result/Fig6/Fig6.tiff",
       width = 6,
       height = 4,
       units = "in",
       dpi = 1000)

############# Legend details ###################################################
t.test(p$PC1 ~ p$cluster)
t.test(p$PC2 ~ p$cluster)
summary(lm(PC1 ~ CENTRE + cluster, data=p))
summary(lm(PC2 ~ CENTRE + cluster, data=p)) 
