######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# HPO categorised as lab or clinical, each term in 10 or more patients
loc <- read.csv("../result/HPO_freq_name_labORclinical.csv")

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters 
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")

# load patient/phenotype matrix, hierarchical clusters
newz <- readRDS("../result/newz.RDS")
thc <- readRDS("../result/thc.RDS")
phc <- readRDS("../result/phc.RDS")

############# Fig4 ############################################################
figname <- "Fig4"

ClinicalGroups = ifelse(colnames(newz) %in%
                           infection_cluster$infectionBronchiectasis,
                         "Infection",
                         ifelse(colnames(newz) %in% 
                                  complex_cluster$complexPatients,
                                "Complex", NA))

# Data that maps to the heatmap should be set at the rownames
annot_col <- data.frame(row.names = colnames(newz),
                        ClinicalGroups)

# 7 key phenotypes; 525 patients
Fig4a <- pheatmap(newz,
         color = c("#f0f0f0", "#636363"),
         show_colnames = FALSE,
         cluster_rows = thc,
         cluster_cols = phc,
         cutree_cols = 2,
         cutree_rows = 1,
         legend = FALSE,
         fontsize = 8,
         fontfamily = "Times",
         labels_row = c("Infection",
                        "Bronchiectasis",
                        "Granuloma",
                        "Splenomegaly",
                        "NRH of the liver",
                        "Lymphoproliferation",
                        "Autoimmunity"),
         annotation_col = annot_col,
         annotation_colors = list(ClinicalGroups = c("Complex" = "#fb9ac7",
                               "Infection" = "#084594")),
         annotation_legend = FALSE)


# Get the column order from the heatmap object
col_order <- Fig4a$tree_col$order  # dendrogram leaf order
col_labels <- colnames(newz)[col_order]  # patient IDs in display order

# Get cluster assignments (cutree_cols = 2)
col_clusters <- cutree(Fig4a$tree_col, k = 2)

# See which patients are in each cluster
cluster1_patients <- names(col_clusters[col_clusters == 1])
cluster2_patients <- names(col_clusters[col_clusters == 2])

# keep infection/complex patients
data <- nbr %>% 
  dplyr::filter(STUDY_ID %in% c(infection_cluster$infectionBronchiectasis, 
                                complex_cluster$complexPatients)) 

# each term in 25 or more patients
loc25 <- loc %>% dplyr::filter(HPO_freq > 25)

# keep HPO categorised as lab or clinical
data$hpo <- sapply(data$hpo, function(x) x[x %in% loc25$HPO_code])

# a column with the HPO terms
data$hpo_long <- sapply(data$hpo,
                        function(y) hpo$name[match(unlist(y), names(hpo$name))])

hpo_long <- unlist(data$hpo_long)

unique_long <- unique(hpo_long)

## create a matching patient id vector
patient <- rep(data$STUDY_ID, 
               sapply(data$hpo, length))

## make empty matrix - default value is 0
patient_hpo_mat <- matrix(data = 0,
                          nrow = length(loc25$HPO_code),
                          ncol = nrow(data),
                          dimnames = list(unique_long, data$STUDY_ID))

## add 1 where the term matches a patient
patient_hpo_mat[cbind(hpo_long, patient)] <- 1
wide_df <- as.data.frame(patient_hpo_mat)

# Assign HPO to Category
LabClinical_info <- data.frame(HPO_term = rownames(patient_hpo_mat),
  LabClinical = loc25$Category[match(rownames(patient_hpo_mat),
                                     loc25$HPO_name)])

# Assign Patient to cluster
InfectionComplex_info <- data.frame(Patient = colnames(patient_hpo_mat),
  InfectionComplex = ifelse(colnames(patient_hpo_mat) %in%
                     infection_cluster$infectionBronchiectasis,
                       "infection", "complex"))

###### plot lab values #########################################################
fn <- "lab"
# make a long df: HPO_code | HPO_term  | Patient | Presence |  
# InfectionComplex |  LabClinical
final_long_df <- wide_df %>%
  rownames_to_column(var = "HPO_term") %>%
  pivot_longer(
    cols = -HPO_term,
    names_to = "Patient",
    values_to = "Presence") %>%
  left_join(LabClinical_info, by = "HPO_term") %>%
  left_join(InfectionComplex_info, by = "Patient") %>%
  dplyr::filter(LabClinical == fn)

# shorten two terms
long_term_t <- final_long_df$HPO_term[grep("Complete or near-complete",
                                           final_long_df$HPO_term)]
names(long_term_t) <- gsub("Complete or near-complete absence of specific",
                           long_term_t, replacement = "Absence* of")

final_long_df$HPO_term[grep("Complete or near-complete",
                            final_long_df$HPO_term)] <- names(long_term_t)

rownames(wide_df)[rownames(wide_df) %in% long_term_t] <-
  unique(names(long_term_t))

HPO_order <- rownames(wide_df)[order(rowSums(wide_df))]

final_long_df <- final_long_df %>% 
  mutate(HPO_term = factor(HPO_term, levels = HPO_order))


# HPO_term, Count, Group (Infection or Complex)
plot_data_right <- final_long_df %>%
  dplyr::filter(Presence == 1) %>%
  group_by(HPO_term, InfectionComplex) %>%
  summarise(Count = n(), .groups = "drop")

Fig4c <- ggplot(plot_data_right, 
                    aes(x = HPO_term,
                        y = Count,
                        fill = InfectionComplex),
                  ) +
    geom_bar(stat = "identity",
             #fill = InfectionComplex,
             show.legend =  TRUE) +
    scale_fill_manual(values = infComp) +
    coord_flip() + 
    labs(x = "", y = "Total Count") +
    scale_y_continuous(breaks = c(0, 250, 500))+
    theme_minimal() + 
    theme(panel.grid         = element_blank(),
          strip.text.y       = element_blank(), 
          strip.background   = element_blank(),
          axis.ticks.y       = element_blank(),
          panel.grid.major.x = element_line(color = "lightgrey",
                                            linewidth = 0.2,
                                            linetype = 2),
          legend.title       = element_blank(),
          text               = element_text(size = 8, family = "Times"),
          axis.text          = element_text(size = 8, family = "Times"))
  




###### plot clinical values #########################################################
fn <- "clinical"

# make a long df: HPO_code | HPO_term  | Patient | Presence |  
# InfectionComplex |  LabClinical
final_long_dfC <- wide_df %>%
  rownames_to_column(var = "HPO_term") %>%
  pivot_longer(
    cols = -HPO_term,
    names_to = "Patient",
    values_to = "Presence") %>%
  left_join(LabClinical_info, by = "HPO_term") %>%
  left_join(InfectionComplex_info, by = "Patient") %>%
  dplyr::filter(LabClinical == fn)

# shorten two terms
long_term_t <- final_long_dfC$HPO_term[grep("Complete or near-complete",
                                           final_long_dfC$HPO_term)]
names(long_term_t) <- gsub("Complete or near-complete absence of specific",
                           long_term_t, replacement = "Absence* of")

final_long_dfC$HPO_term[grep("Complete or near-complete",
                            final_long_dfC$HPO_term)] <- names(long_term_t)


rownames(wide_df)[rownames(wide_df) %in% long_term_t] <-
  unique(names(long_term_t))

HPO_order <- rownames(wide_df)[order(rowSums(wide_df))]

final_long_dfC <- final_long_dfC %>% 
  mutate(HPO_term = factor(HPO_term, levels = HPO_order))


# HPO_term, Count, Group (Infection or Complex)
plot_data_right <- final_long_dfC %>%
  dplyr::filter(Presence == 1) %>%
  group_by(HPO_term, InfectionComplex) %>%
  summarise(Count = n(), .groups = "drop")

Fig4b <- ggplot(plot_data_right, 
                aes(x = HPO_term,
                    y = Count,
                    fill = InfectionComplex),
) +
  geom_bar(stat = "identity",
           #fill = InfectionComplex,
           show.legend = FALSE) +
  scale_fill_manual(values = infComp) +
  coord_flip() + 
  labs(x = "", y = "Total Count") +
  scale_y_continuous(breaks = c(0, 100, 200))+
  theme_minimal() + 
  theme(panel.grid         = element_blank(),
        strip.text.y       = element_blank(), 
        strip.background   = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey",
                                          linewidth = 0.2,
                                          linetype = 2),
        legend.title       = element_blank(),
        text               = element_text(size = 8, family = "Times"),
        axis.text          = element_text(size = 8, family = "Times"))


############# Layout ###########################################################

# align all plots vertically
plots <- cowplot::align_plots(Fig4a$gtable,
                              Fig4b,
                              Fig4c,
                              align = 'v', # vertically
                              axis = 'l') # by the left

# put together the bottom row and then everything
bottom_row <- plot_grid(
  plots[[2]], plots[[3]],
  labels = c("b", "c"),
  rel_widths = c(1, 1, .3),
  nrow = 1,
  label_size = 8,
  label_fontfamily = "Times"
)

plot_grid(plots[[1]],
          bottom_row,
          labels = c("a"),
          ncol = 1,
          rel_heights = c(1, 2),
          label_size = 8,
          label_fontfamily = "Times")


ggsave("../result/Figure4.tiff",
       width  = 10,
       height = 6,
       units  = "in",
       dpi    = 1000,
       bg     = "white") # the annotation background is white

ggsave("../result/Figure4.hamming.jpeg",
       width  = 10,
       height = 6,
       units  = "in",
       dpi    = 1000,
       bg     = "white") # the annotation background is white


############# Legend details ###################################################
# number of patients
nrow(data)
# number of HPO
length(unique(final_long_df$HPO_term)) # lab=33 lab
length(unique(final_long_dfC$HPO_term)) # lab=29 Clinical
