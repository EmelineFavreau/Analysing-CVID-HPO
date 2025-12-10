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

############# Fig4 ############################################################
figname <- "Fig4"
fn <- "lab"

# keep infection/complex patients
data <- nbr %>% 
  dplyr::filter(STUDY_ID %in% c(infection_cluster$infectionBronchiectasis, 
                                complex_cluster$complexPatients)) 

loc25 <- loc %>% dplyr::filter(HPO_freq > 25)

# keep HPO categorised as lab or clinical, each term in 25 or more patients
data$hpo <- sapply(data$hpo, function(x) x[x %in% loc25$HPO_code])

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
                       "InfectionBronchiectasis", "Complex"))

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

p_main <- final_long_df %>%
  ggplot(aes(x = HPO_term,
             y = Patient,
             fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  facet_grid( . ~ InfectionComplex, 
              scales = "free", 
              space = "free") +
  scale_fill_manual(values = presAbsvalues,
                    labels = c("Absent", "Present"),
                    name = "") +
  theme(panel.grid       = element_blank(),
        axis.text.x      = element_blank(), 
        axis.ticks.x     = element_blank(),
        strip.background = element_rect(fill = "grey85", color = NA),
        strip.text       = element_text(size = 6, family = "Times"),
        legend.position  = "bottom",
        text             = element_text(size = 6, family = "Times"),
        axis.text        = element_text(size = 6, family = "Times"),
        legend.key.size  = unit(0.2, "cm"),
        legend.text      = element_text(size=7)) +
  labs(x = "", y = "Patient") + 
  coord_flip()


plot_data_right <- final_long_df %>%
  dplyr::filter(Presence == 1) %>%
  group_by(HPO_term) %>%
  summarise(Count = n(), .groups = "drop")



p_right <- ggplot(plot_data_right, 
                  aes(x = HPO_term,
                      y = Count)) +
  geom_bar(stat = "identity",
           fill = "grey50",
           show.legend = FALSE) +
  coord_flip() + 
  labs(x = "", y = "Total Count") +
  scale_y_continuous(breaks = c(0, 250, 500))+
  theme_minimal() + 
  theme(panel.grid         = element_blank(),
        strip.text.y       = element_blank(), 
        strip.background   = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey",
                                          linewidth = 0.2,
                                          linetype = 2),
        text               = element_text(size = 8, family = "Times"),
        axis.text          = element_text(size = 7, family = "Times"))


Fig4 <- p_main + p_right +
  plot_layout(widths = c(3, 1))

############# Layout ###########################################################
Fig4
ggsave(paste("../result/", figname, "/", figname, ".tiff", sep = ""),
       width = 6,
       height = 4,
       units = "in",
       dpi = 1000)

############# Legend details ###################################################
# number of patients
nrow(data)
# number of HPO
length(HPO_order)
