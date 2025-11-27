######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# cols: name,cluster,patient_number,cluster_proportion,long_name,category
cluster_demographics <- fread("../result/cluster_demographics.csv")

# cols:cluster,biological_category,patient_count,
# total_patient_in_biological_measure,percentage_in_biological_measure,
# biological_measure,BH_adjust_pvalue
cr <- fread( "../result/cluster_biomarker_tests_ratio_pvalue.csv") 


############# Fig5 ############################################################


cr <- cr %>% dplyr::filter(biological_measure %in% c("smb",
                                               "cd21",
                                               "NFKB1",
                                               "anyPathogenic"))

# filter for the significant associations only
# smb, cd21, NFKB1, anyPathogenic
Fig5data <- cluster_demographics %>% 
  dplyr::filter(name %in% c("smb_normal",
                                 "smb_minus",
                                 "cd21_low",
                                 "cd21_plus" ,
                                 "NFKB1",
                                 "anyPathogenic")) 

# add a column matching the p value table (cr$biological_measure)
Fig5data$biological_measure <- gsub("_.*$",
                                    Fig5data$name,
                                    replacement = "")

p_value_data <- Fig5data %>%
  dplyr::filter(name %in% c(
                            "NFKB1",
                            "anyPathogenic")) %>% 
  group_by(long_name, biological_measure) %>%
  summarize(y_pos = max(cluster_proportion), .groups = 'drop') %>%
  left_join(cr, by = "biological_measure") %>%
  mutate(p_label = sprintf("p = %.2g", BH_adjust_pvalue))

psmb <- unique(cr$BH_adjust_pvalue[cr$biological_measure %in% c("smb")]) %>% 
  format(., scientific = FALSE, digits = 1)
pcd21<- unique(cr$BH_adjust_pvalue[cr$biological_measure %in% c("cd21")]) %>% 
  format(., scientific = FALSE, digits = 1)

# Create a data frame for the bracket(s)
bracket_data <- data.frame(
  category = c("phenotype","phenotype"), 
  xmin = c("SmB+", "CD21low normal"),
  xmax = c("SmB-", "CD21low high"),
  y.position = c(45, 40), 
  label = c(paste("p = ",
            psmb,
            sep =""),
            paste("p = ", 
            pcd21,
            sep =""))
)



Fig5 <- ggplot(Fig5data,
       aes(fill = cluster,
           y = cluster_proportion,
           x = long_name)) + 
  geom_bar(position = "dodge",
           stat = "identity") + 
  geom_text(data = p_value_data,
    aes(x = long_name,
      y = y_pos + 0.02,
      label = p_label),
    inherit.aes = FALSE, 
    vjust = 0,       
    hjust = -0.5) +
  geom_bracket(
    data = bracket_data,
    aes(xmin = xmin,
      xmax = xmax,
      y.position = y.position,
      label = label),
    inherit.aes = FALSE,
    tip.length = 0.015,
    vjust = 3,
    hjust = 1) +
  scale_fill_manual(values = infComp) +
  scale_x_discrete(labels = c("SmB+" = "smB+" ,
                              "SmB-" = "smB-",
                              "CD21low high" = expression("CD21"^lo),
                              "CD21low normal" = expression("CD211"^norm),
                              "NFKB1" = "NFKB1",
                              "any Pathogenic" = "Any Pathogenic Variant",
                              "SmB+" = "smB+",
                              "SmB-"  = "smB-" ,
                              "CD21low high" = expression("CD21"^lo),
                              "CD21low normal" = expression("CD211"^norm),
                              "NFKB1" = "NFKB1",
                              "any Pathogenic" = "Any Pathogenic Variant")) +
  ylab("Proportion of patients within the cluster") +
  xlab("") +
  coord_flip() + 
  theme_bw() +
  ylim(c(0, 45)) 

############# Layout ###########################################################
Fig5 
ggsave("../result/Fig5/Fig5.jpeg",
       width = 15,
       height = 10,
       units = "cm")

############# Legend details ###################################################


