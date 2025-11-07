######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

cluster_demographics <- fread("../result/cluster_demographics.csv")
cr <- fread( "../result/bio_cluster_tests.csv") # cluster


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
# add a column matching the p value table
Fig5data$biological_measure <- gsub("_.*$",
                                    Fig5data$name,
                                    replacement = "")

p_value_data <- Fig5data %>%
  dplyr::filter(name %in% c(
                            "NFKB1",
                            "anyPathogenic")) %>% 
  # Group by the x-axis variable and the joining key
  group_by(long_name, biological_measure) %>%
  
  # Find the position of the end of the *longest* bar in the group
  summarize(y_pos = max(cluster_proportion), .groups = 'drop') %>%
  
  # Join the p-value data from 'cr'
  left_join(cr, by = "biological_measure") %>%
  
  # Format the p-value label for plotting
  # "%.2g" is good for scientific notation.
  mutate(p_label = sprintf("p = %.2g", BH_adjust_pvalue))


# Create a data frame for the bracket(s)
bracket_data <- data.frame(
  
  # 1. Specify the facet:
  # match 'Fig5data$category'
  category = c("phenotype","phenotype"), 
  
  # 2. Specify the start and end points (the 'long_name' values):
  xmin = c("SmB+", "CD21low normal"),
  xmax = c("SmB-", "CD21low high"),
  
  # 3. Specify the y-position (on the proportion axis):
  y.position = c(45, 40), 
  
  # 4. Specify the label for the bracket:
  label = c("p = 0.00017" , "p = 0.00002")
)



ggplot(Fig5data,
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
                              "any Pathogenic variant" = "Any Pathogenic Variant",
                              "SmB+" = "smB+",
                              "SmB-"  = "smB-" ,
                              "CD21low high" = expression("CD21"^lo),
                              "CD21low normal" = expression("CD211"^norm),
                              "NFKB1" = "NFKB1",
                              "any Pathogenic variant" = "Any Pathogenic Variant")) +
  ylab("Proportion of patients within the cluster") +
  xlab("") +
  coord_flip() + 
  theme_bw() +
  ylim(c(0, 45)) 
 
 
ggsave("../result/Fig5/Fig5.jpeg",
       width = 15,
       height = 10,
       units = "cm")

