######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
sim_mat <- fread("../result/Fig1/sim_mat.csv")
num_hpo_compL <- fread("../result/Fig1/num_hpo_compL.csv")
num_hpo_comp <- fread("../result/Fig1/num_hpo_comp.csv")
############# Fig1A ############################################################

colnames(sim_mat) <- gsub(x = colnames(sim_mat),
                          pattern = "V",
                          replacement = "Clinician #")
df <- reshape::melt.array(as.matrix(sim_mat))
colnames(df) <- c("x", "y", "Semantic Similarity")
df$`Semantic Similarity` <- as.numeric(substr(df$`Semantic Similarity`, 1, 4))

# colour vector for plotting
RdYlBu_pal <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                  "RdYlBu")))(100)


# save figure
Fig1A <- ggplot(df, aes(x = x,
               y = y,
               fill = `Semantic Similarity`)) +
  geom_tile(color = "white") +
  geom_text(aes(label = `Semantic Similarity`),
            color = "black", size = 3) +
  scale_fill_gradientn(
    colours = RdYlBu_pal,
    na.value = "white",
    breaks = c(min(sort(as.numeric(df$`Semantic Similarity`))),
              max(sort(as.numeric(df$`Semantic Similarity`)))),
    labels = c("Not similar", "Very similar")) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1.5),
        legend.position.inside = c(0.95, 0.05),  # Position legend at the bottom-right
        legend.justification = c("right", "top")) +
  labs(x = NULL, y = NULL,
       fill = "Semantic Similarity")




############# Fig1B ############################################################
# Fig1: violon

num_hpo_compL$Timepoint <- ifelse(num_hpo_compL$Timepoint == "before",
                                  "Pre-Training",
                                  "Post-Training")
Fig1B <- ggplot(num_hpo_compL,
                aes(x = Timepoint, y = Records, fill = Timepoint)) +
  geom_violin(alpha = 0.5, trim = FALSE) +  
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  labs(x = "",
       y = "Number of HPOs per patient") +
  scale_fill_manual(values = training_colours) +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c("Pre-Training","Post-Training"))




############# Fig1C ############################################################
# Fig1C


Fig1C <- ggplot(num_hpo_comp,
                aes(x = before_after_similarity_value,
                    y = differential)) +
  geom_point(size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(
    x = "Pre/Post-Training HPO Semantic Similarity",
    y = "Post - Pre Training HPO Count") +
  geom_density_2d_filled(alpha = 0.5) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1),
                     limits = c(0.2, 1.1)) + 
  ylim(c(-20, 33)) + theme(legend.position = "none")


############# Layout ###########################################################


plots <- cowplot::align_plots(Fig1A, Fig1B,
                              align = 'v', axis = 'l')
bottom_row <- cowplot::plot_grid(plots[[2]], Fig1C,
                                 labels = c('b', 'c'), label_size = 12)
Fig1 <- cowplot::plot_grid(plots[[1]], bottom_row, 
                           labels = c('a', ''), label_size = 12, ncol = 1)


ggsave("../result/Fig1/Fig1.jpeg",
       width = 20,
       height = 20,
       units = "cm")

############# Legend details ###################################################
summary(num_hpo_compL$Records[num_hpo_compL$Timepoint == "Post-Training"])
IQR(num_hpo_compL$Records[num_hpo_compL$Timepoint == "Post-Training"])
