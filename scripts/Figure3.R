######################## source common libraries ###############################
source("common.R")
library(ggsignif)
library(cluster)
library(ggpubr)
library(ggrepel)
library(vegan)
library(ggtext)

######################## import input ##########################################
# import violin data
similarity_df <- readRDS("similarity_df.RDS")
annotations <- readRDS("annotations.RDS")
pairwise_ks <- readRDS("pairwise_ks.RDS")
y_positions <- readRDS("y_positions.RDS")
pooled_label <- readRDS("pooled_label.RDS")
infCompBetween <- readRDS("infCompBetween.RDS")
group_order <- readRDS("group_order.RDS")
step <- readRDS("step.RDS")

# import qq-plot data
qq_df <- readRDS("qq_df.RDS")

# import MDS data
mds_df <- readRDS("mds_df.RDS")
InfComp <- readRDS("InfComp.RDS")

############# plotting ###########################################################

# Plot violin
Fig3a <- ggplot(similarity_df,
              aes(x = Group, y = Similarity, fill = Group)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggsignif::geom_signif(
    comparisons = annotations,
    annotations = pairwise_ks$label,
    y_position  = y_positions,
    tip_length  = 0.01,
    textsize    = 2.5,
    family      = "Times",
    vjust       = -0.2,
    parse = FALSE
  ) +
  labs(
    x       = "",
    y       = "Pairwise Lin's Similarity",
    caption = pooled_label
  ) +
  scale_fill_manual(values = infCompBetween) +
  scale_x_discrete(limits = group_order) +
  coord_cartesian(clip = "off") +
  expand_limits(y = max(y_positions) + step) +
  theme(
    text            = element_text(size = 8, family = "Times"),
    axis.text       = element_text(size = 8),
    legend.position = "none",
    plot.margin     = margin(t = 20, r = 5, b = 5, l = 5),
    plot.caption    = ggtext::element_markdown(size = 7, family = "Times",
                                               hjust = 0, face = "italic")
  )

# qq plot
Fig3b <- ggplot(qq_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey", linewidth = 0.5) +
  labs(
    x = "Between-Groups patient pairwise similarity Quantiles",
    y = "Pooled patient pairwise similarity Quantiles"
  ) +
  #xlim(0.2, 1)+ylim(0.2,1)+
  theme_classic() +
  theme(text = element_text(size = 8, family = "Times"))



# plot MDS colour-coded by groups
Fig3c <- ggscatter(mds_df, x = "V1", y = "V2",
                      color = "groups",
                      palette = InfComp,
                      size = 1,
                      ellipse = TRUE,
                      ellipse.type = "convex",
                      mean.point = TRUE,
                      mean.point.size = 5,
                      shape = "groups",
                      xlab = "MDS Dimension 1",
                      ylab = "MDS Dimension 2") + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.9,0.9)) +
  theme(text = element_text(size = 8, family = "Times")) 

# plot MDS colour-coded by centres
Fig3d <- ggscatter(mds_df, x = "V1", y = "V2",
                   color = "centres",
                   palette = centrescol,
                   size = 2,
                   shape = "centres",
                   xlab = "MDS Dimension 1",
                   ylab = "MDS Dimension 2") + 
  scale_shape_manual(values = c(1:11)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(text = element_text(size = 8, family = "Times")) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(ncol = 1))

############# Layout ###########################################################
top <- Fig3a + Fig3b
bottom <-Fig3c + Fig3d 
top / bottom + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(
    guides = 'collect'
  ) 


ggsave("../result/Figure3.tiff",
       width  = 8.5,
       height = 6,
       units  = "in",
       dpi    = 1000)

ggsave("../result/Figure3.jpeg",
       width  = 8.5,
       height = 6,
       units  = "in",
       dpi    = 1000)

############# Legend details ###################################################
similarity_df %>% 
  group_by(Group) %>% 
  dplyr::summarise(mean = mean(Similarity), n = n())
