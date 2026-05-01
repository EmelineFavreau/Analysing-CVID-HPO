######################## source common libraries ###############################
source("common.R")
library(ggsignif)
library(cluster)
library(ggpubr)
library(ggrepel)
library(vegan)
library(ggtext)

######################## import input ##########################################

# load nbr
nbr <- readRDS("../result/tidy_data")

# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"


######################## Major Comment 3 - Violin Plot #########################
# calculating the pairwise Lin similarity between each patient's set of HPO 
# terms and categorizing these comparisons as (1) complex-to-complex, 
# (2) infection-to-infection, and (3) complex-to-infection. 
# The similarity value vs. categories could be plotted as a 
# histogram/violin plot and an appropriate statistical test could 
# evaluate the differences between categories (e.g. Wilcox or permutation) 
# to assess whether inter-category similarity is lower than intra-category similarity. 


# obtain minimal sets for each patient
minimal_term_sets <- list()
for(i in 1:nrow(nbr)){
  minimal_term_sets[[i]] <- minimal_set(hpo, nbr$hpo[[i]])
}

names(minimal_term_sets) <- nbr$STUDY_ID

######################## Step 1: calculate similarities#########################
# calculate a similarity matrix, containing pairwise term-set similarities:
calc_sim_mat <- function(patient){
  # subset to the group of patients
  minimal_term_sets <- minimal_term_sets[names(minimal_term_sets) %in% patient]
  
  # calculate IC from large corpus
  information_content <- descendants_IC(hpo)
  
  # pairwise Lin similarity 
  sim_mat <- get_sim_grid(ontology = hpo, term_sets = minimal_term_sets)
  
  # keep informative data
  #sim_mat[lower.tri(sim_mat, diag = T)] <- NA
  
  return(sim_mat)
}

# similarity for subgroups
infectionOnly_sim <- calc_sim_mat(infection_cluster$STUDY_ID)
complexOnly_sim <- calc_sim_mat(complex_cluster$STUDY_ID)
infectionVsComplex_sim <- calc_sim_mat(c(infection_cluster$STUDY_ID,
                                         complex_cluster$STUDY_ID))

# subset infectionVsComplex_sim to antagonist pairs
infectionVsComplex_simS <- 
  infectionVsComplex_sim[rownames(infectionVsComplex_sim) %in%
                           complex_cluster$STUDY_ID,
          colnames(infectionVsComplex_sim) %in% infection_cluster$STUDY_ID]

######################## Step 2: prep violin plot###############################
# plot a violin with x is the groups, y is Lin's similarity
inf <- as.vector(infectionOnly_sim[upper.tri(infectionOnly_sim)])
comp <- as.vector(complexOnly_sim[upper.tri(complexOnly_sim)])
cross <- as.vector(infectionVsComplex_simS)

#saveRDS(inf, "inf.RDS")
#saveRDS(comp, "comp.RDS")
#saveRDS(cross, "cross.RDS")

similarity_df <- tibble(
  Similarity = c(inf, comp, cross),
  Group = c(rep("Within-Infection", times = length(inf)),
            rep("Within-Complex", times = length(comp)),
            rep("Between-Groups", times = length(cross))))


# Define the fixed group order
group_order <- c("Within-Infection", "Within-Complex", "Between-Groups")

df <- similarity_df
groups <- unique(df$Group)
pairs <- combn(groups, 2, simplify = FALSE)

######################## Step 3: assess distribution ############################
results <- lapply(pairs, function(pair) {
  x <- df$Similarity[df$Group == pair[1]]
  y <- df$Similarity[df$Group == pair[2]]
  test <- ks.test(x, y)
  data.frame(
    group1    = pair[1],
    group2    = pair[2],
    statistic = test$statistic,
    p_value   = test$p.value,
    pooled    = FALSE
  )
})
ks_results <- bind_rows(results)

# also between and within
pooled_test <- ks.test(cross, c(inf, comp))

pooled_row <- data.frame(
  group1    = "Within-Infection + Within-Complex (pooled)",
  group2    = "Between-Groups",
  statistic = pooled_test$statistic,
  p_value   = pooled_test$p.value,
  pooled    = TRUE
)

ks_results <- bind_rows(ks_results, pooled_row)

ks_results <- ks_results %>%
  mutate(label = case_when(
    p_value < 0.001 ~ "p < 0.001",
    TRUE            ~ paste0("p = ", round(p_value, 3))
  ))

format_p_sci_expr <- function(p) {
  sci   <- formatC(p, format = "e", digits = 2)
  parts <- strsplit(sci, "e")[[1]]
  base  <- parts[1]
  exp   <- as.integer(parts[2])
  # Returns a plotmath expression string for use with parse(text = ...)
  paste0("p == ", base, " %*% 10^", exp)
}

ks_results <- ks_results %>%
  mutate(label = sapply(p_value, format_p_sci_expr))

y_max   <- max(df$Similarity, na.rm = TRUE)
y_range <- diff(range(df$Similarity, na.rm = TRUE))
step    <- y_range * 0.12

get_span <- function(g1, g2, order) {
  abs(match(g1, order) - match(g2, order))
}

pairwise_ks <- ks_results %>%
  dplyr::filter(!pooled) %>%
  dplyr::mutate(span = mapply(get_span, group1, group2,
                       MoreArgs = list(order = group_order))) %>%
  dplyr::arrange(span)

annotations <- mapply(c, pairwise_ks$group1, pairwise_ks$group2, SIMPLIFY = FALSE)
y_positions <- seq(y_max + step, by = step, length.out = nrow(pairwise_ks))


format_p_md <- function(p) {
  sci   <- formatC(p, format = "e", digits = 2)
  parts <- strsplit(sci, "e")[[1]]
  base  <- parts[1]
  exp   <- as.integer(parts[2])
  paste0("p = ", base, " &times; 10<sup>", exp, "</sup>")
}

pooled_label <- paste0(
  "Kolmogorov-Smirnov Test: Between-Groups vs. pooled (Within-Infection + Within-Complex): ",
  format_p_md(ks_results$p_value[ks_results$pooled])
)

# vector of colours
infCompBetween <- c(infComp, "Grey")
names(infCompBetween) <- c("Within-Complex","Within-Infection", "Between-Groups")


######################## Step 4: violon plot####################################
Fig <- ggplot(similarity_df,
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
    parse = TRUE
  ) +
  labs(
    x       = "",
    y       = "Pairwise Lin's Similarity",
    caption = pooled_label          # <-- pooled KS result here
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
  ) + 
ggtitle("Fig7: Lin pairwise similarities within & between groups")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig7.jpg")

Fig
# significant distribution difference:
# ks between-groups and within-groups is D = 0.061308, p-value < 2.2e-16)
# quantile-quantile plot
# each distribution is grouped in intervals with equal size
# most of the plot follows the 1:1 slope: identical similarities
# the most dissimilar pairs between-groups (minimum x=0.27)
# the most dissimilar pairs within-group (inf minimum y= 0.34)
# driver is at the tail end: the most dissimilar pairs are between groups

######################## Step 5: qq plot########################################
pooled <- c(inf, comp)
n      <- min(length(cross), length(pooled))
q_cross  <- quantile(cross,  probs = seq(0, 1, length.out = n))
q_pooled <- quantile(pooled, probs = seq(0, 1, length.out = n))

qq_df <- data.frame(y = q_pooled, x = q_cross)

# Plot
qqFig <- ggplot(qq_df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, colour = "grey", linewidth = 0.5) +
  labs(
    y = "Quantiles — Within-Infection + Within-Complex (pooled)",
    x = "Quantiles — Between-Groups"
  ) +
  xlim(0.2, 1)+ylim(0.2,1)+
  theme_classic() +
  theme(text = element_text(size = 8, family = "Times"))+ 
  ggtitle("Fig8: qqplot of Lin pairwise similarities within VS between groups")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig8.jpg")


qqFig
# the extremes are different, not the distribution in the middle part


######################## Step 6: MDS ########################################
# take the similarity matrix and plot another type of dimensionality reduction
# Deliverable: an answer and accompanying figure to visualise 
# how different the two groups are.

# HPO similarity matrix: infection vs complex patients
ic <- infectionVsComplex_sim

# Perform MDS analysis
mds_ic <- cmdscale(dist(ic))

# Add cluster information to the MDS results
mds_df <- as.data.frame(mds_ic)
mds_df$groups <- ifelse(rownames(mds_df) %in% infection_cluster$STUDY_ID, 
                         "Infection", "Complex")

# colours
InfComp <- infComp
names(InfComp) <- c("Complex", "Infection")

# Plot using ggscatter
mdsFig <- ggscatter(mds_df, x = "V1", y = "V2",
          color = "groups",
          palette = InfComp,
          size = 2,
          ellipse = TRUE,
          shape = "groups",
          ellipse.type = "convex",
          title = "Infection-vs-Complex Lin's Similarity Multidimensionality Scaling ",
          xlab = "MDS Dimension 1",
          ylab = "MDS Dimension 2") + 
  theme_classic() +
  theme(text = element_text(size = 8, family = "Times"))+ 
  ggtitle("Fig9: MDS of between-group Lin pairwise similarities")
ggsave("/Users/emelinefavreau/Library/CloudStorage/OneDrive-UniversityofCambridge/INTREPID/03-HPO-CVID/2026-04-JACI-D-26-00345/plots/Fig9.jpg")


# answer to comment: 
# The HPO set of each Infection patient has a Lin's similarity to
# any given set of Complex patient (in a similarity matrix).
# From this similarity matrix, a Euclidean distance matrix is calculated.
# This distance matrix is the input to a classical multidimensional scaling 
# (MDS) with 2 dimensions.
# The resulting plot shows an overlap between the two groups.

######################## Step 7: permutation ###################################
# take the current mds, permute the groups (ie MDS coordinates do not move) 
# and calculate the variance between each group. 
# repeat 999 times to test the null hypothesis of no variance between groups

set.seed(42)

# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
permanova_result <- adonis2(mds_df[, c("V1", "V2")] ~ groups, 
        data = mds_df, 
        permutations = 999,
        method = "euclidean")
permanova_result
#                  Df            SumOfSqs    R2          F        Pr(>F)
# Model            1             18.32       0.02903     15.636   0.001 ***
# Residual         523           612.81      0.97097                    
# Total            524           631.13      1.00000          

# P value is < 0.05. Reject Null hypothesis; groups are in the different MDS space.
# R2 = 0.02903. group-membership explains 29% of variation.
# F = 15.636. Between-group / within-group variance is above 1, thus not by chance.
# Conclusion: centroids are different between groups.

