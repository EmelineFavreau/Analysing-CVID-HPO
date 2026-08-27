######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# load patient clusters
infection_cluster <- fread("../result/InfectionBronchiectasisPatients.csv")
complex_cluster <- fread("../result/complexPatients.csv")
colnames(infection_cluster) <- "STUDY_ID"
colnames(complex_cluster) <- "STUDY_ID"

# Import data
# Lin's similarities
alllin <- readRDS("../result/all_sim.RDS")

# Deng
all_GICsim <- readRDS("../result/all_GICsim.RDS")

######################## analysis ##############################################

# aim: compare semantic similarities (Lin and Deng)

# infection patients in col
infpat <- infection_cluster$STUDY_ID

# complex patients in row
comppat <- complex_cluster$STUDY_ID

# compare Deng vs Lin in all pairwise comparisons
# combine pairwise similarities 
df <- data.frame()
for (i in c(infpat, comppat)){
  for (c in c(infpat, comppat)){
    lin = alllin[rownames(alllin) == c, colnames(alllin) == i]
    deng = all_GICsim[rownames(all_GICsim) == i,
                      colnames(all_GICsim) == c]
    df  <- rbind(df , c(i, c, lin, deng))
  }
}
# add names
colnames(df) <- c("patient1", "patient2", "lin", "deng")

# change class
df$lin <- as.numeric(df$lin)
df$deng <- as.numeric(df$deng)

# remove self-comparison
df <- df[df$patient1 != df$patient2,]

# deng created NAs, to be removed
df <- df[!is.na(df$deng),]
df <- df[!is.na(df$lin),]

# Correlation test (correlation normal distribution)
dfcorr <- cor.test(df$lin, df$deng, method = "pearson") # p-value < 2.2e-16

# plot
Figure_E3 <- ggplot(df, aes(x = lin, y = deng)) +
  geom_hex(bins = 100) +
  geom_smooth(method = lm,
              se = FALSE,
              col = 'red',linewidth = 0.5
  ) +# a line of best fit 
  labs(
    x = "Lin's similarity",
    y = "Deng similarity"
  ) +
  xlim(0, 1) + ylim(0, 1)+
  theme_classic() +
  theme(text = element_text(size = 8, family = "Times")) + 
  ggtitle(label = "Figure E3: All patient pairwise correlation",
          subtitle = paste("Pearson R =",
                           round(dfcorr$estimate, 3), "; Line of Best Fit in red"
          ))


############# save figure, stats ###############################################
ggsave("../result/SupplFigE3/Figure_E3.tiff",
       width = 3,
       height = 3.3,
       units = "in",
       dpi = 1000)

ggsave("../result/SupplFigE3/Figure_E3.jpg",
       width = 3,
       height = 3.3,
       units = "in",
       dpi = 1000)

