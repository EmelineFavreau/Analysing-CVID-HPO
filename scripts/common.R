# get libraries
library(biomaRt)
library(colorspace) 
library(conflicted)
library(data.table)
library(dplyr)
library(epitools)
library(fastmatch)
library(ggpubr)
library(ggplot2)
library(ontologyIndex)
library(ontologySimilarity)
if(system("whoami",intern=TRUE)!="chrisw")
  library(ontologyPlot)
library(paletteer)
library(patchwork)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape)
library(rstatix)
library(SimReg)
library(stringr)
library(tibble)
library(tidyr)
library(tidytext)
library(tidyverse)
library(uwot)

# Ontology with 18964 terms
# v2024-04-04
hpo <- get_OBO("../input/hp.obo",
               merge_equivalent_terms = FALSE) 


theme_set(theme_bw()) +
  theme(text = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12))

centres_vec <- c(
  CUH      = "Cambridge University Hospitals",
  Glasgow  = "Queen Elizabeth University Hospital Glasgow",
  Midlands = "University Hospital of North Midlands",
  BH       = "Barts Health",
  RF       = "Royal Free",
  FP       = "Frimley Park",
  Papworth = "Papworth",
  SR       = "Salford Royal",
  HE       = "Heart of England",
  LTH      = "Leeds Teaching Hospitals",
  ICH      = "Imperial College Healthcare"
)

# Fig 1
training_colours <- c("Pre-Training" = "#f0f0f0",
                   "Post-Training" = "#636363")
cd21col <- c("expansion of low CD21" = "#E56669FF",
             "CD21 normal" = "#B23539FF")
smbcol <- c("SmB plus" = "#FAB57CFF",
            "SmB minus" = "#FFCBA6FF") 
trcol <- c("Transitional B normal" = "#FFEF98FF", 
           "Transitional B high"  = "#CCAA3E")
Bcellcol <- c(cd21col, smbcol, trcol)
smBPLUSCD21col <- c("#73652DFF", "#9F905FFF")
smBMINUSCD21col <- c("#E79498FF", "#FFAEB1FF")
smBMINUSTrcol <- c("#514289FF", "#8275B8FF")


# Fig 4:
presAbsvalues <-  c("0" = "grey90",
           "1" = "black")
clinLab <-  c("clinical" = "#b2df8a",
                    "lab" = "#1f78b4")
infComp <- c("complex"="#fb9ac7","infection"="#084594")


# Fig 6: PCA
# 11 centres in colours
centrescol <- c("#5F4690FF", "#666666FF","#1D6996FF", 
                "#0F8554FF",
  "#73AF48FF", "#EDAD08FF", "#E17C05FF", "#CC503EFF",
  "#94346EFF", "black", "#994E95FF" )
