# source common.R, which contains common function, data, libraries
source("common.R")


############# prep ############################################################
# load nbr
nbr <- readRDS("../result/tidy_data")

# Each row describes a patient carrying a rare TNFRSF13B variation
trp <- read_excel("../input/Genotype_Final.xlsx", 
                                      sheet = "Summary", range = "B1:B8")

# Each row describes a patient carrying a canonical TNFRSF13B variation
tcp <- read_excel("../input/Genotype_Final.xlsx", 
                                     sheet = "Summary", range = "A1:A31")

# Each row describes a patient carrying NFKB1 variation 
nfp <- read_excel("../input/Genotype_Final.xlsx", 
                             sheet = "Summary", range = "C1:C20")

# Each column lists patients with other diagnostic variants than TACI or NFKB1
pp <- 
  read_excel("../input/Genotype_Final.xlsx", 
             sheet = "Summary", range = "D1:D26")

# Each row describes a patient, gene, variant, protein
pgvp <- read_excel("../input/Table_3_genetics.xlsx")

############# tidy ############################################################
trp <- trp[trp$TNFRsf13B_rare %in% nbr$STUDY_ID,]
tcp <- tcp[tcp$TNFRSF13B_canonical %in% nbr$STUDY_ID,]
nfp <- nfp[nfp$NFKB1 %in% nbr$STUDY_ID,]
pp <- pp[pp$Other %in% nbr$STUDY_ID,]

si <- c(trp$TNFRsf13B_rare,
  tcp$TNFRSF13B_canonical,
  nfp$NFKB1,
  pp$Other)


dg <- c(rep("Rare TNFRSF13B", times = nrow(trp)),
        rep("Canonical TNFRSF13B", times = nrow(tcp)),
        rep("NFKB1", times = nrow(nfp)),
        rep("other", times = nrow(pp)))

tt <- tibble(si = si,
             Patient = paste("Patient ", 1:length(si), sep = ""),
       Sex = nbr$SEX[match(si, nbr$STUDY_ID)],
       Diagnosis_gene = ifelse(dg == "other",
                        nbr$`Diagnosis gene`[match(si, nbr$STUDY_ID)],
                        dg))

tochange <- grep("^HGNC:", tt$Diagnosis_gene, value = TRUE)

# update to the valid code using biomart; get the database 
ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")

# input a vector, get a df with 2 columns
query_genes <- function(x) getBM(
  attributes = c("hgnc_symbol", "hgnc_id"),
  filters = "hgnc_id",
  values = x,
  mart = ensembl)

# two columns: gene symbol, gene id
genes_with_symbol_df <- query_genes(tochange)

for(g in 1:nrow(tt)){
  gg <- tt$Diagnosis_gene[g]
  if(gg %in% genes_with_symbol_df$hgnc_id){
    tt$Diagnosis_gene[g] <- 
      genes_with_symbol_df$hgnc_symbol[genes_with_symbol_df$hgnc_id == gg]
  }
} 

# Patient, sex
# subset pgvp to tt patients
pgvu <- pgvp[pgvp$Patient %in% tt$si, ]
pgvu$sex <- tt$Sex[match(pgvu$Patient, tt$si)] 
pgvf <- rbind(pgvu, c(tt$si[!tt$si %in% pgvu$Patient],
          tt$Diagnosis_gene[!tt$si %in% pgvu$Patient], NA, NA,
          tt$Sex[!tt$si %in% pgvu$Patient]))
pgvf$Patient <- 1:nrow(pgvf)


############# save all #########################################################
fwrite(pgvf, "../result/Table2/Table2.csv")
