######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load presence/absence matrix
patient_hpo_bio_mat <- readRDS("../result/patient_hpo_bio_mat.RDS") 
# 925 HPO+other, 528 patients



############# Part A ###########################################################
# a: take in the HPO data,  create and save embeddings + HPO
# make a edge list: each row contains two columns describing a edge in the HPO
# onthology
# in preparation for https://github.com/eliorc/node2vec/blob/master/README.md


# subset to just the HPO terms and the patients with complete B cell
patient_hpo_mat <- patient_hpo_bio_mat[grep("^HP",
                                            rownames(patient_hpo_bio_mat)), ] 

# keep HPO terms that are present at least once
patient_hpo_mat <- patient_hpo_mat[rowSums(patient_hpo_mat) != 0, ]
# dim(patient_hpo_mat)
# 889 HPO terms, 528 patients

# two colums: each row represents an edge
result <- do.call(rbind, lapply(row.names(patient_hpo_mat), function(term) {
  connectors <- c(hpo$ancestors[[term]], hpo$children[[term]])
  connectors <- unique(na.omit(connectors))  # remove duplicates + NAs
  connectors <- connectors[connectors != term]  # remove self-links
  
  if (length(connectors) > 0) {
    data.frame(term = term,
               connector = connectors,
               stringsAsFactors = FALSE)
  }
}))

# Get all unique HPO codes
all_hpos <- unique(c(result$term, result$connector))

# obsolete terms:
obsolete_terms <- row.names(patient_hpo_mat)[!row.names(patient_hpo_mat) %in%
                                               all_hpos]

# Create a lookup table mapping HPO -> integer ID
hpo_id <- setNames(seq_along(all_hpos), all_hpos)

# Make a copy of result and add numeric IDs
result_with_id <- result
result_with_id$term_id <- hpo_id[result$term]
result_with_id$connector_id <- hpo_id[result$connector]
thisfocus <- "no_Bcell_filtering"
fwrite(result_with_id, 
       paste("../result/reduc_dim_input_with_id_", thisfocus, ".csv", sep = ""))
fwrite(result_with_id[,3:4],
       paste("../result/HPO_id_", thisfocus, ".tsv", sep = ""),
       col.names = F, sep = "\t")

patient_hpo_mat <- patient_hpo_mat[!rownames(patient_hpo_mat) %in% obsolete_terms,]
saveRDS(patient_hpo_mat, 
        paste("../result/patient_hpo_bio_mat_for_reduc_dim_", thisfocus,".rds", sep = ""))

