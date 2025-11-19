######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################

# load presence/absence matrix
phbm <- readRDS("../result/patient_hpo_bio_mat.RDS") 

############# Part A ###########################################################
# take in the HPO data,create and save embeddings + HPO
# make an edge list: each row contains two columns describing a edge in the HPO
# onthology
# in preparation for https://github.com/eliorc/node2vec/blob/master/README.md


# subset to just the HPO terms and the patients with complete B cell
phm <- phbm[grep("^HP", rownames(phbm)), ] 

# keep HPO terms that are present at least once
phm <- phm[rowSums(phm) != 0, ]

# two columns: each row represents an edge
result <- do.call(rbind, lapply(row.names(phm), function(term) {
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
obsolete_terms <- row.names(phm)[!row.names(phm) %in%
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

phm <- phm[!rownames(phm) %in% obsolete_terms,]

############## save all
saveRDS(phm, 
        paste("../result/patient_hpo_bio_mat_for_reduc_dim_",
              thisfocus,".rds", sep = ""))

