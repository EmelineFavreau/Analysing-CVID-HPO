######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# missing terms are designed by Luiza
# the first data frame describe the new terms (term and code) and their parents
mhtp <- 
  read_excel("../input/Missing_HPO_Terms_Proposal_LC_paper.xlsx", 
             sheet = "parents")


############# Suppl Table New terms DAG ####################################
############## creating new HPO codes
# vector of the new terms 
missing_HPO_terms <- mhtp %>% 
  dplyr::filter(ProposedCode == "ProposedCode") %>% 
  dplyr::select(ProposedTerm) %>% 
  unlist() %>% unique()

# create temporary missing codes
nums <- seq(from = 9999901,
            length.out = length(missing_HPO_terms))

missing_HPO_codes <- sprintf("HP:%07d", nums)

# assign them randomly
names(missing_HPO_codes) <- missing_HPO_terms

# input them in the column for Proposed Code and Parent Code
for(i in 1:length(missing_HPO_codes)){
  tt <- missing_HPO_codes[i]
  mhtp$ProposedCode[mhtp$ProposedTerm == names(tt)] <- tt
  mhtp$ProposedParentCode[mhtp$ProposedParentTerm == names(tt)] <- tt
}
# subset of new terms with new terms as parents 
missing_HPO_terms_on_new_branch <- mhtp %>% 
  dplyr::filter(ProposedParentCode %in% missing_HPO_codes) %>% 
  dplyr::select(ProposedParentTerm) %>% 
  unlist() %>% unique()

# add code to terms 
names(missing_HPO_terms_on_new_branch) <- 
  missing_HPO_codes[match(missing_HPO_terms_on_new_branch,
                          names(missing_HPO_codes))]

############## creating hpo database elements with new terms
# 1. missing_HPO_parents: 
# each element is a character vector of HPO codes of the direct parents, (vector is direct parent)
# and is named with the HPO code (name is child)
missing_HPO_parents <- lapply(unname(missing_HPO_codes),
   function(x) unique(mhtp$ProposedParentCode[grep(x, mhtp$ProposedCode)]))
names(missing_HPO_parents) <- missing_HPO_codes

# 2. missing_HPO_children: 
# each element is a character vector of HPO codes of the children, (vector is children)
# and is named with the HPO code (name is the parent)
missing_HPO_children <- lapply(unname(missing_HPO_codes),
                               function(x) character()) # most of new terms are branch-ending, no children
names(missing_HPO_children) <- missing_HPO_codes

# updating missing_HPO_children with missing_HPO_terms_on_new_branch:
for(i in 1:length(missing_HPO_terms_on_new_branch)){
  # find the code for the parent
  cc <- names(missing_HPO_terms_on_new_branch)[i]
  # find the codes for the children
  pc <- mhtp$ProposedCode[mhtp$ProposedParentCode %in% cc]
  # add the code of the child, when the element is the parent
  j <- grep(cc, names(missing_HPO_children))
  missing_HPO_children[[j]] <- pc
}




# 3. missing_HPO_ancestors: 
# list: each element is a character vector of HPO codes of all ancestors, (vector is list of all ancestors)
# including the direct parents, and is named with the HPO code (name is any HPO)
missing_HPO_ancestors <- lapply(missing_HPO_parents,
                                function(x) get_ancestors(hpo, x))

for(i in 1:length(missing_HPO_codes)){
  occ <- grep(missing_HPO_codes[i],
              names(missing_HPO_ancestors),
              fixed = TRUE) 
  pn <- mhtp$ProposedParentCode[grep(missing_HPO_codes[i],
                                     mhtp$ProposedCode,
                                     fixed = TRUE)]
  au <- unique(c(missing_HPO_ancestors[[occ]], pn))
  missing_HPO_ancestors[[occ]] <- unique(c(missing_HPO_ancestors[[occ]], pn))
}


# in case of grand child, add grandparent
for(i in 1:length(missing_HPO_terms_on_new_branch)){
  nn = names(missing_HPO_terms_on_new_branch[i]) # code of grandchild
  np = mhtp$ProposedParentCode[grep(nn, 
                                    mhtp$ProposedCode,
                                    fixed = TRUE)] # codes of parents
  gp = missing_HPO_ancestors[grep(nn,
                                  names(missing_HPO_ancestors),
                                  fixed = TRUE)] %>%
    unlist() %>% unname()# codes of grandchild's ancestors before inputing
  an = unique(c(gp, np))
  ai = match(nn, names(missing_HPO_ancestors))
  missing_HPO_ancestors[[ai]] <- an
}

# in case of grand child, add grandparent
for(i in 1:length(missing_HPO_ancestors)){
  if(sum(missing_HPO_codes %in% missing_HPO_ancestors[[i]]) > 0){
    print(names(missing_HPO_ancestors)[i])
    # add ancestors of missing_HPO_codes
    mm <- missing_HPO_codes[missing_HPO_codes %in% missing_HPO_ancestors[[i]]]
    ai <- match(mm, names(missing_HPO_ancestors))
    ua <- missing_HPO_ancestors[[ai]]
    missing_HPO_ancestors[[i]] <- unique(c(missing_HPO_ancestors[[i]], ua))
  }
}

# logical vector, each entry is false or true, and named with the HPO code
# none of the missing terms are obsolete
missing_HPO_obsolete <- rep(FALSE, length(missing_HPO_codes))
names(missing_HPO_obsolete) <- missing_HPO_codes

# Update the hpo with the new terms
hpo_with_missing_terms <- hpo


# named character: each entry is HPO code, and named with the HPO code
missing_HPO_codes1 <- missing_HPO_codes
names(missing_HPO_codes1) <- missing_HPO_codes1
hpo_with_missing_terms$id <- c(hpo_with_missing_terms$id,
                               missing_HPO_codes1)

# named character: each entry is HPO term, and named with the HPO code
missing_HPO_codes2 <- names(missing_HPO_codes)
names(missing_HPO_codes2) <- missing_HPO_codes
hpo_with_missing_terms$name <- c(hpo_with_missing_terms$name,
                                 missing_HPO_codes2)

# list: each element is a character vector of HPO codes of the direct parents,
# and is named with the HPO code. 
hpo_with_missing_terms$parents <- c(hpo_with_missing_terms$parents,
                                    missing_HPO_parents)

# list: each element is a character vector of HPO codes of the children,
# and is named with the HPO code
hpo_with_missing_terms$children <- c(hpo_with_missing_terms$children,
                                     missing_HPO_children)

# list: each element is a character vector of HPO codes of all ancestors,
# including the direct parents, and is named with the HPO code
hpo_with_missing_terms$ancestors <- c(hpo_with_missing_terms$ancestors, 
                                      missing_HPO_ancestors)

## add parent terms just in case
for(i in seq_along(missing_HPO_ancestors)) {
  missing_HPO_ancestors[[i]] <- unique(c(missing_HPO_ancestors[[i]],
                                         missing_HPO_parents[[i]]))
}

# list: each element is a FALSE/TRUE
# with names of HPO codes of obsolete terms
hpo_with_missing_terms$obsolete <- c(hpo_with_missing_terms$obsolete, 
                                      missing_HPO_obsolete)


# Update the hpo with the new terms
hpo2 <- hpo_with_missing_terms

# plot 1: plot all missing terms
jpeg("../result/SI/missing_terms_dag.jpeg",
     quality = 900,
     width = 2200,
     height = 1400)

# # specify the exact set of terms to appear in the plot
# terms_to_plot <- unique(c(missing_HPO_codes,
#                           unlist(missing_HPO_parents),
#                           unlist(missing_HPO_children)))
# terms_to_plot <- unique(c(missing_HPO_codes,
#                           unlist(missing_HPO_parents)))
terms_to_plot <- remove_links(hpo2, unique(c(missing_HPO_codes, 
                   get_ancestors(hpo2, missing_HPO_codes))))
# Character vector of colours 
colours_for_terms <- ifelse(terms_to_plot %in% missing_HPO_codes,
                            "#af8dc3","powderblue")



onto_plot(
  # ontology_index object
  ontology = hpo2,
  # specify the exact set of terms to appear in the plot
  terms = terms_to_plot,
  # Character vector of colours 
  fillcolor = colours_for_terms,
  # Numeric vector of term frequencies named by term IDs
  frequencies = get_term_frequencies(hpo2,
                                     missing_HPO_codes), 
  # List of character vectors of ontological term IDs
  term_sets = missing_HPO_codes,
  fontsize = 30)

dev.off()
