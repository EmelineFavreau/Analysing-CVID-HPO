######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################


# missing terms are designed by Luiza
# the first data frame describe the new terms (term and code) and their parents
mhtp <- 
  read_excel("../input/Missing_HPO_Terms_Proposal.LC.cleaned.xlsx", 
             sheet = "parents")


############# Suppl Table New terms DAG ####################################

# vector of the terms 
missing_HPO_terms <- unique(mhtp$ProposedTerm[is.na(mhtp$ProposedCode)])

# create temporary missing codes
nums <- seq(from = 9999901,
            length.out = length(missing_HPO_terms))

missing_HPO_codes <- sprintf("HP:%07d", nums)

# assign them randomly
names(missing_HPO_codes) <- missing_HPO_terms

# input them in the df
for(i in 1:length(missing_HPO_codes)){
  tt <- missing_HPO_codes[i]
  mhtp$ProposedCode[mhtp$ProposedTerm == names(tt)] <- tt
}

# update parent code if relationship is known
u <- mhtp$ProposedParentTerm[mhtp$ProposedParentCode == "NA"]
for(i in 1:length(u)){
  tc <- missing_HPO_codes[match(u[i], names(missing_HPO_codes))]
  mhtp$ProposedParentCode[mhtp$ProposedParentTerm == names(tc)] <- tc
}


# each element is a character vector of HPO codes of the direct parents,
# and is named with the HPO code
missing_HPO_parents <- lapply(unname(missing_HPO_codes),
   function(x) mhtp$ProposedParentCode[grep(x, mhtp$ProposedCode)])
names(missing_HPO_parents) <- missing_HPO_codes



# add the following:
# "Tuberculous mycobacterial infection"(#7) is parent to
# "Disseminated tuberculous mycobacterial infection" (#7)
# "Viral meningitis" is parent to "Herpes (simplex) meningitis"
# "Viral meningitis" is parent to  "Recurrent herpes meningitis"
# "Unusual bacterial infection" is parent to "Bacterial meningitis"
# "Non-infectious encephalitis" is parent to "Autoimmune encephalitis"
# vector of parents 
pcv <- c("Tuberculous mycobacterial infection",
         "Viral meningitis",
         "Viral meningitis",
         "Unusual bacterial infection",
         "Non-infectious encephalitis")
# add the vector names, representing the children
# e.g. "Disseminated tuberculous mycobacterial infection" is child of 
# "Tuberculous mycobacterial infection"
names(pcv) <- c("Disseminated tuberculous mycobacterial infection",
                "Herpes (simplex) meningitis",
                "Recurrent herpes meningitis",
                "Bacterial meningitis",
                "Autoimmune encephalitis")

# update relationships
for(i in 1:length(pcv)){
  oc <- grep(pcv[i], names(missing_HPO_codes), fixed = TRUE) # the order of the parent
  occ <- grep(names(pcv[i]), names(missing_HPO_codes), fixed = TRUE) # the order of the child
  missing_HPO_parents[[occ]] <- c(missing_HPO_parents[[occ]], missing_HPO_parents[[oc]]) 
}

# each element is a character vector of HPO codes of the children,
# and is named with the HPO code
missing_HPO_children <- lapply(unname(missing_HPO_codes),
                               function(x) character())
names(missing_HPO_children) <- missing_HPO_codes

for(i in 1:length(pcv)){
  # find the code for the parent (e.g. "HP:999990")
  pc <- missing_HPO_codes[names(missing_HPO_codes) == pcv[i]]
  # find the code for the child (e.g. "HP:9999908" )
  cc <- missing_HPO_codes[names(missing_HPO_codes) == names(pcv)[i]]
  # add the code of the child, when the element is the parent
  missing_HPO_children[names(missing_HPO_children) == pc] <- cc
}





# list: each element is a character vector of HPO codes of all ancestors,
# including the direct parents, and is named with the HPO code
missing_HPO_ancestors <- lapply(missing_HPO_parents,
                                function(x) get_ancestors(hpo, x))


# update relationships
for(i in 1:length(pcv)){
  oc <- grep(pcv[i], names(missing_HPO_codes), fixed = TRUE) # the order of the parent
  occ <- grep(names(pcv[i]), names(missing_HPO_codes), fixed = TRUE) # the order of the child
  missing_HPO_ancestors[[occ]] <- unique(c(missing_HPO_ancestors[[occ]],
                                  unname(missing_HPO_codes[oc])))
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



