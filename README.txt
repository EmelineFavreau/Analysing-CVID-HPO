Scripts used for HPO analyses:

# Set up
common.R: list all dependencies

# Clean data by inputting lab values as HPO, excluding records that do not match analysis design

# For primary analyses:
01-epidemiology_table.R: count patients for age, sex, CVID presentations
02-analysis_training_data.R: HPO similarity between training sets
03-number_patients_per_category.R: count patients within subsets (B cell, immunoglobulin, genotype)
04-clustering-patients.R: using key HPO terms, cluster patients as infectionBronchiectasis or complex
05.1-matrix_design.R: create matrix with presence/absence of phenotype and genotype for all patients
05.2-cluster_demographics.R: count infection/complex patients within subsets (B cell, immunoglobulin, genotype)
05.3-Euroclass_demographics.R: count patients following Euroclass subsets
05.4-Bcell_HPO_associations.R: Cochran-Mantel-Haenszel tests for SmB/CD21/Tr level and each phenotype (HPO, GLILD, key HPO groups)
05.5-immuno_HPO_associations.R: Cochran-Mantel-Haenszel tests for Immunoglobulin groups and each phenotype (HPO, GLILD, key HPO groups)
05.6-genetic_HPO_associations.R: Cochran-Mantel-Haenszel tests for genetic variants groups and each phenotype (HPO, GLILD, key HPO groups)
05.7-clusters_associations.R: Cochran-Mantel-Haenszel tests for infection/complex and each biomarker


# For manuscript figures and tables:
Figure1.R
Figure2.R
Figure3.R
Figure4.R
Figure5.R
Figure6_prep.R
Figure6_generate-embedding.txt
create_embedding.py
Figure6.R
SuppFigure-Missing-terms.R
SuppTable-complete-HPO-list.R
Table2.R
Table3.R
