## Search examples spreadsheet: https://docs.google.com/spreadsheets/d/1kLycqoQAKMgyJcb8orQ6x4Oh39sYg22G/edit#gid=1882786248

setwd("/home/rebecca/omicon/search/search_results/09-12-23")

library(ontologyIndex)

source("/home/rebecca/omicon/search/search_query_fxns.R")

## Ontologies are used for sample attribute searches:
MONDO <- get_ontology("https://projects.nextwaretech.net/omicon/mondo/mondo.obo", extract_tags="everything")
UBERON <- get_ontology("https://projects.nextwaretech.net/omicon/uberon/basic.obo", extract_tags="everything", 
                       propagate_relationships=c("is_a", "part_of"), merge_equivalent_terms=F)

## Path to mapping tables for mapping identifiers:
tables_dir <- "/home/rebecca/omicon/mapping_tables/2022-11-02"

## Use latest data paths (for both normal and glioma samples):
data_dirs1 <- read.csv("/home/rebecca/omicon/search/data_directory_glioma_2023-06-01.csv") 
data_dirs2 <- read.csv("/home/rebecca/omicon/search/data_directory_normal_2023-06-01.csv") 
data_dirs <- rbind(data_dirs1, data_dirs2)
data_dirs[data_dirs == ""] <- NA

## Run searches:
## *Note: working directory is where results will be saved
# search_queries_sheet1(data_dirs)
# search_queries_sheet2(data_dirs)
# search_queries_covariation(data_dirs)
