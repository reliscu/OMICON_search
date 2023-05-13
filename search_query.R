setwd("/home/rebecca/omicon/search")

source("/home/rebecca/omicon/search/search_query_fxns.R")

## Search examples: https://docs.google.com/spreadsheets/d/1kLycqoQAKMgyJcb8orQ6x4Oh39sYg22G/edit#gid=1882786248

MONDO <- get_ontology("https://projects.nextwaretech.net/omicon/mondo/mondo.obo", extract_tags="everything")
UBERON <- get_ontology("https://projects.nextwaretech.net/omicon/uberon/basic.obo", extract_tags="everything", 
                       propagate_relationships=c("is_a", "part_of"), merge_equivalent_terms=F)

tables_dir <- "/home/rebecca/omicon/mapping_tables/2022-11-02"

data_dirs <- read.csv("/home/rebecca/omicon/advanced_search/data_directory_glioma_2023-02-21.csv") %>% na_if("")

# nonbrain <- c("cranial nerve II", "C1 segment of cervical spinal cord")
