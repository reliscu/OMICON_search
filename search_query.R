## Search examples spreadsheet: https://docs.google.com/spreadsheets/d/1kLycqoQAKMgyJcb8orQ6x4Oh39sYg22G/edit#gid=1882786248

source("/home/rebecca/omicon/search/search_query_fxns.R")

MONDO <- get_ontology("https://projects.nextwaretech.net/omicon/mondo/mondo.obo", extract_tags="everything")
UBERON <- get_ontology("https://projects.nextwaretech.net/omicon/uberon/basic.obo", extract_tags="everything", 
                       propagate_relationships=c("is_a", "part_of"), merge_equivalent_terms=F)

tables_dir <- "/home/rebecca/omicon/mapping_tables/2022-11-02"

data_dirs1 <- read.csv("/home/rebecca/omicon/advanced_search/data_directory_glioma_2023-02-21.csv") %>% na_if("")
data_dirs2 <- read.csv("/home/rebecca/omicon/advanced_search/data_directory_normal_2023-02-21.csv") %>% na_if("")
data_dirs <- rbind(data_dirs1, data_dirs2)

# nonbrain <- c("cranial nerve II", "C1 segment of cervical spinal cord")

setwd("/home/rebecca/omicon/search/search_results/05-15-23")

search_queries_table1(data_dirs)

search_queries_table2(data_dirs)