setwd("/home/rebecca/omicon/search")

library(ontologyIndex)

source("/home/rebecca/omicon/search/search_query_fxns.R")

## Ontologies are used for sample attribute searches:
# MONDO <- get_ontology("https://projects.nextwaretech.net/omicon/mondo/mondo.obo", extract_tags="everything")
# UBERON <- get_ontology("https://projects.nextwaretech.net/omicon/uberon/basic.obo", extract_tags="everything", 
#                        propagate_relationships=c("is_a", "part_of"), merge_equivalent_terms=F)

## Path to mapping tables for mapping identifiers:
tables_dir <- "/home/rebecca/omicon/mapping_tables/2022-11-02"

## Use latest data paths (for both normal and glioma samples):
data_dirs1 <- read.csv("/home/rebecca/omicon/search/data_directory_glioma_2023-06-01.csv") 
data_dirs2 <- read.csv("/home/rebecca/omicon/search/data_directory_normal_2023-06-01.csv") 
data_dirs <- rbind(data_dirs1, data_dirs2)
data_dirs[data_dirs == ""] <- NA

#############################################

# Result Type: Covariation Module
# Set ID: M642_v7.4
# Enrichment P-value: < 1e-10 (Module Definition: Bonferroni)
# Feature Input Type: Gene Symbol
# Feature Input ID: MKI67, BUB1, CKAP2L (Module Definition: Seed)
# Feature Match Count: > 1 (OR)

setid <- "M642_v7.4"
pval_cut <- 1e-10
feature_list <- c("MKI67", "BUB1", "CKAP2L")

covariation_queries <- data.frame(Result_Type=c("Covariation modules"),
                                  Module_Definition=c("Bonferroni"),
                                  Set_ID=c("M642_v7.4"),
                                  Enrichment_Pval=c(1e-10),
                                  Sample_Match_Percent=NA,
                                  Disease=NA,
                                  Feature_Input_Type="Gene Symbol",
                                  Feature_Input_ID=c("MKI67, BUB1, CKAP2L"),
                                  RE_Count=NA)

example_list <- lapply(1:nrow(data_dirs), function(i) {
  
  if (!is.na(data_dirs$FM_dir[i])) {
    
    networks <- list.files(path = data_dirs$FM_dir[i], pattern = "signum", full.names = T)
    networks <- networks[unlist(lapply(networks, function(x) length(list.files(path = x)) > 0))]
    
    DS_attr <- read.csv(list.files(path = data_dirs$SN_dir[i], pattern = "DS_attributes", full.names = T)[1])
    unique_id <- DS_attr$Value[DS_attr$Attribute == "Unique Identifier"]
    platform <- DS_attr$Value[DS_attr$Attribute == "Mapping Tables"]
    
    networks_list <- future_lapply(1:length(networks), function(j) {
      
      enrich_dirs <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
      enrich_out <- covariation_enrich_search(enrich_dirs, setid=setid, pval_cut=pval_cut)
      enrich_out <- enrich_out[enrich_out$Mod_Def=="TOPMODPOSBC",]
      
      if (nrow(enrich_out) > 0) {
        
        kME <- fread(list.files(path=networks[j], pattern="kME", full.names=T), data.table=F)
  
        if(unique_id=="SYMBOL"){
          
          kME <- mapAlias2Symbol(features=kME, unique_id_col=2, 
                                 tables_dir, keep_all=T, fill_NAs=T)
          
        } else {
          
          kME <- map2Any(features=kME, unique_id, map_to="SYMBOL", 
                         unique_id_col=2, platform, tables_dir, 
                         keep_all=T)
          
        }
        
        kME <- kME[,!is.element(colnames(kME), "SYMBOL.y")]
        kME[kME == "NA" | kME == ""] <- NA
        kME <- kME |>
          dplyr::filter(ModSeed %in% enrich_out$Module) %>%
          tidyr::separate_rows(SYMBOL, sep=" \\| ") |>
          as.data.frame()
        
        feature_mods <- covariation_feature_search(kME, feature_list, feature_type="SYMBOL", mod_def="Seed", and_or="OR")
        
        if(!is.null(feature_mods)){
          network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
          return(data.frame(Dataset=data_dirs$Title[i], 
                            Network=network, feature_mods))
        }
      }
      
    }, future.seed=T) 
    
    return(do.call(rbind, networks_list))
    
  } 
  
})
example1 <- do.call(rbind, example_list)

write.csv(example1, file=paste0("covariation_query_results_all_siginificant_", Sys.Date(), ".csv"), row.names=F)

covariation_queries$RE_Count[1] <- paste("All significant:", nrow(example1))

example1 <- example1 |>
  dplyr::group_by(Dataset) |>
  dplyr::slice_min(Pval, with_ties=T)

print("Covariation example 1 complete.")

write.csv(example1, file=paste0("covariation_query_results_most_significant_", Sys.Date(), ".csv"), row.names=F)

#############################################