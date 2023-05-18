library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(future.apply)
library(ontologyIndex)

options(future.globals.maxSize=Inf)
plan(multicore, workers=10)

source("/home/rebecca/code/misc/map_identifiers/map_identifiers_function.R")

search_queries_sheet2 <- function(data_dirs, MONDO, UBERON){
  
  ############################################# Sheet 2, example 1 ############################################# 
  
  print("Starting sheet 2 example 1.")
  
  ## Find all human brain samples with any kind of disease (disease DOES NOT EQUAL Normal)
  
  DC_dirs <- unique(sapply(strsplit(na.omit(unique(data_dirs$DC_dir)), "/"), function(x){
    paste(x[-length(x)], collapse="/")
  }))
  
  example1_list <- lapply(1:length(DC_dirs), function(i){
    
    DC_attr <- read.csv(list.files(path=DC_dirs[i], pattern="DC_attributes", full.names=T))
    sampleinfo <- fread(list.files(path=DC_dirs[i], pattern="sample_attributes", 
                                   full.names=T), data.table=F)
    sampleinfo <- sampleinfo[get_brain_samples(uberon_vec=sampleinfo$UBERON_ID),]
    sampleinfo |> na_if("") |>
      dplyr::filter(sampleinfo$Disease!="Normal") |>
      dplyr::mutate(Data_Collection=DC_attr[1,2]) |>
      dplyr::select(Data_Collection, Label, 
                    Organism, Tissue, Disease)
    
  })
  example1 <- do.call(rbind, example1_list)
  
  print("Sheet 2 example 1 complete.")
  
  write.csv(example1, file=paste0("sheet2_example1_search_results_", Sys.Date(), ".csv"), row.names=F)

  ############################################# Sheet 2, example 2 ############################################# 
  
  print("Starting sheet 2 example 2.")
  
  ## Find all covariation networks from human samples without disease where minimum module size >= 10 and # modules ≥ 50
  
  example2_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(sum(is.element(sampleinfo$Disease, "Normal"))>0){
        
        networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
        networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
        minsize <- as.numeric(gsub("minSize", "", sapply(strsplit(networks, "_"), function(x){
          x[grep("minSize", x)]
        })))
        networks <- networks[minsize>=10]
        
        network_list <- lapply(1:length(networks), function(j){
          
          modstats <- fread(list.files(path=networks[j], pattern="Module_statistics", full.names=T)[1], data.table=F)
          
          if(nrow(modstats)>=50){
            DS_attr <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="attributes", full.names=T))
            network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
            return(data.frame(Dataset=data_dirs$Title[i], Network=network, 
                              No.Modules=nrow(modstats)))
          } 
          
        }) 
        
        return(do.call(rbind, network_list))
        
      } 
      
    } 
    
  })
  example2 <- do.call(rbind, example2_list)
  
  print("Sheet 2 example 2 complete.")
  
  write.csv(example2, file=paste0("sheet2_example2_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 2, example 3 ############################################# 
  
  print("Starting sheet 2 example 3.")
  
  ## Find all covariation modules in datasets consisting SOLELY of normal adult human samples that contain 
  ## the following genes: BUB1, MKI67, PBK, and WEE1 (ranked by Jaccard index)
  
  feature_list <- c("BUB1", "MKI67", "PBK", "WEE1")
  
  example3_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      brain_samples <- get_brain_samples(uberon_vec=sampleinfo$UBERON_ID)
      
      if(length(brain_samples)==nrow(sampleinfo)){
       
        if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
          
          if(length(intersect(grep("Normal", sampleinfo$Disease), 
                              grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
            
            networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
            networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
            
            DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
            unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
            platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
            
            networks_list <- future_lapply(1:length(networks), function(j){
              
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
              kME <- kME |> na_if("NA") |> na_if("") |>
                tidyr::separate_rows(SYMBOL, sep=" \\| ") |>
                as.data.frame()
              
              feature_mods <- covariation_feature_search(kME, feature_list, feature_type="SYMBOL", mod_def="Any", and_or="AND")
              
              if(!is.null(feature_mods)){
                network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
                return(data.frame(Dataset=data_dirs$Title[i], 
                                  Network=network, feature_mods))
              }
              
            }, future.seed=T) 
            
            return(do.call(rbind, networks_list))
            
          }
          
        }
        
      }
      
    } 
    
  })
  example3 <- do.call(rbind, example3_list)
  
  if(!is.null(example3)>0){
    
    print("Sheet 2 example 3 complete.")
    write.csv(example3, file=paste0("sheet2_example3_search_results_", Sys.Date(), ".csv"), row.names=F)
    
  } else {
    print("No results found for example query 3.")
  }
 
  ############################################# Sheet 2, example 4 ############################################# 
  
  print("Starting sheet 2 example 4.")
  
  ## Find all covariation modules in datasets consisting SOLELY of human normal brain samples 
  ## that are significantly enriched with microglial markers, ranked by enrichment P-values
    
  example4_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      brain_samples <- get_brain_samples(uberon_vec=sampleinfo$UBERON_ID)

      if(length(intersect(brain_samples, grep("Normal", sampleinfo$Disease, ignore.case=T)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
            enrich_out <- covariation_enrich_search(enrich, setname="microglia")
            network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
            return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out))
            
          }, future.seed=T) 
          
          return(do.call(rbind, networks_list))
            
        }
        
    } 
    
  })
  example4 <- do.call(rbind, example4_list)
  
  example4$Mod_ID <- paste(example4$Dataset, example4$Network, 
                           example4$Module, example4$Mod_Def)
  
  example4 <- example4 |> 
    dplyr::group_by(Mod_ID) |>
    dplyr::slice_min(Pval, with_ties=T) |> 
    dplyr::arrange(Pval)
  
  write.csv(example4, file=paste0("sheet2_example4_search_results_gene_sets_", Sys.Date(), ".csv"), row.names=F)
  
  example4 <- example4 |> 
    dplyr::group_by(Mod_ID) |>
    dplyr::slice_min(Pval, with_ties=F)

  print("Sheet 2 example 4 complete.")
  
  write.csv(example4, file=paste0("sheet2_example4_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 2, example 5 ############################################# 
  
  print("Starting sheet 2 example 5.")
  
  ## Export a list of unique gene symbols for covariation modules identified in datasets that include normal
  ## adult male brain samples and are maximally enriched with markers of radial glia (cell type), ranked by kME / Fidelity
  
  example5_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
    
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
        
        brain_samples <- get_brain_samples(uberon_vec=sampleinfo$UBERON_ID)
        
        if(length(intersect(brain_samples, 
                            intersect(grep("M", sampleinfo$Sex), 
                                      intersect(grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T), 
                                                grep("Normal", sampleinfo$Disease)))))>0){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          networks_list <- future_lapply(1:length(networks), function(j){

            enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
            enrich_out <- covariation_enrich_search(enrich, setname="radial_glia")
            if(nrow(enrich_out)>0){
              DS_attr <- list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1]
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out, 
                                DS_Attr=DS_attr))
            }
            
          }, future.seed=T)
          modules <- do.call(rbind, networks_list)
          
          if(nrow(modules)>0){
            
            modules <- modules |> dplyr::slice_min(Pval, with_ties=T)
            
            ## Get module genes
            
            gene_list <- lapply(1:nrow(modules), function(j){
              
              kME <- fread(list.files(path=file.path(data_dirs$FM_dir[i], modules$Network[j]), 
                                      pattern="kME", full.names=T), data.table=F)
              col <- kME[,grep(toupper(modules$Mod_Def[j]), toupper(colnames(kME)))]
              kME <- kME[is.element(col, as.character(modules$Module[j])),]
              
              ## Map identifiers to SYMBOL:
              
              DS_attr <- read.csv(modules$DS_Attr[j])
              unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
              platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
              
              if(unique_id=="SYMBOL"){
                
                kME <- mapAlias2Symbol(features=kME, unique_id_col=2, 
                                       tables_dir, keep_all=T, fill_NAs=T)
                
              } else {
                
                kME <- map2Any(features=kME, unique_id, map_to="SYMBOL", 
                               unique_id_col=2, platform, tables_dir, 
                               keep_all=T)
                
              }
              
              kME <- kME[,!is.element(colnames(kME), "SYMBOL.y")]
              
              kME <- kME |> 
                tidyr::separate_rows(SYMBOL, sep=" \\| ") 
              
              return(data.frame(SYMBOL=unique(kME$SYMBOL),
                                Dataset=modules$Dataset[j], Network=modules$Network[j],
                                SetName=modules$SetName[j], Module=modules$Module[j],
                                Mod_Def=modules$Mod_Def))
              
            })
            
            return(do.call(rbind, gene_list))
            
          } 
          
        } 
        
      } 
      
    }
    
  })
  example5 <- do.call(rbind, example5_list)
  
  example5 <- example5 |>
    na_if("NA") |>
    dplyr::filter(!is.na(SYMBOL))
  
  write.csv(example5, file=paste0("sheet2_example5_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example5 <- example5 |> 
    dplyr::group_by(SYMBOL) |> 
    dplyr::summarise(No.Datasets=n()) |>
    dplyr::arrange(desc(No.Datasets))
  
  print("Sheet 2 example 5 complete.")
  
  write.csv(example5, file=paste0("sheet2_example5_search_results_", Sys.Date(), ".csv"), row.names=F)
  
}

search_queries_sheet1 <- function(data_dirs, MONDO, UBERON){
  
  ############################################# Sheet 1, example 1 ############################################# 
  
  print("Starting sheet 1 example 1.")
  
  ## Find all RNA-seq datasets with at least 100 samples from adult human gliomas that contain the gene MALAT1

  DS_dirs <- na.omit(c(unique(data_dirs$SN_dir), unique(data_dirs$DC_dir)))
  
  example1_list <- lapply(1:length(DS_dirs), function(i){
    
    DS_files <- list.files(path=DS_dirs[i], full.names=T)
    DS_attr <- read.csv(DS_files[grep("DS[0-9]*_attributes", DS_files)][1])
    
    if(DS_attr$Value[DS_attr$Attribute=="Technology"]=="Sequencer"){
      
      unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
      platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
      
      DS_list <- c()
      
      if(sum(grepl("outliers_removed", DS_files))>0){
        
        OR_sampleinfo <- read.csv(DS_files[grep("outliers_removed_sample", DS_files)])
        
        if(is.element("Developmental_Epoch", colnames(OR_sampleinfo))){
          
          glioma_samples <- get_glioma_samples(mondo_vec=OR_sampleinfo$MONDO_ID)
          
          OR_sampleinfo <- OR_sampleinfo[intersect(glioma_samples, grep("adult", OR_sampleinfo$Developmental_Epoch, ignore.case=T)),]
          
          if(nrow(OR_sampleinfo)>=100){
            
            OR_feature <- fread(DS_files[grep("outliers_removed_feature", DS_files)], data.table=F)
            
            if(unique_id=="SYMBOL"){
              
              OR_feature <- mapAlias2Symbol(features=OR_feature, unique_id_col=2, 
                                            tables_dir, keep_all=T, fill_NAs=T)
              
            } else {
              
              OR_feature <- map2Any(features=OR_feature, unique_id, map_to="SYMBOL", 
                                    unique_id_col=2, platform, tables_dir, keep_all=T)

            }
            
            OR_feature <- OR_feature |> 
              tidyr::separate_rows(SYMBOL, sep=" \\| ") 
            
            if(length(grep("MALAT1", OR_feature$SYMBOL))>0){
              
              OR_DS_attr <- read.csv(DS_files[grep("outliers_removed_DS_attributes", DS_files)])
              dataset <- paste(OR_DS_attr$Value[OR_DS_attr$Attribute=="Title"],
                               OR_DS_attr$Value[OR_DS_attr$Attribute=="Version"])
              DS_list <- dataset
              
            }
            
            if(length(DS_files[grep("Qnorm", DS_files)])>0){
              
              QN_feature <- fread(DS_files[grep("Qnorm_.*feature", DS_files)], data.table=F)
              
              if(unique_id=="SYMBOL"){
                
                QN_feature <- mapAlias2Symbol(features=QN_feature, unique_id_col=2, 
                                              tables_dir, keep_all=T, fill_NAs=T)
                
              } else {
                
                QN_feature <- map2Any(features=QN_feature, unique_id, map_to="SYMBOL", 
                                      unique_id_col=2, platform, tables_dir, keep_all=T)
                
              }
              
              QN_feature <- QN_feature |> 
                tidyr::separate_rows(SYMBOL, sep=" \\| ") 
              
              if(length(grep("MALAT1", QN_feature$SYMBOL))>0){
                
                QN_DS_attr <- read.csv(DS_files[grep("Qnorm_.*DS_attributes", DS_files)])
                dataset <- paste(QN_DS_attr$Value[QN_DS_attr$Attribute=="Title"],
                                 QN_DS_attr$Value[QN_DS_attr$Attribute=="Version"])
                DS_list <- c(DS_list, dataset)
                
              }
              
            }
            
            if(length(DS_files[grep("ComBat", DS_files)])>0){
              
              BC_feature <- fread(DS_files[grep("ComBat_feature", DS_files)], data.table=F)
              
              if(unique_id=="SYMBOL"){
                
                BC_feature <- mapAlias2Symbol(features=BC_feature, unique_id_col=2, 
                                              tables_dir, keep_all=T, fill_NAs=T)
                
              } else {
                
                BC_feature <- map2Any(features=BC_feature, unique_id, map_to="SYMBOL", 
                                      unique_id_col=2, platform, tables_dir, keep_all=T)
                
              }
              
              BC_feature <- BC_feature |> 
                tidyr::separate_rows(SYMBOL, sep=" \\| ") 
              
              if(length(grep("MALAT1", BC_feature$SYMBOL))>0){
                
                BC_DS_attr <- read.csv(DS_files[grep("ComBat_DS_attributes", DS_files)])
                dataset <- paste(BC_DS_attr$Value[BC_DS_attr$Attribute=="Title"],
                                 BC_DS_attr$Value[BC_DS_attr$Attribute=="Version"])
                DS_list <- c(DS_list, dataset)
                
              }
              
            }
            
          }
          
        } 
          
      } else {
        
        sampleinfo <- read.csv(DS_files[grep("sample_attributes", DS_files)])
        
        if(is.element("Developmental_Epoch", colnames(sampleinfo))){
          
          glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
          
          sampleinfo <- sampleinfo[intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)),]
          
          if(nrow(sampleinfo)>=100){
            
            DS_features <- fread(DS_files[grep("feature_metadata", DS_files)], data.table=F)
            
            if(unique_id=="SYMBOL"){
              
              DS_features <- mapAlias2Symbol(features=DS_features, unique_id_col=2, 
                                             tables_dir, keep_all=T, fill_NAs=T)
              
            } else {
              
              DS_features <- map2Any(features=DS_features, unique_id, map_to="SYMBOL", 
                                     unique_id_col=2, platform, tables_dir, keep_all=T)
              
            }
            
            DS_features <- DS_features |> 
              tidyr::separate_rows(SYMBOL, sep=" \\| ") 
            
            if(length(grep("MALAT1", DS_features$SYMBOL))>0){
              
              dataset <- paste(DS_attr$Value[DS_attr$Attribute=="Title"],
                               DS_attr$Value[DS_attr$Attribute=="Version"])
              DS_list <- dataset
              
            }
            
          }
          
        }
        
      }

      return(DS_list)
      
    }
    
  })
  example1 <- data.frame(Dataset=unique(do.call(c, example1_list)))
  
  print("Sheet 1 example 1 complete.")
  
  write.csv(example1, file=paste0("sheet1_example1_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 2 ############################################# 
  
  print("Starting sheet 1 example 2.")
  
  ## Find all data collections that include paired RNA-seq and microarray datasets from the same human 
  ## brain samples, with sample size ≥ 100
  
  DC_dirs <- unique(sapply(strsplit(na.omit(unique(data_dirs$DC_dir)), "/"), function(x){
    paste(x[-length(x)], collapse="/")
  }))
  
  example2_list <- lapply(1:length(DC_dirs), function(i){
    
    sampleinfo <- fread(list.files(path=DC_dirs[i], pattern="sample_attributes", 
                                   full.names=T), data.table=F)
    
    ## Restrict to data collections with at least 100 samples:
    
    if(nrow(sampleinfo>=100)){
      
      DS_dirs <- list.dirs(path=DC_dirs[i], recursive=F)
      
      techs <- sapply(1:length(DS_dirs), function(j){
        DS_attr <- fread(list.files(path=DS_dirs[j], pattern="DS[0-9]+_attributes", 
                                    full.names=T), data.table=F)
        return(DS_attr$Value[DS_attr$Attribute=="Technology"])
      })
      
      ## Restrict to data collections that utilized both technologies:
      
      if(sum(is.element(unique(techs), c("Microarray", "Sequencer")))==2){
        
        DS_samples <- lapply(1:length(DS_dirs), function(j){
          DS_sampleinfo <- fread(list.files(path=DS_dirs[j], pattern="sample_attributes", full.names=T), data.table=F)
          return(data.frame(Label=DS_sampleinfo$Label, Technology=techs[j],
                            Tissue=DS_sampleinfo$Tissue, UBERON_ID=DS_sampleinfo$UBERON_ID))
        })
        
        samples <- do.call(rbind, DS_samples)
        
        ## Restrict to brain samples:
        
        samples <- samples[get_brain_samples(uberon_vec=samples$UBERON_ID),]
        
        ## Identify samples that were profiled with both technologies:
        
        samples <- samples |>
          dplyr::group_by(Label) |>
          dplyr::summarise(
            Technology=paste(sort(unique(Technology)), collapse=", ")
          )
        
        nsamples <- sum(grepl("Microarray, Sequencer", samples$Technology))
        
        if(nsamples>0){
          
          DC_attr <- read.csv(list.files(path=DC_dirs[i], pattern="DC_attributes", full.names=T))
          dataset <- paste(DC_attr$Value[DC_attr$Attribute=="Title"])
          return(data.frame(Data_Collection=dataset, No.Samples=nsamples))
          
        }
        
      }
      
    }
    
  })
  example2 <- do.call(rbind, example2_list)
  
  print("Sheet 1 example 2 complete.")
  
  write.csv(example2, file=paste0("sheet1_example2_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 4 ############################################# 

  print("Starting sheet 1 example 4.")
  
  ## Find all adult human malignant glioma samples with deleterious mutations in IDH1

  DC_dirs <- unique(sapply(strsplit(na.omit(unique(data_dirs$DC_dir)), "/"), function(x){
    paste(x[-length(x)], collapse="/")
  }))
  
  example4_list <- lapply(1:length(DC_dirs), function(i){
    
    DC_attr <- read.csv(list.files(path=DC_dirs[i], pattern="DC_attributes", full.names=T))
    sampleinfo <- fread(list.files(path=DC_dirs[i], pattern="sample_attributes", 
                                   full.names=T), data.table=F)
    
    if(is.element("Developmental_Epoch", colnames(sampleinfo)) & is.element("Tumor_Grade", colnames(sampleinfo))){
      
      glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
      
      sampleinfo <- sampleinfo[intersect(intersect(grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T), 
                                                   which(sampleinfo$Tumor_Grade>=2)), glioma_samples),]
      
      if(nrow(sampleinfo)>0 & sum(grepl("IDH1", colnames(sampleinfo)))>0){

        sampleinfo <- sampleinfo[,!grepl("Source", colnames(sampleinfo))]
        cols <- grep("IDH1", colnames(sampleinfo))
        mut_samples <- apply(sampleinfo[cols], 1, function(x){
          sum(is.element(x, "NCIT:C172343"))>0
        })
        
        if(sum(mut_samples)>0){
          
          sampleinfo <- sampleinfo[mut_samples,]
          sampleinfo$IDH1_Mut <- apply(sampleinfo[cols], 1, function(x){
            x <- x[x!="WT"]
            x[x==""] <- NA
            paste(na.omit(unique(x)), collapse=", ")
          })
          sampleinfo |> 
            dplyr::mutate(Data_Collection=DC_attr[1,2]) |>
            dplyr::select(Data_Collection, Label, Organism,
                          Disease, Developmental_Epoch, 
                          Tumor_Grade, IDH1_Mut)
          
        }
        
      }
      
    }
    
  })
  example4 <- do.call(rbind, example4_list)
  
  print("Sheet 1 example 4 complete.")
  
  write.csv(example4, file=paste0("sheet1_example4_search_results_", Sys.Date(), ".csv"), row.names=F)

  ############################################# Sheet 1, example 5 ############################################# 

  print("Starting sheet 1 example 5.")
  
  ## Find all adult human oligodendroglioma samples with deletion of chr1p OR chr19q 
  ## AND no mutation in IDH1

  DC_dirs <- unique(sapply(strsplit(na.omit(unique(data_dirs$DC_dir)), "/"), function(x){
    paste(x[-length(x)], collapse="/")
  }))
  
  example5_list <- lapply(1:length(DC_dirs), function(i){
    
    DC_attr <- read.csv(list.files(path=DC_dirs[i], pattern="DC_attributes", full.names=T))
    sampleinfo <- fread(list.files(path=DC_dirs[i], pattern="sample_attributes", 
                                   full.names=T), data.table=F)
    
    if(is.element("Developmental_Epoch", colnames(sampleinfo))){
      
      sampleinfo <- sampleinfo[intersect(grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T),
                                         grep("oligodendroglioma", sampleinfo$Disease)),]
      
      if(nrow(sampleinfo)>0 & 
         sum(grepl("IDH1", colnames(sampleinfo)))>0 & 
         sum(grepl("chr1p|chr19q", colnames(sampleinfo)))>0){
        
        sampleinfo <- sampleinfo[,!grepl("Source", colnames(sampleinfo))]
        
        ## Get all IDH1 WT samples:
        
        IDH1_cols <- grep("IDH1", colnames(sampleinfo))
        WT_samples <- which(apply(sampleinfo[IDH1_cols], 1, function(x){
          sum(is.element(x, "WT"))==length(IDH1_cols)
        }))
        
        ## Get all samples with chr1p and/or chr19q deletion:
        
        del_cols <- grep("chr1p|chr19q", colnames(sampleinfo))
        del_samples <- which(apply(sampleinfo[del_cols], 1, function(x){
          sum(is.element(x, "SO:0001743"))>0
        }))

        mut_samples <- intersect(WT_samples, del_samples)
        
        if(length(mut_samples)>0){
          
          sampleinfo <- sampleinfo[mut_samples,]
          
          sampleinfo$IDH1_Mut <- apply(sampleinfo[IDH1_cols], 1, function(x){
            x[x==""] <- NA
            paste(na.omit(unique(x)), collapse=", ")
          })
          
          sampleinfo$Deletion <- apply(sampleinfo[del_cols], 1, function(x){
            x <- colnames(sampleinfo)[del_cols][x=="SO:0001743"]
            paste(gsub("_Mut.*", "", na.omit(unique(x))), collapse=", ")
          })
          
          sampleinfo |> 
            dplyr::mutate(Data_Collection=DC_attr[1,2]) |>
            dplyr::select(Data_Collection, Label, Organism,
                          Disease, Developmental_Epoch, 
                          IDH1_Mut, Deletion)
          
        }
      
      }

    }
    
  })
  example5 <- do.call(rbind, example5_list)
  
  print("Sheet 1 example 5 complete.")
  
  write.csv(example5, file=paste0("sheet1_example5_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 6 ############################################# 
  
  print("Starting sheet 1 example 6.")
  
  ## Find all Analyses involving FM by Rebecca that have been created since 4/1/21
  
  example6_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      proj_attr <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="project_attributes", full.names=T))
      owner <- proj_attr$Value[proj_attr$Attribute=="Owner"]
      if(owner=="Rebecca Eliscu"){
        analysis <- proj_attr$Value[proj_attr$Attribute=="Title"]
        return(data.frame(Analysis=analysis))
      }
    }
    
  })
  example6 <- do.call(rbind, example6_list)
  
  print("Sheet 1 example 6 complete.")
  
  write.csv(example6, file=paste0("sheet1_example6_search_results_", Sys.Date(), ".csv"), row.names=F)

  ############################################# Sheet 1, example 7 ############################################# 
  
  print("Starting sheet 1 example 7.")
  
  ## Find all covariation networks from human gliomas on Affymetrix U133A microarrays where 
  ## min module size >= 10 and # modules >= 50
  
  example7_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      platform <- DS_attr$Value[DS_attr$Attribute=="Platform"]
      
      glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
      
      if(length(glioma_samples)>0 & platform=="Affymetrix U133A"){
        
        networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
        networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
        minsize <- as.numeric(gsub("minSize", "", sapply(strsplit(networks, "_"), function(x){
          x[grep("minSize", x)]
        })))
        networks <- networks[minsize>=10]
        
        network_list <- lapply(1:length(networks), function(j){
          
          modstats <- fread(list.files(path=networks[j], pattern="Module_statistics", full.names=T)[1], data.table=F)
          
          if(nrow(modstats)>=50){
            DS_attr <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="attributes", full.names=T))
            network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
            return(data.frame(Dataset=data_dirs$Title[i], Network=network, 
                              No.Modules=nrow(modstats)))
          } 
          
        }) 
        
        return(do.call(rbind, network_list))
        
      } 
      
    } 
    
  })
  example7 <- do.call(rbind, example7_list)
  
  print("Sheet 1 example 7 complete.")
  
  write.csv(example7, file=paste0("sheet1_example7_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 8 ############################################# 

  print("Starting sheet 1 example 8.")
  
  ## Find all covariation modules in datasets consisting SOLELY of adult human glioma samples that contain 
  ## the following genes: BUB1, MKI67, PBK, and WEE1

  feature_list <- c("BUB1", "MKI67", "PBK", "WEE1")

  example8_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
        
        glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
        
        if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
          unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
          platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
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
            kME <- kME |> na_if("NA") |> na_if("") |>
              tidyr::separate_rows(SYMBOL, sep=" \\| ") |>
              as.data.frame()
            
            feature_mods <- covariation_feature_search(kME, feature_list, feature_type="SYMBOL", mod_def="Any", and_or="OR")
            
            if(!is.null(feature_mods)){
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], Network=network, feature_mods))
            }
            
          }, future.seed=T) 
          
          return(do.call(rbind, networks_list))
          
        }
        
      }
      
    } 
    
  })
  example8 <- do.call(rbind, example8_list)
  
  print("Sheet 1 example 8 complete.")
  
  write.csv(example8, file=paste0("sheet1_example8_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 9 ############################################# 
  
  print("Starting sheet 1 example 9.")
  
  ## Find all covariation modules in datasets consisting SOLELY of adult human oligodendroglioma samples that
  ## are significantly enriched with microglial markers, ranked by enrichment P-values
  
  example9_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){
        
        if(length(intersect(grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T),
                            grep("oligodendroglioma", sampleinfo$Disease)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
            enrich_out <- covariation_enrich_search(enrich, setname="microglia")
            network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
            return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out))
            
          }, future.seed=T) 
          
          return(do.call(rbind, networks_list))
          
        }
        
      }

    } 
    
  })
  example9 <- do.call(rbind, example9_list)
  
  example9$Mod_ID <- paste(example9$Dataset, example9$Network,
                           example9$Module, example9$Mod_Def)

  example9 <- example9 |>
    dplyr::group_by(Mod_ID) |>
    dplyr::slice_min(Pval, with_ties=T) |>
    dplyr::arrange(Pval)
  
  write.csv(example9, file=paste0("sheet1_example9_search_results_gene_sets_", Sys.Date(), ".csv"), row.names=F)
  
  example9 <- example9 |>
    dplyr::group_by(Mod_ID) |>
    dplyr::slice_min(Pval, with_ties=F)

  print("Sheet 1 example 9 complete.")
  
  write.csv(example9, file=paste0("sheet1_example9_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 10 ############################################# 
  
  print("Starting sheet 1 example 10.")
  
  ## Find all Analyte Data files generated EXCLUSIVELY from adult human glioma samples by RNA-seq that 
  ## were used as input for covariation analysis by FindModules
  
  example10_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      
      if(DS_attr$Value[DS_attr$Attribute=="Technology"]=="Sequencer"){
        
        SN_attr_files <- list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)
        proj_attr <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="project_attributes", full.names=T))
        input_dataset <- data_dirs$Title[i]

        ## Get sample attributes for dataset used as input for FindModules:
        
        for(j in 1:length(SN_attr_files)){
          DS_attr <- read.csv(SN_attr_files[j])
          dataset <- paste(DS_attr$Value[DS_attr$Attribute=="Title"],
                           DS_attr$Value[DS_attr$Attribute=="Version"])
          if(dataset==input_dataset) break
        }
        
        sampleinfo <- read.csv(gsub("DS_attributes", "sample_attributes", SN_attr_files[j]))
        
        if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
          
          glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
          
          if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
            
            analyte_file <- gsub("_DS_attributes", "", SN_attr_files[j])
            analyte_file <- sapply(strsplit(analyte_file, "/"), function(x){
              x[length(x)]
            })
            
            return(data.frame(Dataset=input_dataset, 
                              Analyte_File=analyte_file))
            
          }
          
        }
        
      }
      
    }
    
  })
  example10 <- do.call(rbind, example10_list)
  
  print("Sheet 1 example 10 complete.")
  
  write.csv(example10, file=paste0("sheet1_example10_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 11 ############################################# 
  
  print("Starting sheet 1 example 11.")
  
  ## Find all Gene Set Enrichment Results (.pdf) for Covariation Networks produced from adult human glioma 
  ## samples using the top 1% or higher of biweight midcorrelations and a minimum module size ≥ 20
  
  example11_list <- lapply(1:nrow(data_dirs), function(i){

    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){
        
        glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
        
        if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))>0){
          
          ## Identify networks with relsig >= .99 & min size >= 20:
          
          netstats <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="network_statistics", full.names=T))
          netstats <- netstats[,intersect(which(netstats[netstats$Network=="RelSignum",]>=.99), 
                                          which(netstats[netstats$Network=="MinSize",]>=20))]
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum")
          networks <- networks[is.element(make.names(networks), colnames(netstats))]
          networks <- file.path(data_dirs$FM_dir[i], networks)
          
          ## Get signifcantly enriched gene sets:
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            pdf_files <- list.files(path=networks[j], pattern="enrichment.*pdf", 
                                    full.names=T, recursive=T)
            
            if(length(pdf_files)>0){
              pdf_files <- sapply(strsplit(pdf_files, "/"), function(x) x[length(x)])
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], 
                                Network=network, 
                                PDF=pdf_files))
            }
            
          }, future.seed=T) 
          
          return(do.call(rbind, networks_list))
          
        }
        
      }
      
    } 
    
  })
  example11 <- do.call(rbind, example11_list)
  
  print("Sheet 1 example 11 complete.")
  
  write.csv(example11, file=paste0("sheet1_example11_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 12 ############################################# 
  
  print("Starting sheet 1 example 12.")
  
  ## Export a list of unique gene symbols for covariation modules identified in bulk RNA-seq datasets 
  ## that include adult human male gliomas and are maximally enriched with markers of radial glia (cell type)
  
  example12_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      
      if(DS_attr$Value[DS_attr$Attribute=="Technology"]=="Sequencer"){
        
        sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
        
        if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
          
          glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
          
          if(length(intersect(grep("M", sampleinfo$Sex), 
                              intersect(grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T), 
                                        glioma_samples)))>0){
            
            networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
            networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
            
            networks_list <- future_lapply(1:length(networks), function(j){
              
              enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
              enrich_out <- covariation_enrich_search(enrich, setname="radial_glia")
              if(nrow(enrich_out)>0){
                DS_attr <- list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1]
                network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
                return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out, 
                                  DS_Attr=DS_attr))
              }
              
            }, future.seed=T)
            modules <- do.call(rbind, networks_list)
            
            if(nrow(modules)>0){
              
              modules <- modules |> dplyr::slice_min(Pval, with_ties=T)
              
              ## Get module genes
              
              gene_list <- lapply(1:nrow(modules), function(j){
                
                kME <- fread(list.files(path=file.path(data_dirs$FM_dir[i], modules$Network[j]), 
                                        pattern="kME", full.names=T), data.table=F)
                col <- kME[,grep(toupper(modules$Mod_Def[j]), toupper(colnames(kME)))]
                kME <- kME[is.element(col, as.character(modules$Module[j])),]
                
                ## Map identifiers to SYMBOL:
                
                DS_attr <- read.csv(modules$DS_Attr[j])
                unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
                platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
                
                if(unique_id=="SYMBOL"){
                  
                  kME <- mapAlias2Symbol(features=kME, unique_id_col=2, 
                                         tables_dir, keep_all=T, fill_NAs=T)
                  
                } else {
                  
                  kME <- map2Any(features=kME, unique_id, map_to="SYMBOL", 
                                 unique_id_col=2, platform, tables_dir, 
                                 keep_all=T)
                  
                }
                
                kME <- kME[,!is.element(colnames(kME), "SYMBOL.y")]
                
                kME <- kME |> 
                  tidyr::separate_rows(SYMBOL, sep=" \\| ") 
                
                return(data.frame(SYMBOL=unique(kME$SYMBOL),
                                  Dataset=modules$Dataset[j], Network=modules$Network[j],
                                  SetName=modules$SetName[j], Module=modules$Module[j],
                                  Mod_Def=modules$Mod_Def))
                
              })
              
              return(do.call(rbind, gene_list))
              
            } 
            
          } 
          
        } 
        
      }
      
    }
    
  })
  example12 <- do.call(rbind, example12_list)
  
  write.csv(example12, file=paste0("sheet1_example12_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example12 <- example12 |> 
    na_if("NA") |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::group_by(SYMBOL) |> 
    dplyr::summarise(No.Datasets=n()) |>
    dplyr::arrange(desc(No.Datasets))
  
  print("Sheet 1 example 12 complete.")
  
  write.csv(example12, file=paste0("sheet1_example12_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 14 ############################################# 
  
  print("Starting sheet 1 example 14.")
  
  ## Find all unique human Entrez IDs for covariation modules generated EXCLUSIVELY from adult human gliomas 
  ## that were maximally enriched with markers of OPCs (cell type)
  
  example14_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){

      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
        
        glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
        
        if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
            enrich_out <- covariation_enrich_search(enrich, setname="OPC")
            if(nrow(enrich_out)>0){
              DS_attr <- list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1]
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out, 
                                DS_Attr=DS_attr))
            }
            
          }, future.seed=T)
          modules <- do.call(rbind, networks_list)
          
          if(nrow(modules)>0){
            
            modules <- modules |> dplyr::slice_min(Pval, with_ties=T)
            
            ## Get module genes
            
            gene_list <- lapply(1:nrow(modules), function(j){
              
              kME <- fread(list.files(path=file.path(data_dirs$FM_dir[i], modules$Network[j]), 
                                      pattern="kME", full.names=T), data.table=F)
              col <- kME[,grep(toupper(modules$Mod_Def[j]), toupper(colnames(kME)))]
              kME <- kME[is.element(col, as.character(modules$Module[j])),]
              
              ## Map identifiers to SYMBOL:
              
              DS_attr <- read.csv(modules$DS_Attr[j])
              unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
              platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
            
              if(unique_id!="ENTREZID"){
                
                kME <- map2Any(features=kME, unique_id, map_to="ENTREZID", 
                               unique_id_col=2, platform, tables_dir, 
                               keep_all=T)
                
              }
              
              kME <- kME |> 
                tidyr::separate_rows(ENTREZID, sep=" \\| ") 
              
              return(data.frame(ENTREZID=unique(kME$ENTREZID),
                                Dataset=modules$Dataset[j], Network=modules$Network[j],
                                SetName=modules$SetName[j], Module=modules$Module[j],
                                Mod_Def=modules$Mod_Def))
              
            })
            
            return(do.call(rbind, gene_list))
            
          } 
          
        } 
        
      } 
      
    }
    
  })
  example14 <- do.call(rbind, example14_list)
  
  write.csv(example14, file=paste0("sheet1_example14_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example14 <- example14 |> 
    na_if("NA") |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::group_by(ENTREZID) |> 
    dplyr::summarise(No.Datasets=n()) |>
    dplyr::arrange(desc(No.Datasets))
  
  print("Sheet 1 example 14 complete.")
  
  write.csv(example14, file=paste0("sheet1_example14_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 15 ############################################# 
  
  print("Starting sheet 1 example 15.")
  
  ## Return a list of all unique enzymes encoded by genes from covariation modules identified 
  ## in bulk gene expression datasets generated EXCLUSIVELY from adult human gliomas 
  ## that were MAXIMALLY enriched with markers of T cells
  
  example15_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
        
        glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
        
        if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
            enrich_out <- covariation_enrich_search(enrich, setname="t_cell")
            if(nrow(enrich_out)>0){
              DS_attr <- list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1]
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out, 
                                DS_Attr=DS_attr))
            }
            
          }, future.seed=T)
          modules <- do.call(rbind, networks_list)
          
          if(nrow(modules)>0){
            
            modules <- modules |> dplyr::slice_min(Pval, with_ties=T)
            
            ## Get module genes
            
            gene_list <- lapply(1:nrow(modules), function(j){
              
              kME <- fread(list.files(path=file.path(data_dirs$FM_dir[i], modules$Network[j]), 
                                      pattern="kME", full.names=T), data.table=F)
              col <- kME[,grep(toupper(modules$Mod_Def[j]), toupper(colnames(kME)))]
              kME <- kME[is.element(col, as.character(modules$Module[j])),]
              
              ## Map identifiers to SYMBOL:
              
              DS_attr <- read.csv(modules$DS_Attr[j])
              unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
              platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
              
              kME <- map2Any(features=kME, unique_id, map_to="ENZYME", 
                             unique_id_col=2, platform, tables_dir, 
                             keep_all=T)
              
              kME <- kME |> 
                tidyr::separate_rows(ENZYME, sep=" \\| ") 
              
              return(data.frame(ENZYME=unique(kME$ENZYME),
                                Dataset=modules$Dataset[j], Network=modules$Network[j],
                                SetName=modules$SetName[j], Module=modules$Module[j],
                                Mod_Def=modules$Mod_Def))
              
            })
            
            return(do.call(rbind, gene_list))
            
          } 
          
        } 
        
      } 
      
    }
    
  })
  example15 <- do.call(rbind, example15_list)
  
  example15 <- example15 |>
    na_if("NA") |>
    dplyr::filter(!is.na(ENZYME))
    
  write.csv(example15, file=paste0("sheet1_example15_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example15 <- example15 |> 
    dplyr::group_by(ENZYME) |> 
    dplyr::summarise(No.Datasets=n()) |>
    dplyr::arrange(desc(No.Datasets))
  
  print("Sheet 1 example 15 complete.")
  
  write.csv(example15, file=paste0("sheet1_example15_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 16 ############################################# 
  
  print("Starting sheet 1 example 16.")
  
  ## Return a list of unique OMIM IDs associated with seed genes for covariation modules identified
  ## in bulk gene expression datasets generated EXCLUSIVELY from adult human glioma samples that include
  ## GATA1, GATA3, and AIRE as seed genes, ranked by frequency of occurrence
  
  feature_list <- c("GATA1", "GATA3", "AIRE")
  
  example16_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
      
      if(is.element("Developmental_Epoch", colnames(sampleinfo))){ 
        
        glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
        
        if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))==nrow(sampleinfo)){
          
          networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
          networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
          
          DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
          unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
          platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
          
          networks_list <- future_lapply(1:length(networks), function(j){
            
            kME <- fread(list.files(path=networks[j], pattern="kME", full.names=T), data.table=F)
            
            ## Map to SYMBOL to find modules that contain features from feature list:
            
            if(unique_id=="SYMBOL"){
              
              kME1 <- mapAlias2Symbol(features=kME, unique_id_col=2, 
                                      tables_dir, keep_all=T, fill_NAs=T)
              
            } else {
              
              kME1 <- map2Any(features=kME, unique_id, map_to="SYMBOL",
                              unique_id_col=2, platform, tables_dir, 
                              keep_all=T)
              
            }
            
            ## Get names of modules:
            
            feature_mods <- mod2list(kME1, feature_type="SYMBOL", feature_list, mod_def="Seed", and_or="AND")
            
            if(length(feature_mods)>0){
              
              ## Get unique identifiers from matching modules and map to OMIM:
              
              mod_col <- grep(mod_def, colnames(kME))
              kME <- kME[is.element(kME[,mod_col], names(feature_mods)),]
              kME <- map2Any(features=kME, unique_id, map_to="OMIM", 
                             unique_id_col=2, platform, tables_dir, 
                             keep_all=T)
              kME <- kME |> na_if("NA") |> na_if("") |>
                tidyr::separate_rows(OMIM, sep=" \\| ") |>
                dplyr::filter(!is.na(OMIM)) |>
                as.data.frame()
              
              if(nrow(kME)>0){
                
                kME <- kME |>
                  dplyr::select(
                    UNIQUE.ID, OMIM, 
                    !!sym(colnames(kME)[grep(mod_def, colnames(kME))])
                  )
                
                colnames(kME)[grep("Mod", colnames(kME))] <- "Module"
                network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
                return(data.frame(Dataset=data_dirs$Title[i], Network=network, 
                                  kME, Mod_Def="Seed"))
                
              }
              
            }
            
          }, future.seed=T) 
          
          return(do.call(rbind, networks_list))
          
        }
        
      }
      
    } 
    
  })
  example16 <- do.call(rbind, example16_list)
  
  if(!is.null(example16)>0){
    
    write.csv(example16, file=paste0("sheet1_example16_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
    
    example16 <- example16 |> 
      na_if("NA") |>
      dplyr::filter(!is.na(OMIM)) |>
      dplyr::group_by(OMIM) |> 
      dplyr::summarise(No.Datasets=n()) |>
      dplyr::arrange(desc(No.Datasets))
    
    print("Sheet 1 example 16 complete.")
    
    write.csv(example16, file=paste0("sheet1_example16_search_results_", Sys.Date(), ".csv"), row.names=F)
    
   
  } else {
    print("No results found for example query 16.")
  }
  
  ############################################# Sheet 1, example 18 ############################################# 
  
  print("Starting sheet 1 example 18.")
  
  ## Find all covariation modules for which >= 1% of seed genes belong to the Gene Ontology Biological Processs
  ## 'oxidative phosphorylation' (GO:0006119), ranked by percent match
  
  feature_list <- c("GO:0006119")
  
  example18_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
      networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      unique_id <- DS_attr$Value[DS_attr$Attribute=="Unique Identifier"]
      platform <- DS_attr$Value[DS_attr$Attribute=="Mapping Tables"]
      
      networks_list <- future_lapply(1:length(networks), function(j){
        
        kME <- fread(list.files(path=networks[j], pattern="kME", full.names=T), data.table=F)
        
        ## Get # of features in each module prior to mapping:
        
        kME <- kME[kME$ModSeed!="",]
        mod_lengths <- table(kME$ModSeed)
        
        ## Map to GO terms to find modules that contain features from feature list:
        
        kME <- map2Any(features=kME, unique_id, map_to="GO",
                       unique_id_col=2, platform, tables_dir, 
                       keep_all=T)
        kME <- kME |> 
          tidyr::separate_rows(GO, sep=" \\| ") 
        
        ## Calculate fraction of module seed genes that are associated with GO term:
        
        GO_mods <- tapply(kME$GO, kME$ModSeed, "[")
        mod_lengths <- mod_lengths[match(names(GO_mods), names(mod_lengths))]
        GO_percent <- unlist(lapply(1:length(GO_mods), function(k){
          sum(is.element(GO_mods[[k]], feature_list))/mod_lengths[k]*100
        }))
        GO_percent <- GO_percent[GO_percent>=1]
        
        if(length(GO_percent)>0){
          
          network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
          return(data.frame(Dataset=data_dirs$Title[i], Network=network, 
                            Module=names(GO_percent), Mod_Def="Seed", 
                            Percent=GO_percent))
          
        }
        
      }, future.seed=T) 
      
      return(do.call(rbind, networks_list))
      
    } 
    
  })
  example18 <- do.call(rbind, example18_list)
  
  example18 <- example18 |>
    dplyr::arrange(desc(Percent))
  
  print("Sheet 1 example 18 complete.")
  
  write.csv(example18, file=paste0("sheet1_example18_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 19 ############################################# 
  
  print("Starting sheet 1 example 19.")
  
  ## Find all unique gene sets associated with human neuronal cell types that were significantly enriched 
  ## (P < 1e-10) in covariation modules derived from adult human glioma samples analyzed by RNA-seq
  
  example19_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      
      if(DS_attr$Value[DS_attr$Attribute=="Technology"]=="Sequencer"){
        
        sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
        
        if(is.element("Developmental_Epoch", colnames(sampleinfo))){
          
          glioma_samples <- get_glioma_samples(mondo_vec=sampleinfo$MONDO_ID)
          
          if(length(intersect(glioma_samples, grep("adult", sampleinfo$Developmental_Epoch, ignore.case=T)))>0){
            
            networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum", full.names=T)
            networks <- networks[unlist(lapply(networks, function(x) length(list.files(path=x))>0))]
            
            networks_list <- future_lapply(1:length(networks), function(j){

              enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
              enrich_out <- covariation_enrich_search(enrich, setname="neuron", pval_cut=1e-10, 
                                                      organism="Homo sapiens", category="CellType")
              if(nrow(enrich_out)>0){
                enrich_out <- enrich_out |>
                  dplyr::group_by(SetID, SetName) |>
                  dplyr::slice_min(Pval, with_ties=F)
                network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
                return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out))
              }
              
            }, future.seed=T) 
            
            return(do.call(rbind, networks_list))
            
          }
          
        }
        
      }
      
    } 
    
  })
  example19 <- do.call(rbind, example19_list)
  
  write.csv(example19, file=paste0("sheet1_example19_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example19 <- example19 |> 
    dplyr::group_by(SetID, SetName) |>
    dplyr::summarise(SetName=unique(SetName),
                     SetID=unique(SetID))
  
  print("Sheet 1 example 19 complete.")

  write.csv(example19, file=paste0("sheet1_example19_search_results_", Sys.Date(), ".csv"), row.names=F)
  
  ############################################# Sheet 1, example 20 ############################################# 
  
  print("Starting sheet 1 example 20.")
  
  ## Find all gene sets that were significantly enriched (P < 1e-10) in covariation modules derived from the top 1% 
  ## of pairwise cors with a minimum module size >= 12 in RNA-seq datasets produced with at least 10% grade 2 oligodendrogliomas
  
  example20_list <- lapply(1:nrow(data_dirs), function(i){
    
    if(!is.na(data_dirs$FM_dir[i])){
      
      DS_attr <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="DS_attributes", full.names=T)[1])
      
      if(DS_attr$Value[DS_attr$Attribute=="Technology"]=="Sequencer"){
        
        sampleinfo <- read.csv(list.files(path=data_dirs$SN_dir[i], pattern="sample_attributes", full.names=T)[1])
        
        if(is.element("Tumor_Grade", colnames(sampleinfo))){
          
          if(length(intersect(grep(2, sampleinfo$Tumor_Grade), grep("oligodendroglioma", sampleinfo$Disease)))>=.1*nrow(sampleinfo)){
            
            ## Identify networks with relsig >= .99 & min size >= 12:
            
            netstats <- read.csv(list.files(path=data_dirs$FM_dir[i], pattern="network_statistics", full.names=T))
            netstats <- netstats[,intersect(which(netstats[netstats$Network=="RelSignum",]>=.99), 
                                            which(netstats[netstats$Network=="MinSize",]>=12))]
            
            networks <- list.files(path=data_dirs$FM_dir[i], pattern="signum")
            networks <- networks[is.element(make.names(networks), colnames(netstats))]
            networks <- file.path(data_dirs$FM_dir[i], networks)
            
            ## Get signifcantly enriched gene sets:
            
            networks_list <- future_lapply(1:length(networks), function(j){
              
              enrich <- list.files(path=networks[j], pattern="GSHyperG", full.names=T)
              enrich_out <- covariation_enrich_search(enrich, pval_cut=1e-10)
              network <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
              return(data.frame(Dataset=data_dirs$Title[i], Network=network, enrich_out))
              
            }, future.seed=T) 
            
            return(do.call(rbind, networks_list))
            
          }
          
        }
        
      }
      
    } 
    
  })
  example20 <- do.call(rbind, example20_list)
  
  write.csv(example20, file=paste0("sheet1_example20_search_results_modules_", Sys.Date(), ".csv"), row.names=F)
  
  example20 <- example20 |> 
    dplyr::group_by(SetID, SetName) |>
    dplyr::summarise(SetName=unique(SetName),
                     SetID=unique(SetID))
  
  print("Sheet 1 example 20 complete.")
  
  write.csv(example20, file=paste0("sheet1_example20_search_results_", Sys.Date(), ".csv"), row.names=F)
  
}

covariation_enrich_search <- function(enrich, setname=NULL, pval_cut=.05, organism=NULL, category=NULL){
  
  enrich_list <- lapply(1:length(enrich), function(j){
    
    enrich_pvals <- fread(list.files(path=enrich[j], pattern=".csv", full.names=T), data.table=F)
    
    if(!is.null(setname)){
      enrich_pvals <- enrich_pvals[grep(setname, enrich_pvals$SetName, ignore.case=T),]
    }
    
    if(!is.null(organism)){
      enrich_pvals <- enrich_pvals[grep(organism, enrich_pvals$Species, ignore.case=T),]
    }
    
    if(!is.null(category)){
      enrich_pvals <- enrich_pvals[grep(category, enrich_pvals$Category, ignore.case=T),]
    }

    if(nrow(enrich_pvals)>0){
      
      temp <- reshape2::melt(enrich_pvals[,c(1,8:ncol(enrich_pvals))])
      colnames(temp) <- c("SetID", "Module", "Pval")
      temp$SetName <- enrich_pvals$SetName[match(temp$SetID, enrich_pvals$SetID)]
      mod_def <- sapply(strsplit(enrich[j], "_"), function(x) x[length(x)])
      return(data.frame(Mod_Def=mod_def, temp))
      
    }
    
  }) 
  enrich_out <- do.call(rbind, enrich_list)
  enrich_out <- enrich_out[enrich_out$Pval<pval_cut,]
  return(enrich_out)
  
}

covariation_feature_search <- function(kME, feature_list, feature_type, 
                                       mod_def=c("BC", "FDR", "Seed"), 
                                       and_or="AND"){
  
  if(mod_def=="Any"){

    modbc <- mod2list(kME, feature_type, feature_list, mod_def="BC", and_or)
    modfdr <- mod2list(kME, feature_type, feature_list, mod_def="FDR", and_or)
    modseed <- mod2list(kME, feature_type, feature_list, mod_def="Seed", and_or)

    if(sum(length(modbc), length(modfdr), length(modseed))>0){
      modules <- c(names(modbc), names(modfdr), names(modseed))
      mod_def <- c(rep("Bonferroni", length(modbc)),
                   rep("FDR", length(modfdr)),
                   rep("Seed", length(modseed)))
      return(data.frame(Module=modules,
                        Mod_Def=mod_def))
    }

  } else {

    modules <- mod2list(kME, feature_type, feature_list, mod_def, and_or)
    if(length(modules)>0){
      return(data.frame(Module=names(modules), Mod_Def=mod_def))
    }

  }

}

mod2list <- function(kME, feature_type, feature_list=NULL, 
                     mod_def=c("BC", "FDR", "Seed"), 
                     and_or=c("AND", "OR")){
  
  mod_col <- grep(mod_def, colnames(kME))
  kME <- kME |> na_if("")
  kME <- kME[!is.na(kME[,mod_col]),]
  
  module_list <- tapply(kME[,feature_type], kME[,mod_col], "[")
  module_list <- lapply(module_list, na.omit)
  
  if(!is.null(feature_list)){
    
    match_mods <- unlist(lapply(module_list, function(x){
      if(and_or=="AND"){
        length(intersect(x, feature_list))==length(feature_list)
      } else {
        length(intersect(x, feature_list))>0
      }
    }))
    
    module_list <- module_list[match_mods]
    
  }
  
  return(module_list)
  
}

get_brain_samples <- function(uberon_vec){
  
  ## Iterate to get ALL children of 'brain':

  id <- UBERON$id[UBERON$name=="brain"]
  brain_children_current <- unlist(UBERON$children[names(UBERON$children)==id], 
                                   use.names=F) 
  
  brain_children_all <- id
  brain_children_previous <- c()
  while(length(brain_children_previous)<length(brain_children_current)){
    brain_children_previous <- brain_children_current
    brain_children_current <- unlist(
      UBERON$children[is.element(names(UBERON$children), 
                                 brain_children_previous)], 
      use.names=F
    )
    brain_children_all <- c(brain_children_all, brain_children_current)
  }
  brain_children <- unique(c(brain_children_previous, brain_children_all))
  
  uberon_vec[is.na(uberon_vec)] <- "NA"
  unlist(sapply(1:length(uberon_vec), function(i){
    if(sum(is.element(brain_children, uberon_vec[i]))>0) 
      return(i)
  }))
  
}

get_glioma_samples <- function(mondo_vec){
  mondo_vec[is.na(mondo_vec)] <- "NA"
  id <- MONDO$id[MONDO$name=="glioma"]
  unlist(sapply(1:length(mondo_vec), function(i){
    parents <- get_ancestors(MONDO, terms=mondo_vec[i])
    if(sum(is.element(parents, id))>0) return(i)
  }))
}

make_data_directory <- function(root_dir){
  
  DC <- list.dirs(path=root_dir, full.names=T, recursive=F)
  
  suppressWarnings(rm(DS_directory))
  
  for(i in 1:length(DC)){
    
    print(DC[i])
    
    makeDC <- list.files(path=list.files(path=DC[i], pattern="makeDC", full.names=T), 
                         pattern="_makeDC_", full.names=T)
    
    if(length(makeDC)>0){
     
      DS <- list.files(path=makeDC, pattern="_DS[0-9]", full.names=T)
      SN_FM_dirs <- list.dirs(path=file.path(DC[i], "all"), full.names=T, recursive=F)
      
      if(length(SN_FM_dirs)==0){
        SN_FM_dirs <- list.dirs(path=DC[i], full.names=T, recursive=F)
      }
      
      if(sum(grepl("by_region", SN_FM_dirs))>0){
        SN_FM_dirs <- c(SN_FM_dirs, list.dirs(path=file.path(DC[i], "by_region"), 
                                              full.names=T, recursive=F))
      }
      
      SN <- SN_FM_dirs[grep("makeSN", SN_FM_dirs)]
      FM <- SN_FM_dirs[grep("makeFM", SN_FM_dirs)]
      
      DS_list <- lapply(1:length(DS), function(j){
        
        DS_attr <- read.csv(list.files(path=DS[j], pattern="DS[0-9]+_attr", full.names=T))
        title <- paste(DS_attr$Value[DS_attr$Attribute=="Title"], DS_attr$Value[DS_attr$Attribute=="Version"])
        return(data.frame(Title=title, DC_dir=DS[j]))
        
      })
      
      SN_list <- lapply(1:length(SN), function(j){
        
        SN_run <- list.files(path=SN[j], pattern="SampleNetworks", full.names=T)
        
        run_list <- lapply(1:length(SN_run), function(i){
          
          makeSN <- list.files(path=SN_run[i], pattern="_makeSN", full.names=T)
          groups <- list.dirs(path=makeSN, full.names=T, recursive=F)
          makeSN_list <- lapply(1:length(groups), function(k){
            
            OR <- list.files(path=groups[k], full.names=T)[grep("outlier", list.files(path=groups[k]))]
            OR_attr <- read.csv(OR[grep("DS_attributes", OR)])
            title <- paste(OR_attr$Value[OR_attr$Attribute=="Title"], 
                           OR_attr$Value[OR_attr$Attribute=="Version"])
            
            makeSN_directory <- data.frame(Title=title, SN_dir=groups[k])
            
            QN <- list.files(path=groups[k], full.names=T)[grep("Qnorm", list.files(path=groups[k]))]
            
            if(length(QN)>0){
              
              QN_attr <- read.csv(QN[grep("DS_attributes", QN)])
              title <- paste(QN_attr$Value[QN_attr$Attribute=="Title"], 
                             QN_attr$Value[QN_attr$Attribute=="Version"])
              makeSN_directory <- rbind(makeSN_directory, data.frame(Title=title, SN_dir=groups[k]))
              
            }
            
            BC <- list.files(path=groups[k], full.names=T)[grep("ComBat", list.files(path=groups[k]))]
            
            if(length(BC)>0){
              
              BC_attr <- read.csv(BC[grep("DS_attributes", BC)])
              title <- paste(BC_attr$Value[BC_attr$Attribute=="Title"], 
                             BC_attr$Value[BC_attr$Attribute=="Version"])
              makeSN_directory <- rbind(makeSN_directory, data.frame(Title=title, SN_dir=groups[k]))
              
            }
            
            return(makeSN_directory)
            
          })
          return(do.call(rbind, makeSN_list))
          
        })
        
        return(do.call(rbind, run_list))
        
      })
      
      FM_list <- lapply(1:length(FM), function(j){
        
        makeFM <- list.files(path=list.files(path=FM[j], pattern="_Modules", full.names=T), 
                             pattern="_makeFM_", full.names=T)
        
        makeFM_list <- lapply(1:length(makeFM), function(i){
          
          FM_attr <- read.csv(list.files(path=makeFM[i], pattern="project_attributes", full.names=T))
          title <- FM_attr$Value[FM_attr$Attribute=="Input Datasets"]
          return(c(Title=title, FM_dir=makeFM[i]))
          
        })
        
        return(do.call(rbind, makeFM_list))
        
      })
      
      if(exists("DS_directory")){
        
        DS_directory <- rbind(DS_directory, do.call(rbind, DS_list))
        SN_directory <- rbind(SN_directory, do.call(rbind, SN_list))
        FM_directory <- rbind(FM_directory, do.call(rbind, FM_list))
        
      } else {
        
        DS_directory <- do.call(rbind, DS_list)
        SN_directory <- do.call(rbind, SN_list)
        FM_directory <- do.call(rbind, FM_list)
        
      }
      
    }
    
  }
  
  data_directory <- merge(DS_directory, SN_directory, by="Title", all=T)
  data_directory <- merge(data_directory, FM_directory, by="Title", all=T)
  
  return(data_directory)
  
}
