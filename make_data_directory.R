source("/home/rebecca/omicon/search/search_query_fxns.R")

## Glioma samples:

root_dir <- "/mnt/bdata/rebecca/omicon/glioma_meta_datasets_processed_new_attributes_new"

data_dirs <- make_data_directory(root_dir)

write.csv(data_dirs, file=paste0("/home/rebecca/omicon/search/data_directory_glioma_", Sys.Date(), ".csv"), row.names=F)

## Normal samples:

root_dir <- "/mnt/bdata/rebecca/omicon/kevin_data"

data_dirs <- make_data_directory(root_dir)

write.csv(data_dirs, file=paste0("/home/rebecca/omicon/search/data_directory_normal_", Sys.Date(), ".csv"), row.names=F)