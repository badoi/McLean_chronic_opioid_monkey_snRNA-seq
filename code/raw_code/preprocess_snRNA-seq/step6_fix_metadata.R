library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(here)

## load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "OUD_Striatum_refined_all_SeuratObj_N16.h5Seurat") %>% 
  LoadH5Seurat()

DefaultAssay(object = obj_merged) = 'RNA'
head(obj_merged[[]])

## read in the correct metadata
df = here('data/raw_data/tables/RM_metadata_030923_Age.xlsx') %>% 
  readxl::read_xlsx() %>% rename_with(make.names) %>% 
  rename('Age' = 'Age..Y.', 'Weight' = 'Weight..kg.') %>% 
  mutate(Monkey = gsub('ample', '', Sample))

meta = obj_merged[[]] %>% rownames_to_column('cellBarcode') %>% 
  inner_join(df) %>% column_to_rownames('cellBarcode') %>% 
  dplyr::select(-c(Sample, Monkey))
meta = meta[colnames(obj_merged),]

obj_merged2 = obj_merged
for(n in names(meta)){
  obj_merged2 = AddMetaData(obj_merged2, meta %>% pull(all_of(n)), col.name = n)
}

table(obj_merged2$Sex)
table(obj_merged2$Pair)

## Save the Seurat object w/ the correct metadata
obj_merged2 %>% SaveH5Seurat(filename= here('data/tidy_data/Seurat_projects', 
                              "OUD_Striatum_refined_all_SeuratObj_N16.h5Seurat"), 
                             overwrite = T)


## save in the MSN dataset for label refinement
msn_type = grep('^D', obj_merged2$celltype3, value = T) %>% unique()
obj_msn2 = obj_merged2 %>% 
  subset(subset = celltype3 %in% msn_type) %>% 
  RunPCA(verbose = FALSE, assay = 'integrated') %>% 
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

save_subset_msn = here('data/tidy_data/Seurat_projects', 
                       "OUD_Striatum_refined_msn_SeuratObj_N16.h5Seurat")
SaveH5Seurat(obj_msn2, filename = save_subset_msn, overwrite = TRUE)

