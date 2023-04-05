## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(DropletQC)
library(future)

library(Rmagic)
library(phateR)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

#########################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 8)
options(future.globals.maxSize = 40000 * 1024^2)
options(future.rng.onMisuse = 'ignore')

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
save_fn = list.files(here('data/raw_data/Seurat_objects'),
                     pattern = 'STARsolo_SoupX_rawCounts', full.names = T)
names(save_fn) = basename(save_fn) %>% ss('STARsolo_SoupX_rawCounts_|.rds', 2)
num_samples = length(save_fn)
objList = lapply(save_fn, readRDS)

########################################################
# 2) use Seurat reciprocal PCA to join samples together
## find integrating features
features <- SelectIntegrationFeatures(object.list = objList, nfeatures = 3000)
objList <- PrepSCTIntegration(object.list = objList, anchor.features = features)
objList <- lapply(X = objList, FUN = RunPCA, features = features, verbose = FALSE)

## select representative caudate and putamen samples for integration
## both these samples are control subjects/good number cells & depth
## and one of each sex based on:
## https://satijalab.org/seurat/articles/integration_large_datasets.html
ref = which(names(objList) %in% c('S1_1', 'S3_1', 'S5_1', 'S8_1'))

## find pair-wise anchoring cells between samples and each reference
obj.anchors <- FindIntegrationAnchors(
  object.list = objList, normalization.method = "SCT", reference = ref,
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

## merging samples together into a joint space
obj_merged <- obj.anchors %>% 
  IntegrateData(normalization.method = "SCT", dims = 1:30) %>%
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 2, algorithm = 2, verbose = TRUE)


#####################################
# 3) add in monkey/sample metadata
XIST = GetAssayData(object = obj_merged[["RNA"]], slot = "data")['XIST',]
df_sex = cbind(obj_merged[[]], XIST = XIST) %>% group_by(orig.ident) %>% 
  summarise(XIST = max(XIST), 
            Sex = case_when(XIST > 0 ~ 'Female', T ~ 'Male') %>% factor()) %>% 
  dplyr::rename('sample_name' = 'orig.ident') %>% dplyr::select(-XIST)

pheno = readxl::read_xlsx(here('data/raw_data/tables/fenster2020.xlsx')) %>%
  rename_with(make.names) %>%
  rename_with(~ gsub("(\\.){2,}", '\\.', .x)) %>%
  rename_with(~ gsub('\\.$', '', .x)) %>%
  mutate(Monkey = ss(sample_name, '_', 1)) %>% 
  inner_join(df_sex) %>% 
  column_to_rownames('sample_name')

## take a look at the phenotype table
head(pheno)

## look at the per-cell meta data
head(obj_merged@meta.data)

## append pt. phenotypes to single cell meta data 
obj_merged@meta.data = cbind(obj_merged@meta.data, pheno[obj_merged[[]][,'orig.ident'],])
head(obj_merged@meta.data)

###############################################
# 5) filter the lower quality cells
## look at which clusters should be kept by miQC per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(miQC.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by doublet SCDS per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(scds.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by both metrics
(t1 = obj_merged@meta.data %>% group_by(seurat_clusters) %>%
    summarise(num = n(), 
              numKeep = sum(scds.keep == 'keep' & miQC.keep == 'keep'), 
              prop = sum(numKeep) / n() ) %>% arrange(prop))

## keep cells in the clusters that have more than 50% of OK cells
good_clusters <- t1 %>% filter(prop > 0.50) %>% pull(seurat_clusters)

## export unfiltered per-cell QC metric table
dir.create('data/tidy_data/tables')
save_qcTale_fn = here('data/tidy_data/tables', 
                      paste0("McLean_chronic_opioid_monkey_unfiltered_QC_table_N",num_samples,'.txt.gz'))
write_tsv(obj_merged@meta.data, save_qcTale_fn)

## subset cells to those not predicted low QC or doublet
obj_filtered = subset(obj_merged, subset = miQC.keep == "keep" & 
                        scds.keep == "keep" &  seurat_clusters %in% good_clusters)

## recompute PCA and UMAP embedding post-filtering
obj_filtered = obj_filtered %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 2, algorithm = 2, verbose = TRUE)

##############################################################
# 4) smooth out the genes with imputation across full data
DefaultAssay(obj_filtered) <- "RNA"
obj_filtered <- obj_filtered %>%NormalizeData() %>% RunALRA(Assay="RNA")

#########################################################
# 6) save projects, as h5 object for on-disk computation
dir.create('data/tidy_data/Seurat_projects')
save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
                       paste0("McLean_chronic_opioid_monkey_filtered_SCT_SeuratObj_N",num_samples,'.h5Seurat'))
SaveH5Seurat(obj_filtered, filename = save_proj_h5_fn,overwrite = TRUE)

## save normalized, UMAP embedded, object for downstream analyses
save_proj_fn = here('data/tidy_data/Seurat_projects', 
                      paste0("McLean_chronic_opioid_monkey_filtered_SCT_SeuratObj_N",num_samples,'.rds'))
saveRDS(obj_filtered, save_proj_fn)

