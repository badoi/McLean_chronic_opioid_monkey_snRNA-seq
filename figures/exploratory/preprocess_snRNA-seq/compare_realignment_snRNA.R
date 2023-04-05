library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(here)

FIGDIR='figures/exploratory/preprocess_snRNA-seq'
sapply(file.path(FIGDIR, c('plots', 'tables', 'rdas')), dir.create, showWarnings =F)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

pheno = readxl::read_xlsx(here('data/raw_data/tables/fenster2020.xlsx')) %>%
  rename_with(make.names) %>%
  rename_with(~ gsub("(\\.){2,}", '\\.', .x)) %>%
  rename_with(~ gsub('\\.$', '', .x)) %>%
  mutate(Monkey = ss(sample_name, '_', 1))

######################################################
## 1) load in the unfiltered realigned files
STARsolo_fn = 'data/raw_data/STARsolo_out2' %>% 
  list.dirs(full.names = T,recursive = T) %>% str_subset('GeneFull')%>%
  str_subset('filtered$')%>%
  str_subset('/\\.', negate = TRUE)  

names(STARsolo_fn) = ss(STARsolo_fn,"(STARsolo_out2/)|(.Solo.out)",2)

## loops over each file name, runs the Read10X function
## create a SeuratObject, one for each 10x sample
dataList <- lapply(STARsolo_fn, Read10X, strip.suffix = T)
objList <- mapply(CreateSeuratObject, min.cells = 3, min.features = 200,
                  counts = dataList, project = names(dataList))
rm(dataList); gc() ## this file is really big. don't need anymore
df.new_unfiltered = lapply(objList, function(obj) obj[[]] %>% rownames_to_column('CellBarcode')) %>% 
  rbindlist() %>% dplyr::rename("sample_name" = 'orig.ident') %>% 
  inner_join(pheno)

dim(df.new_unfiltered) #  29241      7

## get the unfiltered counts from original analyses
obj.old_unfiltered = readRDS('data/tidy_data/fenster_round1_analyses/seurat.bcbio.RDS')
head(obj.old[[]])
df.old_unfiltered = obj.old_unfiltered[[]] %>% rownames_to_column('CellBarcode') %>% 
  mutate(bcbio_id = ss(CellBarcode, ':'), 
         CellBarcode = CellBarcode %>% ss(':', 2) %>% str_replace('-', '')) %>% 
  dplyr::select(-orig.ident) %>% 
  dplyr::rename('nCount_RNA_old' = 'nCount_RNA', 'nFeature_RNA_old' = 'nFeature_RNA')

dim(df.old_unfiltered) #  187980      4

df_unfiltered = inner_join(df.new_unfiltered, df.old_unfiltered)
dim(df_unfiltered) #  26953    9

######################################################
## 2) plot the per-cell gene w/ the unfiltered alignments
plot_fn1 = here(FIGDIR, 'plots', 'compare_realignment_nCountRNA_unfiltered.pdf')
pdf(plot_fn1, height = 8, width = 6)
ggplot(df_unfiltered, aes(x = nCount_RNA_old, y = nCount_RNA)) +
  geom_point(pch = 20, alpha = .5) + 
  geom_abline(slope = 1, intercept =0, color = 'black') +
  facet_wrap(~sample_name) + 
  theme_bw()+ ggtitle('Raw aligned UMI per cell')
dev.off()
  
plot_fn2 = here(FIGDIR, 'plots', 'compare_realignment_nFeature_RNA_unfiltered.pdf')
pdf(plot_fn2, height = 8, width = 6)
ggplot(df_unfiltered, aes(x = nFeature_RNA_old, y = nFeature_RNA)) +
  geom_point(pch = 20, alpha = .5) + 
  geom_abline(slope = 1, intercept =0, color = 'black') +
  facet_wrap(~sample_name) + 
  theme_bw()+ ggtitle('Raw aligned detected Gene per cell')
dev.off()


######################################################
## 3) load in the filtered new and old Seurat objects
save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
                       paste0("McLean_chronic_opioid_monkey_filtered_SCT_SeuratObj_N16.h5Seurat"))
obj.new = LoadH5Seurat(file = save_proj_h5_fn)
head(obj.new[[]])
df_new = obj.new[[]] %>% rownames_to_column('CellBarcode') %>% 
  mutate(CellBarcode = ss(CellBarcode, '_'))
dim(df_new) #  23310    18

obj.old = readRDS('data/tidy_data/fenster_round1_analyses/seurat.filtered.RDS')
head(obj.old[[]])
df_old = obj.old[[]] %>% rownames_to_column('CellBarcode') %>% 
  mutate(bcbio_id = ss(CellBarcode, ':'), 
         CellBarcode = CellBarcode %>% ss(':', 2) %>% str_replace('-', '')) %>% 
  dplyr::select(-orig.ident) %>% 
  dplyr::rename('nCount_RNA_old' = 'nCount_RNA', 'nFeature_RNA_old' = 'nFeature_RNA')
dim(df_old) #  31443    13

df = inner_join(df_new, df_old)
dim(df) #  19667    28



#############################################################################
## 4) plot the per-cell gene w/ the SoupX corrected counts vs. filtered cells
plot_fn3 = here(FIGDIR, 'plots', 'compare_realignment_nCountRNA_SoupX.pdf')
pdf(plot_fn3, height = 8, width = 6)
ggplot(df, aes(x = nCount_RNA_old, y = nCount_RNA)) +
  geom_point(pch = 20, alpha = .5) + 
  geom_abline(slope = 1, intercept =0, color = 'black') +
  facet_wrap(~sample_name) + 
  theme_bw() + ggtitle('SoupX corrected UMI per cell')
dev.off()



plot_fn4 = here(FIGDIR, 'plots', 'compare_realignment_nFeature_RNA_SoupX.pdf')
pdf(plot_fn4, height = 8, width = 6)
ggplot(df, aes(x = nFeature_RNA_old, y = nFeature_RNA)) +
  geom_point(pch = 20, alpha = .5) + 
  geom_abline(slope = 1, intercept =0, color = 'black') +
  facet_wrap(~sample_name) + 
  theme_bw() + ggtitle('SoupX corrected detected Gene per cell')
dev.off()







