library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(here)
library(rcartocolor)

PLOTDIR='figures/exploratory/preprocess_snRNA-seq'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Shell/OT', 'D2-Shell/OT',  'D1-NUDAP', 'D1-ICj')
subtypes_col = c('#1f78b4',  '#e31a1c',  '#6a3d9a', 'green')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(4, "Safe"), 
                   carto_pal(length(othertypes) -4 , "Vivid"))
names(othertypes_col) = othertypes
typecolors = c(subtypes_col, othertypes_col)

## plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

###################################
# 1) read in Logan snRNA dataset to plot
## load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "OUD_Striatum_refined_all_SeuratObj_N16.h5Seurat") %>% 
  LoadH5Seurat()

names(obj_merged[[]] )
table(obj_merged$celltype3 )

fig1_split_umap_allCells_fn = 
  here(PLOTDIR, 'plots', 'figNN_McLean_monkey_umap_NAc.pdf')

pdf(fig1_split_umap_allCells_fn, width = 100/in2mm, height =  80/in2mm)
DimPlot(object = obj_merged, reduction = "umap", group.by = 'celltype3', label = T,
        label.size = 1.8, pt.size = 1, ncol = 1, cols = typecolors, raster = T) +
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 2), byrow = F)) +
  theme(legend.position = 'right', plot.title= element_blank(),
        legend.spacing.x = unit(2, 'mm'),
        legend.spacing.y = unit(2, 'mm'),
        legend.key.size = unit(3, "mm"), 
        plot.margin = unit(rep(1,4), "mm"))
dev.off()



