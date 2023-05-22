ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## num parallel threads for Seurat
plan(sequential) # no parallel
options(future.globals.maxSize = 80e9)
library(here)

DATADIR='data/raw_data'

## make for this subdirs
PLOTDIR='figures/exploratory/preprocess_snRNA-seq'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Shell/OT', 'D2-Shell/OT',  'D1-NUDAP', 'D1-ICj')
subtypes_col = c('#1f78b4',  '#e31a1c',  '#6a3d9a', 'green')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 'Interneuron', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(5, "Safe"), 
                   carto_pal(length(othertypes) -5 , "Vivid"))
names(othertypes_col) = othertypes
typecolors = c(subtypes_col, othertypes_col)

## plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "OUD_Striatum_refined_all_SeuratObj_N16.h5Seurat") %>% 
  LoadH5Seurat()

# set the groupings by the refined cell type column
## merge interneurons again
obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 
                              'Interneuron',obj_merged$celltype3)
table(obj_merged$celltype3)

Idents(object = obj_merged) <- "celltype3"

# change the ordering of the cell types
levels(obj_merged) <- names(typecolors)[names(typecolors) %in% unique(obj_merged$celltype3)]


###################################
# 2) plot the diagonal matrix of reported marker genes & cell types
markerGenes = c(
  'RBFOX3', # NeuN-neurons'
  'GAD2', # Gaba-ergic neurons
  'PPP1R1B', # DARPP-32 MSNs
  'DRD1', # direct MSNs
  'TAC1', # direct MSNs
  'DRD2',  # indirect MSNs
  'PENK',  # indirect MSNs
  'DRD3', # ICj 
  'FOXP2', # NUDAP/ D1/D2H
  'OPRM1', # NUDAP
  'RXFP1', # NUDAP/ D1/D2H
  'LHX6', # MGE interneurons
  'PTHLH', # MGE interneurons
  'TRH', # MGE interneurons
  'SST', # MGE interneurons
  'AQP4', # Astrocytes
  'C3', # Microglia
  'MOG', # Oligodendrocyte
  'PDGFRA' # OPC
)


fig2_diagonal_matrix_dotplot_fn = 
  here(PLOTDIR, 'plots', 'figNN_McLean_monkey_celltype_diagonal_matrix_dotplot.all.pdf')
pdf(fig2_diagonal_matrix_dotplot_fn, width = 180/in2mm, height =  60/in2mm)
DotPlot( obj_merged, features = markerGenes, cols = c('lightgray', 'navyblue'),
         cluster.idents = F, scale = T, scale.by = "radius", assay = 'alra') +
  my_theme + scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  theme(legend.position = 'bottom', plot.title= element_blank(),
        # axis.text.x = element_text(angle = -30, vjust = 0, hjust=1),
        axis.title = element_blank(), 
        legend.key.size = unit(2, "mm"),
        legend.spacing.x = unit(2, 'mm'))
dev.off()


