## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)
library(writexl)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## main differential gene expression package
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(Matrix.utils)
library(limma)
library(edgeR)
library(sva)
# BiocManager::install("swfdr")
library(swfdr)

## regress out the surrogate variables
# BiocManager::install("LieberInstitute/jaffelab")
library(jaffelab)

xx.dat.dir = "data/tidy_data/Seurat_projects/"
DATADIR = "data/tidy_data/differential_expression"
h5Seurat.dat.fn = "OUD_Striatum_refined_all_SeuratObj_N8.h5Seurat"
dir.create(here(DATADIR, 'rdas'))
dir.create(here(DATADIR, 'tables'))
dir.create(here(DATADIR, 'plots'))

# 1) load in cell type labels for label transfer
## read in Logan BU snRNA dataset to label transfer
save_merged_fn = here(xx.dat.dir, h5Seurat.dat.fn)
## load only the scaled and raw RNA counts
obj = save_merged_fn %>% LoadH5Seurat(assay = 'RNA') 
head(obj[[]]) #12171 samples in all and 6843 in msn

## convert to SingleCellExperiment for GlmGamPoi
sce = as.SingleCellExperiment(obj)

##############################################################
# 2) compute per-sample per cell type pseudobulk DGE profiles
## save to hdf5 file format for quick save/load
h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
saveHDF5SummarizedExperiment(sce, h5Dir, prefix=h5Seurat.dat.fn, replace=TRUE)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

# xx.dat.dir= "/bgfs/ctseng/xix66/MonkeyScRNA/BaDoi/data/tidy_data/Seurat_projects"
h5Seurat.dat.fn = "OUD_Striatum_refined_all_SeuratObj_N8.h5Seurat"
pseudo.fn = "OUD_Striatum_refined_all_PseudoBulk_N8.sce2.rds"
sva.fn = "OUD_Striatum_refined_all_PseudoBulk_N8.sva.rds"
voom.fn = "voomLimma_norm_object_N8pb.rds"
voom.sva.fn = "voomLimma_norm_sva_object_N8pb.rds"
lmfit.fn = "voomLimma_diffGene_bigModelFit.rds"
lmfit.sva.fn = "voomLimma_diffGene_bigModelFitSVA.rds"
design.fn = "bigModelFit_designMatrix.rds"
design.sva.fn = "bigModelFitSVA_designMatrix.rds"
celltype_prop.fn = "rdas/OUD_Striatum_refined_celltype3_proportion.rds"
sumtab.fn = "plots/voomLimma_nDE_pdf"
sumtab.sva.fn = "plots/voomLimma_sva_nDE_pdf"
DElist.fn = paste0("Monkey_voom", c("_nosva_", "_sva_"), "limma_bigModelSVA_N8.celltype.rds")
DEsheet.fn = paste0("Monkey_voom", c("_nosva_", "_sva_"), "limma_bigModelSVA_N8.celltype.xlsx")

#######################################################
# 0) Seurat uses the future package for parallelization
plan("multicore", workers = 14)
options(future.globals.maxSize = 80 * 1024^3)

##################################################
# 1) create or load pseudobulk sce object
save_pseudobulk =here(DATADIR, 'rdas', pseudo.fn)
if(!file.exists(save_pseudobulk)){
  ## load the single cell counts
  h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
  sce = loadHDF5SummarizedExperiment(h5Dir, prefix=h5Seurat.dat.fn)
  
  ## merge interneurons again
  sce$celltype3 = ifelse(grepl('Int', sce$celltype3), 'Interneuron',sce$celltype3)
  table(sce$celltype3, sce$Monkey)
  #calculate cell proportion 
  a.tab = table(sce$celltype3, sce$Monkey)
  cell_prop = a.tab %>% apply(2, function(x){x/sum(x)}) %>% apply(1, mean)
  #(sce$celltype3 %>% table )/length(sce$celltype3 ) actually very close
  sum(cell_prop)==1
  saveRDS(cell_prop, paste0(DATADIR, "/", celltype_prop.fn))

  ## aggregate by cluster-sample to create pseudo-bulk count matrix
  colData(sce)
  # groups <- colData(sce)[, c("celltype3", "Case", 'Region')]
  groups <- colData(sce)[, c("celltype3", "Monkey")]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  dim(pb) #72 22138
  
  ## split by cluster, transform & rename columns
  pb_colData = colData(sce) %>% as.data.frame() %>%
    rownames_to_column('match') %>% 
    mutate(Pair = factor(Pair), match = paste(celltype3, Monkey, sep = '_')) %>%  
    filter(!duplicated(match)) %>% column_to_rownames('match')
  pb_colData = pb_colData[rownames(pb),]
  
  ## make sure this PD is correct
  # with(pb_colData, table(Case, celltype3, Region))
  with(pb_colData, table(celltype3, Monkey))
  with(pb_colData, table(celltype3, Pair))
  
  ## add number of cells per aggregate
  num_cells = groups %>% as.data.frame() %>% 
    # mutate(tmp = paste(celltype3, Case, Region,sep= '_')) %>% 
      mutate(tmp = paste(celltype3, Monkey,sep= '_')) %>% 
    pull(tmp) %>% table()
  num_cells %>% as.numeric() %>% summary() #should we filter some out?
  pb_colData$numCells = num_cells[rownames(pb_colData)]
  
  ## add the gene detection rate
  pb_colData$cdr <- scale(rowMeans(pb > 0)) 
  
  ## create SingleCellExperiment from pseudo bulk counts across all cell types and region
  (pb <- SingleCellExperiment(assays = t(pb), colData = pb_colData))
  
  ## remap case index nested inside OUD tx
  # remap_case_itx = split(pb$Case, pb$DSM.IV.OUD) %>% 
  remap_case_itx = split(pb$Monkey, pb$condition) %>% 
    lapply(function(x){
      x = setNames(LETTERS[as.numeric(factor(x))], x)
      x[!duplicated(x)]
    }) %>% unlist()
  names(remap_case_itx) = ss(names(remap_case_itx), '\\.', 2)
  # pb$CaseItx = remap_case_itx[as.character(pb$Case)]
  pb$CaseItx = remap_case_itx[as.character(pb$Monkey)]
  
  table(pb$condition, pb$CaseItx)
  table(pb$celltype3, pb$CaseItx)
  table(pb$celltype3, pb$Monkey)
  
  saveRDS(pb, save_pseudobulk)
} else {
  pb = readRDS(save_pseudobulk)
}

####################################################
## 2) filter pseudobulk samples that have too few cells
pb = pb[, pb$numCells > 15]# 70 out of 72 
pb = pb[, pb$celltype3 != 'D1-ICj'] # drop D1-ICj b/c too few animals w/ these
pb$celltype3 = make.names(pb$celltype3) %>% as.factor()
pb$Sex = as.factor(pb$Sex)
pb$Pair = as.factor(pb$Pair)
pb$numCells = as.numeric(pb$numCells)

## make an interaction term for all the combinations
pb$celltype_tx = interaction(pb$condition, pb$celltype3) %>% #female samples are too few to consider sex effect
  as.factor() %>% droplevels()
table(pb$celltype_tx)

## interaction term w/o the tx
# pb$celltype_rg_sex = interaction(pb$celltype3,pb$Sex, pb$Region)  %>% 
pb$celltype0 = interaction(pb$celltype3)  %>% 
  as.factor() %>% droplevels()
table(pb$celltype0)

## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ 0 + celltype_tx  + # term capturing the main effects
                         Pair + numCells  + cdr, # co-variates 
                       data = colData(pb))

## construct the null model, used in regressing out the various factors
design0 <- model.matrix(~ 0 +  celltype0+ Pair + cdr + numCells, data = colData(pb))

####################################
# 3) normalization using voom-limma
y <- DGEList(counts = assays(pb)[[1]])
dim(y) # 22138    64

## filter out genes w/ low counts
A <- rowMeans(y$counts)
isexpr <- A > 2
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y) # 10733   210

## filter out ribosomal genes, filter out mitochondria genes
drop.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = F)
drop.genes %>% sort %>% data.frame() %>% 
  write_tsv(here(DATADIR, 'tables', 'dropped_mito_ribo_genes.tsv'))
keep.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = T)
y = y[keep.genes, , keep.lib.size = FALSE]
dim(y) #  10733   210

# normalize counts
y <- calcNormFactors(y)

## voom precision weights and sample-wise quality weights normalization
v <- voomWithQualityWeights(y, design)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Monkey)
cor$consensus # 0.05806561

## recalculate weights after adjusting for correlated samples from same subject
v <- voomWithQualityWeights(y, design, block = colData(pb)$Monkey, 
          correlation = cor$consensus)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Monkey)
cor$consensus # 0.0581864

#XX: also save v without sva components
save_voom = here(DATADIR, 'rdas', voom.fn)
saveRDS(v, file = save_voom)
# v = readRDS(here::here(DATADIR, 'rdas', voom.fn))

############################################################
# 4) Use surrogate variables to estimate unmodeled variation

## estimate the number of SVs from the adjusted
# save_sva =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sva.rds')
save_sva =here(DATADIR, 'rdas', sva.fn)
if(! file.exists(save_sva)){
  (n.sv = num.sv(v$E, design, method="be", seed = set.seed(1))) #23
  svobj = sva(v$E, design, design0, n.sv=n.sv, B = 20)
  saveRDS(svobj, save_sva)
} else {
  svobj = readRDS(save_sva)
}

## add the SVs to the model matrix
designSV = cbind(design, svobj$sv)
design0SV = cbind(design0, svobj$sv)

## recalculate sample quality weights after calculating the SVs
v.sva <- voomWithQualityWeights(y, designSV, block = colData(pb)$Monkey, 
                            correlation = cor$consensus)
cor.sva <- duplicateCorrelation(v.sva, designSV, block = colData(pb)$Monkey)
cor.sva$consensus # 0.0513175

save_voom_sva = here(DATADIR, 'rdas', voom.sva.fn)
# saveRDS(v, file = save_voom)
saveRDS(v.sva, file = save_voom_sva)

####################################################################
input.list = list(v = list(v, v.sva), cor = list(cor, cor.sva),
                  design = list(design, designSV))
fn.list = list(lmfit.fn = list(lmfit.fn, lmfit.sva.fn), 
               design.fn = list(design.fn, design.sva.fn), 
               sumtab.fn = list(sumtab.fn, sumtab.sva.fn), 
               DElist.fn = DElist.fn, 
               DEsheet.fn = DEsheet.fn, 
               DEup.fn = DEup.fn, 
               DEdown.fn = DEdown.fn)
p.all = list()

for(i in 1:2){

  ##XX: 6.1 fit the model without SVA components
  ## 6.2 fit the model with SVA components
  fit <- lmFit(input.list$v[[i]], input.list$design[[i]], 
               block = colData(pb)$Monkey, 
               correlation = input.list$cor[[i]]$consensus)
  fit <- eBayes(fit, robust = TRUE)
  
  save_fit = here(DATADIR, 'rdas', fn.list$lmfit.fn[[i]])
  saveRDS(fit, file = save_fit)

  save_design = here(DATADIR, 'rdas', fn.list$design.fn[[i]])
  saveRDS(input.list$design[[i]], file = save_design)
  # design = readRDS(here::here(DATADIR, 'rdas', fn.list$design.fn[[1]]))
  # designSV = readRDS(here::here(DATADIR, 'rdas', fn.list$design.fn[[2]]))

  ###########################################################
  ## 7) compute the differences b/t tx within each cell type
  celltypes=levels(factor(pb$celltype3 %>% make.names()))
  designSV = input.list$design[[i]]
  designSV2 =designSV
  colnames(designSV2) = make.names(colnames(designSV2))

  ## make the cell type contrasts
  con_celltypes = sapply(setNames(celltypes, celltypes),function(cell) {
  cell = colnames(designSV2) %>% make.names() %>% str_subset(paste0('\\.', cell, '$'))
  OUD = cell %>% str_subset('morphine'); CTL = cell %>% str_subset('control')

  N_OUD = OUD %>% length(); OUD = OUD %>% paste(collapse = ' + ')
  N_CTL = CTL %>% length(); CTL = CTL %>% paste(collapse = ' + ')
  paste('(',OUD,')/',N_OUD, '-(',CTL,')/',N_CTL)
  })


  ## proportion of each cell type
  # df_prop = 'data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt' %>% 
  #   read_tsv() %>% deframe()
  df_prop = readRDS(paste0(DATADIR, "/", celltype_prop.fn))
  names(df_prop) = names(df_prop) %>% make.names()
  df_prop = df_prop[celltypes]

  ind_neur = grepl('^D|^Int',celltypes)
  ind_glia = !grepl('^D|^Int',celltypes)

  ## create the contrasts for OUD effect Between all cells or major classes
  con_groups = c('All' = paste0('(', con_celltypes,')*', df_prop) %>% paste(collapse = ' + '), 
                 'Neuron' = paste0('(', con_celltypes[ind_neur],')*', 
                                   df_prop[ind_neur]/sum(df_prop[ind_neur])) %>% paste(collapse = ' + '), 
                 'Glia' =  paste0('(', con_celltypes[ind_glia],')*', 
                                df_prop[ind_glia]/sum(df_prop[ind_glia])) %>% paste(collapse = ' + '))

  ## refit the model based on these contrasts
  cont.matrix <- makeContrasts(contrasts= c(con_groups, con_celltypes), levels=designSV2)
  rownames(cont.matrix) = colnames(designSV)
  fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()

  ## compute the DEGs from these contrasts
  deg_list = lapply(setNames(colnames(cont.matrix),  names(c(con_groups, con_celltypes))), 
                    function(coef){
  topTable(coef = coef, fit =fit2, n=Inf) %>% arrange(P.Value) %>% 

  ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
  ## https://pubmed.ncbi.nlm.nih.gov/30581661/
  mutate(adj.P.Val.Within =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
  dplyr::select(-adj.P.Val) %>% rownames_to_column('gene') 
                      })

  ## FDR correction Between all tests
  deg_list = deg_list %>% data.table::rbindlist(idcol = 'celltype') %>% 
    mutate(adj.P.Val.Between =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
    split(by = 'celltype')

  # FDR cutoff
  sapply(deg_list, function(x) x[x$adj.P.Val.Within < 0.05,] %>% nrow())
  sapply(deg_list, function(x) x[x$adj.P.Val.Between < 0.05,] %>% nrow())

  sum.tab = rbind.data.frame(sapply(deg_list, function(x) x[x$adj.P.Val.Between < 0.05,] %>% nrow()), 
                             sapply(deg_list, function(x) x[x$adj.P.Val.Within < 0.05,] %>% nrow()), 
                             sapply(deg_list, function(x) x[x$P.Value < 0.01,] %>% nrow()))
  colnames(sum.tab) = names(deg_list)
  rownames(sum.tab) = c( "adj.P.Val.Between < 0.05", "adj.P.Val.Within < 0.05", "P.Value < 0.01")

  pdf(paste0(DATADIR, "/", fn.list$sumtab.fn[[i]]), height=3, width=15)
  gridExtra::grid.table(sum.tab)
  dev.off()

  ####################################################################
  ## 8) save the output of voom_limma differential state analyses
  rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
  save_res_fn = here(rdasDir, fn.list$DElist.fn[i])
  saveRDS(deg_list, save_res_fn)

  tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
  save_res_fn2 = here(tablesDir, fn.list$DEsheet.fn[i])
  deg_list %>% lapply(function(x) x %>% arrange(P.Value)) %>% writexl::write_xlsx(save_res_fn2)

  p.all[[i]] = do.call(rbind.data.frame, lapply(deg_list, function(a.list){
    data.frame(celltype = a.list[, 1], gene = a.list[, 2], P.value = a.list[, "P.Value"])
  }))
}

saveRDS(p.all, paste0(DATADIR, "/rdas/", "pvals_all.rds"))
p.all = readRDS(paste0(DATADIR, "/rdas/", "pvals_all.rds"))
p.all[[1]]$model = "limma"
p.all[[2]]$model = "limma_sva"
p.all = rbind.data.frame(p.all[[1]], p.all[[2]])
library(ggplot2)
#make histogram
pdf(paste0(DATADIR, "/plots/", "pval_dist.pdf"),height = 12*4, width = 2*4)
ggplot(p.all, aes(x = P.Value))+
  geom_histogram()+
  facet_grid(celltype~model)
dev.off()

#make pair-wise scatter plots
p.all = readRDS(paste0(DATADIR, "/rdas/", "pvals_all.rds"))
colnames(p.all[[1]])[3] = "p_limma"
colnames(p.all[[2]])[3] = "p_limma_sva"
p.all2 = inner_join(p.all[[1]], p.all[[2]])
library(ggplot2)
pdf(paste0(DATADIR, "/plots/", "pval_scatter.pdf"),height = 3*4, width = 4*4)
ggplot(p.all2, aes(x = -log10(p_limma), y = -log10(p_limma_sva)))+
  geom_point(alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gold")+
  facet_wrap(~celltype)
dev.off()

pdf(paste0(DATADIR, "/plots/", "pval_scatter2.pdf"),height = 3*4, width = 4*4)
ggplot(p.all2, aes(x = -(log10(p_limma)+log10(p_limma_sva))/2, 
                   y = -log10(p_limma_sva)+log10(p_limma)))+
  geom_point(alpha = 0.5)+
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "gold")+
  facet_wrap(~celltype)
dev.off()

# mean((p.all2 %>% mutate(x = p_limma_sva-p_limma) %>% filter(celltype=="All") %>% pull(x))>0)

## ------------------------------------------------------------------------------------------------------------
## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
# install.packages("rcartocolor")
library(rcartocolor)
library(ggpubr)

library(data.table)
library(fgsea)
# BiocManager::install("swfdr")
library(swfdr)
# install.packages("msigdbr")
library(msigdbr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DElist.fn = paste0(DATADIR,"/rdas/Monkey_voom", "_sva_", "limma_bigModelSVA_N8.celltype.rds")

## make for this subdirs
PLOTDIR = "figures/exploratory/differential_expression"
paste0(PLOTDIR, "/", c('plots', 'tables', 'rdas')) %>% 
  sapply(dir.create, showWarnings = F, recursive = TRUE)

############################################################
# 1) read in the DEG lists per comparison of OUD vs. Control
res.celltype = DElist.fn %>% readRDS()
names(res.celltype) = paste0(names(res.celltype), '#All')

deg_rank_list1 = lapply(res.celltype, function(deg){
  deg %>% mutate(tmp = -log10(P.Value) * sign(logFC)) %>% 
    dplyr::select(gene, tmp) %>% arrange(tmp) %>% deframe()
})

deg_rank_list = deg_rank_list1

#############################################
## 2) get gene ontologies, use human genes
## grab the H, Hallmark set of gene pathways
## grab the C2, which is the curated canonical pathway sets
## grab the C5, which is the Gene Ontology sets
## https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
pathways_df =  bind_rows(msigdbr("human", category="H"), 
                         msigdbr("human", category="C2"), 
                         msigdbr("human", category="C5"))

## get the SynGO gene ontologies, use human genes
syngo_df = readxl::read_xlsx("/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/SynGO_bulk_download_release_20210225/tables/syngo_annotations.xlsx")
pathways_df = rbindlist(list(pathways_df, syngo_df), fill = T) 

## reshape/label for pathway naming purposes
pathways <-pathways_df %>% 
  mutate(
    gs_subcat = ifelse(is.na(gs_subcat) | gs_subcat == '', gs_cat, gs_subcat),
    gs_name = paste(gs_subcat, gs_name, sep ='#')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)

## exclude the really really big gene sets
lengths(pathways) %>% summary()
pathways = pathways[lengths(pathways)<500]
length(pathways) # 21058

table(pathways_df$gs_cat)

pathways_df2 = pathways_df %>% 
  dplyr::select(gs_subcat, gs_name, gs_description) %>% distinct() %>% 
  dplyr::rename('pathway_group' = 'gs_subcat', 'pathway' = 'gs_name',
                'description' = 'gs_description')

#############################
## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = pathways,
                   minSize=15, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

alpha = 0.05
gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% 
  arrange(pval) %>% filter(!is.na(pval)) %>% 
  mutate(
    MSigDb_Group = ss(pathway, '#', 1), 
    celltype = group %>% ss('#', 1), 
    MORPHINE.v.CTL.in = group %>% ss('#', 1),
    pathway = ss(pathway, '#', 2),
    padj = lm_qvalue(pval, X=size)$q, 
    celltype = group %>% ss('#', 1), 
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  inner_join(pathways_df2) %>% dplyr::select(-group) %>% 
  relocate(MORPHINE.v.CTL.in, celltype, MSigDb_Group, description, .before= everything()) %>% 
  split(f = .$MORPHINE.v.CTL.in)

## save the enrichment w/ all the pathways, significant and otherwise
gsea_df %>% saveRDS(paste0(DATADIR, "/rdas/Monkey_sva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.unfiltered.rds"))
## filter out just the significant pathways
gsea_df = lapply(gsea_df, function(x){
  x %>% filter(padj < alpha)
})
sapply(gsea_df, nrow)

out_fn = paste0(DATADIR,'/tables/Monkey_sva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)
gsea_df %>%saveRDS(paste0(DATADIR, "/rdas/Monkey_sva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.rds"))

## take a look at the enrichments
gsea_df2 = gsea_df %>% rbindlist()
a.tab = table(gsea_df2$MORPHINE.v.CTL.in, gsea_df2$MSigDb_Group)

pdf(paste0(PLOTDIR, "/plots/Monkey_sva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.pdf"), height=6, width=10)
gridExtra::grid.table(a.tab, rows = rownames(a.tab), cols = colnames(a.tab))
dev.off()




####################################################
## do the pathway enrichment from DEGs w/o the SVA 
DElist.fn = paste0(DATADIR,"/rdas/Monkey_voom", "_nosva_", "limma_bigModelSVA_N8.celltype.rds")

############################################################
# 1) read in the DEG lists per comparison of OUD vs. Control
res.celltype = DElist.fn %>% readRDS()
names(res.celltype) = paste0(names(res.celltype), '#All')

deg_rank_list1 = lapply(res.celltype, function(deg){
  deg %>% mutate(tmp = -log10(P.Value) * sign(logFC)) %>% 
    dplyr::select(gene, tmp) %>% arrange(tmp) %>% deframe()
})

deg_rank_list = deg_rank_list1

## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = pathways,
                   minSize=15, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

alpha = 0.05
gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% 
  arrange(pval) %>% filter(!is.na(pval)) %>% 
  mutate(
    MSigDb_Group = ss(pathway, '#', 1), 
    celltype = group %>% ss('#', 1), 
    MORPHINE.v.CTL.in = group %>% ss('#', 1),
    pathway = ss(pathway, '#', 2),
    padj = lm_qvalue(pval, X=size)$q, 
    celltype = group %>% ss('#', 1), 
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  inner_join(pathways_df2) %>% dplyr::select(-group) %>% 
  relocate(MORPHINE.v.CTL.in, celltype, MSigDb_Group, description, .before= everything()) %>% 
  split(f = .$MORPHINE.v.CTL.in)

## save the enrichment w/ all the pathways, significant and otherwise
gsea_df %>% saveRDS(paste0(DATADIR, "/rdas/Monkey_nosva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.unfiltered.rds"))
# gsea_df = readRDS(paste0(DATADIR, "rdas/Monkey_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.unfiltered.rds"))
## filter out just the significant pathways
gsea_df = lapply(gsea_df, function(x){
  x %>% filter(padj < alpha)
})
sapply(gsea_df, nrow)

out_fn = paste0(DATADIR,'/tables/Monkey_nosva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)
gsea_df %>% saveRDS(paste0(DATADIR, "/rdas/Monkey_nosva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.rds"))

## take a look at the enrichments
gsea_df2 = gsea_df %>% rbindlist()
a.tab = table(gsea_df2$MORPHINE.v.CTL.in, gsea_df2$MSigDb_Group)

pdf(paste0(PLOTDIR, "/plots/Monkey_nosva_GSEA_enrichment_in_human_msigdb_H_C2_C5_SynGO.pdf"), height=6, width=10)
gridExtra::grid.table(a.tab, rows = rownames(a.tab), cols = colnames(a.tab))
dev.off()


