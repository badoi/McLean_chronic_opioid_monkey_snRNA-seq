## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

## AUCell object for DNA damage estimates
library(AUCell)

## statistical analyses
library(tidymodels)
library(lme4)
library(broom.mixed)
library(lmerTest)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/compare_dna_damage_welch'
PLOTDIR='figures/exploratory/compare_dna_damage_welch'
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F)

#########################################
# 1) load in the DNA damage estimates

## load in the glia only AUCell scores
save_meta_fn = here('data/tidy_data/Seurat_projects',
                    'OUD_Striatum_refined_all_SeuratObj_N16.meta.rds')
meta = readRDS(file=save_meta_fn) %>% rownames_to_column('cellBarcode') %>% 
  mutate(condition = case_when(condition=='control' ~ 'CTL', T ~'Morphine'))
table(meta$celltype3)
table(meta$condition)

save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_gliatype_N16.rds')
glia_AUC = readRDS(save_fn)

## calculate adaptive thresholds to call a "damaged" glial cell
AUCell_thr_fn = here(PLOTDIR, 'plots', 'AUCell_welch_dna_damage_scores_thresholds.glia.pdf')
pdf(AUCell_thr_fn, onefile = T)
glia_assignment <- AUCell_exploreThresholds(glia_AUC, plotHist=T, assign=T, nCores = 1) 
dev.off()

lapply(glia_assignment, '[[', 'aucThr')

glia_AUC_df = getAUC(glia_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .) %>% 
  filter(!grepl('^D', celltype3), celltype3!= 'Interneuron' ) %>% 
  mutate(celltype3 = droplevels(celltype3), 
         Stage2_class = cellBarcode %in% glia_assignment$Stage2$assignment)
table(glia_AUC_df$celltype3)
table(glia_AUC_df$Stage2_class)

## load in the neuron-only AUCell scores
save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_neurontype_N16.rds')
neuron_AUC = readRDS(save_fn)

## calculate adaptive thresholds to call a "damaged" neuronal cell
AUCell_thr_fn2 = here(PLOTDIR, 'plots', 'AUCell_welch_dna_damage_scores_thresholds.neuron.pdf')
pdf(AUCell_thr_fn2, onefile = T)
neuron_assignment <- AUCell_exploreThresholds(neuron_AUC, plotHist=T, assign=T, nCores = 1) 
dev.off()

lapply(neuron_assignment, '[[', 'aucThr')

neuron_AUC_df =getAUC(neuron_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .) %>% 
  mutate(celltype3 = droplevels(celltype3), 
         Stage1_class = cellBarcode %in% neuron_assignment$Stage1$assignment)
table(neuron_AUC_df$celltype3)
table(neuron_AUC_df$Stage1_class)

##############################################
# 2) average neuron scores per sample neurons
neuron_dam_per_sample = neuron_AUC_df %>% group_by(Monkey) %>% 
  mutate(Stage1= mean(Stage1), 
         numCell = n()) %>% ungroup() %>% 
  distinct(condition, Sex, Monkey, Age, Pair, numCell, .keep_all = T)

## stats
neuron_modBySample_s1 = lm(Stage1 ~ condition + Pair + numCell,
                    data = neuron_dam_per_sample) 
summary(neuron_modBySample_s1) 

#   lm(formula = Stage1 ~ condition + Pair + numCell, data = neuron_dam_per_sample)
# 
# Residuals:
#   1          2          3          4          5          6          7 
# -2.479e-04  9.383e-05  2.479e-04 -7.122e-05 -9.383e-05  8.288e-05 -8.288e-05 
# 8 
# 7.122e-05 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)        4.313e-02  9.453e-04  45.626  0.00048 ***
#   conditionmorphine  2.021e-03  4.898e-04   4.126  0.05404 .  
#   PairB             -2.2512e-03  3.133e-04  -8.654  0.01309 *  
#   PairC              2.368e-03  3.247e-04   7.291  0.01830 *  
#   PairD             -4.262e-03  3.856e-04 -11.054  0.00808 ** 
#   numCell           -6.609e-06  6.687e-07  -9.883  0.01008 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0002867 on 2 degrees of freedom
# Multiple R-squared:  0.9987,	Adjusted R-squared:  0.9956 
# F-statistic: 317.7 on 5 and 2 DF,  p-value: 0.003141

## plots
pdf(here(PLOTDIR, 'plots', 'McLean_Striatum.neuron.dnaDamStage1.perSample.pdf'), width = .85, height = 1)
ggplot(neuron_dam_per_sample, aes(x =condition, y = Stage1, fill = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = .5, size = 1, aes(x = condition, y = Stage1)) + 
  geom_line(aes(group = Pair), alpha = 0.5) +
  scale_fill_manual(values =c('white', '#C61B8C')) +
  facet_wrap(~'Neuron') + 
  theme_bw(base_size = 5) + ylim(c(.025, .05)) +
  ylab('NeuN+/gH2AX+\nDamage Score') +
  theme(legend.position = 'none', axis.title.x=element_blank()) 
dev.off()



####################################################
# 4) average glia DNA damage scores per sample glia
glia_dam_per_sample = glia_AUC_df %>% group_by(Monkey) %>% 
  mutate(Stage2 = mean(Stage2), 
         numCell = n()) %>% ungroup() %>% 
  distinct(condition, Sex, Monkey, Age, Pair, numCell, .keep_all = T)

glia_modBySample_s2 = lm(Stage2 ~ condition + Pair + numCell,
                    data = glia_dam_per_sample)
summary(glia_modBySample_s2) 

# lm(formula = Stage2 ~ condition + Pair + numCell, data = glia_dam_per_sample)
# 
# Residuals:
#   1         2         3         4         5         6         7         8 
# 0.003021 -0.004518 -0.003021 -0.009414  0.004518 -0.007918  0.007918  0.009414 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
#   (Intercept)        2.471e-01  4.242e-02   5.824   0.0282 *
#   conditionmorphine -8.726e-04  1.122e-02  -0.078   0.9451  
#   PairB             -1.218e-02  1.382e-02  -0.881   0.4711  
#   PairC             -5.479e-03  1.345e-02  -0.407   0.7233  
#   PairD              5.551e-03  1.546e-02   0.359   0.7539  
#   numCell           -2.964e-05  6.354e-05  -0.466   0.6868  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01345 on 2 degrees of freedom
# Multiple R-squared:  0.5672,	Adjusted R-squared:  -0.5148 
# F-statistic: 0.5242 on 5 and 2 DF,  p-value: 0.7577


pdf(here(PLOTDIR, 'plots', 'McLean_Striatum.glia.dnaDamStage2.perSample.pdf'), width = .85, height = 1)
ggplot(neuron_dam_per_sample, aes(x =condition, y = Stage2, fill = condition)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = .5, size = 1, aes(x = condition, y = Stage2)) + 
  geom_line(aes(group = Pair), alpha = 0.5) +
  scale_fill_manual(values =c('white', '#C61B8C')) +
  scale_shape_manual(values =c(21, 22)) +
  facet_wrap(~'Glia') + 
  theme_bw(base_size = 5) +
  ylab('NeuN-/gH2AX+\nDamage Score') +
  theme(legend.position = 'none', axis.title.x=element_blank()) 
dev.off()


######################################################################
# 5) compute neuron type average DNA damage ratios for neuronal damage

## stratified by cell type and sample
dam_per_neuronxSample = neuron_AUC_df %>% 
  arrange(celltype3) %>% 
  group_by(Monkey, celltype3) %>% 
  mutate_if(is.numeric, mean) %>% 
  mutate(numCell = n()) %>% 
  ungroup() %>% 
  distinct(condition, Sex, Monkey, Age, Pair, numCell, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'numCell')), ~ (. - mean(.))/sd(.)) %>% 
  group_by(celltype3) %>%   
  mutate(celltype3 = factor(celltype3, unique(celltype3)))

lm(Stage1 ~ celltype3 + celltype3:condition + Pair + numCell, 
     data = dam_per_neuronxSample) %>% summary()

## stage 1 in glia
modBySampleAndcellStage1 = lm(Stage1 ~ celltype3 + celltype3:condition + Pair + numCell, 
                          data = dam_per_neuronxSample) %>% tidy() %>% 
  filter(grepl('^celltype3', term) & grepl('morphine$', term)) %>% arrange(p.value) %>% 
  mutate(celltype = ss(term, '3|:', 2))

# lm(formula = Stage1 ~ celltype3 + celltype3:condition + Pair + 
#      numCell, data = dam_per_neuronxSample)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -2.550e-03 -1.169e-03  4.181e-05  1.001e-03  2.496e-03 
# 
# Coefficients:                            Estimate Std. Error t value Pr(>|t|)
# (Intercept)                             0.0311909  0.0011491  27.143  < 2e-16
# celltype3D1-NUDAP                       0.0057143  0.0011926   4.791 0.000111
# celltype3D2-Shell/OT                    0.0027073  0.0011802   2.294 0.032761
# celltype3Interneuron                    0.0060512  0.0013444   4.501 0.000218
# PairB                                  -0.0015707  0.0008440  -1.861 0.077524
# PairC                                   0.0031578  0.0008501   3.715 0.001369
# PairD                                  -0.0020076  0.0008851  -2.268 0.034533
# numCell                                -0.0001265  0.0004743  -0.267 0.792379
# celltype3D1-Shell/OT:conditionmorphine  0.0074769  0.0012978   5.761 1.23e-05
# celltype3D1-NUDAP:conditionmorphine     0.0052331  0.0012985   4.030 0.000655
# celltype3D2-Shell/OT:conditionmorphine  0.0061345  0.0014136   4.340 0.000318
# celltype3Interneuron:conditionmorphine  0.0068387  0.0011956   5.720 1.35e-05
# 
# (Intercept)                            ***
#   celltype3D1-NUDAP                      ***
#   celltype3D2-Shell/OT                   *  
#   celltype3Interneuron                   ***
#   PairB                                  .  
#   PairC                                  ** 
#   PairD                                  *  
#   numCell                                   
#   celltype3D1-Shell/OT:conditionmorphine ***
#   celltype3D1-NUDAP:conditionmorphine    ***
#   celltype3D2-Shell/OT:conditionmorphine ***
#   celltype3Interneuron:conditionmorphine ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001661 on 20 degrees of freedom
# Multiple R-squared:  0.9214,	Adjusted R-squared:  0.8781 
# F-statistic:  21.3 on 11 and 20 DF,  p-value: 1.143e-08


## plots
pdf(here(PLOTDIR, 'plots', 'McLean_Striatum_dnaDamValStage1.perSampleByCell.pdf'), width = 2.15, height = 1)
ggplot(dam_per_neuronxSample, aes(x =condition, y = Stage1, fill = condition)) +
  geom_violin(size = .25, bw = .002) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  geom_point(size = 1, alpha = .5)+
  geom_line(aes(group = Pair), alpha = 0.5) +
  scale_fill_manual(values =c('white', '#C61B8C')) +
  scale_color_manual(values =c('white', '#C61B8C')) +
  facet_wrap(~celltype3, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) + ylim(c(.025, .05)) +
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()



################################################################
# 6) compute cell type average DNA damage ratios for glial damage

## stratified by cell type and sample
dam_per_gliaxSample = glia_AUC_df %>% 
  arrange(celltype3) %>% 
  group_by(Monkey, celltype3) %>% 
  mutate_if(is.numeric, mean) %>% 
  mutate(numCell = n()) %>% 
  ungroup() %>% 
  distinct(condition, Sex, Monkey, Age, Pair, numCell, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'numCell')), ~ (. - mean(.))/sd(.)) %>% 
  group_by(celltype3) %>%   
  mutate(celltype3 = factor(celltype3, unique(celltype3)))

## take a look
dam_per_gliaxSample %>% count(Monkey)

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lm(Stage2 ~ celltype3 + celltype3:condition + Pair + numCell, 
   data = dam_per_gliaxSample) %>% summary()

# lm(formula = Stage2 ~ celltype3 + celltype3:condition + Pair + 
#      numCell, data = dam_per_gliaxSample)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.014984 -0.004387  0.000073  0.005410  0.012741 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                            0.2418975  0.0066193  36.544  < 2e-16
# celltype3Microglia                    -0.0466047  0.0075444  -6.177 4.92e-06
# celltype3Oligos                       -0.0015150  0.0068467  -0.221 0.827125
# celltype3Oligos_Pre                   -0.0392734  0.0090110  -4.358 0.000304
# PairB                                 -0.0048278  0.0048130  -1.003 0.327810
# PairC                                 -0.0002582  0.0047859  -0.054 0.957506
# PairD                                  0.0084593  0.0049414   1.712 0.102372
# numCell                               -0.0026974  0.0030822  -0.875 0.391875
# celltype3Astrocytes:conditionmorphine -0.0012833  0.0073283  -0.175 0.862750
# celltype3Microglia:conditionmorphine   0.0061090  0.0069316   0.881 0.388612
# celltype3Oligos:conditionmorphine     -0.0054663  0.0068753  -0.795 0.435910
# celltype3Oligos_Pre:conditionmorphine  0.0079684  0.0068092   1.170 0.255660
# 
# (Intercept)                           ***
#   celltype3Microglia                    ***
#   celltype3Oligos                          
#   celltype3Oligos_Pre                   ***
#   PairB                                    
#   PairC                                    
#   PairD                                    
#   numCell                                  
#   celltype3Astrocytes:conditionmorphine    
#   celltype3Microglia:conditionmorphine     
#   celltype3Oligos:conditionmorphine        
#   celltype3Oligos_Pre:conditionmorphine    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## stage 2 in glia
modBySampleAndcellStage2 = lm(Stage2 ~ celltype3 + celltype3:condition+ Age + Sex + 
                                PMI + RIN + numCell + Region + (1|Monkey), 
                              data = dam_per_gliaxSample) %>% tidy() %>% 
  filter(grepl('^celltype3', term) & grepl('morphine$', term)) %>% arrange(p.value) %>% 
  mutate(celltype = ss(term, '3|:', 2))


## plots
pdf(here(PLOTDIR, 'plots', 'McLean_Striatum_dnaDamValStage2.perSampleByCell.pdf'), width = 2.25, height = 1)
ggplot(dam_per_gliaxSample, aes(x =condition, y = Stage1, fill = condition)) +
  geom_violin(size = .25, bw = .002) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = .5, 
              position = position_jitterdodge(), alpha = .5)+
  scale_fill_manual(values =c('white', '#2A62AD')) +
  scale_color_manual(values =c('white', '#2A62AD')) +
  scale_shape_manual(values =c(21, 22)) +
  facet_wrap(~celltype3, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('NeuN-/gH2AX+\nDamage Score') +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), axis.title.x=element_blank(),
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()

#############################
# 7) export to spreadsheet
dfList = list(
  'Stage1_NeuN+_gH2AX+_neuron' = 
    modBySampleAndcellStage1 %>% as.data.frame() %>% relocate('celltype', .after = 'term'), 
  'Stage2_NeuN-_gH2AX+_glia' = 
    modBySampleAndcellStage2 %>% as.data.frame() %>% relocate('celltype', .after = 'term'))

stats_fn = here(PLOTDIR, 'tables', 'McLean_Striatum_dnaDam_statistics.perSampleByCell.xlsx')
writexl::write_xlsx(dfList, stats_fn)

