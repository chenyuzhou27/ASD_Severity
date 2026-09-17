# Dissecting Gut Microbiome Variations in Autism Spectrum Disorder Clinical Symptoms and Severity  

## Overview

Autism spectrum disorder (ASD) is a heterogeneous neurodevelopmental condition with individualised needs. Although evidence implicates the role of gut microbiome in ASD pathophysiology, prior studies have overlooked its clinical heterogeneity, leaving unclear whether microbial variations contribute to clinical severity and specific symptoms. Here, we performed shotgun metagenomic sequencing on stool samples from 720 children with ASD (median age 8 years, 87% male) and 151 matched children without ASD. Children with ASD were stratified into three severity groups based on 14 symptom subscales covering core ASD symptoms and co-occurring psychopathologies. Both overall severity and individual symptoms were significantly associated with gut microbiome composition (R2 = 0.0017 - 0.0065, p < 0.05). Children with more severe ASD exhibited enrichment of opportunistic pathogens, such as Citrobacter freundii and Enterococcus faecalis, and depleted anti-inflammatory Pseudoflavonifractor capillosus and GABA synthesis. Symptom-specific analyses revealed that sensory atypicalities had the largest number of associated microbial markers and specifically linked to excess excitatory neurotransmitters (q < 0.2). These findings demonstrate that gut microbiome composition varies across the ASD severity spectrum. Imbalance in microbial excitatory/inhibitory neurotransmitter metabolism emerges as a plausible mechanism underlying sensory atypicalities. Such insights may inform the development of targeted microbial therapeutics for ASD. 

**Data flow overview**

```
raw abundance tables (species / pathway / ko) + metadata
        │
        └─ 06.batch_effect_adjust.R  →  batch-adjusted abundance tables (*.adjust_bacth_cluster.txt)
                │
                ├─ 07.datacleaning.microbiomeprecess.R  →  MPSE object, alpha/beta diversity, species table
                │        │                                  (supFig1a, supFig1b, Manuscript.RData, metadata.carb.rds)
                │        └─ 07_Figure2.final.R  →  Fig.2 panels (alpha/beta diversity, PERMANOVA, distances)
                │
                ├─ 08.taxa.MaAsLin.batch01confounder.R     →  taxa/Log_lm_*      (species MaAsLin2)
                ├─ 09.pathway.batch01confounder.R           →  pathway/Log_lm_*   (pathway MaAsLin2)
                └─ 10.module.batch01confounder.R            →  module/Log_lm_*    (GMM/GBM module MaAsLin2)
                        │
                        ├─ 11.MaAsLin.cluster.batch01confonder.R  →  Fig.3 (severity-gradient forest plot)
                        └─ 12.MaAsLin.subscore.pdf.batch.01confounder.R  →  Fig.4-6 + SFigure6-9
                                │
                                └─ 13.asd.index.module.R  →  Fig.1d, Fig.2b, Fig.3b (ASD index analysis)
```

---

## 06.batch_effect_adjust.R

Uses `MMUPHin::adjust_batch` with `Batch` as the batch variable and `Cluster` as a covariate to correct batch effects in three abundance profiles (species, pathway, and KO gene families), normalizing each to relative abundance.

Output files:

- `species.profile.adjust_bacth_cluster.txt`
- `pathway.profile.adjust_bacth_cluster.txt`
- `ko.profile.adjust_bacth_cluster.txt`

## 07.datacleaning.microbiomeprecess.R

Cleans the data and builds the analysis object: merges the batch-adjusted species profile with a species tree (`sp.tree.rds`) and diet metadata (`metadata_diet_per_cal.xlsx`) into a phyloseq / MicrobiotaProcess (MPSE) object, then computes alpha diversity and PCA and derives the relative-abundance species table filtered at >5% prevalence.

Output files (with figure numbers):

- `supFig1a.alpha_asd_control.pdf` — Supplementary Figure 1a (alpha diversity, ASD vs Control)
- `supFig1b.pca.2groups.pdf` — Supplementary Figure 1b (PCA, ASD vs Control)
- `metadata.carb.rds` — cleaned sample metadata used by scripts 08-13
- `Manuscript.RData` — intermediate workspace loaded by scripts 08 and 13

## 07_Figure2.final.R

Draws the microbiome-overview panels: alpha diversity across severity clusters with Wilcoxon tests, PCA by cluster, PERMANOVA (R² bar plot with BH-adjusted q values), and within- and between-group Bray-Curtis distance comparisons.

Output files (with figure numbers):

- `alpha.div.test.result.txt` — Fig 2a/2b (alpha diversity Wilcoxon tests, BH-adjusted)
- `Fig2c.permonova.withcontrol.txt` — Fig 2c (PERMANOVA R² and q values)
- `Fig2d.distance.within.group.txt` — Fig 2d (within-group Bray-Curtis distance tests)
- `Fig2e.distance.between.group.txt` — Fig 2e (between-group Bray-Curtis distance tests)


## 08.taxa.MaAsLin.batch01confounder.R

Runs MaAsLin2 on the filtered species table, adjusting for the confounder set (Age, Gender, BMI, Medication, atopic_disease, protein and fibre intake per 1000 kcal) to test Cohort and Cluster in the full cohort, and then tests each questionnaire score in ASD samples only.

Output files (data for taxonomy panels of Fig.3 / Fig.4 / Fig.5 and downstream supplementary figures/tables):

- `taxa/Log_lm_Cohort/` — ASD vs Control
- `taxa/Log_lm_Cluster/` — severity gradient, full cohort
- `taxa/Log_lm_asd_<variable>/` — ASD only: T_SRS_total, T_SRS_SCI, T_SRS_RRB, CBCL_AP_T, CBCL_Internalizing_T, CBCL_Externalizing_T, ASC_total, M_SEQ_total, M_SEQ_hypo, M_SEQ_hyper, M_SEQ_seeking

Each output folder contains the standard MaAsLin2 files.

## 09.pathway.batch01confounder.R

Mirrors script 08 at the pathway level: runs MaAsLin2 on the batch-adjusted pathway profile with the same confounder set to test Cohort and Cluster in the full cohort and each questionnaire score in ASD samples only.

Output files (data for pathway panels of Fig.3 / Fig.4 / Fig.5 and downstream supplementary figures/tables):

- `pathway/Log_lm_Cohort/`
- `pathway/Log_lm_Cluster/`
- `pathway/Log_lm_asd_<variable>/` — same variable list as script 08

## 10.module.batch01confounder.R

Reconstructs GMM and GBM functional module abundances from the batch-adjusted KO profile with `omixerRpm` and runs MaAsLin2 on module abundance with the same confounder set against Cohort, Cluster, and each questionnaire score in ASD samples only.

Output files (data for module panels of Fig.3 / Fig.4 / Fig.5 and downstream supplementary figures/tables):

- `module/Log_lm_Cohort/`
- `module/Log_lm_Cluster/`
- `module/Log_lm_asd_<variable>/` — same variable list as script 08

## 11.MaAsLin.cluster.batch01confonder.R

Combines the MaAsLin2 results from scripts 08-10 (q < 0.2, |coef| > 0.3), and draws the severity-gradient forest plot for the ordinal Cluster linear term, marking features that are also associated with ASD diagnosis (Cohort) using different point shapes.

Output files (with figure numbers):

- `Figure3.cluster_ord_forest_plots.L.q.0.2.withcontrol_coef.0.3.withASD.pdf` — Fig 3a 

## 12.MaAsLin.subscore.pdf.batch.01confounder.R

Groups the ASD-only MaAsLin2 results from scripts 08-10 by questionnaire domain (SRS, CBCL, ASC-ASD, SEQ) and draws forest plots per domain at q < 0.2, plus UpSet plots of the overlap of significant features between domains and a parallel set of p < 0.05 supplementary figures and tables.

Output files (with figure numbers):

- `Figure4.SRS_forest_plot.qv.pdf` — Fig 4 (SRS)
- `Figure4.seq_forest_plot.qv.pdf` — Fig 4 (SEQ)
- `Figure5.cbcl_forest_plot.qv.pdf` — Fig 5 (CBCL)
- `Figure5.asc_asd_forest_plot.qv.pdf` — Fig 5 (ASC-ASD)
- `Figure6.UpSet_Pathway.pdf`, `Figure6.UpSet_Taxonomy.pdf`, `Figure6.UpSet_Module.pdf` — Fig 6 (feature overlap across domains)
- `SFigure6.SRS_forest_plot.pv.pdf`, `SFigure7.seq_forest_plot.pv.pdf`, `SFigure8.cbcl_forest_plot.pv.pdf`, `SFigure9.asc_asd_forest_plot.pv.pdf` — Supplementary Figures 5-9 at p < 0.05
- `symptom.qvalue.xlsx` — q < 0.2 marker table per domain
- `SFigure5_8.txt` — p < 0.05 marker table behind SFigure6-9

## 13.asd.index.module.R

Analyses the continuous ASD index : draws a ridge plot of the index across severity clusters and runs MaAsLin2 of modules, pathways, and species against the ASD index with the same confounder set, then merges the hits (q < 0.2, |coef| > 0.3), flags those that also appear among ASD-diagnosis (Cohort) markers, and draws the resulting forest plot.

Output files (with figure numbers):

- `Figure1d.pdf` — Fig 1d (ASD index distribution across severity clusters)
- `Figure2b.asd_index_equal_distance_forest_plots.L.q.0.2_withcontrol.0.3.txt` — Fig 2b (ASD index marker table)
- `FigureS2.qv.txt` — Supplementary S2 
- `Figure3b.asd_index_forest_plots.L.q0.2.coef0.3.withcontrol.withASD.pdf` — Fig 3b (ASD index forest plot)

---


## Environment & Dependencies

- **R Version:** 4.5.3 (2026-03-11)
- **Platform:** x86_64-redhat-linux-gnu
- **OS Platform:** AlmaLinux 9.7 (Moss Jungle Cat)

<details>
<summary><b>Click to expand full sessionInfo()</b></summary>

```R
R version 4.5.3 (2026-03-11)
Platform: x86_64-redhat-linux-gnu
Running under: AlmaLinux 9.7 (Moss Jungle Cat)

Matrix products: default
BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_HK.UTF-8          LC_NUMERIC=C                  LC_TIME=en_HK.UTF-8           LC_COLLATE=en_HK.UTF-8       
 [5] LC_MONETARY=en_HK.UTF-8       LC_MESSAGES=en_HK.UTF-8       LC_PAPER=en_HK.UTF-8          LC_NAME=en_HK.UTF-8          
 [9] LC_ADDRESS=en_HK.UTF-8        LC_TELEPHONE=en_HK.UTF-8      LC_MEASUREMENT=en_HK.UTF-8    LC_IDENTIFICATION=en_HK.UTF-8

time zone: Asia/Hong_Kong
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] omixerRpm_0.3.3          ggforestplot_0.1.0       survminer_0.5.2          Maaslin2_1.24.1          xlsx_0.6.5              
 [6] clusterSim_0.51-6        MASS_7.3-65              cluster_2.1.8.2          ggstatsplot_1.0.0        corrr_0.4.5             
[11] vegan_2.7-3              permute_0.9-10           phyloseq_1.54.2          purrr_1.2.2              tidyr_1.3.2             
[16] broom_1.0.12             ggpubr_0.6.3             rstatix_0.7.3            MicrobiotaProcess_1.22.1 stringr_1.6.0           
[21] Hmisc_5.2-5              ggplot2_4.0.3            readxl_1.4.5             dplyr_1.2.1              tibble_3.3.1
[26] MMUPHin_1.18.1   

loaded via a namespace (and not attached):
  [1] splines_4.5.3               ggplotify_0.1.3             cellranger_1.1.0            datawizard_1.3.1           
  [5] rpart_4.1.24                lifecycle_1.0.5             lattice_0.22-9              insight_1.5.0              
  [9] backports_1.5.1             magrittr_2.0.5              rmarkdown_2.31              otel_0.2.0                 
 [13] DBI_1.3.0                   RColorBrewer_1.1-3          ade4_1.7-24                 multcomp_1.4-30            
 [17] abind_1.4-8                 GenomicRanges_1.62.1        BiocGenerics_0.56.0         yulab.utils_0.2.4          
 [21] nnet_7.3-20                 TH.data_1.1-5               xlsxjars_0.9.0              rappdirs_0.3.4             
 [25] sandwich_3.1-1              gdtools_0.5.0               IRanges_2.44.0              S4Vectors_0.48.1           
 [29] ggrepel_0.9.8               correlation_0.8.8           tidytree_0.4.7              codetools_0.2-20           
 [33] coin_1.4-3                  DelayedArray_0.36.1         tidyselect_1.2.1            aplot_0.2.9                
 [37] farver_2.1.2                effectsize_1.0.2            matrixStats_1.5.0           stats4_4.5.3               
 [41] base64enc_0.1-6             Seqinfo_1.0.0               jsonlite_2.0.0              ggtreeExtra_1.20.1         
 [45] multtest_2.66.0             e1071_1.7-17                Formula_1.2-5               survival_3.8-6             
 [49] iterators_1.0.14            systemfonts_1.3.2           foreach_1.5.2               tools_4.5.3                
 [53] ggnewscale_0.5.2            treeio_1.34.0               ggstar_1.0.6                Rcpp_1.1.1-1.1             
 [57] glue_1.8.1                  gridExtra_2.3               SparseArray_1.10.10         xfun_0.57                  
 [61] mgcv_1.9-4                  MatrixGenerics_1.22.0       withr_3.0.2                 fastmap_1.2.0              
 [65] ggh4x_0.3.1                 digest_0.6.39               R6_2.6.1                    gridGraphics_0.5-1         
 [69] colorspace_2.1-2            utf8_1.2.6                  generics_0.1.4              fontLiberation_0.1.0       
 [73] data.table_1.18.4           robustbase_0.99-7           class_7.3-23                htmlwidgets_1.6.4          
 [77] S4Arrays_1.10.1             parameters_0.28.3           pkgconfig_2.0.3             rJava_1.0-18               
 [81] gtable_0.3.6                modeltools_0.2-24           statsExpressions_2.0.0      S7_0.2.2                   
 [85] XVector_0.50.0              pcaPP_2.0-5                 htmltools_0.5.9             fontBitstreamVera_0.1.1    
 [89] carData_3.0-6               biomformat_1.38.3           scales_1.4.0                Biobase_2.70.0             
 [93] optparse_1.8.2              ggfun_0.2.0                 knitr_1.51                  rstudioapi_0.18.0          
 [97] reshape2_1.4.5              checkmate_2.3.4             nlme_3.1-168                proxy_0.4-29               
[101] zoo_1.8-15                  parallel_4.5.3              libcoin_1.0-12              foreign_0.8-91             
[105] pillar_1.11.1               grid_4.5.3                  vctrs_0.7.3                 car_3.1-5                  
[109] htmlTable_2.5.0             paletteer_1.7.0             evaluate_1.0.5              mvtnorm_1.3-7              
[113] cli_3.6.6                   compiler_4.5.3              rlang_1.2.0                 crayon_1.5.3               
[117] rstantools_2.6.0            ggsignif_0.6.4              labeling_0.4.3              rematch2_2.1.2             
[121] plyr_1.8.9                  fs_2.1.0                    ggiraph_0.9.6               stringi_1.8.7              
[125] Biostrings_2.78.0           lazyeval_0.2.3              bayestestR_0.17.0           fontquiver_0.2.1           
[129] Matrix_1.7-4                patchwork_1.3.2             SummarizedExperiment_1.40.0 igraph_2.3.0               
[133] RcppParallel_5.1.11-2       biglm_0.9-3                 ggtree_4.0.5                DEoptimR_1.1-4             
[137] ape_5.8-1
