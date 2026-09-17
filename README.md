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


## R packages used

`MMUPHin`, `Maaslin2`, `vegan`, `rstatix`, `usedist`, `broom`, `phyloseq`, `MicrobiotaProcess`, `omixerRpm`, `clusterSim`, `dplyr`, `tibble`, `tidyr`, `purrr`, `stringr`, `Hmisc`, `readxl`, `xlsx`, `ggplot2`, `ggpubr`, `ggforestplot`, `ggridges`, `UpSetR`, `ComplexHeatmap`, `survminer`, `corrr`, `ggstatsplot`, `stats`, `utils`
