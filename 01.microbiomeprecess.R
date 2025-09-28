rm(list=ls())

library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(ggplot2)
library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(gghalves)
library(corrr)
library(ggstatsplot)
library(ggpubr)
library(Hmisc)
library(stringr)

species.batch.cluster <- read.table("species.profile.adjust_bacth_cluster.txt")
pathway <- read_xlsx("function_aug28_filter_norm.xlsx",sheet = 1)
species <- species.batch.cluster %>%  t() %>% as.data.frame() %>%  rownames_to_column(var = "Subject.ID")

metadata <- species.batch.cluster[,1:115]

sp.tree <-species.batch.cluster %>% dplyr::select("clade_name", "clade_taxid") %>%
  filter(grepl("s__", clade_name) & !grepl("t__", clade_name))%>% 
  mutate(sp.name = strsplit(clade_name,"\\|") %>% sapply(., "[[", 7) %>% gsub("s__","",.))

species.trans <- species.batch.cluster %>%
  rownames_to_column(var = "sp.name") %>%
  mutate(sp.name = gsub("^s__","",sp.name)) %>%
  merge(sp.tree,.,by = "sp.name",all.x = T) %>%
  dplyr::select(-sp.name)

sample.id <- colnames(species.trans)[-c(1:2)]

variables_list <- c(
  "T_SRS_total","T_SRS_AWR", "T_SRS_COG", "T_SRS_COMM", "T_SRS_MOT", "T_SRS_RRB","T_SRS_SCI",
  "CBCL_Total_T","CBCL_AP_T", "CBCL_Internalizing_T", "CBCL_Externalizing_T",
  "ASC_total","ASC_PA", "ASC_AA", "ASC_SA", "ASC_uncertainty",
  "M_SEQ_total","M_SEQ_hypo", "M_SEQ_hyper", "M_SEQ_seeking",
  "FNVD", "FAPD", "FDD", "Bristol_stool_chart"
)

diet.list <- c("Carb" , "Prot" , "Fat" ,"TotFib" )
medication.list <- c("Medication","ADHD.med.1","ADHD.med.2", "antipsychotic.med")

diet <- read.xlsx2("diet_metadata.xlsx",sheetIndex = 1)
diet.filter <- diet %>% 
  filter(ID %in% metadata.select$sample_id) 

metadata.select.updata.diet <- merge(metadata.select %>% dplyr::select(-all_of(diet.list)),
                                     diet.filter %>% dplyr::select(ID,all_of(diet.list)),
                                     by.x = "sample_id",by.y="ID",all.x = T)

##### MicrobiotaProcess
source("perm_univari.R")
sp.otu <- species.trans %>% dplyr::select(-c(clade_name,clade_taxid)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(OTU = paste("OTU",OTU,sep="")) %>%
  column_to_rownames(var = "OTU") %>%
  mutate_all(~ as.numeric(.))
sp.tax.info <- species.trans %>%
  dplyr::select(clade_name) %>%
  lapply(function(x) strsplit(as.character(x), split = "\\|")) 
sp.tax <- do.call(rbind,sp.tax.info$clade_name) %>% as.data.frame() %>%
  setNames(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  mutate_all(~ gsub(".*__", "", .)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(OTU = paste("OTU",OTU,sep="")) %>%
  column_to_rownames(var = "OTU") %>%
  as.matrix()

metadata.select.updata.diet <- metadata.select.updata.diet %>% remove_rownames() %>%
  column_to_rownames("sample_id") 

OTU = otu_table(sp.otu, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(sp.tax)
samples = sample_data(metadata.select.updata.diet)

phy.data <- phyloseq(OTU, TAX, samples)
phy.data
phy.mpse <- phy.data %>% as.MPSE() 

###alpha plot
phy.mpse %<>%
  mp_cal_alpha(.abundance=Abundance,force=TRUE )


mycomp <- list(c("ASD", "Control"))
alpha_2g <- phy.mpse %>%
  mp_plot_alpha(
    .group=Cohort,
    .alpha=c(Observe, Shannon, Simpson),
    #comparisons=mycomp,map_signif_level=T
  ) +
  scale_fill_manual(values=c( "#467897","#E7CD79"),guide="none") +
  scale_color_manual(values=c("#467897","#E7CD79"), guide="none")

pdf("supFig1a.alpha_asd_control.pdf",width = 5,height = 3)
alpha_2g
dev.off()


alpha_match.sev <- phy.mpse %>% filter(Cohort != "Control") %>%
  mp_plot_alpha(
    .group=Cluster,
    .alpha=c(Observe, Shannon, Simpson),
    comparisons=list(c("Mildest", "Moderate"),c("Mildest", "Severe"),c("Moderate", "Severe")),
    map_signif_level=T
  ) +
  scale_fill_manual(values=c( "#F9DB6D","#36827F","#464D77"),guide="none") +
  scale_color_manual(values=c( "#F9DB6D","#36827F","#464D77"), guide="none")

pdf("Fig.1d.sever_sample_alpha_div_asd.pdf",width = 5,height = 4)
ggarrange(alpha_match.sev)
dev.off()

###beta diversity
phy.mpse %<>%
  mp_decostand(.abundance=Abundance)
###PCA
p.pca.sev <- phy.mpse %>% filter(Cohort != "Control") %>%
  mp_cal_pca(.abundance=hellinger, action="add") %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cluster,
    .color = Cluster,
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("#F9DB6D","#36827F","#464D77"),guide="none") +
  scale_color_manual(values=c("#F9DB6D","#36827F","#464D77"), guide="none")
dev.off()
pdf("Fig.1e.match_sample.pca.3groups.pdf",height = 5,width = 5)
ggarrange(p.pca.sev)
dev.off()

p.pca <- phy.mpse %>%
  mp_cal_pca(.abundance=hellinger, action="add") %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cohort,
    .color = Cohort,
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("#467897","#E7CD79"),guide="none") +
  scale_color_manual(values=c("#467897","#E7CD79"), guide="none")


pdf("supFig1b.pca.2groups.pdf",height = 5,width = 5)
ggarrange(p.pca)
dev.off()

##Adonise
#univariable
phy.mpse_cohort <- phy.mpse %>%
  mp_adonis(.abundance=hellinger, .formula=~Cohort, distmethod="bray", permutations=999, action="add")
adonis.asd <- phy.mpse_cohort %>% mp_extract_internal_attr(name="adonis")
adonis.asd


phy.mpse_cluster <- phy.mpse  %>%
  mp_adonis(.abundance=hellinger, .formula=~Cluster , distmethod="bray", permutations=999, action="add")
adonis.asd.cluster <- phy.mpse_cluster %>% mp_extract_internal_attr(name=adonis)
adonis.asd.cluster

###univariable permanova (with control sample)
phy.mpse_bray <- phyloseq::distance(phy.data, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phy.data))
na_counts <- sampledf %>%
  summarise(across(everything(), ~ sum(is.na(.))))
perm.variable_list <- c("Cohort","edu","Cluster", "Age", "Gender", "BMI", "Siblings", "T_SRS_SCI", "T_SRS_RRB", 
                        "M_SEQ_seeking", "M_SEQ_hypo", "M_SEQ_hyper", "CBCL_AP_T", "CBCL_Externalizing_T", 
                        "ASC_total", "Bristol_stool_chart", "Medication", "atopic_disease", "Carb", "Prot", "Fat","TotFib")

permanova.all <- purrr::map_dfr(perm.variable_list, process_variable)

permanova.table <- permanova.all %>% as.data.frame() %>%
  mutate(Variable = perm.variable_list) %>%
  dplyr::select(Variable, R2 = R2  , p.value ='Pr(>F)')

permanova.table <- permanova.table %>%
  arrange(desc(R2))

permanova.table <- permanova.table  %>% 
  mutate(Variable = case_when(
  Variable == "Prot" ~ "Protein",
  Variable == "TotFib" ~ "Fibre",
  Variable == "Gender" ~ "Sex",
  Variable == "atopic_disease" ~ "Atopic Diseases",
  Variable == "Bristol_stool_chart" ~ "Bristol Stool Form",
  Variable == "ASC_total" ~ "ASC-ASD Total",
  Variable == "T_SRS_SCI" ~ "SRS-2 SCI",
  Variable == "T_SRS_RRB" ~ "SRS-2 RRB",
  Variable == "M_SEQ_hypo" ~ "SEQ Hypo",
  Variable == "M_SEQ_hyper" ~ "SEQ Hyper",
  Variable == "M_SEQ_seeking" ~ "SEQ Seeking",
  Variable == "CBCL_AP_T" ~ "CBCL AP",
  Variable == "CBCL_Externalizing_T" ~ "CBCL Ext",
  Variable == "Cluster" ~ "ASD Severity",
  TRUE ~ Variable
))

per.R2<-ggplot(permanova.table, aes(x = reorder(Variable, -R2), y = R2, fill = p.value < 0.05)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "grey")) +
  labs(x = "Variable", y = "R2", fill = "p < 0.05") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
write.table(permanova.table,"Fig1f.permonova.withcontrol.txt",sep="\t")
pdf("Fig.1f.permonova.withcontorl.pdf",width = 7,height = 5)
per.R2
dev.off()


###bray distance 
phy.mpse %<>%
  mp_cal_dist(
    .abundance = hellinger,
    distmethod = "bray"
  )

tbl <- phy.mpse %>%
  mp_extract_dist(distmethod="bray", .group=Cohort)

boxplot.braydis <- tbl %>% 
  ggplot(aes(x=GroupsComparison, y=bray)) + 
  geom_boxplot(aes(fill=GroupsComparison)) + 
  #geom_jitter(width=0.1) + 
  xlab(NULL) +
  theme(legend.position="none") +
  theme_classic2() +
  ggsignif::geom_signif(comparisons = list(c("ASD-vs-ASD", "ASD-vs-Control"),c("ASD-vs-ASD", "Control-vs-Control"),c("ASD-vs-Control", "Control-vs-Control")),
                        map_signif_level = TRUE,y_position =c(1,1.1,1.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

medians <- tbl %>%
  group_by(GroupsComparison) %>%
  summarise(Median = median(bray, na.rm = TRUE))

boxplot.braydis <- boxplot.braydis + 
  geom_text(data = medians, aes(x = GroupsComparison, y = Median, label = round(Median, 3)), vjust = -0.5) + 
  ylab("Bray−Curtis Distance")

pdf("SupFig2.bray.distance.asd.control.pdf",height = 6,width = 5)
boxplot.braydis
dev.off()

bray_dist_df <- as.data.frame(as.matrix(phy.mpse_bray))

bray_dist_df$Sample <- rownames(bray_dist_df)
bray_dist_df <- bray_dist_df %>%
  gather(key = "Sample2", value = "Distance", -Sample) %>%
  left_join(sampledf %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Cluster), by = c("Sample" = "SampleID")) %>%
  left_join(sampledf %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Cluster), by = c("Sample2" = "SampleID"))

bray_dist_df <- bray_dist_df %>%
  rename(Cluster1 = Cluster.x, Cluster2 = Cluster.y)

bray_dist_df.within.cluster <- bray_dist_df %>%
  filter(Sample != Sample2) %>%
  filter(Sample < Sample2) %>%
  filter(Cluster1 ==Cluster2) %>%
  mutate(Cohort = ifelse(Cluster1 == "Control", "Control", "ASD"))

p.cluster <- ggplot(bray_dist_df.within.cluster, aes(x = Cluster1, y = Distance, fill = Cluster1)) +
  geom_violin(alpha = 1) +
  geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("gray","#F9DB6D", "#36827F", "#464D77")) +
  labs(x = "Cluster", y = "Bray-Curtis Distance", title = "Bray-Curtis Distance by Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.cohort <-  ggplot(bray_dist_df.within.cluster, aes(x = Cohort, y = Distance, fill = Cohort)) +
  geom_violin(alpha = 1) +
  geom_boxplot(width=0.1) +
  scale_fill_manual(values = c("#467897","#E7CD79")) +
  labs(x = "Cluster", y = "Bray-Curtis Distance", title = "Bray-Curtis Distance by Cohort") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

means <- bray_dist_df.within.cluster %>%
  group_by(Cluster2) %>%
  summarise(Mean = mean(Distance, na.rm = TRUE),
            Count = sum(!is.na(Distance)))

medians <- bray_dist_df.within.cluster %>%
  group_by(Cluster2) %>%
  summarise(Median = median(Distance, na.rm = TRUE),
            Count = sum(!is.na(Distance)))

p<- ggplot(bray_dist_df.within.cluster, aes(x = Cluster1, y = Distance, fill = Cluster1)) +
  geom_boxplot() +
  scale_fill_manual(values = c("gray","#F9DB6D", "#36827F", "#464D77")) +
  labs(x = "Cluster", y = "Bray-Curtis Distance", title = "Bray-Curtis Distance by Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bray.dis.sev <-  p + stat_compare_means(
   comparisons = list(c("Control", "Mildest"),c("Control", "Moderate"),c("Control", "Severe"),
                     c("Mildest", "Moderate"),c("Mildest", "Severe"),c("Moderate", "Severe")),
   method = "wilcox.test",label = "p.signif") + 
   stat_compare_means(label.y = 1.5) 


pdf("Figure1.pdf",height = 10,width = 10)
ggarrange(alpha_match.sev,p.pca.sev,per.R2,bray.dis.sev,labels = "auto",heights = c(2,1.5))
dev.off()

pdf("FigureS2_bary_betweencluster.pdf",height = 8,width = 12)
ggarrange(boxplot.braydis,bray_between.cluster.control,bray_between.cluster.sever,legend = "none",align = "hv",
          labels = "auto",ncol = 3,nrow = 1)
dev.off()

###taxonomy
sp.prof.long <- phy.mpse %>%
  mp_extract_abundance(taxa.class=Species) %>%
  tidyr::unnest(cols=AbundanceBySample) %>% dplyr::rename(species="label")
sp.prof <- sp.prof.long[,-c(2,4)] %>% pivot_wider(names_from = species, values_from = RelAbundanceBySample)

zero.num <-apply(sp.prof[,-c(1:114)], 2, function(x){sum(x!=0)})
pick.species <- names(which(zero.num/dim(sp.prof)[1] > 0.05))

save.image("Manuscript.RData")

