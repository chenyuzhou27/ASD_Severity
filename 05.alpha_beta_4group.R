rm(list=ls())
library(tibble)
library(dplyr)
library(readxl)
library(tibble)
library(dplyr)
library(ggplot2)      
library(readxl)       
library(dplyr)        
library(tibble)      
library(Hmisc )
library(stringr)
library(MicrobiotaProcess)
library(rstatix)
library(ggpubr)
load("Manuscript.RData")
 
alpha.all.sev <- phy.mpse %>%
  mp_plot_alpha(
  .group=Cluster,
  .alpha=c(Observe, Shannon, Simpson),
  map_signif_level=T,
  
  ) +
  scale_fill_manual(values=c( "gray","#F9DB6D","#36827F","#464D77"),guide="none") +
  scale_color_manual(values=c("gray", "#F9DB6D","#36827F","#464D77"), guide="none")
alpha.all.sev

alph.index <- phy.mpse %>% mp_extract_sample() 
alph.index %>%
  wilcox_test(Observe ~ Cluster) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

alph.index %>%
  wilcox_test(Simpson ~ Cluster) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()

p.pca.all.sev <- phy.mpse %>% 
  mp_cal_pca(.abundance=hellinger, action="add") %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cluster,
    .color = Cluster,
    .size = 1,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE 
  ) +
  scale_fill_manual(values=c( "gray","#F9DB6D","#36827F","#464D77"),guide="none") +
  scale_color_manual(values=c( "gray","#F9DB6D","#36827F","#464D77"), guide="none")
dev.off()



###asd index
index <- read.xlsx2("index_asd_equal_distance.xlsx",sheetIndex = 1)
sampledf.asd.index <- merge(index %>% select(Subject.ID,combined_tau_1.0), sampledf,by.x="Subject.ID",by.y = "row.names",all.y = T) %>%
  rename("ASD_index" = "combined_tau_1.0") %>%
  mutate(ASD_index = as.numeric(ASD_index))

alph.index.asd.index <- merge(alph.index,sampledf.asd.index,by.x = "Sample",by.y = "Subject.ID",all.y = T)
alph.index.asd.index %>%
  select(Observe, Shannon, Simpson, ASD_index) %>%
  cor_test(method = "spearman") %>%
  filter(var1 == "ASD_index")

library(ggplot2)
library(dplyr)
library(broom)
library(tidyr)
library(purrr)

spearman_test <- cor.test(
  ~ Simpson + ASD_index, 
  data = alph.index.asd.index, 
  method = "spearman"
)

p <- ggplot(data = alph.index.asd.index) +
  geom_point(mapping = aes(x = Simpson, y = ASD_index),
    alpha = 0.6,
    size = 3,
    color = "steelblue",
    shape = 21,
    fill = "steelblue",
    stroke = 0.5
  ) +
  geom_smooth(mapping = aes(x = Simpson, y = ASD_index))

###permanova
perm.variable_list <- c("Batch","Cohort","Cluster", "Age", "Gender", "BMI", "Siblings", "T_SRS_SCI", "T_SRS_RRB", 
                        "M_SEQ_seeking", "M_SEQ_hypo", "M_SEQ_hyper", "CBCL_AP_T", "CBCL_Externalizing_T", 
                        "ASC_total", "Bristol_stool_chart", "Medication", "atopic_disease", "Carb", "Prot", "Fat","TotFib","ASD_index")

###
process_variable <- function(variable) {
  na_samples <- which(is.na(sampledf[[variable]]))
  
  if (length(na_samples) == 0) {
    sampledf_clean <- sampledf.asd.index
    phy.mpse_bray_clean <- phy.mpse_bray
  } else {
    sampledf_clean <- sampledf[-na_samples, ]
    phy.mpse_bray_matrix <- as.matrix(phy.mpse_bray)
    phy.mpse_bray_clean <- as.dist(phy.mpse_bray_matrix[-na_samples, -na_samples])
  }
  
  formula <- as.formula(paste("phy.mpse_bray_clean ~", variable))
  
  univ.adonis <- vegan::adonis2(formula, data = sampledf_clean)
  
  tmp <- univ.adonis[1, c( 3, 5)]
  rownames(tmp) <- variable
  return(tmp)
}

permanova.all <- purrr::map_dfr(c, process_variable)

permanova.table <- permanova.all %>% as.data.frame() %>%
  mutate(Variable = perm.variable_list) %>%
  dplyr::select(Variable, R2 = R2  , p.value ='Pr(>F)')

permanova.table <- permanova.table %>%
  arrange(desc(R2))

permanova.table <- permanova.table  %>% 
  dplyr::filter(Variable != "Batch") %>%
  mutate(Variable = case_when(
    Variable == "Prot" ~ "Protein",
    Variable == "TotFib" ~ "Fibre",
    Variable == "Gender" ~ "Sex",
    Variable == "Cohort" ~ "ASD",
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

###asd index
library(ggridges)

stat.test <- sampledf.asd.index %>%
  wilcox_test(ASD_index ~ Cluster) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test


asd.index.ridges <- ggplot(sampledf.asd.index, aes(x = ASD_index, y = Cluster, fill = Cluster)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = c("gray","#F9DB6D", "#36827F", "#464D77")) +
  labs(x = "ASD Index", y = "Cluster") +
  theme_classic()


pdf("Fig1d_2a_2b_2c.pdf",height = 10,width = 12)
ggarrange(asd.index.ridges,alpha.all.sev,p.pca.all.sev,per.R2,labels = "auto",ncol = 2,nrow = 2,align = "hv")
dev.off()

##########
###bray distance violine plot
phy.mpse_bray <- phyloseq::distance(phy.data, method = "bray")
phy.mpse %>%
  mp_cal_dist(
    .abundance = hellinger,
    distmethod = "bray"
  )
pdist.cluster <- phy.mpse %>%
  mp_plot_dist(
    .distmethod = bray,
    .group = Cluster,
    group.test = TRUE
  ) 

###cluster
tbl <- phy.mpse %>%
  mp_extract_dist(distmethod="bray", .group=Cluster)

bray.within.group.data <- tbl %>% filter(GroupsComparison %in% c("Control-vs-Control","Mildest-vs-Mildest","Moderate-vs-Moderate","Severe-vs-Severe")) 

summary_stats <- bray.within.group.data%>%
  group_by(GroupsComparison) %>%
  summarise(
    median_val = median(bray, na.rm = TRUE),
    mean_val = mean(bray, na.rm = TRUE),
    .groups = 'drop'
  )

cluster.violin <- bray.within.group.data %>%
  ggplot(., aes(x = GroupsComparison, y = bray, fill = GroupsComparison)) +
  geom_violin(alpha = 1) +
  geom_boxplot(width=0.1,fill = "white") +
  scale_fill_manual(values = c("gray","#F9DB6D", "#36827F", "#464D77")) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, fatten = 1, color = "black", lwd = 0.5) +
  geom_text(
    data = summary_stats,
    aes(x = GroupsComparison, y = median_val, label = sprintf("%.3f", median_val)),
    color = "black",
    size = 3,
    vjust = -0.5,  
    hjust = 0.5
  ) +
  labs(x = "Cluster", y = "Bray-Curtis Distance", title = "Bray-Curtis Distance by Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.35,1.2) 

cluster.violin <- cluster.violin + stat_compare_means(
  method = "wilcox.test",
  comparisons = combn(levels(factor(bray.within.group.data$GroupsComparison)), 2, simplify = FALSE),
  label = "p.signif",
  tip.length = 0.01,
  step.increase = 0.05,
  size = 3
)
####control vs other groups
plot_data <- tbl %>%
  filter(GroupsComparison %in% c(
    "Control-vs-Mildest", "Control-vs-Moderate", "Control-vs-Severe"
  ))

summary_stats <- plot_data %>%
  group_by(GroupsComparison) %>%
  summarise(
    median_val = median(bray, na.rm = TRUE),
    mean_val = mean(bray, na.rm = TRUE),
    .groups = 'drop'
  )


p <- ggplot(plot_data, aes(x = GroupsComparison, y = bray, fill = GroupsComparison)) +
  geom_violin(alpha = 1) +
  geom_boxplot(width = 0.1, fill = "white") +
  
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.3, fatten = 1, color = "black", lwd = 0.5) +
  
  geom_text(
    data = summary_stats,
    aes(x = GroupsComparison, y = median_val, label = sprintf("%.3f", median_val)),
    color = "black",
    size = 3,
    vjust = -0.5,  
    hjust = 0.5
  ) +
  
  scale_fill_manual(values = c(
    "Control-vs-Mildest" = "#F9DB6D",
    "Control-vs-Moderate" = "#36827F",
    "Control-vs-Severe" = "#464D77"
  )) +
  
  labs(
    x = "Group Comparisons", 
    y = "Bray-Curtis Distance", 
    title = "Bray-Curtis Distance by Group Comparisons"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.35,1.2) 

p_control <- p + stat_compare_means(
  method = "wilcox.test",
  comparisons = combn(levels(factor(plot_data$GroupsComparison)), 2, simplify = FALSE),
  label = "p.signif",
  tip.length = 0.01,
  step.increase = 0.05,
  size = 3
)

pdf("Fig2d_Fig2e.pdf", width = 10, height = 5)
ggarrange(cluster.violin,p_control,nrow = 1,ncol = 2,widths = c(1.2,1))
dev.off()

