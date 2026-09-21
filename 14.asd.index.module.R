rm(list=ls())
###module
library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)
library(stringr)
library(ggplot2)

load("module.massiline.RData")
metadata <- readRDS("metadata.carb.rds")

asd.index.ridges <- ggplot(metadata, aes(x = combined_tau_1.0, y = Cluster, fill = Cluster)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = c("gray","#F9DB6D", "#36827F", "#464D77")) +
  labs(x = "ASD Index", y = "Cluster") +
  theme_classic()

pdf("Figure1d.pdf",height = 5,width = 6)
ggarrange(asd.index.ridges)
dev.off()

library(Maaslin2)
metadata.sp.tmp<- metadata %>%
  dplyr::select(sample_id, combined_tau_1.0,Cohort, Cluster, Age, Gender, BMI, Siblings, Medication, atopic_disease,  Prot.g.1000kcal, TotFib.g.1000kcal) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(combined_tau_1.0 = as.numeric(combined_tau_1.0),
         Gender = as.factor(Gender),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease))) 

###Maaslin2
abundance.module.tmp<-  module %>% remove_rownames() %>% dplyr::select(Sample,all_of(pick.module)) %>% 
  filter(Sample %in% rownames(metadata.sp.tmp)) %>%
  tibble::column_to_rownames(var="Sample") 

all(rownames(abundance.module.tmp) == rownames(metadata.sp.tmp))

fit_data <- Maaslin2(
    input_data     = abundance.module.tmp,        
    input_metadata = metadata.sp.tmp,           
    min_prevalence = 0.05,                     
    normalization  = "none",                   
    analysis_method = "LM",                    
    transform      = "LOG",                      
    output         = "module/Log_lm_ASD_index",                
    fixed_effects  = c("combined_tau_1.0","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal"),  
    reference      = c("Gender,Female","Medication,No","atopic_disease,No")  
  )



###pathway
pathway <- read.table("pathway.profile.adjust_bacth_cluster.txt",sep="\t")
pathway <- t(pathway) %>% as.data.frame() %>% rownames_to_column( var = "Subject.ID")

zero.num <-apply(pathway, 2, function(x){sum(x!=0)})
pick.pathway <- names(which(zero.num/dim(pathway)[1] > 0.05))

pid <- pmatch(metadata.index$sample_id,pathway$Subject.ID)
pathway <- pathway[pid ,]
all(pathway$Subject.ID==metadata.index$sample_id)

abundance.pathway.tmp<-  pathway %>% remove_rownames() %>% 
  dplyr::select(Subject.ID, all_of(pick.pathway)) %>%
  filter(Subject.ID %in% rownames(metadata.sp.tmp)) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "Subject.ID")
all(rownames(abundance.pathway.tmp) == rownames(metadata.sp.tmp))

fit_data <- Maaslin2(
    input_data     = abundance.pathway.tmp,          
    input_metadata = metadata.sp.tmp,           
    min_prevalence = 0.05,                      
    normalization  = "none",                    
    analysis_method = "LM",                     
    transform      = "LOG",                      
    output         = "pathway/Log_lm_ASD_index",                
    
    fixed_effects  = c("combined_tau_1.0","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal"),  
    reference      = c("Gender,Female","Medication,No","atopic_disease,No")  
  )

###species
load("Manuscript.RData")
abundance.sp.tmp<- sp.prof %>% 
  dplyr::select(Sample,all_of(pick.species)) %>% 
  filter(Sample %in% rownames(metadata.sp.tmp)) %>%
  tibble::column_to_rownames(var="Sample") 
all(rownames(abundance.sp.tmp) == rownames(metadata.sp.tmp))

fit_data <- Maaslin2(
    input_data     = abundance.sp.tmp,         
    input_metadata = metadata.sp.tmp,           
    min_prevalence = 0.05,                     
    normalization  = "none",                   
    analysis_method = "LM",                    
    transform      = "LOG",                      
    output         = "taxa/Log_lm_ASD_index",                 
    
    fixed_effects  = c("combined_tau_1.0","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal"),  
    reference      = c("Gender,Female","Medication,No","atopic_disease,No")  
  )
  

####merge asd index result for plot
#q-value
rm(list=ls())

read_and_filter.qvalue <- function(file_path, metadata_value,qval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(qval < qval_cutoff)  %>%
    filter(feature != "MF0030" & feature != "MF0032")
}

module.des <- read.table("../../../module.descript.txt",sep="\t")
###MF0030 and MGB050 is the same moduel, only keep 1
###MGB051	and MF0032 is the same module, only keep 1
module.des <- module.des %>% filter(Module != "MF0030" & Module != "MF0032")
head(module.des)

module_select<- read_and_filter.qvalue("module/Log_lm_ASD_index/all_results.tsv", "combined_tau_1.0", 0.2)
pathway_select<- read_and_filter.qvalue("pathway/Log_lm_ASD_index/all_results.tsv", "combined_tau_1.0", 0.2)
species_select<- read_and_filter.qvalue("taxa//Log_lm_ASD_index/all_results.tsv", "combined_tau_1.0", 0.2)

combine_data <- function(data_list, category) {
  do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
}

with_control <- combine_data(list(species_select), "Species") %>%
  bind_rows(combine_data(list(pathway_select), "Metabolic pathways")) %>%
  bind_rows(combine_data(list(module_select), "Gut-brain modules"))

process_data <- function(data) {
  lowercase_words <- c("To", "Of", "And", "Sp")
  
  data <- data %>% mutate(
    feature_rn = case_when(
      grepl("^s__", feature) ~ gsub("s__", "", feature) %>% gsub("_"," ",.),
      grepl("\\.\\.", feature) ~ sapply(strsplit(feature, "\\.\\."), function(x) paste(x[-1], collapse = "..")),
      TRUE ~ module.des$Description[match(feature, module.des$Module)]
    ),
    feature_rn = as.character(feature_rn),  
    feature_rn = str_replace_all(feature_rn, "\\b([a-z])", function(x) toupper(x)), 
    feature_rn = gsub("\\.\\.", " (", feature_rn),  
    feature_rn = gsub("\\(([^\\)]*)\\.", "(\\1)", feature_rn),  
    feature_rn = gsub("(\\d+)\\.(\\d+)", "\\1,\\2", feature_rn),  
    feature_rn = gsub("(\\d+)\\.([a-zA-Z])", "\\1-\\2", feature_rn),  
    feature_rn = gsub("\\.", " ", feature_rn),  
    feature_rn = gsub("([A-Z])\\s+([a-zA-Z])", "\\1-\\2", feature_rn)  
  )
  
  for (word in lowercase_words) {
    data <- data %>% mutate(
      feature_rn = gsub(paste0("\\b", word, "\\b"), tolower(word), feature_rn)
    )
  }
  
  data <- data %>% rename("Questionnaire" = value)
  
  ###manualy rename table
  label.match <- read.csv("cluster.marker.nodes.mathc_lable.csv")
  data <- data %>%
    mutate(feature_rn = if_else(
      feature %in% label.match$ID,
      label.match$Lable[match(feature, label.match$ID)],
      feature_rn
    ))
  
  data
}


count_significant_markers <- function(data) {
  data %>% filter(qval < 0.2) %>% nrow()
}

with_control <- process_data(with_control)
with_control <- with_control %>%
  filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
  mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
         feature_rn = gsub("\\."," ",feature_rn)) %>%
  mutate(sig_level = ifelse(qval < 0.1,"qv0.1","qv0.2"))


with_control.0.3 <- with_control %>% filter(abs(coef) > 0.3)
write.table(with_control.0.3,"Figure2b.asd_index_equal_distance_forest_plots.L.q.0.2_withcontrol.0.3.txt",sep="\t",quote = F)

###add diagnosis marker
taxa_cohort <- read_and_filter.qvalue(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
path_cohort  <- read_and_filter.qvalue(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
module_cohort <- read_and_filter.qvalue(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)

cohort_related <-   combine_data(list(taxa_cohort), "Species") %>%
  bind_rows(combine_data(list(path_cohort), "Metabolic pathways")) %>%
  bind_rows(combine_data(list(module_cohort), "Gut-brain modules"))

cohort_related <- process_data(cohort_related) 
cohort_related <- cohort_related %>%  filter(!grepl("GGB", feature)) 

write.table(cohort_related,"FigureS2.qv.txt",sep="\t",quote = F)

cohort_related.0.3 <- cohort_related %>% filter(abs(coef) > 0.3)
with_control$feature[!with_control$feature %in% cohort_related$feature]
intersect(cohort_related$feature,with_control$feature)
intersect(cohort_related.0.3$feature,with_control$feature)

create_forest_plot_node <- function(data, title) {
  data <- data %>% arrange(category, desc(coef))
  
  data$category <- factor(data$category, levels = c("Species", "Gut-brain modules", "Metabolic pathways"))
  
  data <- data %>% mutate(
    xmin = coef - 1.96 * stderr,
    xmax = coef + 1.96 * stderr
  )
  
  data <- data %>% mutate(
    q_value = ifelse(qval < 0.2, "q_sig", "q_notsig")
  )
  
  x_min <- min(data$xmin, na.rm = TRUE)
  x_max <- max(data$xmax, na.rm = TRUE)
  x_range <- x_max - x_min
  x_buffer <- 0.1 * x_range 
  
  category_colors <- c(
    "Species" = "#F8766D", 
    "Gut-brain modules" = "#619CFF",    
    "Metabolic pathways" = "#00BA38"    
  )
  
  node_shape_palette <- c(
    "also_in_ASD_diag_with_q_value0.2" = 24,  
    "only_in_Cluster" = 21     
  )
  
  sig_leve <- c(
    "qv0.1" = "black", 
    "qv0.2" = "gray"  
  )
  
  p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), 
                        xmin = xmin, xmax = xmax)) +
    geom_linerange(aes(color = category), linewidth = 3, 
                   position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(aes(fill = sig_level, shape = node_shape), size = 2.5,
               colour = "white", stroke = 0.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = category_colors, name = "Category") +
    scale_fill_manual(values = sig_leve, name = "Significance") +
    
    scale_shape_manual(values = node_shape_palette, name = "Feature Type",
                       labels = c("also_in_ASD_diag_with_q_value0.2" = "Also in ASD Diagnosis (q-value < 0.2)", 
                                  "only_in_Cluster" = "Only in Cluster")) +
    
    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, color = "black", size = 3)
      ),
      shape = guide_legend(
        override.aes = list(fill = "black", color = "black", size = 3)
      ),
      color = guide_legend(override.aes = list(linewidth = 3))
    ) +
    scale_x_continuous(name = "Coefficient", limits = c(x_min - x_buffer, x_max + x_buffer)) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    theme(
      legend.position = "left",
      legend.box.just = "center",
      legend.margin = margin(t = 5, b = 5),      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", size = 12),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    ggtitle(title)
  
  return(p)
}


with_control.0.3_node <- with_control.0.3 %>%
  mutate(node_shape = case_when(
    feature %in% cohort_related$feature ~ "also_in_ASD_diag_with_q_value0.2",
    TRUE ~ "only_in_Cluster"
  ))


cluster.qv.plots.add.asd <- create_forest_plot_node(with_control.0.3_node,"Withcontrol.asd.index")
pdf("Figure3b.asd_index_forest_plots.L.q0.2.coef0.3.withcontrol.withASD.pdf", width = 15, height = 7)
cluster.qv.plots.add.asd
dev.off()

