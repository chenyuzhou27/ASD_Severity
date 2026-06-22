rm(list=ls())

library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)
library(stringr)
library(ggplot2)

load("module.massiline.RData")

index <- read.xlsx2("index_asd_equal_distance.xlsx",sheetIndex = 1)
metadata.index <- merge(index %>% select(Subject.ID,starts_with("combined_tau")), metadata,by.x="Subject.ID",by.y = "sample_id") %>%
  rename("sample_id" = "Subject.ID")

module <- module[pmatch(metadata.index$sample_id,module$Sample),]

library(Maaslin2)
metadata.sp.tmp<- metadata.index %>%
  dplyr::select(sample_id,starts_with("combined_tau"), Age, Gender, BMI, Siblings, Medication, atopic_disease,  TotFib) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(across(starts_with("combined_tau"), as.numeric),
         Gender = as.factor(Gender),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease))) 

###Maaslin2
abundance.module.tmp<-  module %>% remove_rownames() %>% dplyr::select(Sample,all_of(pick.module)) %>% 
  filter(Sample %in% rownames(metadata.sp.tmp)) %>%
  tibble::column_to_rownames(var="Sample") 

tau_var <- grep("^combined_tau_1.0", colnames(metadata.sp.tmp), value = TRUE)

output_dir <- file.path("asd_value", tau_var, "module")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  fit_data <- Maaslin2(
    input_data     = abundance.module.tmp,        
    input_metadata = metadata.sp.tmp,           
    min_prevalence = 0.05,                     
    normalization  = "none",                    
    analysis_method = "LM",                     
    transform      = "LOG",                     
    output         = output_dir,                

    fixed_effects  = c(tau_var,"Age", "Medication", "atopic_disease", "TotFib"),
    reference      = c("Medication,No", "atopic_disease,No")  
  )

read_and_filter.qvalue <- function(file_path, metadata_value,qval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(qval < qval_cutoff)  
}

read_and_filter.pvalue <- function(file_path, metadata_value,pval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(pval < pval_cutoff)  
}

module_score.result <- c()
in_dir <- file.path("asd_value", tau_var, "module")
module_tmp<- read_and_filter.qvalue(paste0(in_dir,"/all_results.tsv"), tau_var, 0.2)
module_tmp <- module_tmp %>% mutate(score = tau_var)
module_score.result <- rbind.data.frame(module_score.result,module_tmp)

module.des <- read.table("module.descript.txt",sep="\t")
module_score.result.des <- merge(module.des,module_score.result,by.x = "Module",by.y= "feature")

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

output_dir <- file.path("asd_value", tau_var, "pathway")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
  fit_data <- Maaslin2(
    input_data     = abundance.pathway.tmp,          
    input_metadata = metadata.sp.tmp,           
    min_prevalence = 0.05,                     
    normalization  = "none",                   
    analysis_method = "LM",                     
    transform      = "LOG",                     
    output         = output_dir,                 
    
    fixed_effects  = c(tau_var,"Age", "Medication", "atopic_disease", "TotFib"), 
    reference      = c("Medication,No", "atopic_disease,No")  
  )

###species
load("Manuscript.RData")
abundance.sp.tmp<- sp.prof %>% 
  dplyr::select(Sample,all_of(pick.species)) %>% 
  filter(Sample %in% rownames(metadata.sp.tmp)) %>%
  tibble::column_to_rownames(var="Sample") 
all(rownames(abundance.sp.tmp) == rownames(metadata.sp.tmp))

dim(abundance.sp.tmp)
output_dir <- file.path("asd_value", tau_var, "species")
  if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
  fit_data <- Maaslin2(
    input_data     = abundance.sp.tmp,         
    input_metadata = metadata.sp.tmp,          
    min_prevalence = 0.05,                     
    normalization  = "none",                   
    analysis_method = "LM",                   
    transform      = "LOG",                     
    output         = output_dir,                 
    
    fixed_effects  = c(tau_var,"Age", "Medication", "atopic_disease", "TotFib"),  
    reference      = c("Medication,No", "atopic_disease,No")  
  )


####merge tau1 result for plot
#q-value
module_select<- read_and_filter.qvalue("asd_value/combined_tau_1.0/module/all_results.tsv", "combined_tau_1.0", 0.2)
pathway_select<- read_and_filter.qvalue("asd_value/combined_tau_1.0/pathway/all_results.tsv", "combined_tau_1.0", 0.2)
species_select<- read_and_filter.qvalue("asd_value/combined_tau_1.0/species/all_results.tsv", "combined_tau_1.0", 0.2)

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
with_control <- with_control %>% filter(!grepl("GGB", feature))
with_control <- with_control %>%
  filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
  mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
         feature_rn = gsub("\\."," ",feature_rn))

with_control <- with_control %>% mutate(consis = "consis") 
with_control.0.3 <- with_control %>% mutate(consis = "consis") %>% filter(abs(coef) > 0.3)

###add diagnosis marker
taxa_cohort <- read_and_filter.qvalue(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
path_cohort  <- read_and_filter.qvalue(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
module_cohort <- read_and_filter.qvalue(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)

cohort_related <-   combine_data(list(taxa_cohort), "Taxonomy") %>%
  bind_rows(combine_data(list(path_cohort), "Pathway")) %>%
  bind_rows(combine_data(list(module_cohort), "Module"))

cohort_related <- process_data(cohort_related) 
cohort_related <- cohort_related %>%  filter(!grepl("GGB", feature)) 

cohort_related.0.3 <- cohort_related %>% filter(abs(coef) > 0.3)
with_control$feature[!with_control$feature %in% cohort_related$feature]
intersect(cohort_related$feature,with_control$feature)
intersect(cohort_related.0.3$feature,with_control$feature)


taxa_cohort_pv <- read_and_filter.pvalue(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
path_cohort_pv  <- read_and_filter.pvalue(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
module_cohort_pv <- read_and_filter.pvalue(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)

cohort_pv.related <-   combine_data(list(taxa_cohort_pv), "Taxonomy") %>%
  bind_rows(combine_data(list(path_cohort_pv), "Pathway")) %>%
  bind_rows(combine_data(list(module_cohort_pv), "Module")) %>%
  process_data(.) %>% filter(!grepl("GGB", feature)) 

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
  
  # 为 category 创建颜色映射
  category_colors <- c(
    "Species" = "#F8766D", 
    "Gut-brain modules" = "#619CFF",    
    "Metabolic pathways" = "#00BA38"    
  )
  
  node_color_palette <- c(
    "also_in_ASD_diag_with_q_value0.2" = "black",  
    "also_in_ASD_diag_with_p_value0.05" = "darkgray",
    "only_in_Cluster" = "white"     
  )
  
  p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), 
                        xmin = xmin, xmax = xmax)) +
    geom_linerange(aes(color = category), linewidth = 3, 
                   position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(aes(fill = node_color), size = 2.5, shape = 21, 
               colour = "white", stroke = 0.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = category_colors, name = "Category") +
    scale_fill_manual(values = node_color_palette, name = "Feature Type",
                      labels = c("also_in_ASD_diag" = "Also in ASD Diagnosis", 
                                 "only_in_Cluster" = "Only in Cluster")) +
    scale_x_continuous(name = "Coefficient", limits = c(x_min - x_buffer, x_max + x_buffer)) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    theme(
      legend.position = "bottom",
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
  mutate(node_color = case_when(
    feature %in% cohort_related$feature ~ "also_in_ASD_diag_with_q_value0.2",
    feature %in% cohort_pv.related$feature ~ "also_in_ASD_diag_with_p_value0.05",
    TRUE ~ "only_in_Cluster"
  ))

cluster.qv.plots.add.asd <- create_forest_plot_node(with_control.0.3_node,"Withcontrol.asd.index")
pdf("Figure3b.pdf", width = 10, height = 7)
cluster.qv.plots.add.asd
dev.off()
saveRDS(with_control.0.3_node,"Figure.asd.index_forest_plots.L.p.0.05.withcontrol.withASD.rds")

