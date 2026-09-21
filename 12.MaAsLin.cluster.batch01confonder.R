# 清除工作区
rm(list=ls())

# 加载必要的库
library(tibble)
library(dplyr)
library(xlsx)
library(ggplot2)
library(readxl)
library(Hmisc)
library(survminer)
library(ggforestplot)
library(stringr)

module.des <- read.table("module.descript.txt",sep="\t")

# read and filter
read_and_filter.pvalue <- function(file_path, metadata_value,pval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(pval < pval_cutoff)  %>%
    filter(feature != "MF0030" & feature != "MF0032")
}

read_and_filter.qvalue <- function(file_path, metadata_value,qval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(qval < qval_cutoff)  %>%
    filter(feature != "MF0030" & feature != "MF0032")
}

# reformat the name
process_data <- function(data) {
  # 定义需要替换为小写的单词列表
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
  data %>% filter(pval < 0.05) %>% nrow()
}


###q-value 0.2
create_forest_plot <- function(data, title) {
  data <- data %>% arrange(category, desc(coef))
  data$category <- factor(data$category, levels = c("Taxonomy", "Module", "Pathway"))
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
    "Taxonomy" = "#F8766D",  
    "Module" = "#619CFF",    
    "Pathway" = "#00BA38"    
  )
  sig_leve <- c(
    "qv0.1" = "black", 
    "qv0.2" = "gray"  
  )
  
  p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), xmin = xmin, xmax = xmax)) +
    geom_linerange(aes(color = category), linewidth = 3, 
                   position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(aes(fill = sig_level), size = 2.5, shape = 21, color = "gray", stroke = 0.3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = category_colors, name = "Category") +
    scale_fill_manual(values = sig_leve, name = "Significance") +
    scale_x_continuous(name = "Coefficient", limits = c(x_min - x_buffer, x_max + x_buffer)) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    theme(legend.position = "bottom") 
  
  return(p)
}


####cluster ord
var <- "Cluster"
class.ass <- c(".L")

plots <- list()
feature_table <- data.frame()
###severity with control
taxa_with_control <- read_and_filter.qvalue(paste0("taxa/Log_lm_", var, "/all_results.tsv"), var, 0.2)
path_with_control <- read_and_filter.qvalue(paste0("pathway/Log_lm_", var, "/all_results.tsv"), var, 0.2)
module_with_control <- read_and_filter.qvalue(paste0("module/Log_lm_", var, "/all_results.tsv"), var, 0.2)


combine_data <- function(data_list, category) {
    do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
  }
  
with_control <- combine_data(list(taxa_with_control), "Taxonomy") %>%
    bind_rows(combine_data(list(path_with_control), "Pathway")) %>%
    bind_rows(combine_data(list(module_with_control), "Module"))
  

with_control <- process_data(with_control)

##remove unknown sp
with_control <- with_control %>% filter(!grepl("GGB", feature)) %>% filter(Questionnaire == class.ass) %>%  
   mutate(sig_level = ifelse(qval < 0.1,"qv0.1","qv0.2"))

# coef 0.3
with_control.0.3 <- with_control %>%
  filter(abs(coef) >0.3)

####compared with ASD diagnosis biomarker
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
intersect(cohort_related.0.3$feature,with_control.0.3$feature)

create_forest_plot_node <- function(data, title) {
  data <- data %>% arrange(category, desc(coef))
  
  data$category <- factor(data$category, levels = c("Taxonomy", "Module", "Pathway"))
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
    "Taxonomy" = "#F8766D",  
    "Module" = "#619CFF",    
    "Pathway" = "#00BA38"    
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

cluster.qv.plots.add.asd <- create_forest_plot_node(with_control.0.3_node,"Withcontrol.cluster")
pdf("Figure3.cluster_ord_forest_plots.L.q.0.2.withcontrol_coef.0.3.withASD.pdf", width = 15, height = 13)
cluster.qv.plots.add.asd
dev.off()

