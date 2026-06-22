rm(list=ls())
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

read_and_filter.qvalue <- function(file_path, metadata_value,qval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(qval < qval_cutoff) 
}

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


count_significant_markers.qv <- function(data) {
  data %>% filter(qval < 0.2) %>% nrow()
}
count_significant_markers.pv <- function(data) {
  data %>% filter(pval < 0.05) %>% nrow()
}

severity_colors <- c(
  "Mildest"  = "#1B9E77",  
  "Moderate" = "#D95F02",  
  "Severe"   = "#7570B3"   
)

create_forest_plot.q0.2 <- function(data, title) {
  data <- data %>% arrange(category, desc(coef))
  data$category <- factor(data$category, levels = c("Taxonomy", "Module", "Pathway"))
  
  max_coef_data <- data %>%
    group_by(feature_rn) %>%
    summarise(max_coef = max(coef))
  
  data <- data %>%
    left_join(max_coef_data, by = "feature_rn")
  
  unique_features <- unique(data$feature_rn)
  levels.feat <- unique_features[order(max_coef_data$max_coef[match(unique_features, max_coef_data$feature_rn)], decreasing = TRUE)]
  data$feature_rn_sorted <- factor(data$feature_rn, levels = levels.feat)
  
  p_taxonomy <- ggforestplot::forestplot(
    name = feature_rn_sorted,
    df = data %>% filter(category == "Taxonomy"),
    se = stderr,
    estimate = coef,
    pvalue = qval,
    psignif = 0.2,
    xlab = "Coefficient",
    title = paste("a.", title, "- Taxonomy"),
    colour = Questionnaire
  ) +
    scale_color_manual(
      values = severity_colors, 
      breaks = c("Mildest", "Moderate", "Severe"),  
      name = "Severity"  
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12)
    )
  
  p_module <- ggforestplot::forestplot(
    name = feature_rn_sorted,
    df = data %>% filter(category == "Module"),
    se = stderr,
    estimate = coef,
    pvalue = qval,
    psignif = 0.2,
    xlab = "Coefficient",
    title = paste("b.", title, "- Module"),
    colour = Questionnaire
  ) +
    scale_color_manual(
      values = severity_colors,  
      breaks = c("Mildest", "Moderate", "Severe"),  
      name = "Severity"  
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12)
    )
  
  p_pathway <- ggforestplot::forestplot(
    name = feature_rn_sorted,
    df = data %>% filter(category == "Pathway"),
    se = stderr,
    estimate = coef,
    pvalue = qval,
    psignif = 0.2,
    xlab = "Coefficient",
    title = paste("c.", title, "- Pathway"),
    colour = Questionnaire
  ) +
    scale_color_manual(
      values = severity_colors,  
      breaks = c("Mildest", "Moderate", "Severe"),  
      name = "Severity"  
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12)
    )
  
  return(list(taxonomy = p_taxonomy, module = p_module, pathway = p_pathway))
}

subgroup <- c("Mildest",  "Moderate", "Severe")

read_data <- function(variables, prefix) {
  taxa<- lapply(variables, function(var) read_and_filter.qvalue(paste0(var,"/taxa/all_results.tsv"), "Cluster",0.2))
  path<- lapply(variables, function(var) read_and_filter.qvalue(paste0(var,"/pathway/all_results.tsv"),  "Cluster",0.2))
  module<- lapply(variables, function(var) read_and_filter.qvalue(paste0(var,"/module/all_results.tsv"),  "Cluster",0.2))
  
  
  combine_data <- function(data_list, category) {
    do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
  }
  
  compated_control <- combine_data(taxa, "Taxonomy") %>%
    bind_rows(combine_data(path, "Pathway")) %>%
    bind_rows(combine_data(module, "Module"))
  
  compated_control <- process_data(compated_control)
  
  compated_control <- compated_control %>%
    filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
    mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
           feature_rn = gsub("\\."," ",feature_rn)) %>%
    mutate(feature_rn = trimws(feature_rn))
  
  significant_markers_compated_control<- count_significant_markers.qv(compated_control)
  
  cat(paste("Significant markers compared control with subgroup:", significant_markers_compated_control, "\n"))
  
  return(compated_control)
}

# 处理每个问卷
subgroup.result <- create_forest_plot.q0.2(read_data(subgroup, "Subgroup"), "Subgroup")
layout_matrix <- matrix(c(1, 3, 2,3), ncol = 2, nrow = 2, byrow = TRUE)

p_taxonomy <- subgroup.result$taxonomy
p_module <- subgroup.result$module
p_pathway <- subgroup.result$pathway

pdf("subgroup_forest_plot.qv.pdf",height = 7,width = 14)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",heights = c(1.3, 1),nrow = 2,common.legend = F),
          p_pathway,widths = c(1 , 1.3),ncol = 2,common.legend = F,legend = "right")
dev.off()

