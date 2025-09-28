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

read_and_filter <- function(file_path, metadata_value,pval_cutoff) {
  read.table(file_path, header = TRUE) %>%
    filter(metadata == metadata_value) %>%
    filter(pval < pval_cutoff) 
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

create_forest_plot <- function(data, title) {
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
    theme_minimal()
  
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
    theme_minimal()
  
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
    theme_minimal()
  
  return(list(taxonomy = p_taxonomy, module = p_module, pathway = p_pathway))
}

srs_variables <- c("T_SRS_total",  "T_SRS_RRB", "T_SRS_SCI")
cbcl_variables <- c( "CBCL_AP_T", "CBCL_Externalizing_T")
asc_asd_variables <- c("ASC_total")
seq_variables <- c( "M_SEQ_hypo", "M_SEQ_hyper","M_SEQ_seeking")


read_data <- function(variables, prefix) {
  taxa_without_control <- lapply(variables, function(var) read_and_filter(paste0("taxa/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  path_without_control <- lapply(variables, function(var) read_and_filter(paste0("pathway/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  module_without_control <- lapply(variables, function(var) read_and_filter(paste0("module/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  
  
  combine_data <- function(data_list, category) {
    do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
  }
  
  without_control <- combine_data(taxa_without_control, "Taxonomy") %>%
    bind_rows(combine_data(path_without_control, "Pathway")) %>%
    bind_rows(combine_data(module_without_control, "Module"))
  
  without_control <- process_data(without_control)
  
  without_control <- without_control %>%
    filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
    mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
           feature_rn = gsub("\\."," ",feature_rn)) %>%
    mutate(feature_rn = trimws(feature_rn))
  
  significant_markers_without_control <- count_significant_markers(without_control)
  
  cat(paste("Significant markers without control for", prefix, ":", significant_markers_without_control, "\n"))
  
  return(without_control)
}

plot_list <- list()

plot_list$srs <- create_forest_plot(read_data(srs_variables, "SRS"), "SRS")
plot_list$cbcl <- create_forest_plot(read_data(cbcl_variables, "CBCL"), "CBCL")
plot_list$asc_asd <- create_forest_plot(read_data(asc_asd_variables, "ASC-ASD"), "ASC-ASD")
plot_list$seq <- create_forest_plot(read_data(seq_variables, "SEQ"), "SEQ")

###srs
p_taxonomy <- plot_list[["srs"]]$taxonomy
p_module <- plot_list[["srs"]]$module
p_pathway <- plot_list[["srs"]]$pathway
  
pdf("Figure5.SRS_forest_plot.pdf",height = 7,width = 14)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(2.3, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.3),ncol = 2,common.legend = T,legend = "right")
dev.off()
  
###cbcl
p_taxonomy <- plot_list[["cbcl"]]$taxonomy
p_module <- plot_list[["cbcl"]]$module
p_pathway <- plot_list[["cbcl"]]$pathway

pdf("Figure4.cbcl_forest_plot.pdf",height = 11,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.6),ncol = 2,common.legend = T,legend = "right")
dev.off()

###seq
p_taxonomy <- plot_list[["seq"]]$taxonomy
p_module <- plot_list[["seq"]]$module
p_pathway <- plot_list[["seq"]]$pathway

pdf("Figure3.seq_forest_plot.pdf",height = 12,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = T,legend = "right")
dev.off()

###srs
p_taxonomy <- plot_list[["asc_asd"]]$taxonomy
p_module <- plot_list[["asc_asd"]]$module
p_pathway <- plot_list[["asc_asd"]]$pathway

pdf("Figure6.asc_asd_forest_plot.pdf",height = 10,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = T,legend = "right")
dev.off()


