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

barCOLS <- c(
  "consis_notsig" = "#f9b282",
  "consis_sig" = "#f9b282",
  "noconsis_notsig" = "#a6d8f0",
  "noconsis_sig" = "#a6d8f0"
)
dotCOLS <- c(
  "consis_notsig" = "#de6b35",
  "consis_sig" = "black",
  "noconsis_notsig" = "#008fd5",
  "noconsis_sig" = "black"
)


create_forest_plot <- function(data, title) {
  data <- data %>% arrange(category, desc(coef))
  data$category <- factor(data$category, levels = c("Taxonomy", "Module", "Pathway"))
  
  data <- data %>% mutate(
    xmin = coef - 1.96 * stderr,
    xmax = coef + 1.96 * stderr
  )
  
  data <- data %>% mutate(
    consis = ifelse(qval < 0.2, paste(consis, "_sig",sep = ""), paste(consis, "_notsig",sep = ""))
  )
  
  x_min <- min(data$xmin, na.rm = TRUE)
  x_max <- max(data$xmax, na.rm = TRUE)
  x_range <- x_max - x_min
  x_buffer <- 0.1 * x_range  
  

  p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), xmin = xmin, xmax = xmax, col = consis, fill = consis)) +
    geom_linerange(linewidth = 3, position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, lty = 2, color = "black") +
    geom_point(size = 2.5, shape = 21, colour = "white", stroke = 0.5, position = position_dodge(width = 0.5)) +
    scale_fill_manual(breaks = names(dotCOLS), values = dotCOLS, labels = c("Present in ASD/Control & q-value > 0.2", "Present in ASD/Control & q-value < 0.2","Not Present in ASD/Control & q-value > 0.2", "Not Present in ASD/Control & q-value < 0.2")) +
    scale_color_manual(breaks = names(barCOLS), values = barCOLS, labels = c("Present in ASD/Control & q-value > 0.2", "Present in ASD/Control & q-value < 0.2","Not Present in ASD/Control & q-value > 0.2", "Not Present in ASD/Control & q-value < 0.2"))  +
    scale_x_continuous(name = "Coefficient", limits = c(x_min - x_buffer, x_max + x_buffer)) +
    scale_y_discrete(name = "") +
    theme_minimal() +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    theme(legend.position = "bottom") +
    ggtitle(title)
  
  return(p)
}


####cluster 
total_variables <- c("Cluster")
class.ass <- c(".L")

plots <- list()
feature_table <- data.frame()
for (var in total_variables) {
  taxa_without_control <- read_and_filter(paste0("taxa/Log_lm_asd_", var, "/all_results.tsv"), var, 0.05)
  path_without_control <- read_and_filter(paste0("pathway/Log_lm_asd_", var, "/all_results.tsv"), var, 0.05)
  module_without_control <- read_and_filter(paste0("module/Log_lm_asd_", var, "/all_results.tsv"), var, 0.05)
 
  taxa_with_control <- read_and_filter(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
  path_with_control <- read_and_filter(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
  module_with_control <- read_and_filter(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
  
  combine_data <- function(data_list, category) {
    do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
  }
  
  without_control <- combine_data(list(taxa_without_control), "Taxonomy") %>%
    bind_rows(combine_data(list(path_without_control), "Pathway")) %>%
    bind_rows(combine_data(list(module_without_control), "Module"))
  
  with_control <- combine_data(list(taxa_with_control), "Taxonomy") %>%
    bind_rows(combine_data(list(path_with_control), "Pathway")) %>%
    bind_rows(combine_data(list(module_with_control), "Module"))
  
  without_control <- process_data(without_control)
  with_control <- process_data(with_control)
  
  without_control <- without_control %>% filter(!grepl("GGB", feature)) %>% filter(Questionnaire == class.ass)
  with_control <- with_control %>% filter(!grepl("GGB", feature))
  
  without_control <- without_control %>%
    mutate(consis = ifelse(feature %in% with_control$feature &
                             coef * with_control$coef[match(feature, with_control$feature)] > 0,
                           "consis", "noconsis")) %>%
    filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
    mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
           feature_rn = gsub("\\."," ",feature_rn))
  
  without_control.conflit <- without_control %>%
    mutate(consis = ifelse(feature %in% with_control$feature &
                             coef * with_control$coef[match(feature, with_control$feature)] < 0,
                           "confilit", "noconsis")) %>%
    filter(!grepl("GGB", feature) & !grepl("un_f__",feature)) %>%
    mutate(feature_rn = gsub("\\.\\.","-",feature_rn),
           feature_rn = gsub("\\."," ",feature_rn))
  
  significant_markers_without_control <- count_significant_markers(without_control)
  
  cat(paste("Significant markers without control for", var, ":", significant_markers_without_control, "\n"))
  
  plots[[var]] <- create_forest_plot(without_control, var)
  
  feature_table <- rbind(feature_table, data.frame(Variable = var, Feature = without_control$feature,coef = without_control$coef))
}

table(without_control$category)
table(without_control %>% filter(qval<0.2) %>% select(category))

pdf("Figure2.cluster_ord_forest_plots.L.pdf", width = 10, height = 13)
plots$Cluster
dev.off()
write.table(feature_table,"cluster.feature_ord.txt",sep="\t",quote = F )

###control vs asd
var = "Cohort"
  taxa_with_control <- read_and_filter(paste0("taxa/Log_lm_", var, "/all_results.tsv"), var, 0.05)
  path_with_control <- read_and_filter(paste0("pathway/Log_lm_", var, "/all_results.tsv"), var, 0.05)
  module_with_control <- read_and_filter(paste0("module/Log_lm_", var, "/all_results.tsv"), var, 0.05)
  
  combine_data <- function(data_list, category) {
    do.call(rbind, lapply(data_list, function(df) df %>% mutate(category = category)))
  }
  
  with_control <- combine_data(list(taxa_with_control), "Taxonomy") %>%
    bind_rows(combine_data(list(path_with_control), "Pathway")) %>%
    bind_rows(combine_data(list(module_with_control), "Module"))
  
  with_control <- process_data(with_control)
  
  with_control <- with_control %>% filter(!grepl("GGB", feature))
  
  barCOLS <- c(
    "q_sig" = "#f9b282",
    "q_notsig" = "#a6d8f0"
  )
  dotCOLS <- c(
    "q_sig" = "#de6b35",
    "q_notsig" = "#008fd5"
  )
  
  create_one_dataset <- forest_plot <- function(data, title) {
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
    
    p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), xmin = xmin, xmax = xmax, col = q_value, fill = q_value)) +
      geom_linerange(linewidth = 3, position = position_dodge(width = 0.5)) +
      geom_vline(xintercept = 0, lty = 2, color = "black") +
      geom_point(size = 2.5, shape = 21, colour = "white", stroke = 0.5, position = position_dodge(width = 0.5)) +
      scale_fill_manual(breaks = names(dotCOLS), values = dotCOLS) +
      scale_color_manual(breaks = names(barCOLS), values = barCOLS)  +
      scale_x_continuous(name = "Coefficient", limits = c(x_min - x_buffer, x_max + x_buffer)) +
      scale_y_discrete(name = "") +
      theme_minimal() +
      facet_grid(category ~ ., scales = "free_y", space = "free_y") +
      theme(legend.position = "bottom") +
      ggtitle(title)
    
    return(p)
  }
  
plots.corhot<- create_one_dataset(with_control, var)
  
write.table(with_control,"cohort.biomarker.txt",sep="\t",quote = F,row.names = F)

pdf("FigureS3.cohort_forest_plots.L.pdf", width = 10, height = 22)
plots.corhot
dev.off()
