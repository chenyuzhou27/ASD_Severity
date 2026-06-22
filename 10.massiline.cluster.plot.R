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


count_significant_markers <- function(data) {
  data %>% filter(pval < 0.05) %>% nrow()
}
barCOLS <- c(
  "q_sig" = "#f9b282",
  "q_notsig" = "#a6d8f0"
)
dotCOLS <- c(
  "q_sig" = "#de6b35",
  "q_notsig" = "#008fd5"
)


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
  
  p <- ggplot(data, aes(x = coef, y = reorder(feature_rn, coef), xmin = xmin, xmax = xmax, col = consis, fill = consis)) +
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
with_control <- with_control %>% filter(!grepl("GGB", feature)) %>% filter(Questionnaire == class.ass) %>% mutate(consis = "q_sig") 

cluster.plots.all<- create_forest_plot(without_control, var)
pdf("FigureS3", width = 10, height = 19)
cluster.plots.all
dev.off()

# coef 0.3
with_control.0.3 <- with_control %>%
  filter(abs(coef) >0.3)

####compared with ASD diagnosis biomarker
taxa_cohort <- read_and_filter.qvalue(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
path_cohort  <- read_and_filter.qvalue(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)
module_cohort <- read_and_filter.qvalue(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.2)

taxa_cohort_pv <- read_and_filter.pvalue(paste0("taxa/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
path_cohort_pv  <- read_and_filter.pvalue(paste0("pathway/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)
module_cohort_pv <- read_and_filter.pvalue(paste0("module/Log_lm_Cohort/all_results.tsv"), "Cohort", 0.05)

cohort_related <-   combine_data(list(taxa_cohort), "Taxonomy") %>%
  bind_rows(combine_data(list(path_cohort), "Pathway")) %>%
  bind_rows(combine_data(list(module_cohort), "Module"))


cohort_related <- process_data(cohort_related) 
cohort_related <- cohort_related %>%  filter(!grepl("GGB", feature)) 

cohort_related.0.3 <- cohort_related %>% filter(abs(coef) > 0.3)

cohort_pv.related <-   combine_data(list(taxa_cohort_pv), "Taxonomy") %>%
  bind_rows(combine_data(list(path_cohort_pv), "Pathway")) %>%
  bind_rows(combine_data(list(module_cohort_pv), "Module")) %>%
  process_data(.) %>% filter(!grepl("GGB", feature)) 

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
                      labels = c("also_in_ASD_diag_with_q_value0.2" = "Also in ASD Diagnosis (q-value < 0.2)", 
                                 "also_in_ASD_diag_with_p_value0.05" = "Also in ASD Diagnosis (p-value < 0.05)",
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

# 使用新数据创建图表
with_control.0.3_node <- with_control.0.3 %>%
  mutate(node_color = case_when(
    feature %in% cohort_related$feature ~ "also_in_ASD_diag_with_q_value0.2",
    feature %in% cohort_pv.related$feature ~ "also_in_ASD_diag_with_p_value0.05",
    TRUE ~ "only_in_Cluster"
  ))


cluster.qv.plots.add.asd <- create_forest_plot_node(with_control.0.3_node,"Withcontrol.cluster")
pdf("Figure3a", width = 10, height = 13)
cluster.qv.plots.add.asd
dev.off()

