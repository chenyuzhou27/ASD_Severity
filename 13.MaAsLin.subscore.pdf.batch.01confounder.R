# 清除工作区
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
# 模块描述数据
module.des <- read.table("module.descript.txt",sep="\t")
###MF0030 and MGB050 is the same moduel, only keep 1
###MGB051	and MF0032 is the same module, only keep 1
module.des <- module.des %>% filter(Module != "MF0030" & Module != "MF0032")

# 读取数据并过滤
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
    psignif = 0.1,
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
    psignif = 0.1,
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
    psignif = 0.1,
    xlab = "Coefficient",
    title = paste("c.", title, "- Pathway"),
    colour = Questionnaire
  ) +
    theme_minimal()
  
  sumary.plot <- data %>% select(feature, category) %>% unique() %>% count(category)
  return(list(taxonomy = p_taxonomy, module = p_module, pathway = p_pathway,summary.plot = sumary.plot))
}

# 变量列表
srs_variables <- c("T_SRS_total",  "T_SRS_RRB", "T_SRS_SCI")
cbcl_variables <- c( "CBCL_AP_T", "CBCL_Externalizing_T")
asc_asd_variables <- c("ASC_total")
seq_variables <- c( "M_SEQ_hypo", "M_SEQ_hyper","M_SEQ_seeking")


read_data <- function(variables, prefix) {
  taxa_without_control <- lapply(variables, function(var) read_and_filter.qvalue(paste0("taxa/Log_lm_asd_", var, "/all_results.tsv"), var,0.2))
  path_without_control <- lapply(variables, function(var) read_and_filter.qvalue(paste0("pathway/Log_lm_asd_", var, "/all_results.tsv"), var,0.2))
  module_without_control <- lapply(variables, function(var) read_and_filter.qvalue(paste0("module/Log_lm_asd_", var, "/all_results.tsv"), var,0.2))
  
  
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
  
  significant_markers_without_control <- count_significant_markers.qv(without_control)
  cat(paste("Significant markers without control for", prefix, ":", significant_markers_without_control, "\n"))
  
  return(without_control)
}

plot_list <- list()

plot_list$srs <- create_forest_plot.q0.2(read_data(srs_variables, "SRS"), "SRS")
plot_list$cbcl <- create_forest_plot.q0.2(read_data(cbcl_variables, "CBCL"), "CBCL")
plot_list$asc_asd <- create_forest_plot.q0.2(read_data(asc_asd_variables, "ASC-ASD"), "ASC-ASD")
plot_list$seq <- create_forest_plot.q0.2(read_data(seq_variables, "SEQ"), "SEQ")



symptom.qv <- rbind.data.frame(read_data(srs_variables, "SRS"),
                               read_data(cbcl_variables, "CBCL"),
                               read_data(asc_asd_variables, "ASC-ASD"),
                               read_data(seq_variables, "SEQ"))

write.xlsx(symptom.qv, "symptom.qvalue.xlsx", sheetName = "symptom.qvalue.marker")

layout_matrix <- matrix(c(1, 3, 2,3), ncol = 2, nrow = 2, byrow = TRUE)

###srs
p_taxonomy <- plot_list[["srs"]]$taxonomy
p_module <- plot_list[["srs"]]$module
p_pathway <- plot_list[["srs"]]$pathway
ggarrange(p_taxonomy,ggarrange(p_module,p_pathway,align = "hv",labels = "auto",heights = c(1, 2),nrow = 2,common.legend = F),
          widths = c(1 , 1.6),ncol = 2,common.legend = F,legend = "right")

pdf("Figure4.SRS_forest_plot.qv.pdf",height = 5,width = 16)
ggarrange(p_taxonomy,ggarrange(p_module,p_pathway,align = "hv",labels = "auto",heights = c(1, 2),nrow = 2,common.legend = T),
          widths = c(1 , 1.6),ncol = 2,common.legend = T,legend = "right")
dev.off()
  
###cbcl
p_taxonomy <- plot_list[["cbcl"]]$taxonomy
p_module <- plot_list[["cbcl"]]$module
p_pathway <- plot_list[["cbcl"]]$pathway

my_colors <- c(
  "CBCL_Externalizing_T" = "#E41A1C", 
  "CBCL_AP_T"            = "black"
)

p_module <- p_module + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)


ggarrange(p_taxonomy,ggarrange(p_module,p_pathway,align = "hv",labels = "auto",heights = c(1, 3.5),nrow = 2,common.legend = F),
          widths = c(1 , 1.6),ncol = 2,common.legend = F,legend = "right")


pdf("Figure5.cbcl_forest_plot.qv.pdf",height = 8,width = 16)
ggarrange(p_taxonomy,ggarrange(p_module,p_pathway,align = "hv",labels = "auto",heights = c(1, 3.5),nrow = 2,common.legend = T),
          widths = c(1 , 1.6),ncol = 2,common.legend = T,legend = "right")
dev.off()

###seq
p_taxonomy <- plot_list[["seq"]]$taxonomy
p_module <- plot_list[["seq"]]$module
p_pathway <- plot_list[["seq"]]$pathway

ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = F),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = F,legend = "right")

pdf("Figure4.seq_forest_plot.qv.pdf",height = 10,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = T,legend = "right")
dev.off()

###asc
p_taxonomy <- plot_list[["asc_asd"]]$taxonomy
p_module <- plot_list[["asc_asd"]]$module
p_pathway <- plot_list[["asc_asd"]]$pathway

ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(2.5, 1),nrow = 2,common.legend = F),
          p_pathway,widths = c(1 , 1.8),ncol = 2,common.legend = F,legend = "right")

pdf("Figure5.asc_asd_forest_plot.qv.pdf",height = 6,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(2.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.8),ncol = 2,common.legend = T,legend = "right")
dev.off()

####data summary
result <- list()
result$srs <- read_data(srs_variables, "SRS")
result$cbcl <- read_data(cbcl_variables, "CBCL")
result$asc_asd <- read_data(asc_asd_variables, "ASC-ASD")
result$seq <- read_data(seq_variables, "SEQ")

colnames(result$srs)
srs_features <- result$srs %>% distinct(feature,category)
cbcl_features <- result$cbcl %>% distinct(feature,category)
asc_asd_features <- result$asc_asd %>% distinct(feature,category)
seq_features <- result$seq %>% distinct(feature,category)

###plot
library(ComplexHeatmap)
library(dplyr)

srs_features_vec <- srs_features$feature %>% unique()
cbcl_features_vec <- cbcl_features$feature %>% unique()
asc_asd_features_vec <- asc_asd_features$feature %>% unique()
seq_features_vec <- seq_features$feature %>% unique()
set_colors <- c("SRS" = "#F9DB6D", "CBCL" = "#36827F", "ASC-ASD" = "#464D77", "SEQ" = "#A8D08D")
library(UpSetR)
library(dplyr)
feature_sets <- list(
  "SRS" = srs_features_vec,
  "CBCL" = cbcl_features_vec,
  "ASC-ASD" = asc_asd_features_vec,
  "SEQ" = seq_features_vec
)

###3 catalogs
srs_features$data_source <- "SRS"
cbcl_features$data_source <- "CBCL"
asc_asd_features$data_source <- "ASC-ASD"
seq_features$data_source <- "SEQ"

# 2. 合并所有数据
all_features <- bind_rows(
  srs_features,
  cbcl_features,
  asc_asd_features,
  seq_features
)

# 3. 对每个 category 分别绘图
categories <- c( "Pathway", "Taxonomy","Module")
for (cat in categories) {
  
  # 筛选当前 category
  sub <- all_features %>% filter(category == cat)
  
  # 构建列表
  feature_sets <- list(
    "SRS"     = unique(sub$feature[sub$data_source == "SRS"]),
    "CBCL"    = unique(sub$feature[sub$data_source == "CBCL"]),
    "ASC-ASD" = unique(sub$feature[sub$data_source == "ASC-ASD"]),
    "SEQ"     = unique(sub$feature[sub$data_source == "SEQ"])
  )
  
  # 去掉空集合
  feature_sets <- feature_sets[sapply(feature_sets, length) > 0]

  # 保存
  p.tmp <- upset(
    fromList(feature_sets),
    order.by = "freq",
    decreasing = TRUE,
    number.angles = 45,
    point.size = 3.5,
    line.size = 1.5,
    mainbar.y.label = paste0("Feature Intersections-",cat),
    sets.x.label = "Number of Features",
    sets.bar.color = set_colors,
    main.bar.color = "lightblue",
    
  )
  assign(paste0("p.", cat), p.tmp)
  

}


pdf("Figure6.UpSet_Pathway.pdf",width = 8,height = 7)
p.Pathway
dev.off()
pdf("Figure6.UpSet_Taxonomy.pdf",width = 8,height = 7)
p.Taxonomy
dev.off()
pdf("Figure6.UpSet_Module.pdf",width = 4,height = 7)
p.Module
dev.off()


###pvalue 0.05
create_forest_plot.p0.05 <- function(data, title) {
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
  taxa_without_control <- lapply(variables, function(var) read_and_filter.pvalue(paste0("taxa/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  path_without_control <- lapply(variables, function(var) read_and_filter.pvalue(paste0("pathway/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  module_without_control <- lapply(variables, function(var) read_and_filter.pvalue(paste0("module/Log_lm_asd_", var, "/all_results.tsv"), var,0.05))
  
  
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
  
  significant_markers_without_control <- count_significant_markers.qv(without_control)
  
  cat(paste("Significant markers without control for", prefix, ":", significant_markers_without_control, "\n"))
  
  return(without_control)
}

sym.pv.df <- rbind.data.frame(read_data(srs_variables, "SRS"),read_data(cbcl_variables, "CBCL"),read_data(asc_asd_variables, "ASC-ASD"),read_data(seq_variables, "SEQ"))
write.table(sym.pv.df,"SFigure5_8.txt",sep="\t",quote = F)

plot_list.pv <- list()
plot_list.pv$srs <- create_forest_plot.p0.05(read_data(srs_variables, "SRS"), "SRS")
plot_list.pv$cbcl <- create_forest_plot.p0.05(read_data(cbcl_variables, "CBCL"), "CBCL")
plot_list.pv$asc_asd <- create_forest_plot.p0.05(read_data(asc_asd_variables, "ASC-ASD"), "ASC-ASD")
plot_list.pv$seq <- create_forest_plot.p0.05(read_data(seq_variables, "SEQ"), "SEQ")

layout_matrix <- matrix(c(1, 3, 2,3), ncol = 2, nrow = 2, byrow = TRUE)


###srs
p_taxonomy <- plot_list.pv[["srs"]]$taxonomy
p_module <- plot_list.pv[["srs"]]$module
p_pathway <- plot_list.pv[["srs"]]$pathway

pdf("SFigure6.SRS_forest_plot.pv.pdf",height = 7,width = 14)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(2.3, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.3),ncol = 2,common.legend = T,legend = "right")
dev.off()

###cbcl
p_taxonomy <- plot_list.pv[["cbcl"]]$taxonomy
p_module <- plot_list.pv[["cbcl"]]$module
p_pathway <- plot_list.pv[["cbcl"]]$pathway

pdf("SFigure8.cbcl_forest_plot.pv.pdf",height = 11,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.6),ncol = 2,common.legend = T,legend = "right")
dev.off()

###seq
p_taxonomy <- plot_list.pv[["seq"]]$taxonomy
p_module <- plot_list.pv[["seq"]]$module
p_pathway <- plot_list.pv[["seq"]]$pathway

pdf("SFigure7.seq_forest_plot.pv.pdf",height = 12,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = T,legend = "right")
dev.off()

###ASD
###srs
p_taxonomy <- plot_list.pv[["asc_asd"]]$taxonomy
p_module <- plot_list.pv[["asc_asd"]]$module
p_pathway <- plot_list.pv[["asc_asd"]]$pathway

pdf("SFigure9.asc_asd_forest_plot.pv.pdf",height = 10,width = 16)
ggarrange(ggarrange(p_taxonomy, p_module,align = "hv",labels = "auto",heights = c(1.5, 1),nrow = 2,common.legend = T),
          p_pathway,widths = c(1 , 1.4),ncol = 2,common.legend = T,legend = "right")
dev.off()

