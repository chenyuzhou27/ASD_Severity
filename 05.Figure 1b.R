library(MicrobiotaProcess)
library(readxl)
library(dplyr)
library(vegan)
library(ggplot2)


df <- read_excel("Transformed.xlsx")
df <- df[!df$`Subject ID` %in% c("A0122", "A0132", "A0153", "A0176", "A0351", "A0494", "A1240"), ]

sample_id_col = 'Subject ID'


feature_cols <- setdiff(colnames(df), c(sample_id_col, "Cluster"))
df_meta <- df[, c(sample_id_col, "Cluster")] %>% mutate_all(as.character)
df_t <- as.data.frame(t(df[, feature_cols]))
colnames(df_t) <- df_meta[[sample_id_col]]

mpse_obj <- MPSE(assays = df_t, colData = df_meta)

abundance_matrix <- as.matrix(mpse_obj@assays@data@listData[["Abundance"]])
storage.mode(abundance_matrix) <- "double"

mpse_obj@assays@data@listData[["Abundance"]] <- abundance_matrix
mpse_obj@colData@listData[[sample_id_col]] <- df_meta[[sample_id_col]]
mpse_obj@colData@listData[["Cluster"]] <- df_meta$Cluster

mpse_asd <- mpse_obj %>%
  mp_cal_pca(.abundance = Abundance, action = "add")

pdf("Figure 1b.pdf", height = 5, width = 13)

cluster_colors <- c(
  Control = "#BEBEBE",
  Mild = "#F9DB6D",
  Moderate = "#36827F",
  Severe = "#464D77"
)

p.pca <- mpse_asd %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cluster,
    .color = Cluster,
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_color_manual(values = cluster_colors, guide = "none")

print(p.pca)
dev.off()

abundance_for_dist <- t(abundance_matrix)

set.seed(123)
dist_matrix_all <- vegdist(abundance_for_dist, method = "euclidean")
permanova_with_control <- adonis2(
  dist_matrix_all ~ Cluster,
  data = df_meta,
  permutations = 999
)

asd_keep <- df_meta$Cluster %in% c("Mild", "Moderate", "Severe")
df_meta_asd <- df_meta[asd_keep, , drop = FALSE]
abundance_asd <- abundance_for_dist[asd_keep, , drop = FALSE]
dist_matrix_asd <- vegdist(abundance_asd, method = "euclidean")
permanova_without_control <- adonis2(
  dist_matrix_asd ~ Cluster,
  data = df_meta_asd,
  permutations = 999
)

print("PERMANOVA including Control:")
print(permanova_with_control)

print("PERMANOVA excluding Control (ASD severity only):")
print(permanova_without_control)

