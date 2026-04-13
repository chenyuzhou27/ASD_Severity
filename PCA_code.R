# 加载必要的包
library(MicrobiotaProcess)
library(readxl)
library(dplyr)
library(vegan)
library(ggplot2)

# Step 1: Read and clean data
df <- read_excel("./kmeans_filled.xlsx") %>%
  filter(!is.na(Cluster), Cluster != "Control")

# Step 2: Z-score normalization for numeric columns
numeric_cols <- sapply(df, is.numeric)
df[, numeric_cols] <- scale(df[, numeric_cols])

# Step 3: Prepare data for MPSE object
df_meta <- df[, 1:2] %>% mutate_all(as.character)
df_t <- as.data.frame(t(df[, 3:ncol(df)]))
colnames(df_t) <- df_meta$`Subject ID`

# Step 4: Create MPSE object
mpse_obj <- MPSE(assays = df_t, colData = df_meta)

# Ensure abundance data is numeric matrix
abundance_matrix <- as.matrix(mpse_obj@assays@data@listData[["Abundance"]])
storage.mode(abundance_matrix) <- "double"

mpse_obj@assays@data@listData[["Abundance"]] <- abundance_matrix
mpse_obj@colData@listData[["Subject ID"]] <- df_meta$`Subject ID`
mpse_obj@colData@listData[["Cluster"]] <- df_meta$Cluster

# Step 5: Perform PCA
mpse_asd <- mpse_obj %>%
  mp_cal_pca(.abundance = Abundance, action = "add")

# Step 6: Create PCA plot with PERMANOVA results
pdf("./HATCH/draft/pcoa.pca.2groups.pdf", height = 5, width = 13)

# Calculate PERMANOVA on original abundance data
abundance_for_dist <- t(abundance_matrix)
dist_matrix <- vegdist(abundance_for_dist, method = "euclidean")

set.seed(123)
adonis_result <- adonis2(dist_matrix ~ Cluster, 
                         data = df_meta, 
                         permutations = 999)

# Create plot
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
  scale_fill_manual(values = c("#F9DB6D", "#36827F", "#464D77"), guide = "none") +
  scale_color_manual(values = c("#F9DB6D", "#36827F", "#464D77"), guide = "none") +
  labs(caption = sprintf("PERMANOVA: R²=%.2f, p=%.3f", 
                         adonis_result$R2[1], 
                         adonis_result$`Pr(>F)`[1])) +
  theme(plot.caption = element_text(size = 10, hjust = 0.5))

print(p.pca)
dev.off()

# Save results
saveRDS(p.pca, file = "./HATCH/draft/pca.rds")
print(adonis_result)