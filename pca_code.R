library(MicrobiotaProcess)
library(tidyr)
library(gghalves)
library(corrr)
library(ggstatsplot)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(ggside)
library(tibble)
library(flextable)
library(readxl)  
library(dplyr)          

# Step 1: Read the Excel file
df <- read_excel("./kmeans_filled.xlsx")

# Step 2: Data cleaning
df <- df %>%
  filter(!is.na(Cluster)) %>%   # Remove rows with NA in the "Cluster" column
  filter(Cluster != "Control") # Remove rows where "Cluster" is "Control"

# Step 3: Scale numeric columns (Z-score normalization)
numeric_cols <- sapply(df, is.numeric) # Identify numeric columns
df[, numeric_cols] <- scale(df[, numeric_cols])

df_num = df[,3:ncol(df)]
df_meta = df[,1:2]

df_t <- as.data.frame(t(df))

colnames(df_t) <- df_t[1, ]  
df_t <- df_t[-1, ]           
df_t <- df_t[-1, ]  
df_t <- cbind(feature = rownames(df_t), df_t)

rownames(df_t) <- NULL

df_meta <- df_meta %>% mutate_all(as.character)
mpse_obj <- MPSE(assays=df_t[,2:ncol(df_t)],colData = df_meta)

abundance_df <- mpse_obj@assays@data@listData[["Abundance"]]

if (is.data.frame(abundance_df)) {

  abundance_double <- as.matrix(abundance_df)
  storage.mode(abundance_double) <- "double"  
} else {
  stop("Abundance is not a data.frame!")
}


print(dim(abundance_double))
print(typeof(abundance_double)) 

mpse_obj@assays@data@listData[["Abundance"]] = abundance_double

mpse_obj@colData@listData[["Subject ID"]]=df_meta$'Subject ID'
mpse_obj@colData@listData[["Cluster"]]=df_meta$Cluster

###PCA
mpse_asd <- mpse_obj %>%
  mp_cal_pca(.abundance=Abundance, action="add")

pdf("./pca.pdf", height = 5, width = 13)

p.pca <- mpse_asd %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cluster,
    .color = Cluster,
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values = c("#F9DB6D", "#36827F", "#464D77"), guide = "none") +
  scale_color_manual(values = c("#F9DB6D", "#36827F", "#464D77"), guide = "none")

saveRDS(p.pca, file = "./pca.rds")
pca <- readRDS("./pca.rds") 
library(vegan)
library(MicrobiotaProcess)

permanova_data <- p.pca[["data"]] %>% 
  select(`Subject ID`, 'PC1 (47.58%)', 'PC2 (12.75%)', Cluster) 

dist_matrix <- vegdist(permanova_data[, c('PC1 (47.58%)', 'PC2 (12.75%)')], method = "euclidean")

set.seed(123)
adonis_result <- adonis2(
  dist_matrix ~ Cluster,
  data = permanova_data,
  permutations = 999
)
print(adonis_result)

p.pca +
  labs(caption = sprintf(
    "PERMANOVA: R²=%.2f, p=%.3f", 
    adonis_result$R2[1], 
    adonis_result$`Pr(>F)`[1]
  )) +
  theme(plot.caption = element_text(size=10, hjust=0.5))


dev.off()

