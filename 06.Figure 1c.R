library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)

df_feature <- read_excel("Transformed.xlsx")
df_feature <- df_feature[!df_feature$`Subject ID` %in% c("A0122", "A0132", "A0153", "A0176", "A0351", "A0494", "A1240"), ]
feature_columns <- setdiff(colnames(df_feature), c("Cluster", "Subject ID"))


euclidean_dist_df.within.cluster <- data.frame()

clusters <- unique(df_feature$Cluster)
for(cluster in clusters) {
  group_data <- df_feature %>% 
    filter(Cluster == cluster) %>% 
    select(all_of(feature_columns))
  
  if(nrow(group_data) > 1) {
    distances <- as.vector(dist(group_data, method = "euclidean"))
    distances <- round(distances, digits = 4)
    temp_df <- data.frame(Cluster = cluster, Distance = distances)
    euclidean_dist_df.within.cluster <- rbind(euclidean_dist_df.within.cluster, temp_df)
  }
}

medians <- euclidean_dist_df.within.cluster %>%
  group_by(Cluster) %>%
  summarise(
    Median = round(median(Distance, na.rm = TRUE), digits = 3),
    Count = sum(!is.na(Distance))
  )

p <- ggplot(euclidean_dist_df.within.cluster, aes(x = Cluster, y = Distance, fill = Cluster)) +
  geom_boxplot(
    outlier.size = 1.5,
    lwd = 1.4,
    fatten = 2
  ) +
stat_summary(
  fun = median,
  geom = "text",
  aes(label = sprintf("%.3f", after_stat(y))),
  vjust = -0.5,
  size = 4,
  fontface = "bold",
  color = "black"
) +
scale_fill_manual(values = c("Control" = "gray", 
                             "Mild" = "#F9DB6D", 
                             "Moderate" = "#36827F", 
                             "Severe" = "#464D77")) +
  labs(x = "ASD Severity", y = "Euclidean Distance", title = "Euclidean Distance by ASD Severity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- p + stat_compare_means(
  method = "kruskal.test", 
  label.y = max(euclidean_dist_df.within.cluster$Distance, na.rm = TRUE) * 1.15,
  label = "p.format"
)

comparisons <- list(
  c("Control", "Mild"),
  c("Control", "Moderate"),
  c("Control", "Severe"),
  c("Mild", "Moderate"),
  c("Mild", "Severe"),
  c("Moderate", "Severe")
)

pairwise_results <- euclidean_dist_df.within.cluster %>%
  rstatix::wilcox_test(
    Distance ~ Cluster,
    comparisons = comparisons,
    p.adjust.method = "BH"
  ) %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(
    x = "Cluster",
    fun = "max",
    step.increase = 0.1
  )

p <- p + stat_pvalue_manual(
  pairwise_results,
  label = "p.adj.signif",
  size = 6,
  bracket.size = 0.8,
  tip.length = 0.01,
  inherit.aes = FALSE
)

print(pairwise_results)


print(medians)
print(p)

ggsave("boxplot.pdf")

kruskal_test <- kruskal.test(Distance ~ Cluster, data = euclidean_dist_df.within.cluster)

cat(sprintf(
  "Kruskal-Wallis H statistic = %.3f, p-value = %.3e\n",
  kruskal_test$statistic, 
  kruskal_test$p.value
))

library(openxlsx)

write.xlsx(pairwise_results, "Figure 1c.xlsx")
