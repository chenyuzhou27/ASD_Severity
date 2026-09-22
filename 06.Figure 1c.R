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

# Retain the same asymptotic Wilcoxon test, ties and continuity correction.
# Log probabilities avoid underflow when the ordinary p-value is zero.
wilcox_log_p <- function(group1, group2) {
  x <- euclidean_dist_df.within.cluster$Distance[
    euclidean_dist_df.within.cluster$Cluster == group1]
  y <- euclidean_dist_df.within.cluster$Distance[
    euclidean_dist_df.within.cluster$Cluster == group2]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- as.double(length(x))
  n2 <- as.double(length(y))
  stopifnot(n1 >= 50L, n2 >= 50L)
  ranks <- rank(c(x, y))
  tie_sizes <- as.numeric(table(ranks))
  u <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  sigma <- sqrt((n1 * n2 / 12) * ((n1 + n2 + 1) -
    sum(tie_sizes^3 - tie_sizes) / ((n1 + n2) * (n1 + n2 - 1))))
  stopifnot(is.finite(sigma), sigma > 0)
  centered_u <- u - n1 * n2 / 2
  z <- (centered_u - sign(centered_u) * 0.5) / sigma
  log_p <- min(0, log(2) + pnorm(-abs(z), log.p = TRUE))
  ordinary_p <- wilcox.test(x, y, alternative = "two.sided",
                           exact = FALSE, correct = TRUE)$p.value
  if (ordinary_p > .Machine$double.xmin) {
    stopifnot(isTRUE(all.equal(log_p, log(ordinary_p), tolerance = 1e-10)))
  }
  log_p
}

bh_log_adjust <- function(log_p) {
  n <- length(log_p)
  order_desc <- order(log_p, decreasing = TRUE)
  adjusted <- pmin(0, cummin(log_p[order_desc] + log(n / seq.int(n, 1))))
  adjusted[order(order_desc)]
}

scientific_log <- function(log_value, digits = 3L) {
  vapply(log_value, function(value) {
    if (is.na(value)) return(NA_character_)
    if (value == -Inf) return("0")
    exponent <- floor(value / log(10))
    mantissa <- exp(value - exponent * log(10))
    if (round(mantissa, digits) >= 10) {
      mantissa <- 1
      exponent <- exponent + 1
    }
    sprintf(paste0("%.", digits, "fe%+d"), mantissa, as.integer(exponent))
  }, character(1))
}

pairwise_results$log.p <- mapply(wilcox_log_p,
                                pairwise_results$group1, pairwise_results$group2)
pairwise_results$log.q <- bh_log_adjust(pairwise_results$log.p)
pairwise_results$p <- exp(pairwise_results$log.p)
pairwise_results$p.adj <- exp(pairwise_results$log.q)
pairwise_results$p.scientific <- scientific_log(pairwise_results$log.p)
pairwise_results$q.scientific <- scientific_log(pairwise_results$log.q)
pairwise_results$q.label <- paste0("q = ", pairwise_results$q.scientific)
pairwise_results <- rstatix::add_significance(pairwise_results, "p.adj")

p <- p + stat_pvalue_manual(
  pairwise_results,
  label = "q.label",
  size = 3.5,
  bracket.size = 0.8,
  tip.length = 0.01,
  inherit.aes = FALSE
)

print(as.data.frame(pairwise_results[, c("group1", "group2", "p.scientific",
                                        "q.scientific", "log.p", "log.q")]),
      digits = 16, row.names = FALSE)


print(medians)
print(p)

ggsave("Figure 1c.pdf")

kruskal_test <- kruskal.test(Distance ~ Cluster, data = euclidean_dist_df.within.cluster)

cat(sprintf(
  "Kruskal-Wallis H statistic = %.3f, p-value = %.3e\n",
  kruskal_test$statistic, 
  kruskal_test$p.value
))

library(openxlsx)

# Numeric p/p.adj may underflow to zero; scientific-text and log columns retain
# the computed asymptotic probabilities. q denotes BH-adjusted p, not Storey q.
write.xlsx(pairwise_results, "Figure 1c.xlsx")
