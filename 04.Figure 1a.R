library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)



data = read_excel("Transformed.xlsx")
data <- data[!data$`Subject ID` %in% c("A0122", "A0132", "A0153", "A0176", "A0351", "A0494", "A1240"), ]
clinical_feature = data[c("Subject ID",'Cluster','T_SRS_AWR', 'T_SRS_COG','T_SRS_COMM','T_SRS_MOT','T_SRS_RRB', 'CBCL_AP_T',
          'CBCL_Externalizing_T', 'ASC_PA', 'ASC_AA', 'ASC_SA',
          'ASC_uncertainty', 'M_SEQ_hypo','M_SEQ_hyper',
          'M_SEQ_seeking')]


df = clinical_feature
df_long <- df %>%
  pivot_longer(cols = 3:ncol(df), 
               names_to = "feature", 
               values_to = "value")

summary_stats <- df_long %>%
  group_by(Cluster, feature) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    lower_ci = mean(value, na.rm = TRUE) - qt(0.975, df = n() - 1) * sd(value, na.rm = TRUE) / sqrt(n()),
    upper_ci = mean(value, na.rm = TRUE) + qt(0.975, df = n() - 1) * sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
custom_colors <- c("Control" ='#BEBEBE' ,"Mild" = "#F9DB6D", "Moderate" = "#36827F", "Severe" = "#464D77")

p_line <- ggplot(summary_stats, aes(x = feature, y = mean_value, group = Cluster, color = Cluster)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Cluster), alpha = 0.5, color = NA) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Mean Line Plot with 95% Confidence Interval",
       x = "Feature",
       y = "Z-Score Mean Value",
       color = "Cluster",
       fill = "Cluster") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.y = element_line(size = 1),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_x_discrete(breaks = summary_stats$feature) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5))

ggsave("Figure 1a.pdf", plot = p_line, height = 5, width = 6)


