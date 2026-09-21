library(openxlsx)
library(MatchIt)

ct = read.xlsx("Raw_control.xlsx", sheet = "Sheet1")
metadata = read.xlsx("HATCH_metadata.xlsx", sheet = "metadata")
cluster_result = read.xlsx("Cluster_result.xlsx", sheet = "Sheet1")

asd = metadata[match(cluster_result$Subject.ID, metadata$Subject.ID),c("Subject.ID", "Age", "Gender", "Height.(cm)", "Weight.(kg)", "BMI")]
asd$Cluster = cluster_result$Cluster

ct = ct[, c('Subject.ID', 'Cluster','Age','Gender','Height.(cm)','Weight.(kg)','BMI')]
ct$ASD = 0
asd = asd[, c('Subject.ID', 'Cluster','Age','Gender','Height.(cm)','Weight.(kg)','BMI')]
asd$ASD = 1

run_psm = function(gender_value, severe_ratio) {
  merged0 <- rbind(asd, ct)
  merged0$Age <- as.numeric(merged0$Age)
  merged0 = merged0[merged0$Gender == gender_value,]

  merged <- as.data.frame(lapply(merged0[, c(3, 5:7)], scale))
  merged$Gender = merged0$Gender
  merged$Gender <- factor(merged$Gender)
  merged$Subject.ID = merged0$Subject.ID
  merged$Cluster = merged0$Cluster
  merged$ASD = merged0$ASD

  names(merged)[which(names(merged) == "Height..cm.")] <- "Height"
  names(merged)[which(names(merged) == "Weight..kg.")] <- "Weight"

  ct_mild = merged[merged$Cluster == "Mild" | merged$Cluster == "control", ]
  matching_ct_mild = matchit(
    ASD ~ Age + Height + Weight + BMI,
    data = ct_mild,
    method = "nearest",
    ratio = 1,
    replace = TRUE
  )
  matched_ct_mild = match.data(matching_ct_mild)

  ct_moder = merged[merged$Cluster == "Moderate" | merged$Cluster == "control", ]
  matching_ct_moder = matchit(
    ASD ~ Age + Height + Weight + BMI,
    data = ct_moder,
    method = "nearest",
    ratio = 1,
    replace = TRUE
  )
  matched_ct_moder = match.data(matching_ct_moder)

  ct_seve = merged[merged$Cluster == "Severe" | merged$Cluster == "control", ]
  matching_ct_seve = matchit(
    ASD ~ Age + Height + Weight + BMI,
    data = ct_seve,
    method = "nearest",
    ratio = severe_ratio,
    replace = TRUE
  )
  matched_ct_seve = match.data(matching_ct_seve)

  matched_ct_mild_ct = matched_ct_mild[matched_ct_mild$Cluster == "control", ]
  matched_ct_moder_ct = matched_ct_moder[matched_ct_moder$Cluster == "control", ]
  matched_ct_seve_ct = matched_ct_seve[matched_ct_seve$Cluster == "control", ]

  id_mild <- matched_ct_mild_ct$Subject.ID
  id_moder <- matched_ct_moder_ct$Subject.ID
  id_seve <- matched_ct_seve_ct$Subject.ID

  inter = intersect(intersect(id_mild, id_moder), id_seve)
  all_ids <- c(id_mild, id_moder, id_seve)
  at_least_twice <- unique(all_ids[duplicated(all_ids) | duplicated(all_ids, fromLast = TRUE)])
  union_ids <- union(union(id_mild, id_moder), id_seve)

  print(union_ids)
  print(at_least_twice)

  ct2 <- subset(ct, Subject.ID %in% at_least_twice)
  invisible(ct2)
}

fe_result = run_psm(
  gender_value = '2',
  severe_ratio = 5
)

ma_result = run_psm(
  gender_value = '1',
  severe_ratio = 2
)

combined_result = rbind(fe_result, ma_result)

combined_result$Cluster[grepl("^control", combined_result$Cluster, ignore.case = TRUE)] <- "Control"

write.xlsx(
  combined_result,
  "Control_selected.xlsx"
)

