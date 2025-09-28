rm(list=ls())
###module
library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)
library(clusterSim)
library(omixerRpm)
ko.prof <- read.csv("../final_profile/ko_aug28_filter.csv")
metadata <- read.table("batch_effect/metadata.rename.questionnaire.txt",sep="\t") 

ko.prof <- read.table("batch_effect/ko.profile.adjust_bacth_cluster.txt",sep="\t")
summary(colSums(ko.prof))

ko.prof.nor <- ko.prof %>% as.data.frame() %>% rownames_to_column(var = "entry")

mods <- rpm(ko.prof.nor, minimum.coverage=0.3, annotation = 1,normalize.by.length = T, distribute = T,threads = 4)

# Load the default mapping database
db <- loadDefaultDB()
# get the name of the first predicted module
getNames(db, mods@annotation[1,])
coverage.gmm <- asDataFrame(mods, "coverage")
abundance.gmm <- asDataFrame(mods,"abundance")
write.table(abundance.gmm,"batch_effect/GMM.871s.profile.txt",sep="\t",quote = F)

listDB()
gbm.db <- loadDB("GBMs.v1.0")
mods.gbm <- rpm(ko.prof.nor, minimum.coverage=0.3, annotation = 1,module.db=gbm.db,normalize.by.length = T, distribute = T, threads = 4)
getNames(gbm.db, mods.gbm@annotation[1,])
coverage.gbm <- asDataFrame(mods.gbm, "coverage")
abundance.gbm <- asDataFrame(mods.gbm,"abundance")
write.table(abundance.gbm,"batch_effect/GBM.871s.profile.txt",sep="\t",quote = F)

gbm.des <- abundance.gbm[,1:2]
rownames(abundance.gbm) <- abundance.gbm[,1]
abundance.gbm <- abundance.gbm[,-c(1:2)]
abundance.gbm.ori <- abundance.gbm
summary(abundance.gbm[abundance.gbm>0])

abundance.gbm <- abundance.gbm %>% t() %>% as.data.frame %>% rownames_to_column(var="Sample")
###GMM
gmm.des <- abundance.gmm[,1:2]
rownames(abundance.gmm) <- abundance.gmm[,1]
abundance.gmm <- abundance.gmm[,-c(1:2)]
abundance.gmm.ori <- abundance.gmm
summary(abundance.gmm[abundance.gmm>0])

abundance.gmm <- abundance.gmm %>% t() %>% as.data.frame %>% rownames_to_column(var="Sample")

module <- merge(abundance.gmm,abundance.gbm)
pick.module <- colnames(module)
module.des <- rbind(gmm.des,gbm.des)

pid <- pmatch(metadata$sample_id,module$Sample)
all(pid == c(1:871))
module <- module[pid,]

save.image("batch_effect/module.massiline.RData")


library(Maaslin2)
abundance.sp.extract1 <- module %>% remove_rownames() %>% dplyr::select(Sample,all_of(pick.module)) %>% 
  tibble::column_to_rownames(var="Sample")

# metadata
abundance.sp.extract2 <- metadata %>% 
  dplyr::select(sample_id, Cohort, Cluster, Age, Gender, BMI, Siblings, Medication, atopic_disease, Prot, TotFib) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))

#####
fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "module/Log_lm_Cohort",
                    fixed_effects  = c("Cohort","Age","Medication","atopic_disease","TotFib"),
                    reference      = c("Cohort,Control","Medication,No","atopic_disease,No")
)
###cluster ord
fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "module/Log_lm_Cluster",
                    fixed_effects  = c("Cluster","Age","Medication","atopic_disease","TotFib"),
                    reference      = c("Medication,No","atopic_disease,No")
)


###cluster num
fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "module/Log_lm_Cluster_num",
                    fixed_effects  = c("Cluster_num","Age","Medication","atopic_disease","TotFib"),
                    reference      = c("Medication,No","atopic_disease,No")
)

run_maaslin2 <- function(fixed_effects, output_path) {
  Maaslin2(
    input_data = abundance.sp.extract1,
    input_metadata = abundance.sp.extract2,
    min_prevalence = 0.05,
    normalization = "none",
    analysis_method = "LM",
    transform = "LOG",
    output = output_path,
    fixed_effects = fixed_effects,
    reference = c("Medication,No","atopic_disease,No")
  )
}

# Define the base components for the fixed effects and output paths
base_fixed_effects <- c("Age","Medication","atopic_disease","TotFib")
base_output_path <- "module/Log_lm_"

# Define the variable parts for the fixed effects and output paths
variable_parts <- c(
  "T_SRS_total", "T_SRS_SCI", "T_SRS_RRB",
  "CBCL_AP_T", "CBCL_Internalizing_T", "CBCL_Externalizing_T",
  "ASC_total",
  "M_SEQ_total", "M_SEQ_hypo", "M_SEQ_hyper", "M_SEQ_seeking",
  "FNVD", "FAPD", "FDD", "Bristol_stool_chart"
)

abundance.sp.extract2 <- metadata %>%
  dplyr::select("sample_id","Cohort", "Cluster","Age","Gender","BMI","Siblings","Medication","atopic_disease","TotFib",all_of(variable_parts)) %>% 
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id")  %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))


fixed_effects_list <- list()
output_paths <- c()

for (var in variable_parts) {
  fixed_effects_list[[length(fixed_effects_list) + 1]] <- c(var, base_fixed_effects)
  output_paths <- c(output_paths, paste0(base_output_path, var))
}

for (i in seq_along(fixed_effects_list)) {
  run_maaslin2(fixed_effects_list[[i]], output_paths[i])
}


###without control
# Extract the metadata
abundance.sp.extract2 <- metadata %>% 
  filter(Cohort != "Control") %>%
  dplyr::select("sample_id","Cohort", "Cluster","Age","Gender","BMI","Siblings","Medication","atopic_disease","TotFib",all_of(variable_parts)) %>% 
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id")  %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))


abundance.sp.extract1 <- module %>% 
  filter(Sample %in% rownames(abundance.sp.extract2)) %>%
  remove_rownames() %>% 
  dplyr::select(Sample,all_of(pick.module)) %>% 
  tibble::column_to_rownames(var="Sample")


base_fixed_effects <- c("Age","Medication","atopic_disease","TotFib")
base_output_path <- "module/Log_lm_asd_"

variable_parts <- c(
  "Cluster","Cluster_num",
  "T_SRS_total", "T_SRS_SCI", "T_SRS_RRB",
  "CBCL_AP_T", "CBCL_Internalizing_T", "CBCL_Externalizing_T",
  "ASC_total",
  "M_SEQ_total", "M_SEQ_hypo", "M_SEQ_hyper", "M_SEQ_seeking",
  "FNVD", "FAPD", "FDD", "Bristol_stool_chart"
)

fixed_effects_list <- list()
output_paths <- c()

for (var in variable_parts) {
  fixed_effects_list[[length(fixed_effects_list) + 1]] <- c(var, base_fixed_effects)
  output_paths <- c(output_paths, paste0(base_output_path, var))
}

for (i in seq_along(fixed_effects_list)) {
  run_maaslin2(fixed_effects_list[[i]], output_paths[i])
}
