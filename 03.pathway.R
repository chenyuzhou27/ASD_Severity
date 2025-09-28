rm(list=ls())
###pathway
library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)

###used the batch effect adjust profile
pathway <- read.table("pathway.profile.adjust_bacth_cluster.txt",sep="\t")
pathway <- t(pathway) %>% as.data.frame() %>% rownames_to_column( var = "Subject.ID")

zero.num <-apply(pathway, 2, function(x){sum(x!=0)})
pick.pathway <- names(which(zero.num/dim(pathway)[1] > 0.05))

pid <- pmatch(pathway$Subject.ID,metadata$sample_id)
all(pid == c(1:871))

# Extract the abundance data
abundance.sp.extract1 <- pathway %>%
  dplyr::select(Subject.ID, all_of(pick.pathway)) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "Subject.ID")

# Extract the metadata
abundance.sp.extract2 <- metadata %>% 
  dplyr::select("sample_id","Cohort", "Cluster","Age","Gender","BMI","Siblings","Medication","atopic_disease","TotFib") %>% 
  remove_rownames() %>%
  tibble::column_to_rownames(var = "sample_id")  %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))


fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "pathway/Log_lm_Cohort",
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
                    output         = "pathway/Log_lm_Cluster",
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
                    output         = "pathway/Log_lm_Cluster_num",
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

base_fixed_effects <- c("Age","Medication","atopic_disease","TotFib")
base_output_path <- "pathway/Log_lm_"

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


abundance.sp.extract1 <- pathway %>%
  dplyr::select(Subject.ID, all_of(pick.pathway)) %>% 
  filter(Subject.ID %in% rownames(abundance.sp.extract2)) %>%
  remove_rownames() %>%
  tibble::column_to_rownames(var = "Subject.ID")

base_fixed_effects <- c("Age","Medication","atopic_disease","TotFib")
base_output_path <- "pathway/Log_lm_asd_"


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

