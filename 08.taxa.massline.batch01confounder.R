rm(list = ls())
###
load("Manuscript.RData")

####Maaslin2
###with control

library(dplyr)
abundance.sp.extract1 <- sp.prof %>% 
  dplyr::select(Sample,all_of(pick.species)) %>% 
  tibble::column_to_rownames(var="Sample")

metadata <- readRDS("metadata.carb.rds")
abundance.sp.extract2 <- metadata %>%
  dplyr::select(Subject.ID, Cohort, Cluster, Age, Gender, BMI, Siblings, Medication, atopic_disease,  Prot.g.1000kcal, TotFib.g.1000kcal) %>%
  tibble::column_to_rownames(var = "Subject.ID") %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))

all(rownames(abundance.sp.extract2) == rownames(abundance.sp.extract1))

###cohort
library(Maaslin2)
fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "taxa/Log_lm_Cohort",
                    fixed_effects  = c("Cohort","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal"),
                    reference      = c("Cohort,Control","Gender,Female","Medication,No","atopic_disease,No")
)

###cluster ord
fit_data = Maaslin2(input_data     = abundance.sp.extract1,
                    input_metadata = abundance.sp.extract2,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = "taxa/Log_lm_Cluster",

                    fixed_effects  = c("Cluster","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal"),
                    reference      = c("Gender,Female","Medication,No","atopic_disease,No")
                    
)


# Define a function to run Maaslin2
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
    reference = c("Gender,Female","Medication,No","atopic_disease,No")
  )
}

# Define the base components for the fixed effects and output paths
base_fixed_effects <- c("Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal")
base_output_path <- "taxa/Log_lm_"

###without control for score
sp.prof.asd <- sp.prof %>% filter(Cohort != "Control")

abundance.sp.extract1 <- sp.prof.asd %>% 
  dplyr::select(Sample,all_of(pick.species)) %>% 
  tibble::column_to_rownames(var="Sample")


# Define the output paths
base_output_path <- "taxa/Log_lm_asd_"

# Define the variable parts for the fixed effects and output paths
variable_parts <- c(
  "T_SRS_total", "T_SRS_SCI", "T_SRS_RRB",
  "CBCL_AP_T", "CBCL_Internalizing_T", "CBCL_Externalizing_T",
  "ASC_total",
  "M_SEQ_total", "M_SEQ_hypo", "M_SEQ_hyper", "M_SEQ_seeking"
)

abundance.sp.extract2 <- metadata %>%
  dplyr::select("Subject.ID","Cohort", "Cluster","Age","Gender", "BMI","Medication","atopic_disease","Prot.g.1000kcal","TotFib.g.1000kcal",all_of(variable_parts[-c(1,2)])) %>% 
  filter(Subject.ID %in% rownames(abundance.sp.extract1)) %>%
  tibble::column_to_rownames(var = "Subject.ID") %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         Cluster = factor(Cluster,levels = c("Control",  "Mildest", "Moderate",   "Severe" ),ordered = T),
         Cluster_num = as.numeric(Cluster),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease)))
str(abundance.sp.extract2)
all(rownames(abundance.sp.extract2) == rownames(abundance.sp.extract1))
# Generate the fixed_effects_list and output_paths using loops
fixed_effects_list <- list()
output_paths <- c()

for (var in variable_parts) {
  fixed_effects_list[[length(fixed_effects_list) + 1]] <- c(var, base_fixed_effects)
  output_paths <- c(output_paths, paste0(base_output_path, var))
}

rm(var)
# Run Maaslin2 for each set of fixed effects and output paths
for (i in seq_along(fixed_effects_list)) {
  run_maaslin2(fixed_effects_list[[i]], output_paths[i])
}

