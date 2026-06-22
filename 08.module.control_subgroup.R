rm(list=ls())
###module
library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)

load("module.massiline.RData")
library(Maaslin2)

subgroup <- unique(metadata$Cluster)[-4]
for(i in subgroup){
  metadata.sp.tmp<- metadata %>%
    dplyr::select(sample_id, Cohort, Cluster, Age, Gender, BMI, Siblings, Medication, atopic_disease,  TotFib) %>%
    remove_rownames() %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
           Gender = as.factor(Gender),
           BMI = as.numeric(BMI),
           Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
           atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease))) %>%
    filter(Cluster %in% c("Control",  i)) 
  
  abundance.sp.tmp<-  module %>% remove_rownames() %>% dplyr::select(Sample,all_of(pick.module)) %>% 
    filter(Sample %in% rownames(metadata.sp.tmp)) %>%
    tibble::column_to_rownames(var="Sample") 
  
  ###cohort
  fit_data = Maaslin2(input_data     = abundance.sp.tmp,
                      input_metadata = metadata.sp.tmp,
                      min_prevalence = 0.05,
                      normalization  = "none",
                      analysis_method = "LM",
                      transform = "LOG",
                      output         = paste(i,"/module",sep=""),
                      fixed_effects  = c("Cluster","Age","Medication","atopic_disease","TotFib"),
                      reference      = c("Cluster,Control","Medication,No","atopic_disease,No")
  )
}

