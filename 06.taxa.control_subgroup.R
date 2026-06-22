rm(list = ls())
###
load("Manuscript.RData")

####Maaslin2

library(Maaslin2)

folderlist <- file.path(rep(unique(sp.prof$Cluster)[-4], each=3), c("taxa","pathway","module"))
sapply(folderlist, dir.create, recursive = TRUE, showWarnings = FALSE)

subgroup <- unique(sp.prof$Cluster)[-4]
for(i in subgroup){
metadata.sp.tmp<- sp.prof %>%
  dplyr::select(Sample, Cohort, Cluster, Age, Gender, BMI, Siblings, Medication, atopic_disease,  TotFib) %>%
  tibble::column_to_rownames(var = "Sample") %>%
  mutate(Cohort = factor(Cohort,levels = c("Control","ASD")),
         Gender = as.factor(Gender),
         BMI = as.numeric(BMI),
         Medication = as.factor(ifelse(is.na(Medication), "No", as.character(Medication))),
         atopic_disease =  as.factor(ifelse(is.na(atopic_disease),"No",atopic_disease))) %>%
  filter(Cluster %in% c("Control",  i))

abundance.sp.tmp<- sp.prof %>% 
  dplyr::select(Sample,all_of(pick.species)) %>% 
  filter(Sample %in% rownames(metadata.sp.tmp)) %>%
  tibble::column_to_rownames(var="Sample") 


###cohort
fit_data = Maaslin2(input_data     = abundance.sp.tmp,
                    input_metadata = metadata.sp.tmp,
                    min_prevalence = 0.05,
                    normalization  = "none",
                    analysis_method = "LM",
                    transform = "LOG",
                    output         = paste(i,"/taxa",sep=""),
                    fixed_effects  = c("Cluster","Age","Medication","atopic_disease","TotFib"),
                    #fixed_effects  = c("Cohort","Age","Medication","TotFib"),
                    reference      = c("Cluster,Control","Medication,No","atopic_disease,No")
)
}


