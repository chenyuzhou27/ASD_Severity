rm(list=ls())
###pathway
library(dplyr)
library(tibble)
library(Maaslin2)
library(tibble)

# Filter the data
metadata <- read.table("metadata.rename.questionnaire.txt",sep="\t") 
pathway <- read.table("pathway.profile.adjust_bacth_cluster.txt",sep="\t")
pathway <- t(pathway) %>% as.data.frame() %>% rownames_to_column( var = "Subject.ID")

zero.num <-apply(pathway, 2, function(x){sum(x!=0)})
pick.pathway <- names(which(zero.num/dim(pathway)[1] > 0.05))

pid <- pmatch(pathway$Subject.ID,metadata$sample_id)
all(pid == c(1:871))

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
  
  abundance.sp.tmp<-  pathway %>% remove_rownames() %>% 
    dplyr::select(Subject.ID, all_of(pick.pathway)) %>%
    filter(Subject.ID %in% rownames(metadata.sp.tmp)) %>%
    remove_rownames() %>%
    tibble::column_to_rownames(var = "Subject.ID")
    
  
  ###cohort
  fit_data = Maaslin2(input_data     = abundance.sp.tmp,
                      input_metadata = metadata.sp.tmp,
                      min_prevalence = 0.05,
                      normalization  = "none",
                      analysis_method = "LM",
                      transform = "LOG",
                      output         = paste(i,"/pathway",sep=""),
                      fixed_effects  = c("Cluster","Age","Medication","atopic_disease","TotFib"),
                      #fixed_effects  = c("Cohort","Age","Medication","TotFib"),
                      reference      = c("Cluster,Control","Medication,No","atopic_disease,No")
  )
}


