library("MMUPHin")
library(dplyr)
library(tibble)
library(clusterSim)
library(readxl)
sp.prof <- read_xlsx("species.profile.txt",sheet = 1)
metadata <- read.table("metadata.txt",sep="\t",head=T)

abundance.sp <- sp.prof %>% 
  column_to_rownames(var = "Subject.ID") %>% 
  dplyr::select(-Cluster) %>%
  t() %>% as.data.frame() 

abundance.sp <- abundance.sp / colSums(abundance.sp)
colSums(abundance.sp) %>%  round(0) %>% sum()

metadata.batch <- metadata %>%
  dplyr::select("Subject.ID","Batch","Cohort", "Cluster") %>%
  column_to_rownames(var = "Subject.ID") %>%
  filter(!is.na(Batch)) 

abundance.sp.batch <- abundance.sp[,rownames(metadata.batch)]

abundance.sp.adj.cluster <- adjust_batch(feature_abd = abundance.sp.batch,
                                         batch = "Batch",
                                         covariates = "Cluster",
                                         data = metadata.batch)$feature_abd_adj

write.table(abundance.sp.adj.cluster,"species.profile.adjust_bacth_cluster.txt",sep="\t",quote = F)

#########pathway batch effect
pathway <- read.table("Pathway.profile.txt",sep="\t")
all.equal(rowSums(pathway[,-c(1:2)]), 100)

pathway.t <- pathway %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Subject.ID") %>% 
  dplyr::select(-Cluster) %>% 
  t() %>%
  as.data.frame() %>%
  dplyr::select(all_of(rownames(metadata.batch)))

pathway.t <- pathway.t / colSums(pathway.t)
colSums(pathway.t) %>%  round(0) %>% sum()

abundance.pathway.adj.cluster <- adjust_batch(feature_abd = pathway.t,
                                              batch = "Batch",
                                              covariates = "Cluster",
                                              data = metadata.batch)$feature_abd_adj

write.table(abundance.pathway.adj.cluster,"pathway.profile.adjust_bacth_cluster.txt",sep="\t",quote = F)

####ko batch effect
ko.prof <- read.csv("ko_filter.csv")

ko.prof <- ko.prof %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "X..Gene.Family") %>% 
  t() 
dim(ko.prof)

ko.prof.row.0 <- apply(ko.prof, 1, function(x){sum(x==0)})
sum(ko.prof.row.0==8092)

ko.prof.tmp <- apply(ko.prof[,-c(1:2)], 1, function(row) row / sum(row))
colSums(ko.prof.tmp) [which(colSums(ko.prof.tmp) != 1)]


pid <- pmatch(rownames(metadata.batch),colnames(ko.prof.tmp))
ko.prof.tmp <- ko.prof.tmp[,pid]


abundance.ko.adj.cluster <- adjust_batch(feature_abd = ko.prof.tmp,
                                         batch = "Batch",
                                         covariates = "Cluster",
                                         data = metadata.batch)$feature_abd_adj

write.table(abundance.ko.adj.cluster,"ko.profile.adjust_bacth_cluster.txt",sep="\t",quote = F)
