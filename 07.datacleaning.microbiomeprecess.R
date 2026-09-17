rm(list=ls())

library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(ggplot2)
library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(corrr)
library(ggstatsplot)
library(ggpubr)
library(Hmisc)
library(stringr)
library(xlsx)

species.batch.cluster <- read.table("species.profile.adjust_bacth_cluster.txt")
species <- species.batch.cluster %>%  t() %>% as.data.frame() %>%  rownames_to_column(var = "Subject.ID")

sp.tree <-readRDS("sp.tree.rds")

species.trans <- species.batch.cluster %>%
  rownames_to_column(var = "sp.name") %>%
  #mutate(sp.name = gsub("^s__","",sp.name)) %>%
  merge(sp.tree,.,by = "sp.name",all.x = T) %>%
  dplyr::select(-sp.name)

sample.id <- colnames(species.trans)[-c(1:2)]

metadata.select.updata.diet <- read.xlsx2("metadata_diet_per_cal.xlsx",sheetIndex = 1)

##### MicrobiotaProcess
sp.otu <- species.trans %>% dplyr::select(-c(clade_name,clade_taxid)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(OTU = paste("OTU",OTU,sep="")) %>%
  column_to_rownames(var = "OTU") %>%
  mutate_all(~ as.numeric(.))
sp.tax.info <- species.trans %>%
  dplyr::select(clade_name) %>%
  lapply(function(x) strsplit(as.character(x), split = "\\|")) 
sp.tax <- do.call(rbind,sp.tax.info$clade_name) %>% as.data.frame() %>%
  setNames(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  mutate_all(~ gsub(".*__", "", .)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(OTU = paste("OTU",OTU,sep="")) %>%
  column_to_rownames(var = "OTU") %>%
  as.matrix()

metadata.select.updata.diet <- metadata.select.updata.diet %>% remove_rownames() %>%
  column_to_rownames("sample_id") 

OTU = otu_table(sp.otu, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(sp.tax)
samples = sample_data(metadata.select.updata.diet)

phy.data <- phyloseq(OTU, TAX, samples)
phy.data
phy.mpse <- phy.data %>% as.MPSE() 

###alpha plot
phy.mpse %<>%
  mp_cal_alpha(.abundance=Abundance,force=TRUE )


mycomp <- list(c("ASD", "Control"))
alpha_2g <- phy.mpse %>%
  mp_plot_alpha(
    .group=Cohort,
    .alpha=c(Observe, Shannon, Simpson),
    #comparisons=mycomp,map_signif_level=T
  ) +
  scale_fill_manual(values=c( "#467897","#E7CD79"),guide="none") +
  scale_color_manual(values=c("#467897","#E7CD79"), guide="none")

pdf("supFig1a.alpha_asd_control.pdf",width = 5,height = 3)
alpha_2g
dev.off()

###beta diversity
phy.mpse %<>%
  mp_decostand(.abundance=Abundance)

p.pca <- phy.mpse %>%
  mp_cal_pca(.abundance=hellinger, action="add") %>%
  mp_plot_ord(
    .ord = PCA,
    .group = Cohort,
    .color = Cohort,
    .size = 1.2,
    .alpha = 1,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse
  ) +
  scale_fill_manual(values=c("#467897","#E7CD79"),guide="none") +
  scale_color_manual(values=c("#467897","#E7CD79"), guide="none")

pdf("supFig1b.pca.2groups.pdf",height = 5,width = 5)
ggarrange(p.pca)
dev.off()

###taxonomy
phy.mpse %<>%
  mp_cal_abundance(.abundance=Abundance, add=TRUE, force=TRUE, relative=TRUE)
sp.prof.long <- phy.mpse %>%
  mp_extract_abundance(taxa.class=Species) %>%
  tidyr::unnest(cols=AbundanceBySample) %>% dplyr::rename(species="label")
sp.prof <- sp.prof.long[,-c(2,4)] %>% pivot_wider(names_from = species, values_from = RelAbundanceBySample)

zero.num <-apply(sp.prof[,-c(1:114)], 2, function(x){sum(x!=0)})
pick.species <- names(which(zero.num/dim(sp.prof)[1] > 0.05))

saveRDS(metadata.select.updata.diet,"metadata.carb.rds")
save.image("Manuscript.RData")

