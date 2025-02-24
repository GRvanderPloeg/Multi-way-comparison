library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(ggplot2)
library(parafac4microbiome)
library(CMTFtoolbox)

# Load data
cytokines_data = read.csv("./Data/AP/input_deduplicated_RvdP.csv", sep=" ", header=FALSE) %>% as_tibble()
colnames(cytokines_data) = c("VEGF", "CRP", "GM-CSF", "IL1alpha", "IL1beta", "IL4", "IL6", "IL8", "IL10", "IL12p70", "IL17A", "IFNgamma", "MIP1alpha", "OPG", "TNFalpha", "RANKL")

cytokines_meta_data = read.csv("./Data/AP/input_deduplicated_metadata_RvdP.csv", sep=" ", header=FALSE) %>% as_tibble()
colnames(cytokines_meta_data) = c("SubjectID", "Visit", "Gender", "Age", "Pain_noPain", "case_control")

otherMeta = read.csv("./Data/AP/Root_meta_data_parafac.txt", sep="\t") %>% as_tibble() %>% select(-Gender)

# Put into cube
I = 52
J = 16
K = 6
cube = array(0L, dim=c(I,J,K))

for(k in 1:K){
  temp = cbind(cytokines_data, cytokines_meta_data) %>% as_tibble()
  cube[,,k] = temp %>%
    filter(Visit == k) %>%
    right_join(cytokines_meta_data %>% select(SubjectID) %>% unique()) %>%
    arrange(SubjectID) %>%
    select(-all_of(colnames(cytokines_meta_data))) %>%
    as.matrix()
}

# Take log + pseudovalue
cube = log(cube+0.0050)

# Remove outlier A11-8
temp = cytokines_meta_data %>% select(-Visit) %>% unique() %>% arrange(SubjectID)

# Remove extra samples
remove = c("A11-8") # first two are from paper, last is new
mask2 = !(processedGeorgiou$mode1$SubjectID %in% remove)
processedGeorgiou$data = processedGeorgiou$data[mask2,,]
processedGeorgiou$mode1 = processedGeorgiou$mode1[mask2,]

processedGeorgiou = list()
processedGeorgiou$mode1 = temp[mask2,] %>% left_join(otherMeta %>% select(SubjectID, PainS_NopainA) %>% unique())
processedGeorgiou$mode2 = colnames(cytokines_data) %>% as_tibble()
processedGeorgiou$mode3 = cytokines_meta_data %>% select(Visit) %>% arrange(Visit) %>% unique() %>% mutate(extraction = c(rep("Before extraction",3), rep("After extraction",3)))

processedGeorgiou$data = cube[mask2,,]

# Center and scale
processedGeorgiou$data = multiwayCenter(processedGeorgiou$data, mode=1)
processedGeorgiou$data = multiwayScale(processedGeorgiou$data, mode=2)

saveRDS(processedGeorgiou, "./Data/AP/cytokines_processed.RDS")
