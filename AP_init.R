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

# Filter to case only
temp = cytokines_meta_data %>% select(-Visit) %>% unique() %>% arrange(SubjectID)
caseMask = temp$case_control=="case"

processedGeorgiou = list()
processedGeorgiou$mode1 = temp[caseMask,] %>% left_join(otherMeta %>% select(SubjectID, PainS_NopainA) %>% unique())
processedGeorgiou$mode2 = colnames(cytokines_data) %>% as_tibble()
processedGeorgiou$mode3 = cytokines_meta_data %>% select(Visit) %>% arrange(Visit) %>% unique() %>% mutate(extraction = c(rep("Before extraction",3), rep("After extraction",3)))

processedGeorgiou$data = cube[caseMask,,]

# Remove extra samples
remove = c("A11-18", "A11-3", "A11-8") # first two are from paper, last is new
mask2 = !(processedGeorgiou$mode1$SubjectID %in% remove)
processedGeorgiou$data = processedGeorgiou$data[mask2,,]
processedGeorgiou$mode1 = processedGeorgiou$mode1[mask2,]

# Center and scale
processedGeorgiou$data = multiwayCenter(processedGeorgiou$data, mode=1)
processedGeorgiou$data = multiwayScale(processedGeorgiou$data, mode=2)

saveRDS(processedGeorgiou, "./Data/AP/cytokines_processed.RDS")

# Load microbiome data - these are case only by default
microbiome_raw = read.csv("./Data/AP/20240429_microbiome_counts.csv", sep=" ", header=FALSE)
taxonomy = read.csv("./Data/AP/20240429_taxonomy.csv", sep=" ", header=FALSE)
subjectMeta2 = read.csv("./Data/AP/20240429_microbiome_sampleMeta.csv", sep=" ", header=FALSE)

# Remove extra samples
remove = c("A11-18", "A11-3", "A11-8 36", "A11-10 17", "A11-15 17", "A11-8 46") # last is new
mask = !(subjectMeta2[,3] %in% remove)

microbiome_raw = microbiome_raw[mask,]
subjectMeta2 = subjectMeta2[mask,]

# Select ASVs based on sparsity per group
sparsityThreshold = 0.5

maskA = subjectMeta2[,4] == "A"
maskS = subjectMeta2[,4] == "S"

microbiomeA = microbiome_raw[maskA,]
microbiomeS = microbiome_raw[maskS,]

sparsityA = colSums(microbiomeA==0) / nrow(microbiomeA)
sparsityS = colSums(microbiomeS==0) / nrow(microbiomeS)

featureMask = (sparsityA <= sparsityThreshold) | (sparsityS <= sparsityThreshold)

# CLR transformation
geomeans = pracma::geomean(as.matrix(microbiome_raw+1), dim=2)
microbiome_clr = log(sweep(microbiome_raw+1, 1, geomeans, FUN="/"))

# Reduce to previously selected ASVs
microbiome_selected = microbiome_clr[,featureMask]
taxonomy_selected = taxonomy[featureMask,]

# Center and scale
microbiome_cnt = sweep(microbiome_selected, 2, colMeans(microbiome_selected), FUN="-")
microbiome_cnt_scl = sweep(microbiome_cnt, 2, apply(microbiome_cnt, 2, sd), FUN="/")

saveRDS(microbiome_cnt_scl, "./Data/AP/microbiome_processed.RDS")
saveRDS(taxonomy_selected, "./Data/AP/taxonomy_processed.RDS")

# all.equal(processedGeorgiou$mode1$SubjectID, subjectMeta2$V2)
# TRUE
