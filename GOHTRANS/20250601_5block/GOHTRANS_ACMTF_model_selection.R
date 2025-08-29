library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)
library(vegan)

homogenizedSamples1 = intersect(NPLStoolbox::Cornejo2025$Tongue_microbiome$mode1 %>% filter(!is.na(GenderID)) %>% select(subject) %>% pull(), NPLStoolbox::Cornejo2025$Salivary_microbiome$mode1 %>% filter(!is.na(GenderID)) %>% select(subject) %>% pull())
homogenizedSamples2 = intersect(NPLStoolbox::Cornejo2025$Salivary_cytokines$mode1 %>% filter(!is.na(GenderID)) %>% select(subject) %>% pull(), NPLStoolbox::Cornejo2025$Salivary_biochemistry$mode1 %>% filter(!is.na(GenderID)) %>% select(subject) %>% pull())
homogenizedSamples3 = NPLStoolbox::Cornejo2025$Blood_hormones$mode1 %>% filter(!is.na(GenderID)) %>% select(subject) %>% pull()

homogenizedSamples = intersect(intersect(homogenizedSamples1, homogenizedSamples2), homogenizedSamples3)

# Tongue
newTongue = NPLStoolbox::Cornejo2025$Tongue

mask = newTongue$mode1$subject %in% homogenizedSamples
newTongue$data = newTongue$data[mask,,]
newTongue$mode1 = newTongue$mode1[mask,]

processedTongue = processDataCube(newTongue, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)

# Saliva
newSaliva = NPLStoolbox::Cornejo2025$Salivary_microbiome

mask = newSaliva$mode1$subject %in% homogenizedSamples
newSaliva$data = newSaliva$data[mask,,]
newSaliva$mode1 = newSaliva$mode1[mask,]

processedSaliva = processDataCube(newSaliva, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)

# Cytokines
newCytokines = NPLStoolbox::Cornejo2025$Salivary_cytokines

mask = newCytokines$mode1$subject %in% homogenizedSamples
newCytokines$data = log(newCytokines$data[mask,,] + 0.1300)
newCytokines$mode1 = newCytokines$mode1[mask,]

processedCytokines = processDataCube(newCytokines, CLR=FALSE, centerMode=1, scaleMode=2)

# Biochemistry
newBio = NPLStoolbox::Cornejo2025$Salivary_biochemistry

mask = newBio$mode1$subject %in% homogenizedSamples
newBio$data = log(newBio$data[mask,,])
newBio$mode1 = newBio$mode1[mask,]

processedBiochemistry = processDataCube(newBio, CLR=FALSE, centerMode=1, scaleMode=2)

# Hormones
newHormones = NPLStoolbox::Cornejo2025$Blood_hormones

mask = newHormones$mode1$subject %in% homogenizedSamples
newHormones$data = log(newHormones$data[mask,,])
newHormones$mode1 = newHormones$mode1[mask,]

processedHormones = processDataCube(newHormones, CLR=FALSE, centerMode=1, scaleMode=2)

# Prep data
datasets = list(processedTongue$data, processedSaliva$data, processedCytokines$data, processedBiochemistry$data, processedHormones$data)
modes = list(c(1,2,3),c(1,4,5),c(1,6,7),c(1,8,9),c(1,10,11))
Z = setupCMTFdata(datasets, modes, normalize=TRUE)

CV_ACMTF = ACMTF_modelSelection(datasets, modes, maxNumComponents=10, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(CV_ACMTF, "./GOHTRANS_CV_ACMTF.RDS")
