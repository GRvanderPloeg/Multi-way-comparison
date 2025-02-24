library(tidyverse)
library(vegan)
library(ggpubr)
library(ape)
library(Polychrome)
library(parafac4microbiome)

# Tongue
tongue = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/tongueCounts.csv", header=FALSE) %>% as_tibble()
tongue_taxonomy = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/taxonomyTongue_fixed.csv", sep=" ") %>% as_tibble()
tongue_sampleMeta = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/tongueSampleMeta.csv",header=FALSE) %>% as_tibble()

colnames(tongue) = tongue_taxonomy$zOTU
temp = read.csv("./Data/GOHTRANS/sampleInfo_fixed.csv", sep=" ") %>% as_tibble()
colnames(tongue_sampleMeta) = c(temp %>% select(-Description,-subject,-newTimepoint) %>% colnames, "Description", "subject", "newTimepoint")

# Load other metadata and ph
ph_BOMP = read_delim("./Data/GOHTRANS/GOH-TRANS_csv_export_20240205114955/GOH-TRANS_export_20240205.csv",
                     delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% as_tibble()

df1 = ph_BOMP %>% select(`Participant Id`, starts_with("5.")) %>% mutate(subject = 1:42, numTeeth = `5.1|Number of teeth`, DMFT = `5.2|DMFT`, numBleedingSites = `5.3|Bleeding sites`, boppercent = `5.4|BOP%`, DPSI = `5.5|DPSI`, pH = `5.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df2 = ph_BOMP %>% select(`Participant Id`, starts_with("12.")) %>% mutate(subject = 1:42, numTeeth = `12.1|Number of teeth`, DMFT = `12.2|DMFT`, numBleedingSites = `12.3|Bleeding sites`, boppercent = `12.4|BOP%`, DPSI = `12.5|DPSI`, pH = `12.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df3 = ph_BOMP %>% select(`Participant Id`, starts_with("19.")) %>% mutate(subject = 1:42, numTeeth = `19.1|Number of teeth`, DMFT = `19.2|DMFT`, numBleedingSites = `19.3|Bleeding sites`, boppercent = `19.4|BOP%`, DPSI = `19.5|DPSI`, pH = `19.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)
df4 = ph_BOMP %>% select(`Participant Id`, starts_with("26.")) %>% mutate(subject = 1:42, numTeeth = `26.1|Number of teeth`, DMFT = `26.2|DMFT`, numBleedingSites = `26.3|Bleeding sites`, boppercent = `26.4|BOP%`, DPSI = `26.5|DPSI`, pH = `26.8|pH`) %>% select(subject, numTeeth, DMFT, numBleedingSites, boppercent, DPSI, pH)

otherMeta = rbind(df1, df2, df3, df4) %>% as_tibble() %>% mutate(newTimepoint = rep(c(0,3,6,12), each=42))

# Remove pH measurements of lower than 0
otherMeta = otherMeta[!otherMeta$pH < 1,] %>% as_tibble()

# Put into cube
I = tongue_sampleMeta$subject %>% unique() %>% length()
J = ncol(tongue)
K = 4
timepoints = c(0, 3, 6, 12)
tongueCube = array(0L, dim=c(I,J,K))

for(k in 1:K){
  temp = cbind(tongue, tongue_sampleMeta) %>% as_tibble()
  tongueCube[,,k] = temp %>%
    filter(newTimepoint == timepoints[k]) %>%
    right_join(tongue_sampleMeta %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(tongue_sampleMeta))) %>%
    as.matrix()
}

tongueCube_mode1 = tongue_sampleMeta %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
tongueCube_mode2 = tongue_taxonomy
tongueCube_mode3 = tongue_sampleMeta %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()

tongueData = list("data"=tongueCube, "mode1"=tongueCube_mode1, "mode2"=tongueCube_mode2, "mode3"=tongueCube_mode3)

# Saliva
saliva = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/salivaCounts.csv", header=FALSE) %>% as_tibble()
saliva_taxonomy = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/taxonomysaliva_fixed.csv", sep=" ") %>% as_tibble()
saliva_sampleMeta = read.csv("./Data/GOHTRANS/20240503_UNOISE_new/salivaSampleMeta.csv",header=FALSE) %>% as_tibble()

colnames(saliva) = saliva_taxonomy$zOTU
temp = read.csv("./Data/GOHTRANS/sampleInfo_fixed.csv", sep=" ") %>% as_tibble()
colnames(saliva_sampleMeta) = c(temp %>% select(-Description,-subject,-newTimepoint) %>% colnames, "Description", "subject", "newTimepoint")

# Put into cube
I = saliva_sampleMeta$subject %>% unique() %>% length()
J = ncol(saliva)
K = 4
timepoints = c(0, 3, 6, 12)
salivaCube = array(0L, dim=c(I,J,K))

for(k in 1:K){
  temp = cbind(saliva, saliva_sampleMeta) %>% as_tibble()
  salivaCube[,,k] = temp %>%
    filter(newTimepoint == timepoints[k]) %>%
    right_join(saliva_sampleMeta %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(saliva_sampleMeta))) %>%
    as.matrix()
}

salivaCube_mode1 = saliva_sampleMeta %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
salivaCube_mode2 = saliva_taxonomy
salivaCube_mode3 = saliva_sampleMeta %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()

salivaData = list("data"=salivaCube, "mode1"=salivaCube_mode1, "mode2"=salivaCube_mode2, "mode3"=salivaCube_mode3)

# Process microbiome data
processedTongue = processDataCube(tongueData, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)
processedSaliva = processDataCube(salivaData, sparsityThreshold = 0.5, considerGroups=TRUE, groupVariable="GenderID", CLR=TRUE, centerMode=1, scaleMode=2)

saveRDS(processedTongue, "./Data/GOHTRANS/tongue.RDS")
saveRDS(processedSaliva, "./Data/GOHTRANS/saliva.RDS")

# Cytokines
df = read.csv("./Data/GOHTRANS/20241209_cytokines.csv", header=FALSE, sep=" ") %>% as_tibble()
featureMeta = read.csv("./Data/GOHTRANS/20241209_cytokines_featureMeta.csv", header=FALSE) %>% as_tibble()
sampleInfo = read.csv("./Data/GOHTRANS/20241209_cytokines_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("subject", "GenderID", "newTimepoint", "unknown", "unknown2")

temp = sampleInfo %>% select(subject, GenderID) %>% unique()
Y = as.numeric(as.factor(temp$GenderID))
Ycnt = Y - mean(Y)

# Put into cube
I = 27
J = 22
K = 3
timepoints = c(0, 3, 12)
cytokineCube = array(0L, dim=c(I,J,K))

for(k in 1:K){
  temp = cbind(df, sampleInfo) %>% as_tibble()
  cytokineCube[,,k] = temp %>%
    filter(newTimepoint == timepoints[k]) %>%
    right_join(sampleInfo %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(sampleInfo))) %>%
    as.matrix()
}

# Log transform
pseudocount = 0.1300
cytokineCube_log = log(cytokineCube + pseudocount)

# Center and scale
cytokineCube_cnt = multiwayCenter(cytokineCube_log, mode=1)
cytokineCube_cnt_scl = multiwayScale(cytokineCube_cnt, mode=2)

# Prep metadata
cytokineCube_mode1 = sampleInfo %>% select(subject, GenderID) %>% arrange(subject) %>% unique()
cytokineCube_mode2 = featureMeta
cytokineCube_mode3 = sampleInfo %>% filter(newTimepoint %in% timepoints) %>% select(newTimepoint) %>% unique()

cytokineData = list("data"=cytokineCube_cnt_scl, "mode1"=cytokineCube_mode1, "mode2"=cytokineCube_mode2, "mode3"=cytokineCube_mode3)
saveRDS(cytokineData, "./Data/GOHTRANS/cytokine.RDS")
