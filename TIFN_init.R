library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

# Load count data
raw_data = read.csv("./Data/TIFN/20221005_wp2/count-table.tsv", sep="\t") %>% as_tibble()
metadata = raw_data %>% select(sample, subject, visit, group, niche)
counts   = raw_data %>% select(-sample, -subject, -visit, -group, -niche)
taxa = read.csv("./Data/TIFN/20221005_wp2/taxonomic-classification.tsv", sep="\t") %>% as_tibble()

# Load RF data
rf_data = read.csv("./Data/TIFN/RFdata.csv")
colnames(rf_data) = c("subject", "id", "fotonr", "day", "group", "RFgroup", "MQH", "SPS(tm)", "Area_delta_R30", "Area_delta_Rmax", "Area_delta_R30_x_Rmax", "gingiva_mean_R_over_G", "gingiva_mean_R_over_G_upper_jaw", "gingiva_mean_R_over_G_lower_jaw")
rf_data = rf_data %>% as_tibble()

rf_data[rf_data$subject == "VSTPHZ", 1] = "VSTPH2"
rf_data[rf_data$subject == "D2VZH0", 1] = "DZVZH0"
rf_data[rf_data$subject == "DLODNN", 1] = "DLODDN"
rf_data[rf_data$subject == "O3VQFX", 1] = "O3VQFQ"
rf_data[rf_data$subject == "F80LGT", 1] = "F80LGF"
rf_data[rf_data$subject == "26QQR0", 1] = "26QQrO"

# rf_data2 = read.csv("./data-raw/vanderPloeg2024_RFdata2.csv") %>% as_tibble()
# rf_data2 = rf_data2[,c(2,4,181:192)]
# rf_data = rf_data %>% left_join(rf_data2)

rf = rf_data %>% select(subject, RFgroup) %>% unique()

# Attach RF data to metadata
metadata = metadata %>% left_join(rf)

# Remove test
mask = metadata$group == "control"
counts = counts[mask,]
metadata = metadata[mask,]

# Prepare export of metadata
mode1 = metadata %>% select(subject, RFgroup) %>% unique() %>% arrange(subject)
mode2 = taxa
mode3 = metadata %>% select(visit) %>% unique() %>% arrange(visit) %>% mutate(status=c("Baseline", "EG", "EG", "EG", "EG", "EG", "Resolution"))

# Tongue
tongueMask = metadata$niche == "tongue"
df_tongue = counts[tongueMask,]
metadata_tongue = metadata[tongueMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_tongue, metadata_tongue) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

tongue = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Lowling
lowlingMask = metadata$niche == "lower jaw, lingual"
df_lowling = counts[lowlingMask,]
metadata_lowling = metadata[lowlingMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_lowling, metadata_lowling) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

lowling = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Lowinter
lowinterMask = metadata$niche == "lower jaw, interproximal"
df_lowinter = counts[lowinterMask,]
metadata_lowinter = metadata[lowinterMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_lowinter, metadata_lowinter) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

lowinter = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Upling
uplingMask = metadata$niche == "upper jaw, lingual"
df_upling = counts[uplingMask,]
metadata_upling = metadata[uplingMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_upling, metadata_upling) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

upling = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Upinter
upinterMask = metadata$niche == "upper jaw, interproximal"
df_upinter = counts[upinterMask,]
metadata_upinter = metadata[upinterMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_upinter, metadata_upinter) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

upinter = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Saliva
salivaMask = metadata$niche == "saliva"
df_saliva = counts[salivaMask,]
metadata_saliva = metadata[salivaMask,]

I = length(unique(metadata$subject))
J = ncol(counts)
K = max(metadata$visit)
X = array(0L, c(I,J,K))

for(k in 1:K){
  X[,,k] = cbind(df_saliva, metadata_saliva) %>%
    as_tibble() %>%
    filter(visit == k) %>%
    right_join(metadata %>% select(subject) %>% unique()) %>%
    arrange(subject) %>%
    select(-all_of(colnames(metadata))) %>%
    as.matrix()
}

saliva = list("data"=X, "mode1"=mode1, "mode2"=mode2, "mode3"=mode3)

# Process
processedTongue = processDataCube(tongue, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)
processedLowling = processDataCube(lowling, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)
processedLowinter = processDataCube(lowinter, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)
processedUpling = processDataCube(upling, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)
processedUpinter = processDataCube(upinter, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)
processedSaliva = processDataCube(saliva, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Save
saveRDS(processedTongue, "./Data/TIFN/tongue.RDS")
saveRDS(processedLowling, "./Data/TIFN/lowling.RDS")
saveRDS(processedLowinter, "./Data/TIFN/lowinter.RDS")
saveRDS(processedUpling, "./Data/TIFN/upling.RDS")
saveRDS(processedUpinter, "./Data/TIFN/upinter.RDS")
saveRDS(processedSaliva, "./Data/TIFN/saliva.RDS")
