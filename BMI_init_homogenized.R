library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

# Load subject info
infant_anthropometrics = read.csv("./Data/BMI/infant_anthropometrics.csv") %>% as_tibble()

# Faeces subject info
sampleInfo = read.csv("./Data/BMI/faeces_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
sampleInfo = sampleInfo %>% left_join(infant_anthropometrics %>% select(RCID, whz.6m))
subjectMeta_faeces = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis, whz.6m) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_faeces = subjectMeta_faeces[-90,]

# Milk subject info
sampleInfo = read.csv("./Data/BMI/milk_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_milk = subjectMeta %>% left_join(subjectMeta_faeces)

# Milk metab subject info
sampleInfo = read.csv("./Data/BMI/milkMetab_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)
subjectMeta_milkMetab = subjectMeta %>% left_join(subjectMeta_faeces %>% select(subject, whz.6m))

# Identify shared subjects
sharedSubjects = intersect(intersect(subjectMeta_faeces$subject, subjectMeta_milk$subject), subjectMeta_milkMetab$subject)
homogenized_subjectMeta_bmi = subjectMeta_faeces %>% filter(subject %in% sharedSubjects) %>% arrange(subject)
homogenized_subjectMeta_whz = homogenized_subjectMeta_bmi %>% filter(whz.6m != "NA")

# Faecal microbiome
df = read.csv("./Data/BMI/faecesCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./Data/BMI/newTaxonomy_faeces.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./Data/BMI/faeces_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# Filter based on sparsity
threshold = 0.75
sparsity = colSums(df==0) / nrow(df)
featureSelection = sparsity <= threshold

# CLR
df_clr = t(apply(df+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()

# Apply feature filter
df_clr = df_clr[,featureSelection]
taxonomy_filtered = taxonomy[featureSelection,]

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_clr)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_clr, sampleInfo) %>%
    as_tibble() %>%
    mutate(subject=as.character(subject)) %>%
    filter(Days == Day) %>%
    select(c(colnames(df_clr),subject)) %>%
    right_join(subjectMeta) %>%
    arrange(subject) %>%
    select(-colnames(subjectMeta)) %>%
    as.matrix()
}

# Mask based on shared subjects for BMI and WHZ
X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]

# Center and scale
X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)

X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)

# Save
faeces_df_bmi = X_bmi_cnt_scl
faeces_df_whz = X_whz_cnt_scl
faeces_subjectMeta = subjectMeta
faeces_taxonomy = taxonomy_filtered
faeces_timepoints = timepoints

# Milk microbiome
df = read.csv("./Data/BMI/milkCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./Data/BMI/newTaxonomy_milk.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./Data/BMI/milk_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# Filter based on sparsity
threshold = 0.85
sparsity = colSums(df==0) / nrow(df)
featureSelection = sparsity <= threshold

# CLR
df_clr = t(apply(df+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()

# Apply feature filter
df_clr = df_clr[,featureSelection]
taxonomy_filtered = taxonomy[featureSelection,]

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_clr)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_clr, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_clr),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
}

# Mask based on shared subjects for BMI and WHZ
X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]

# Center and scale
X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)

X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)

milk_df_bmi = X_bmi_cnt_scl
milk_df_whz = X_whz_cnt_scl
milk_subjectMeta = subjectMeta
milk_taxonomy = taxonomy_filtered
milk_timepoints = timepoints

# Milk metabolomics
df = read.csv("./Data/BMI/milkMetabNumeric.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./Data/BMI/milk_metab_CAS_numbers.csv", header=TRUE, sep=",") %>% as_tibble()
sampleInfo = read.csv("./Data/BMI/milkMetab_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")
metaboliteCategories = read.csv("./Data/BMI/Metabolite_categories.txt", sep="\t") %>% as_tibble()
metaboliteCategories = metaboliteCategories %>% mutate(Metabolite = vctrs::vec_as_names(Metabolite, repair="universal_quiet"))

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Remove duplicates
drop = c(61,81,146)
subjectMeta = subjectMeta %>% mutate(index=1:nrow(.)) %>% filter(!index %in% drop) %>% select(-index) %>% arrange(subject)

# Log transform
df_log = log(df)

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_log)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_log, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_log),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
}

# Mask based on shared subjects for BMI and WHZ
X_bmi = X[subjectMeta$subject %in% homogenized_subjectMeta_bmi$subject,,]
X_whz = X[subjectMeta$subject %in% homogenized_subjectMeta_whz$subject,,]

# Center and scale
X_bmi_cnt = parafac4microbiome::multiwayCenter(X_bmi, mode=1)
X_bmi_cnt_scl = parafac4microbiome::multiwayScale(X_bmi_cnt, mode=2)

X_whz_cnt = parafac4microbiome::multiwayCenter(X_whz, mode=1)
X_whz_cnt_scl = parafac4microbiome::multiwayScale(X_whz_cnt, mode=2)

milkMetab_df_bmi = X_bmi_cnt_scl
milkMetab_df_whz = X_whz_cnt_scl
milkMetab_subjectMeta = subjectMeta
milkMetab_taxonomy = taxonomy
milkMetab_timepoints = timepoints

# Fix milkMetab feature annotations
df = milkMetab_taxonomy %>% mutate(Metabolite = make.names(X), CAS.Registry = make.names(CAS.Registry)) %>% left_join(metaboliteCategories)

# Fix mismatches by hand
df[df$Metabolite == "X2.Aminobutyrate","Class"] = "Amino acids and derivatives"
df[df$Metabolite == "X2.Fucosyllactose","Class"] = "Oligosaccharides"
df[df$Metabolite == "X2.Hydroxybutyrate","Class"] = "Amino acids and derivatives"
df[df$Metabolite == "X2.Oxoglutarate","Class"] = "Energy related"
df[df$Metabolite == "X3.Fucosyllactose","Class"] = "Oligosaccharides"
df[df$Metabolite == "X3SL.partial","Class"] = "Oligosaccharides"
df[df$Metabolite == "X6SL.partial","Class"] = "Oligosaccharides"
df[df$Metabolite == "Methionine","Class"] = "Amino acids and derivatives"
df[70,"Metabolite"] = "tau.Methylhistidine" # Fix non-ascii tau character

milkMetab_featureMeta = df %>% select(-X)

# Save
saveRDS(milkMetab_df_bmi, "./Data/BMI/milkMetab_homogenized.RDS")
saveRDS(milk_df_bmi, "./Data/BMI/milk_homogenized.RDS")
saveRDS(faeces_df_bmi, "./Data/BMI/faeces_homogenized.RDS")
saveRDS(faeces_subjectMeta, "./Data/BMI/homogenized_subjectMeta.RDS")
saveRDS(faeces_taxonomy, "./Data/BMI/faeces_taxonomy.RDS")
saveRDS(milk_taxonomy, "./Data/BMI/milk_taxonomy.RDS")
saveRDS(milkMetab_featureMeta, "./Data/BMI/milk_metabolites.RDS")
