library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)
library(vegan)

# Import data
rf_data = read.csv("./RFdata.csv")
colnames(rf_data) = c("subject", "id", "fotonr", "day", "group", "RFgroup", "MQH", "SPS(tm)", "Area_delta_R30", "Area_delta_Rmax", "Area_delta_R30_x_Rmax", "gingiva_mean_R_over_G", "gingiva_mean_R_over_G_upper_jaw", "gingiva_mean_R_over_G_lower_jaw")
rf_data = rf_data %>% as_tibble()

rf_data[rf_data$subject == "VSTPHZ", 1] = "VSTPH2"
rf_data[rf_data$subject == "D2VZH0", 1] = "DZVZH0"
rf_data[rf_data$subject == "DLODNN", 1] = "DLODDN"
rf_data[rf_data$subject == "O3VQFX", 1] = "O3VQFQ"
rf_data[rf_data$subject == "F80LGT", 1] = "F80LGF"
rf_data[rf_data$subject == "26QQR0", 1] = "26QQrO"

rf_data2 = read.csv("./red_fluorescence_data.csv") %>% as_tibble()
rf_data2 = rf_data2[,c(2,4,181:192)]
rf_data = rf_data %>% left_join(rf_data2)

rf = rf_data %>% select(subject, RFgroup) %>% unique()

age_gender = read.csv("./Ploeg_subjectMetadata.csv", sep=";")
age_gender = age_gender[2:nrow(age_gender),2:ncol(age_gender)]
age_gender = age_gender %>% as_tibble() %>% filter(onderzoeksgroep == 0) %>% select(naam, leeftijd, geslacht)
colnames(age_gender) = c("subject", "age", "gender")

# Correction for incorrect subject ids
age_gender[age_gender$subject == "VSTPHZ", 1] = "VSTPH2"
age_gender[age_gender$subject == "D2VZH0", 1] = "DZVZH0"
age_gender[age_gender$subject == "DLODNN", 1] = "DLODDN"
age_gender[age_gender$subject == "O3VQFX", 1] = "O3VQFQ"
age_gender[age_gender$subject == "F80LGT", 1] = "F80LGF"
age_gender[age_gender$subject == "26QQR0", 1] = "26QQrO"

age_gender = age_gender %>% arrange(subject)


# Functions

mapping = c(-14,0,2,5,9,14,21)

testMetadata = function(model, obj, componentNum){

  normalSubjectLoadings = cbind(model$Fac[[1]][,componentNum], obj$mode1) %>% as_tibble()
  colnames(normalSubjectLoadings) = c("Loading", "subject", "RFgroup")

  transformedSubjectLoadings = transformPARAFACloadings(model$Fac, 2, moreOutput=TRUE)$Ftilde[,componentNum] %>% as_tibble() %>% mutate(subject = rep(obj$mode1$subject, each=nrow(obj$mode3)), day = rep(mapping[obj$mode3$visit],nrow(obj$mode1)))
  colnames(transformedSubjectLoadings) = c("Loading", "subject", "day")

  uncorrectedP = rep(0, 5)

  # Plaque%
  temp=transformedSubjectLoadings %>% left_join(rf_data %>% select(subject, day, plaquepercent),by=c("subject","day"))
  uncorrectedP[1] = cor.test(temp$Loading, temp$plaquepercent, method="pearson")$p.value

  # Bleeding%
  temp=transformedSubjectLoadings %>% left_join(rf_data %>% select(subject, day, bomppercent), by=c("subject","day"))
  uncorrectedP[2] = cor.test(temp$Loading, temp$bomppercent, method="pearson")$p.value

  # RF%
  temp=transformedSubjectLoadings %>% left_join(rf_data %>% select(subject, day, Area_delta_R30), by=c("subject","day"))
  uncorrectedP[3] = cor.test(temp$Loading, temp$Area_delta_R30, method="pearson")$p.value

  # Gender
  partA = normalSubjectLoadings %>% left_join(age_gender, by="subject") %>% select(1,2,gender) %>% filter(gender==1)
  partB = normalSubjectLoadings %>% left_join(age_gender, by="subject") %>% select(1,2,gender) %>% filter(gender==2)

  uncorrectedP[4] = wilcox.test(partA$Loading, partB$Loading)$p.value # Not what I would do nowadays, but this is what is in the TIFN paper

  # Age
  temp = normalSubjectLoadings %>% left_join(age_gender, by="subject") %>% select(1,2,age)
  uncorrectedP[5] = cor.test(temp$Loading, temp$age, method="pearson")$p.value

  return(uncorrectedP)
}

testMetadata2 = function(model, mode1metadata){
  transformedSubjectLoadings = transformPARAFACloadings(model$Fac, 1)

  metadata = mode1metadata %>% left_join(age_gender,by="subject") %>% left_join(rf_data %>% select(subject,day,Area_delta_R30,plaquepercent,bomppercent) %>% filter(day==14), by="subject")

  result = adonis2(transformedSubjectLoadings ~ plaquepercent + bomppercent + Area_delta_R30 + gender + age, data=metadata, method="euclidean", by="margin", permutations=9999, na.action=na.omit)

  return(result)
}

testFeatures = function(model, metadata, componentNum, metadataVar){
  df = metadata$mode1 %>% left_join(rf_data %>% filter(day==14)) %>% left_join(age_gender)
  topIndices = metadata$mode2 %>% mutate(index=1:nrow(.), Comp = model$Fac[[2]][,componentNum]) %>% arrange(desc(Comp)) %>% head() %>% select(index) %>% pull()
  bottomIndices = metadata$mode2 %>% mutate(index=1:nrow(.), Comp = model$Fac[[2]][,componentNum]) %>% arrange(desc(Comp)) %>% tail() %>% select(index) %>% pull()

  timepoint = which(abs(model$Fac[[3]][,componentNum]) == max(abs(model$Fac[[3]][,componentNum])))

  Xhat = parafac4microbiome::reinflateTensor(model$Fac[[1]][,componentNum], model$Fac[[2]][,componentNum], model$Fac[[3]][,componentNum])
  y = df[metadataVar]

  print("Positive loadings:")
  print(cor(Xhat[,topIndices[1],timepoint], y))
  print("Negative loadings:")
  print(cor(Xhat[,bottomIndices[1],timepoint], y))
}


# ACMTF
homogenizedSubjects = parafac4microbiome::vanderPloeg2024$metabolomics$mode1$subject
mask = parafac4microbiome::vanderPloeg2024$tongue$mode1$subject %in% homogenizedSubjects

# Tongue
temp = parafac4microbiome::vanderPloeg2024$tongue
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedTongue = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Lowling
temp = parafac4microbiome::vanderPloeg2024$lower_jaw_lingual
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedLowling = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Lowinter
temp = parafac4microbiome::vanderPloeg2024$lower_jaw_interproximal
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedLowinter = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Upling
temp = parafac4microbiome::vanderPloeg2024$upper_jaw_lingual
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedUpling = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Upinter
temp = parafac4microbiome::vanderPloeg2024$upper_jaw_interproximal
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedUpinter = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Saliva
temp = parafac4microbiome::vanderPloeg2024$saliva
temp$data = temp$data[mask,,]
temp$mode1 = temp$mode1[mask,]
processedSaliva = processDataCube(temp, sparsityThreshold=0.50, considerGroups=TRUE, groupVariable="RFgroup", CLR=TRUE, centerMode=1, scaleMode=2)

# Metabolomics stays the same
processedMetabolomics = parafac4microbiome::vanderPloeg2024$metabolomics

datasets = list(processedTongue$data, processedLowling$data, processedLowinter$data, processedUpling$data, processedUpinter$data, processedSaliva$data, processedMetabolomics$data)
modes = list(c(1,2,3),c(1,4,5),c(1,6,7),c(1,8,9),c(1,10,11),c(1,12,13),c(1,14,15))

CV_ACMTF = ACMTF_modelSelection(datasets, modes, maxNumComponents=5, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())

saveRDS(CV_ACMTF, "./TIFN2_CV_ACMTF.RDS")
