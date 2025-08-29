library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)
library(vegan)

# Microbiome
processedMicrobiome = CMTFtoolbox::Georgiou2025$Tooth_microbiome

# Remove samples due to low library size
mask = rowSums(processedMicrobiome$data) > 6500
processedMicrobiome$data = processedMicrobiome$data[mask,]
processedMicrobiome$mode1 = processedMicrobiome$mode1[mask,]

# Remove duplicate samples
mask = !(processedMicrobiome$mode1$SampleID %in% c("A11-8 36", "A11-10 17", "A11-15 17"))
processedMicrobiome$data = processedMicrobiome$data[mask,]
processedMicrobiome$mode1 = processedMicrobiome$mode1[mask,]

# Also remove subject A11-8 due to being an outlier
processedMicrobiome$data = processedMicrobiome$data[-23,,]
processedMicrobiome$mode1 = processedMicrobiome$mode1[-23,]

# CLR transformation
df = processedMicrobiome$data + 1
geomeans = pracma::geomean(as.matrix(df), dim=2)
df_clr = log(sweep(df, 1, geomeans, FUN="/"))

# Feature filtering
sparsityThreshold = 0.5
maskA = processedMicrobiome$mode1$PainS_NopainA == "A"
maskS = processedMicrobiome$mode1$PainS_NopainA == "S"

dfA = processedMicrobiome$data[maskA,]
dfS = processedMicrobiome$data[maskS,]

sparsityA = colSums(dfA == 0) / nrow(dfA)
sparsityS = colSums(dfS == 0) / nrow(dfS)

mask = (sparsityA <= sparsityThreshold) | (sparsityS <= sparsityThreshold)

processedMicrobiome$data = df_clr[,mask]
processedMicrobiome$mode2 = processedMicrobiome$mode2[mask,]

# Center and scale
processedMicrobiome$data = sweep(processedMicrobiome$data, 2, colMeans(processedMicrobiome$data), FUN="-")
processedMicrobiome$data = sweep(processedMicrobiome$data, 2, apply(processedMicrobiome$data, 2, sd), FUN="/")

# Cytokines
processedCytokines_case = CMTFtoolbox::Georgiou2025$Inflammatory_mediators

# Select only case subjects
mask = processedCytokines_case$mode1$case_control == "case"
processedCytokines_case$data = processedCytokines_case$data[mask,,]
processedCytokines_case$mode1 = processedCytokines_case$mode1[mask,]

# Select only samples with corresponding microbiome data
mask = processedCytokines_case$mode1$SubjectID %in% processedMicrobiome$mode1$SubjectID
processedCytokines_case$data = processedCytokines_case$data[mask,,]
processedCytokines_case$mode1 = processedCytokines_case$mode1[mask,]

processedCytokines_case$data = log(processedCytokines_case$data + 0.005)
processedCytokines_case$data = multiwayCenter(processedCytokines_case$data, 1)
processedCytokines_case$data = multiwayScale(processedCytokines_case$data, 2)

# Prep data
datasets = list(processedCytokines_case$data, as.matrix(processedMicrobiome$data))
modes = list(c(1,2,3),c(1,4))
Z = setupCMTFdata(datasets, modes, normalize=TRUE)

Y = as.numeric(as.factor(processedCytokines_case$mode1$PainS_NopainA))
Ycnt = Y - mean(Y)
Ynorm = Ycnt / norm(Ycnt, "2")
Ynorm = as.matrix(Ynorm)

result = ACMTFR_modelSelection(datasets, modes, Ynorm, maxNumComponents=10, pi=0.50, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(result, "./AP_CV_ACMTFR_50.RDS")

result = ACMTFR_modelSelection(datasets, modes, Ynorm, maxNumComponents=10, pi=0.60, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(result, "./AP_CV_ACMTFR_60.RDS")

result = ACMTFR_modelSelection(datasets, modes, Ynorm, maxNumComponents=10, pi=0.70, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(result, "./AP_CV_ACMTFR_70.RDS")

result = ACMTFR_modelSelection(datasets, modes, Ynorm, maxNumComponents=10, pi=0.80, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(result, "./AP_CV_ACMTFR_80.RDS")

result = ACMTFR_modelSelection(datasets, modes, Ynorm, maxNumComponents=10, pi=0.90, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
saveRDS(result, "./AP_CV_ACMTFR_90.RDS")
