# KT004_transform_and_normalize_data

```r
rm(list=ls())
setwd("C:\\Users\\helbd\\Documents\\Magic Briefcase\\LSHTM\\Kevin's array\\151118 data")
```

##########################################################################################
## Load package libraries
```r
library(gplots)
library(heatmap.plus)
library(limma)
library(marray)
library(multcomp)
library(robust)
library(samr)
library(vsn)
library(mgcv)
```

##########################################################################################
##load source code
```r
source(".\\R Code\\Raw Scripts\\pArrayFunctions.R")
```

##########################################################################################
##Import data
```r
print(load(".\\Data\\Processed Data\\0003.KTarray.Apac.explored.RData"))

pheSera$Slide2 <- pheSera$Slide
pheSera$Slide <- paste(gsub("Apac_", "", pheSera$study), pheSera$Slide, sep="_")
```

##########################################################################################
##Remove slides where blanks are too high
```r
remove <- pheSera[!is.na(pheSera$Flag),]
pheSera <- pheSera[is.na(pheSera$Flag),]
MNA <- MNA[,colnames(MNA) %in% remove$SampleID==F]
MedA <- MedA[,colnames(MedA) %in% remove$SampleID==F]
```

##########################################################################################
##get dataframes into the same order
```r
annAnti <- annAnti[order(annAnti$Unique.ID),]
pheSera <- pheSera[order(pheSera$SampleID),]
MNA.raw <- MNA[order(rownames(MNA)),order(colnames(MNA))]
MedA.raw <- MedA[order(rownames(MedA)),order(colnames(MedA))]
```

##########################################################################################
##Define probe sets
```r
blank.probes     = grep("blank",annAnti$ID)     #blank spots: left intentionally empty, representing the background noise on the slide itself (should be as close to zero as possible)
pctrl.probes     = c(grep("std",annAnti$ID))    #positive control spots:
ahIgG.probes     = NA                           #anti-human IgG spots: will bind antibodies (not malaria specific) from samples to verify that the sample was correctly hybridized to the array
hIgG.probes      = grep("std",annAnti$ID)       #human IgG spots: verify that the fluorophore-bound anti-human IgG is actually binding to the antibodies in the hybridized samples
ttbs.probes      = grep("PBS",annAnti$ID)       #TTBS spots: buffer only, verify no contamination
ctrl.probes      = c(grep("blank",annAnti$ID), grep("PBS",annAnti$ID), grep("std",annAnti$ID), grep("GST",annAnti$ID)) #control spots (positive or negative ctrls):
ntar.probes      = c(grep("blank",annAnti$ID), grep("PBS",annAnti$ID), grep("std",annAnti$ID), grep("GST",annAnti$ID))  #non-target spots (all spots that are not PF antigens of interest):
GST.probes       = c(1:length(annAnti$ID))[c(1:length(annAnti$ID)) %in% c(ntar.probes, grep("PFAMA", annAnti$ID),grep("PFSE", annAnti$ID))==F] #GST-tagged spots
```

##########################################################################################
##Script settings
```r
transform       = "VSN"     # transormation method to make the data normally-distributed (must choose "LOG2", "ASINH" or "VSN")
Normalization   = "RLMc"    # method for normalizing variation between slides, subarrays, batches (must choose "RLMct", "RLMc", "Median", "Quantile" or "none")
```

##########################################################################################
#transform/normalize MedA (median)

##select dataframe
```r
intensity.raw   = MedA.raw  # select which dataset to use (must choose "MNA.raw" or "MedA.raw")
raw             = "MedA"    # must select("MedA" or "MNA" and must match intensity.raw)
```

##########################################################################################
##Remove Slide Background
```r
expr <- intensity.raw

bkgrd <- intensity.raw[grep("blank", rownames(intensity.raw)),]
bkgrd[bkgrd<0 & is.na(bkgrd)==F] <- 0
bkgrd.lvl <- colMeans(bkgrd, na.rm=T) #for each subject, determine their mean slide background

for(i in 1:ncol(expr)) expr[,i] = expr[,i] - bkgrd.lvl[i] 
intensity.BS <- expr
```

##########################################################################################
##Remove Cross-reactivity to GST
```r
expr <- intensity.BS

gst <- intensity.BS[grep("GST", rownames(intensity.BS)),]
gst[gst<0 & is.na(gst)==F] <- 0
gst.lvl <- colMeans(gst, na.rm=T) #for each subject, determine their level of response to GST

for(i in 1:ncol(expr)) expr[GST.probes,i] = expr[GST.probes,i] - gst.lvl[i] #only subtract GST reactivity from GST-tagged probes
intensity.S <- expr
```

##########################################################################################
##Transform data by VSN
```r
expr <- intensity.S
pheSera <- pheSera[match(colnames(expr),pheSera$SampleID),]

fit <- vsn2(expr[ntar.probes,])
expr <- predict(fit,expr) 
intensity.S_T <- expr
```

##########################################################################################
##Normalization (RLMC: ROBUST LINEAR MODEL NORMALIZATION, BASED ONLY ON CONTROL PROBES)
```r
expr <- intensity.S_T
pheSera <- pheSera[match(colnames(expr),pheSera$SampleID),]

pdf(".\\Figures\\Exploratory\\temp.pdf")
expr = norm.RLM(expr=expr,Method=Normalization,ctrl.probes=ctrl.probes)
dev.off()

intensity.S_T_N <- expr
```


##########################################################################################
##Set minimum intensity reading
```r
intensity.S_T_N [intensity.S_T_N <= 0] <- 0
```

##########################################################################################
##clean up data
```r
do.call("<-",list(paste("intensity", raw, sep="."), intensity.S_T_N)) 
```

##########################################################################################
#transform/normalize MNA (mean)

##select dataframe
```r
intensity.raw   = MNA.raw  # select which dataset to use (must choose "MNA.raw" or "MedA.raw")
raw             = "MNA"    # must select("MedA" or "MNA" and must match intensity.raw)
```

##########################################################################################
##Remove Slide Background
```r
expr <- intensity.raw

bkgrd <- intensity.raw[grep("blank", rownames(intensity.raw)),]
bkgrd[bkgrd<0 & is.na(bkgrd)==F] <- 0
bkgrd.lvl <- colMeans(bkgrd, na.rm=T) #for each subject, determine their mean slide background

for(i in 1:ncol(expr)) expr[,i] = expr[,i] - bkgrd.lvl[i] 
intensity.BS <- expr
```

##########################################################################################
##Remove Cross-reactivity to GST
```r
expr <- intensity.BS

gst <- intensity.BS[grep("GST", rownames(intensity.BS)),]
gst[gst<0 & is.na(gst)==F] <- 0
gst.lvl <- colMeans(gst, na.rm=T) #for each subject, determine their level of response to GST

for(i in 1:ncol(expr)) expr[GST.probes,i] = expr[GST.probes,i] - gst.lvl[i] #only subtract GST reactivity from GST-tagged probes
intensity.S <- expr
```

##########################################################################################
##Transform data by VSN
```r
expr <- intensity.S
pheSera <- pheSera[match(colnames(expr),pheSera$SampleID),]

fit <- vsn2(expr[ntar.probes,])
expr <- predict(fit,expr) 
intensity.S_T <- expr
```

##########################################################################################
##Normalization (RLMC: ROBUST LINEAR MODEL NORMALIZATION, BASED ONLY ON CONTROL PROBES)
```r
expr <- intensity.S_T
pheSera <- pheSera[match(colnames(expr),pheSera$SampleID),]

pdf(".\\Figures\\Exploratory\\temp.pdf")
expr = norm.RLM(expr=expr,Method=Normalization,ctrl.probes=ctrl.probes)
dev.off()

intensity.S_T_N <- expr
```

##########################################################################################
##Set minimum intensity reading
```r
intensity.S_T_N [intensity.S_T_N <= 0] <- 0
```

##########################################################################################
##clean up data
```r
do.call("<-",list(paste("intensity", raw, sep="."), intensity.S_T_N)) 
```

##########################################################################################
#clean up workspace & save
```r
save(annAnti, pheSera, intensity.MedA, intensity.MNA, MedA.raw, MNA.raw, file=".\\Data\\Processed Data\\0004.KTarray.Apac.transformed_normalized.RData")

MedA <- data.frame(intensity.MedA)
MedA$antigen <- gsub("_.+$", "", rownames(MedA))
MedA <- MedA[,c(991,1:990)]
write.csv(MedA, file=".\\Data\\Processed Data\\0004.KTarray.Apac.transformed_normalized.MedA.csv", row.names=T)

MNA <- data.frame(intensity.MNA)
MNA$antigen <- gsub("_.+$", "", rownames(MNA))
MNA <- MNA[,c(991,1:990)]
write.csv(MNA, file=".\\Data\\Processed Data\\0004.KTarray.Apac.transformed_normalized.MNA.csv", row.names=T)
```
