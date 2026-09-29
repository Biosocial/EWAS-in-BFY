################################################### 
####  DNAm Quality Control general and EWAS #######     
################################################### 

#Pre-registration can be found here: https://osf.io/8grhj

#### Preparing the data ####

#Data:
#--/data/3_IDATS

#Processed data:
#--/data/4_Processed Data

#Rscripts:
#--/data analysis


#### Install Packages ####

### `pacman` package
install.packages(setdiff("pacman", rownames(installed.packages()))) #this allows to install and load packages in one go, instead of install() library() them one by one

### CRAN packages
pacman::p_load(char = c("tidyverse", "rio", "data.table", "here", "abind", "RColorBrewer", "MatrixEQTL", "qqman",
                        "nlme", "BiasedUrn", "Hmisc", "BiocManager", "devtools", "MASS", "qvalue", "sva", "isva", "readxl","openxlsx", "dplyr", "scales", "VennDiagram"), install = T)

### `BiocManager` packages
pacman::p_load(char = c(
  "bumphunter", # Differentially methylated region analysis
  "minfi", # Many functions for methylation analysis
  "sva", # Batch effect correction
  "limma", # Single-site association analysis
  "IlluminaHumanMethylationEPICv2manifest", # Manifest files for EPIC V2
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38", # Annotation files for EPIC V2
  "FDb.InfiniumMethylation.hg19", "FlowSorted.Blood.EPIC", "BeadSorted.Saliva.EPIC", 
  "GO.db", # Gene ontology queries
  "missMethyl", "GEOquery", "wateRmelon", "ENmix"),
  install = T, update = T)

### GitHub package for saliva-based cell type estimation
pacman::p_load_gh("hhhh5/ewastools")

### 'RefFreeEWAS' package (outdated)
#install.packages("https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/RefFreeEWAS_2.1.tar.gz", repos=NULL, type="source", dependencies = T);library(RefFreeEWAS)

# install 
BiocManager::install("ENmix")
library(ENmix)


# When packages are already installed, just run code below:
### Load CRAN packages
cran_packages <- c(
  "tidyverse", "rio", "data.table", "here", "abind", "RColorBrewer", 
  "MatrixEQTL", "qqman", "nlme", "BiasedUrn", "Hmisc", "BiocManager", 
  "devtools", "MASS", "qvalue", "sva", "isva", "readxl", "openxlsx", 
  "dplyr", "scales", "VennDiagram"
)
lapply(cran_packages, library, character.only = TRUE)

### Load BiocManager packages
bioc_packages <- c(
  "bumphunter", "minfi", "sva", "limma", "IlluminaHumanMethylationEPICv2manifest", 
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "FDb.InfiniumMethylation.hg19", 
  "FlowSorted.Blood.EPIC", "BeadSorted.Saliva.EPIC", "GO.db", "missMethyl", 
  "GEOquery", "wateRmelon", "ENmix"
)
lapply(bioc_packages, library, character.only = TRUE)

### Load GitHub packages
if ("ewastools" %in% rownames(installed.packages())) library(ewastools)


#### Read in Phenotypic data ####

# this contains phenotypic info such (e.g., sex, treatment/control(blinded)) and sample quality (for later sensitivity checks)
PhenoExcel <- read_excel("--/data/3_IDATS/BFY_plates_layout.xlsx", sheet = "Sample List", skip = 1)

# Check classes of all variables
classes_samplesheet <- sapply(PhenoExcel, class)
classes_samplesheet

### Rename unique_Sample_ID to Sample_Name for merging with samplesheet (see below)
## Sample_ID is a composite of sampleid (e.g. 12345) and whether it is child or mother sample (_C or _M). 
## So Sample_ID looks like p12345_C or p12345_M (p is added by lab for person, and same id means same family)
## Block is a variable that groups mother and child from the same family.
colnames(PhenoExcel)[colnames(PhenoExcel) == "unique_Sample_ID"] <- "Sample_Name"

# Add variable if a sample is a technical replicate
PhenoExcel <- PhenoExcel %>%
  mutate(Techreplicate = ifelse(grepl("_repl", Sample_Name), 1, 0))

table(PhenoExcel$Techreplicate) #1518 samples and 18 replicates

comment(PhenoExcel$Techreplicate) <- "0=Sample, 1=Technical Replicate" # Using comment to add a label
comment(PhenoExcel$Techreplicate)

# Create new variable to distinguish mother and child samples
## there is a "Child.is.Female" variable, which gives YES for both mothers and child themselves
## To create a variable "sex" we give all mothers and all daughters 1, and all sons 0
PhenoExcel <- PhenoExcel %>%
  mutate(Mum0Child1 = ifelse(Code.C..child..M..mother. == "C", 1, 0))

PhenoExcel <- PhenoExcel %>%
  relocate(Mum0Child1, .after = Code.C..child..M..mother.)

# Check if all went right
## so all 1=Children, and all 0=Mothers
table(PhenoExcel$Mum0Child1, PhenoExcel$Code.C..child..M..mother.)

comment(PhenoExcel$Mum0Child1) <- "0=Mothers, 1=Children" # Using comment to add a label
comment(PhenoExcel$Mum0Child1)

# Create sex, with 0=boys, 1=girls (kids) and women (mothers)
PhenoExcel <- PhenoExcel %>%
  mutate(sex = case_when(
    Child.is.Female == "Yes" & Code.C..child..M..mother. == "C" ~ 1,
    Code.C..child..M..mother. == "M" ~ 1,
    Child.is.Female == "No" & Code.C..child..M..mother. == "C" ~ 0,
    TRUE ~ NA_real_  # Optional: sets any other combination to NA
  ))

PhenoExcel <- PhenoExcel %>%
  relocate(sex, .after = Child.is.Female)


comment(PhenoExcel$sex) <- "0=boys, 1=girls/women" # Using comment to add a label
comment(PhenoExcel$sex)

# Check if this went well
## all PhenoExcel$Mum0Child1= 0 (Moms) should score 1 on sex
table(Sex = PhenoExcel$sex,
      Mum0Child1 = PhenoExcel$Mum0Child1)

## all PhenoExcel$Child.is.Female = YES should be PhenoExcel$sex =1
## Note, not all PhenoExcel$Child.is.Female = NO should be PhenoExcel$sex = 0 cause some are mothers who have sons so they get sex =1
table(Sex = PhenoExcel$sex,
      Child.is.Female = PhenoExcel$Child.is.Female)



# Check descriptives in the sample

dim(PhenoExcel) #n = 1536, in total we have 1536 Samples
table(PhenoExcel$Techreplicate) #1518 samples and 18 replicates

#### Read in Samplesheet ####

# This contains sentrixID and sentrix position that correspond with later IDATS
# Read in the Excel file, skipping the first 7 rows and making the 8th row the Column names 
samplesheet <- read_excel("--/data/3_IDATS/M01191_01-16_samplesheet.xlsx", skip = 7)

# check classes
classes_samplesheet <- sapply(samplesheet, class)


#### Merge Samplesheet and Phenodata ####

# Here we merge the samplesheet and the phenoytype data so that we have all the information in one file
samplesheet_pheno <- samplesheet %>%
  left_join(PhenoExcel, by = "Sample_Name")

# save it as .csv file and safe it in same folder as IDAT files, cause they need to be in the same folder in order for them to be linked
write.csv(samplesheet_pheno, "--/data/3_IDATS/samplesheet.csv", row.names = FALSE)

#read samplesheet
samplesheet_pheno <- read.csv("--/data/3_IDATS/samplesheet.csv", stringsAsFactors = FALSE)


# read in IDATS
baseDir<- ("--/data/3_IDATS/")
targets <- read.metharray.sheet(baseDir)
anyDuplicated(targets$Sample_Name) # check for duplicated Oragene IDs, won't work for epi preprocessing otherwise

anyDuplicated(targets$Basename) # check for duplicated Oragene IDs, won't work for epi preprocessing otherwise

# check targets
# paths look good
targets_check <- targets[,c("Sample_Name", "sex", "Basename" )]
head(targets_check)

#Sometimes it copies c(" in the path name. If this is the case, remove it otherwise it does not read it in properly
#targets$Basename <- sub("^c\\(\"", "", targets$Basename)

# This is how we found out some IDATS were missing when they sent us the files in Octobre
# This is a list of missing IDATS, which were send to us later by the Lab
# missing_idats <- targets  %>% 
#  filter(Basename == "character(0)")
#write.xlsx(missing_idats, "--/data/3_IDATS/missingIDATS.xlsx")

#Here we select only 8 samples so we can run through the script without running it across all 1536 samples
#This we did before running the full script on all samples to check if the script worked
#targets <- targets[1:8, ]

#### Create RGSet ####

# RGSet object is a convenient format that simplifies the analysis of raw methylation data from IDAT files, providing a more efficient way to manipulate and process the data for further analysis
# RGSet is a kind of multidimensional object with a multitude of datasheets with info (green/red sheets, pheno data, etc)
RGSet <- read.metharray.exp(targets = targets, extended = T, force=TRUE) # read methylation array IDATS using sample sheet


# save RGChannelSet object so we can load that directly rather than loading in IDATS again (computationally more expensive)
#save(RGSet, file="--/data/4_Processed Data/RGset_full.rda")
load(file.path("--/data/4_Processed Data/RGset_full.rda"))

# N=1536 samples
dim(RGSet)

# save RGset only for 8 samples in case you want to run some checks on part of the scripts, without having to run it for 1536 samples
#save(RGSet, file="--/data/4_Processed Data/RGset_8samples.rda")
load(file.path("--/data/4_Processed Data/RGset_8samples.rda"))
# 8 samples
#dim(RGSet)

#Generate manifest and annotation for EPIC V2
#see https://jokergoo.github.io/IlluminaHumanMethylationEPICv2manifest/articles/IlluminaHumanMethylationEPICv2manifest.html

annotation(RGSet)
annotation(RGSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGSet)["annotation"] = "20a1.hg38"
annotation(RGSet)

### RGChannelSet stores also a manifest object that contains the probe design information of the array
manifest <- getManifest(RGSet)       
annotation <- getAnnotation(RGSet)

dim(annotation)

### extract methylated and unmethylated signals MSet (MethylSet)
MSet <- preprocessRaw(RGSet) # do preprocessNoob on RGSet.drop later

#N=1536, with CpG=936990
dim(MSet)

#get info we need for later analyses
GRset <-mapToGenome(MSet)
pheno <- pData(GRset)                          # phenotypic data

table(pheno$sex)                               # table of reported sex
sum(table(pheno$sex))==nrow(samplesheet_pheno) # check to make sure there is reported sex for everyone; output should be "TRUE" 
rm(GRset)



#### Quality control ####

Meth <-getMeth(MSet) # M signal per probe, per sample
Unmeth <-getUnmeth(MSet) # U signal per probe, per sample

#### Log Median Intensity ####
# `minfi` package provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels
# NOTE: When plotting these two medians against each other, it has been observed that good samples cluster together, while failed samples tend to separate and have lower median intensities
# This is not often used for final selection, but good to check anyways. 

# generate the plot
pdf("--/data analysis/00_QC/Plots_full/LogMedianIntensity.pdf") #open PDF device to save the plot

qc <- getQC(MSet)
plotQC(qc)

dev.off() #closes the PDF device, ensuring that the file is properly saved.


#### Density plots ####
# Look at the densities, and group by sex and by child/mother samples
# Here we do it based on the Mvalues, later you can also do this for the beta values

# Generate the density bean plot, with sex as grouping
pdf("--/data analysis/00_QC/Plots_full/BeanPlot_sex.pdf") #open PDF device to save the plot

minfi::densityBeanPlot(MSet, sampGroups = pheno$sex) # color by sex
legend(x = "topright", inset = c(0, -0.2), legend = c("Female", "Male"), 
       fill = c("darkorange", "#1B9E77"), title = "Sex", xpd = TRUE, horiz = TRUE, bty = "n")

dev.off() #closes the PDF device, ensuring that the file is properly saved.

# Generate the densityBeanPlot, with child/mother grouping
pdf("--/data analysis/00_QC/Plots/BeanPlot_ChildMother.pdf") #open PDF device to save the plot

minfi::densityBeanPlot(MSet, sampGroups = pheno$Mum0Child1) # Colors based on child/mother
legend(x = "topright", inset = c(0, -0.2), legend = c("Child", "Mother"), 
       fill = c("darkorange", "#1B9E77"), title = "Group", xpd = TRUE, horiz = TRUE, bty = "n")

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#### PCA plots ####
# Raw PCA plots
beta.raw = getBeta(MSet)
beta.raw2=beta.raw;
beta.raw2[ is.na(beta.raw)] =0 
rm(beta.raw)

#makes the principal component object
#very memory intensive, only run when you have plenty of time
#pulls out proportion of variance explained by each PC
pdf("--/data analysis/00_QC/Plots_full/PropVarExplained.pdf") #open PDF device to save the plot

prin = prcomp (t(beta.raw2), center=T, scale.=F)
rm(beta.raw2)
screeplot(prin, col="dodgerblue", xlab="Principal Components of Raw Beta Values", main=" ") #, cex.label=1.3

dev.off() #closes the PDF device, ensuring that the file is properly saved.

### create a PC plot matrix -> colored by sex
pdf("--/data analysis/00_QC/Plots/PCplot_sex.pdf") #open PDF device to save the plot

par(mar=c(0,0,0,0))
plot.new()
legend("bottom", levels (as.factor(pheno$sex)), fill = as.factor(pheno$sex), title="Principal Components by Sex")
pairs(prin$x[,1:6], col=as.factor(pheno$sex), labels= c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", och=1, cex=0.5))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

### create a PC plot matrix -> colored by batch
pdf("--/data analysis/00_QC/Plots/PCplot_Array.pdf") #open PDF device to save the plot

par(mar=c(0,0,0,0))
plot.new()
legend("bottom", levels (as.factor(pheno$Array)), fill = as.factor(pheno$Array), title="Principal Components by Array Position")
pairs(prin$x[,1:6], col=as.factor(pheno$Array), labels= c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", och=1, cex=0.5))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

### create a PC plot matrix -> colored by slide
pdf("--/data analysis/00_QC/Plots/PCplot_Slide.pdf") #open PDF device to save the plot

par(mar=c(0,0,0,0))
plot.new()
legend("bottom", levels (as.factor(pheno$Slide)), fill = as.factor(pheno$Slide), title="Principal Components by Slide")
pairs(prin$x[,1:6], col=as.factor(pheno$Slide), labels= c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", och=1, cex=0.5))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#### Control probes plot ####

#several internal control probes that can be used to assess the quality control of different sample preparation steps (bisulfite conversion, hybridization, etc.).
#The values of these control probes are stored in the initial RGChannelSet and can be plotted by using the function controlStripPlot and by specifying the control probe type:
pdf("--/data analysis/00_QC/Plots/ControlProbes.pdf") #open PDF device to save the plot

controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
dev.off() #closes the PDF device, ensuring that the file is properly saved.

#### Look at Sex Mismatch ####

#This is done later, see below

#### Low Intensity Exclusion ####
pdf("--/data analysis/00_QC/Plots_full/IntensityPlots.pdf")

MQC <- log2(colMedians(Meth))
UQC <- log2(colMedians(Unmeth))
rm(Meth, Unmeth)

MSet$Array <- factor(MSet$Array)

#Note: in preregistration QC we state that probes will be excluded log2 <9 

### plot M vs U intensity -> color by position
plot(UQC, MQC, col=MSet$Array, main ="M vs. U QC by Position", pch=16,xlab="Log2 Median Unmethylated Intensity",
     ylab="Log2 Median Methylated Intensity", cex.lab=1.2,cex.main=2)
legend("bottomright", levels(MSet$Array), fill = MSet$Array)

### plot M vs U intensity -> color by slide
plot(UQC, MQC, col=as.factor(MSet$Slide), main ="M vs. U QC by Slide", pch=16,xlab="Log2 Median Unmethylated Intensity",
     ylab="Log2 Median Methylated Intensity", cex.lab=1.2,cex.main=2)
legend("topleft", levels(as.factor(MSet$Slide)), fill = as.factor(MSet$Slide))

dev.off() #closes the PDF device, ensuring that the file is properly saved.


#### Mark exclusion probes ####

#### Based on Methylation intenstiy ####
#Drop (or if really small sample: watch out for): Samples with UQC<9 & MQC<9

length(which(UQC<9))  #here you check how many have <9 #N=8
length(which(MQC<9))  #here you check how many have <9 #N=4

xUQC= MSet$Basename[which(UQC<9)]
yMQC= MSet$Basename[which(MQC<9)]

xUQC= MSet@colData@rownames[which(UQC<9)]
yMQC= MSet@colData@rownames[which(MQC<9)]

sample.exclusion = c(xUQC, yMQC) # n=8 samples in total
save(sample.exclusion, file="--/data/4_Processed Data/SampleExclusionUQCMQC.rda")

load("--/data/4_Processed Data/SampleExclusionUQCMQC.rda")
print(sample.exclusion)

# These 8 samples should be excluded (Sample_Name)

#  8 samples should be excluded (rownames)

#### Detection P ####
#DetectionP gives us a metric for assessing signal/noise ratios. Both samples and probes with higher detection p levels (indicative of unacceptable levels of noise) are cut or flagged. How close to 0 are the values? How confident are we that we are detecing anything that isn't background?

#Note: 
#Preregistration : we exclude CpG probes if they have a detection of p>0.01
#Preregistration : we exclude samples if detection p is >.01 in >10% of their probes

##detection P represents the confidence that a given transcript is expressed above the background defined by negative control probes
#Here we flag probes that have a detection of p>.01
hp = minfi::detectionP(RGSet) > 0.01 
save(hp, file=("--/data/4_Processed Data/hp_detectionp.rda"))
exclude.hpv = rownames(hp)[rowMeans(hp) > 0.01]
length(exclude.hpv) #23581  
save(exclude.hpv, file=("--/data/4_Processed Data/exclude.hpv.rda")) #need this later
load("--/data/4_Processed Data/exclude.hpv.rda")
rm(hp)

#here we flag samples if detection p is >.01 in >10% of their probes, and probes that have a success rate of <.95 across all samples
detP <- minfi::detectionP (RGSet)

### Explore number of detection P failed probes/samples
failed <- detP > 0.01  #hp ==failed
per.samp<-colMeans( failed) # Fraction of failed positions per sample
per.probe<-rowMeans( failed) # Fraction of failed samples per position

#summary (per.samp)
sum(per.samp>0.10)  # How many samples had more than 10% of sites fail? #nr of people to be excluded: 5
sum(per.probe>0.05) # How many positions failed in >5% of samples? #nr of probes  to be excluded: 9104

sample.fail<- per.samp[per.samp>0.10] # How many samples had more than 10% of sites fail? #nr of people to be excluded: 5
probe.fail<- failed[per.probe>0.05,]  # How many positions failed in >5% of samples? #nr of probes  to be excluded: 9104

#Samples that had more than 10% of sites fail

### Exclusion in sample sheet for samples/participants
###Drop any samples that failed by detection P by downstream exclusion
save(sample.fail, file= "--/data/4_Processed Data/SampleExclusionDetectionP.rda") # 5 samples
save(probe.fail, file= "--/data/4_Processed Data/probe.fail.rda") # 9104 probes, but all of these probes are part of the exclude.hpv probes, see code below


#Check if probes that fail in more than 10% of the sample are part of the exclude.hpv, where we flag all probes with detection p>.01
pf <- rownames(probe.fail)

# Extract the overlapping CpG names from pf
pf_overlapping <- pf[pf %in% exclude.hpv]

# If you want it as a dataframe
pf_overlapping_df <- data.frame(CpG_Sites = pf_overlapping)

# Check if all entries in pf overlap with exclude.probes
all_overlap <- length(pf) == nrow(pf_overlapping_df)

# if answer is TRUE, this means all CpG sites in probe.fail are in exclude.hpv
print(all_overlap)

### Plot Detection P
# examine mean detection p-values across all samples to identify any failed samples
# In this plot, check if any values of participants exceed mean value of .01

pdf("--/data analysis/00_QC/Plots/DetectionP-bySampleEx.pdf")

pal <- brewer_pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP),  las=2, #col=pal[factor(targets$sex)],
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.01,col="red")
#legend("topleft", legend=levels(factor(targets$sex)), fill=pal,
#bg="white") 

dev.off() #closes the PDF device, ensuring that the file is properly saved.


# Calculate column means
MeanDetp <- colMeans(detP)

# Convert to a data frame
MeanDetp_df <- data.frame(ID = names(MeanDetp), MeanDetectionP = MeanDetp, row.names = NULL)

# Create a new dataframe with IDs where MeanDetectionP > 0.01
MeanDetp_df_filtered <- MeanDetp_df %>%
  filter(MeanDetectionP > 0.01)

# View the new dataframe
print(MeanDetp_df_filtered)


#Histogram of Detection P per sample
pdf("--/data analysis/00_QC/Plots/HistDetectionP-bySample.pdf")

hist (per.samp , breaks = 40 , col = "dodgerblue" , cex.lab =1.3, xlab  = "Fraction of failed positions per sample",xlim=c(0,0.011))
abline(v=0.01, col="red", lwd=2)

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#Histogram of Detection P per probe
pdf("--/data analysis/00_QC/Plots/HistDetectionP-perProbe.pdf")

hist (per.probe , breaks = 50 , col = "dodgerblue" , cex.lab =1.3, xlab  = "Fraction of failed positions per probe", ylim=c(0,10000))
abline(v=0.01, col="red", lwd=2)
dev.off()

#### Cross-reactive probes ####
### probes with SNPs and in cross-reactive regions
### cross-reactive probes are shown to 'co-hybridize' onto sex chromosomes and may show spurious sex differences in methylation that are merely due to technical artifacts

#The List of Cross-probes is retrieved from Peters et al., (2024)
#Consisting of probes that are cross-reactive for Epic V2
#See Additional file 4, "12864_2024_10027_MOESM4_ESM.csv"
#https://doi.org/10.1186/s12864-024-10027-5
#Based on this, Trey Smith created cross reactive probes file for Epic V2 "EPICv2_cross_reactive_probes.csv"

#note, we created beta file with and without cross-reactive probes
#Q: cross-reactive probes: 30627
cross.probes.info <- read.csv("--/data/4_Processed Data/EPICv2_cross_reactive_probes.csv")


#### NBeads ####
### preregistration: exclude CpG probes if they have fewer than 4 beads in more than 5% of the samples

# Get information for Type I and Type II probes
pi1 = getProbeInfo(RGSet, type = "I")
pi2 = getProbeInfo(RGSet, type = "II")

# Identify probes with fewer than 4 beads
lb = getNBeads(RGSet) < 4

# Find Type I probes with fewer than 4 beads in more than 5% of samples
ex1 = pi1$Name[rowMeans(lb[pi1$AddressA,] | lb[pi1$AddressB,]) > 0.05]

# Find Type II probes with fewer than 4 beads in more than 5% of samples
ex2 = pi2$Name[rowMeans(lb[pi2$AddressA,]) > 0.05]

# Combine and duplicate probe names for exclusion
exclude.bds = unique(c(ex1, ex2))

# Print the number of probes that meet the criteria
# N=11079
length(exclude.bds)
head(exclude.bds)

save(exclude.bds, file=("--/data/4_Processed Data/exclude.bds.rda")) #need this later
load("--/data/4_Processed Data/exclude.bds.rda")

##### Exclusion of SNP Affected Probes  #####
##Morgan Levine uses minfis dropLociWithSnps, Bakulski/Mitchel use gaphunter

### from #https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
#Genetic variants Because the presence of SNPs inside the probe body or at the nucleotide extension can have important consequences on the downstream analysis,
#minfi offers the possibility to remove such probes. The function getSnpInfo, applied to a GenomicRatioSet, returns a data frame with 6 columns containing the SNP information of the probes:
#The return object is a matrix with the columns being the samples and the rows being the different SNP probes:
#These SNP probes are intended to be used for sample tracking and sample mixups.
#Each SNP probe ought to have values clustered around 3 distinct values corresponding to homo-, and hetero-zygotes.
snps <- getSnpInfo(GRset)
#Probe, CpG and SBE correspond the SNPs present inside the probe body, at the CpG interrogation and at the single nucleotide extension respectively.
#The columns with rs give the names of the SNPs while the columns with maf gives the minor allele frequency of the SNPs based on the dbSnp database.
#The function addSnpInfo will add to the GenomicRanges of the GenomicRatioSet the 6 columns:
GRsetsnp <- addSnpInfo(GRset)
#Here is an example where we drop the probes containing a SNP at the CpG interrogation and/or at the single nucleotide extension, for any minor allele frequency:
GRset.snp <- dropLociWithSnps(GRsetsnp, snps=c("SBE","CpG"), maf=0)

# Load necessary package
library(minfi)

# Drop CpG loci with SNPs
GRset.snp <- dropLociWithSnps(GRsetsnp, snps = c("SBE", "CpG"), maf = 0)

# Get the list of CpG loci before and after filtering
cpg_nosnp <- featureNames(GRsetsnp)
cpg_aftersnp <- featureNames(GRset.snp)

# Find the CpG loci that were dropped
droppedSNP_cpg <- setdiff(cpg_before, cpg_after)

# Inspect the dropped CpG list
length(droppedSNP_cpg) # 14623 CpGs exclude based on SNP level 

save(droppedSNP_cpg, file="--/data/4_Processed Data/droppedSNP_cpg.rda") 
load("--/data/4_Processed Data/droppedSNP_cpg.rda")


#### Excluding probes ####

# Build RGSet with dropped probes
RGset.drop <- subsetByLoci(rgSet = RGSet,
                           excludeLoci = c(cross.probes.info$IlmnID,   # 30627 probes
                                           exclude.hpv,                # 23581 probes
                                           exclude.bds,                # 11079 probes
                                           droppedSNP_cpg))            # 14623 probes

# Build RGSet with dropped probes, but keeping cross-reactive probes in cause those are needed for DNAmps creation
RGset.drop <- subsetByLoci(rgSet = RGSet,
                           excludeLoci = c(exclude.hpv,                # 23581 probes
                                           exclude.bds,                # 11079 probes
                                           droppedSNP_cpg))            # 14623 probes

RGset.drop.inclCrossReac <- RGset.drop
save(RGset.drop.inclCrossReac, file="--/data/4_Processed Data/RGset.drop.inclCrossReac.rda") 
load("--/data/4_Processed Data/RGset.drop.inclCrossReac.rda")


save(RGset.drop, file="--/data/4_Processed Data/RGset.drop.rda") 
load("--/data/4_Processed Data/RGset.drop.rda")

# Check N excluded probes
excl_cpg_sites <- unique(c(cross.probes.info$IlmnID, exclude.hpv, exclude.bds, droppedSNP_cpg)) #Excluding in total 73643 probes
save(excl_cpg_sites, file="--/data/4_Processed Data/excl_cpg_sites.rda") 
load("--/data/4_Processed Data/excl_cpg_sites.rda")


# Explore number of excluded probes per exclusion critaria , and their overlap
unique_to_Ncrossprobes <- setdiff(cross.probes.info$IlmnID, union(exclude.hpv, union(exclude.bds, droppedSNP_cpg)))
unique_to_NdetP <- setdiff(exclude.hpv, union(cross.probes.info$IlmnID, union(exclude.bds, droppedSNP_cpg)))
unique_to_Nbeads <- setdiff(exclude.bds, union(cross.probes.info$IlmnID, union(exclude.hpv, droppedSNP_cpg)))
unique_to_droppedSNP_cpg <- setdiff(droppedSNP_cpg, union(cross.probes.info$IlmnID, union(exclude.hpv, exclude.bds)))

# Calculate overlaps
overlap_Ncrossprobes_NdetP <- intersect(cross.probes.info$IlmnID, exclude.hpv)
overlap_Ncrossprobes_Nbeads <- intersect(cross.probes.info$IlmnID, exclude.bds)
overlap_Ncrossprobes_droppedSNP <- intersect(cross.probes.info$IlmnID, droppedSNP_cpg)
overlap_NdetP_Nbeads <- intersect(exclude.hpv, exclude.bds)
overlap_NdetP_droppedSNP <- intersect(exclude.hpv, droppedSNP_cpg)
overlap_Nbeads_droppedSNP <- intersect(exclude.bds, droppedSNP_cpg)
overlap_all <- intersect(intersect(cross.probes.info$IlmnID, exclude.hpv), intersect(exclude.bds, droppedSNP_cpg))

# Calculate total CpG sites to exclude
total_to_exclude <- unique(c(cross.probes.info$IlmnID, exclude.hpv, exclude.bds, droppedSNP_cpg))
length(total_to_exclude) #Excluding in total 73643 probes

# Print results
cat("Unique to Ncrossprobes (cross.probes.info$IlmnID):", length(unique_to_Ncrossprobes), "\n")
cat("Unique to NdetP (exclude.hpv):", length(unique_to_NdetP), "\n")
cat("Unique to Nbeads (exclude.bds):", length(unique_to_Nbeads), "\n")
cat("Unique to droppedSNP_cpg:", length(unique_to_droppedSNP_cpg), "\n")
cat("Overlap between Ncrossprobes and NdetP:", length(overlap_Ncrossprobes_NdetP), "\n")
cat("Overlap between Ncrossprobes and Nbeads:", length(overlap_Ncrossprobes_Nbeads), "\n")
cat("Overlap between Ncrossprobes and droppedSNP_cpg:", length(overlap_Ncrossprobes_droppedSNP), "\n")
cat("Overlap between NdetP and Nbeads:", length(overlap_NdetP_Nbeads), "\n")
cat("Overlap between NdetP and droppedSNP_cpg:", length(overlap_NdetP_droppedSNP), "\n")
cat("Overlap between Nbeads and droppedSNP_cpg:", length(overlap_Nbeads_droppedSNP), "\n")
cat("Overlap among all four:", length(overlap_all), "\n")        # 2 CpG sites overlapp all exclusion criteria
cat("Total CpG sites to exclude:", length(total_to_exclude), "\n")

#### Display the Venn diagram excluded probes ####
venn.plot <- venn.diagram(
  x = list(
    Ncrossprobes = cross.probes.info$IlmnID,
    NdetP = exclude.hpv,
    Nbeads = exclude.bds,
    droppedSNP_cpg = droppedSNP_cpg
  ),
  filename = NULL,
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(0, 0, 0, 0),  
  cat.dist = c(0.20, 0.15, 0.10, 0.10),  # Adjust the distance of labels
  main = "Venn Diagram of CpG Sites"
)

grid.draw(venn.plot)


### Check overlap deleted probes and clocks of interest ####

# read in list with probes included in clocks
# retrieved from this website: https://github.com/bio-learn/biolearn/blob/master/biolearn/data/DunedinPACE.csv
# Note: for Grimage CpGs are not publically available, so not sure whether these are final probes, see https://pmc.ncbi.nlm.nih.gov/articles/PMC11245009/

DunedinProbes <- read.csv("--/data/4_Processed Data/DunedinPACE_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
GrimAgeProbes  <- read.csv("--/data/4_Processed Data/GrimAgev1_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
PhenoAgeProbes  <- read.csv("--/data/4_Processed Data/PhenoAge_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
EpigeneticGProbes<- fread("--/data/4_Processed Data/g_processed_Mean_Beta_PIP.txt")

# View the data
colnames(PhenoAgeProbes)
colnames(DunedinProbes)
colnames(GrimAgeProbes)
colnames(EpigeneticGProbes)

#Get CpG list of those removed due to QC, remove _XXX EpicV2 add on and remove duplicates
excl_cpg_sites_cleaned <- sub("_.*", "", excl_cpg_sites)
excl_cpg_sites_unique <- unique(excl_cpg_sites_cleaned)
head(excl_cpg_sites_unique)

# Count how many of excl_cpg_sites_unique are in DunedinProbes$CpGmarker
num_in_dunedin <- sum(excl_cpg_sites_unique %in% DunedinProbes$CpGmarker)
num_in_dunedin # 2 probes of 174 probes

# Count how many of excl_cpg_sites_unique are in PhenoAgeProbes$CpGmarker
num_in_phenoage <- sum(excl_cpg_sites_unique %in% PhenoAgeProbes$CpGmarker)
num_in_phenoage # 0 probes of 514 probes

# Count how many of excl_cpg_sites_unique are in GrimAgeProbes$var 
num_in_grimage <- sum(excl_cpg_sites_unique %in% GrimAgeProbes$var)
num_in_grimage # 16 of 1141 probes

# Count how many of excl_cpg_sites_unique are in EpigeneticGProbes$CpG
num_in_EpiG <- sum(excl_cpg_sites_unique %in% EpigeneticGProbes$CpG)
num_in_EpiG # 13694 of 764525 probes

# Count how many of the cross-reactive probes are in DunedinProbes
common_probesPaceCR <- intersect(cross.probes.info$Name, DunedinProbes$CpGmarker)

# Count how many are in both
num_common_probes <- length(common_probesPaceCR)



#### Read in RGsets in case you do not rerun the whole script above

#RGset Raw files
load(file.path("--/data/4_Processed Data/RGset_full.rda"))

annotation(RGSet)
annotation(RGSet)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGSet)["annotation"] = "20a1.hg38"
annotation(RGSet)

#RGset After dropped probes
load(file.path("--/data/4_Processed Data/RGset.drop.rda"))

annotation(RGset.drop)
annotation(RGset.drop)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGset.drop)["annotation"] = "20a1.hg38"
annotation(RGset.drop)


#### Get Beta values based on raw RGset  ####

beta_raw <-getBeta(RGSet)

dim(beta_raw) #N CpGs = 936990

save(beta_raw, file="--/data/4_Processed Data/beta_raw.rda") 
load("--/data/4_Processed Data/beta_raw.rda")

#### Get Beta values based on RGset with dropped probes ####

beta_dropCpG <-getBeta(RGset.drop) 
dim(beta_dropCpG) # N CpGs = 863347 (which corresponds to 936990 (CpGs raw data, see code line 233) - 73643 (CpGs to exclude after exclusion criteria, see code line 798))


save(beta_dropCpG, file="--/data/4_Processed Data/beta_dropCpG.rda") 
load("--/data/4_Processed Data/beta_dropCpG.rda")


#### Noob normalization ####

#Noob normalize the sample where you excluded probes based on probes detection >.01 and cross reactive probes
noob <- preprocessNoob(RGset.drop, offset = 15, dyeCorr = TRUE, verbose = TRUE)

save(noob, file = file.path("--/data/4_Processed Data/noob_RGset.drop.rda"))

#Noob normalize the sample where you excluded probes but included cross reactive probes
noob.CrossReac <- preprocessNoob(RGset.drop.inclCrossReac, offset = 15, dyeCorr = TRUE, verbose = TRUE)


#Pull phenotype data from noob object
pheno <- pData(noob)
#Extract beta values from the noob object
beta.dropCpG.noob <-getBeta(noob)
save(beta.dropCpG.noob, file = file.path("--/data/4_Processed Data/beta.dropCpG.noob.rda"))

#create beta file including cross reactive probes
beta.dropCpG.IncCrossReac.noob <-getBeta(noob.CrossReac)
dim(beta.dropCpG.IncCrossReac.noob) # N= 891136 probes, which is higher then beta.dropCpG.noob cause it includes cross reactive probes

#### Deal with duplicates in EpicV2 probes ####

#EPIC V2 beadchips contain some probes that are duplicated multiple times across the array (mostly probes important in cancer biology), see https://moffittcaspi.trinity.duke.edu/dunedinpace/validation-dunedinpace-epic-v2-data
#Two probes in the DunedinPACE algorithm are duplicated (cg26180383 and cg06230206), and if left uncorrected might lead to potentially incorrect DunedinPACE estimates by ‘doubling up’ the values from these probes. 
#An additional complication is that probe names (IlmnID) have suffixes to cg identifiers reflecting design information. 
#If suffixes are not removed from data in pipelines using IlmnID as probe identifiers, calculations will fail. To deal with these issues we did the following:
#we used the ‘rm.cgsuffix’ function in the ‘EN.mix’ package (v1.36.08) to remove suffixes and average methylation values from duplicates into one value.


#for betas where CpGs are dropped based on QC and are noob normalized
beta.dropCpG.noob.nodupl <- rm.cgsuffix(beta.dropCpG.noob)
dim(beta.dropCpG.noob.nodupl) #N CpGs 857807 (So 5540 less than beta.dropCpG.noob. )

save(beta.dropCpG.noob.nodupl, file="--/data/4_Processed Data/beta.dropCpG.noob.nodupl.rda") 
load("--/data/4_Processed Data/beta.dropCpG.noob.nodupl.rda")

#for betas where CpGs are dropped, but includes Cross-reactive probes
beta.dropCpG.IncCrossReac.noob.nodupl <- rm.cgsuffix(beta.dropCpG.IncCrossReac.noob)
dim(beta.dropCpG.IncCrossReac.noob.nodupl) #885195 CpGs , which is higher then beta.dropCpG.noob.nodupl cause it includes cross-reactivee probes


# For reference:
# Karen Sudgen shows the following N of probes after QC and removing duplicates: https://moffittcaspi.trinity.duke.edu/dunedinpace/validation-dunedinpace-epic-v2-data
# After QC and removal of duplicate probes the datasets contained 917,069 and 930,659 probes in the wateRmelon/methyumi and SeSAMe datasets, respectively.

#### Excluding samples ####

# We generate two beta datasets:

# For EWAS: 1) beta's with dropped cpgs that did not pass QC that are noob normalized and excluding samples that do not pass QC
# For Clocks: 1) beta's with dropped cpgs that did not pass QC that are noob normalized and excluding samples that do not pass QC 2) same as first, but also collapsing duplicates and taking the mean (see Karen Sugden link above)
# For Clocks we create these two sets as for some clocks cannot be computed with EpicV2 duplicate probes. 


#For EWAS and Clocks

# 1)  beta's with dropped cpgs that are noob normalized and excluding samples that do not pass QC

beta.dropCpG.noob.dropSamples <- beta.dropCpG.noob[, !colnames(beta.dropCpG.noob) %in% c("208226080076_R07C01", "208226080046_R03C01",     # 8 samples excluded based 
                                                                                         "208226080024_R06C01", "208226080085_R08C01",     # based on low intensity log <9 and detection p >.01     
                                                                                         "208290890120_R08C01", "208290890038_R04C01",
                                                                                         "208291390053_R08C01", "208291390095_R01C01")]


save(beta.dropCpG.noob.dropSamples, file="--/data/4_Processed Data/beta.dropCpG.noob.dropSamples.rda") 
load("--/data/4_Processed Data/beta.dropCpG.noob.dropSamples.rda")

dim(beta.dropCpG.noob.dropSamples) # 863347 probes, 1528 samples

# For Clocks
# 2)  beta's with dropped cpgs that are noob normalized and excluding samples that do not pass QC and no duplicates

beta.dropCpG.noob.nodupl.dropSamples <- beta.dropCpG.noob.nodupl[, !colnames(beta.dropCpG.noob.nodupl) %in% c("208226080076_R07C01", "208226080046_R03C01",     # 8 samples excluded based 
                                                                                                              "208226080024_R06C01", "208226080085_R08C01",     # based on low intensity log <9 and detection p >.01     
                                                                                                              "208290890120_R08C01", "208290890038_R04C01",
                                                                                                              "208291390053_R08C01", "208291390095_R01C01")]


save(beta.dropCpG.noob.nodupl.dropSamples, file="--/data/4_Processed Data/beta.dropCpG.noob.nodupl.dropSamples.rda") 
load("--/data/4_Processed Data/beta.dropCpG.noob.nodupl.dropSamples.rda")

dim(beta.dropCpG.noob.nodupl.dropSamples) # 857807 probes, 1528 samples

# 3)  beta's with dropped cpgs that are noob normalized and excluding samples that do not pass QC and no duplicates, this one Includes cross-reactive probes

beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples <- beta.dropCpG.IncCrossReac.noob.nodupl[, !colnames(beta.dropCpG.IncCrossReac.noob.nodupl) %in% c("208226080076_R07C01", "208226080046_R03C01",     # 8 samples excluded based 
                                                                                                                                                     "208226080024_R06C01", "208226080085_R08C01",     # based on low intensity log <9 and detection p >.01     
                                                                                                                                                     "208290890120_R08C01", "208290890038_R04C01",
                                                                                                                                                     "208291390053_R08C01", "208291390095_R01C01")]
dim(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples) # 885195 probes, 1528 samples

save(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples, file="--/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.rda") 
load("--/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.rda")

# 4)  beta's with dropped cpgs that are noob normalized and excluding samples that do not pass QC and but this includes duplicates, this one also Includes cross-reactive probes and duplicates

beta.dropCpG.IncCrossReac.noob.dropSamples <- beta.dropCpG.IncCrossReac.noob[, !colnames(beta.dropCpG.IncCrossReac.noob) %in% c("208226080076_R07C01", "208226080046_R03C01",     # 8 samples excluded based 
                                                                                                                                "208226080024_R06C01", "208226080085_R08C01",     # based on low intensity log <9 and detection p >.01     
                                                                                                                                "208290890120_R08C01", "208290890038_R04C01",
                                                                                                                                "208291390053_R08C01", "208291390095_R01C01")]

dim(beta.dropCpG.IncCrossReac.noob.dropSamples) #CPGs = 891136   Sample = 1528

save(beta.dropCpG.IncCrossReac.noob.dropSamples, file="--/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.dropSamples.rda") 
load("--/data/4_Processed Data/beta.dropCpG.IncCrossReac.noob.dropSamples.rda")

#Double check if this indeed excludes all probes of QC 
beta_values <- rownames(beta.dropCpG.IncCrossReac.noob.dropSamples)  # Extract CpG IDs

# Find overlap between beta values and exclude.hpv
overlap <- intersect(beta_values, exclude.hpv)
overlap <- intersect(beta_values, exclude.bds)
overlap <- intersect(beta_values, droppedSNP_cpg)

print(overlap) #there is no overlap with QC probes, so this shows that all went well! 

# Extract CpG IDs from beta.dropCpG.IncCrossReac.noob.dropSamples
beta_values <- rownames(beta.dropCpG.IncCrossReac.noob.dropSamples)

# Extract IlmnID from cross.probes.info
cross_probes <- cross.probes.info$IlmnID

# Find overlap between beta_values and cross_probes
overlap_with_cross <- intersect(beta_values, cross_probes) #27789 probes of cross reactive probes are in there, which makes sense that its not all 30627 cause some of those do not meet QC

#### Create phenodata file for final dataset ####

#Pull phenotype data from noob object
pheno <- pData(noob)

#Select samples to exclude
ids_to_exclude <- c("208226080076_R07C01", "208226080046_R03C01",
                    "208226080024_R06C01", "208226080085_R08C01",         
                    "208290890120_R08C01", "208290890038_R04C01",
                    "208291390053_R08C01", "208291390095_R01C01")  

# Filter out rows with the specified IDs
pheno.dropSamples <- pheno[!rownames(pheno) %in% ids_to_exclude, ]
save(pheno.dropSamples, file = file.path("--/data/4_Processed Data/pheno.dropSamples.rda"))
load("--/data/4_Processed Data/pheno.dropSamples.rda")



# Check descriptives in the sample, without the technical replicates

dim(pheno.dropSamples) #n = 1528, in total we have 1528 Samples
table(pheno.dropSamples$Techreplicate) #1512 samples and 16 replicates


#### Ensure that sample order matches of Phenodata and Beta data ####

# Extract the second dimension (column names or sample IDs) from beta matrix
beta_sample_ids <- colnames(beta.dropCpG.noob.nodupl.dropSamples)  # dimnames[[2]] corresponds to colnames

# Check if the order matches rownames(pheno)
is_same_order <- all.equal(beta_sample_ids, rownames(pheno.dropSamples))

# This shows that the order is the same
# The order of sample IDs in the beta matrix matches the row names of pheno.
if (is_same_order == TRUE) {
  cat("The order of sample IDs in the beta matrix matches the row names of pheno.\n")
} else {
  cat("The order of sample IDs in the beta matrix does NOT match the row names of pheno.\n")
  print(is_same_order)  # Show any mismatch details
}



#### Density Plots ####

# Generate the densityBeanPlot, with child/mother grouping
pdf("--/data analysis/00_QC/Plots_full/DensityPlot_beta.dropCpG.noob.pdf") #open PDF device to save the plot


#Density plot on noob normalized betas with dropped CpG, but without excluding samples
plot <- densityPlot(beta.dropCpG.noob, sampGroups = pheno$Mum0Child1, main = "BFY", xlab = "Beta",
                    legend = TRUE)

dev.off() #closes the PDF device, ensuring that the file is properly saved.

# Generate the densityBeanPlot, with child/mother grouping
pdf("--/data analysis/00_QC/Plots_full/DensityPlot_beta.dropCpG.noob.dropSamples.pdf") #open PDF device to save the plot


#Density plot on noob normalized betas, with dropped CpGs and with excluding samples that did not pass QC
plot <- densityPlot(beta.dropCpG.noob.dropSamples,sampGroups = pheno.dropSamples$Mum0Child1,  main = "BFY", xlab = "Beta",
                    legend = TRUE)

dev.off() #closes the PDF device, ensuring that the file is properly saved.

# Generate the densityBeanPlot, with child/mother grouping
pdf("--/data analysis/00_QC/Plots_full/DensityPlot_beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples.pdf") #open PDF device to save the plot

#Density plot on noob normalized betas, with dropped CpGs and with excluding samples that did not pass QC
plot <- densityPlot(beta.dropCpG.IncCrossReac.noob.nodupl.dropSamples,sampGroups = pheno.dropSamples$Mum0Child1,  main = "BFY", xlab = "Beta",
                    legend = TRUE)

dev.off() #closes the PDF device, ensuring that the file is properly saved.
dev.off() #closes the PDF device, ensuring that the file is properly saved.

# Generate the densityBeanPlot, with child/mother grouping
pdf("--/data analysis/00_QC/Plots_full/DensityPlot_beta.dropCpG.IncCrossReac.noob.dropSamples.pdf") #open PDF device to save the plot

#Density plot on noob normalized betas, with dropped CpGs and with excluding samples that did not pass QC
plot <- densityPlot(beta.dropCpG.IncCrossReac.noob.dropSamples,sampGroups = pheno.dropSamples$Mum0Child1,  main = "BFY", xlab = "Beta",
                    legend = TRUE)

dev.off() #closes the PDF device, ensuring that the file is properly saved.


#### check overlap CpGs with clocks we are interested in ####

#For DunedinPace
# Find CpG markers that are in the final beta dataset
common_CpGs <- intersect(DunedinProbes$CpGmarker, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many are in the final beta dataset
num_common <- length(common_CpGs)

# Find CpG markers that are NOT in the final beta dataset
not_in_CpGs <- setdiff(DunedinProbes$CpGmarker, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many are not in the final beta dataset
num_not_in <- length(not_in_CpGs)

# Print the results
cat("Number of CpG markers IN the final beta dataset:", num_common, "\n") # 142 CpGs
cat("Number of CpG markers NOT in the final beta dataset:", num_not_in, "\n") # 32 CpGs are missing = 18%

# For PhenoAge
common_PhenoAge <- intersect(PhenoAgeProbes$CpGmarker, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many PhenoAgeProbes are in the final beta dataset
num_common_PhenoAge <- length(common_PhenoAge)

# Find PhenoAgeProbes that are NOT in the final beta dataset
not_in_PhenoAge <- setdiff(PhenoAgeProbes$CpGmarker, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many PhenoAgeProbes are not in the final beta dataset
num_not_in_PhenoAge <- length(not_in_PhenoAge)

# Print the results
cat("Number of PhenoAge CpG markers IN the final beta dataset:", num_common_PhenoAge, "\n") # 495 CpGs
cat("Number of PhenoAge CpG markers NOT in the final beta dataset:", num_not_in_PhenoAge, "\n") # 19 CpGs are missing = 4%


# GrimAge
common_GrimAge <- intersect(GrimAgeProbes$var, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many GrimAgeProbes are in the final beta dataset
num_common_GrimAge <- length(common_GrimAge)

# Find GrimAgeProbes that are NOT in the final beta dataset
not_in_GrimAge <- setdiff(GrimAgeProbes$var, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many GrimAgeProbes are not in the final beta dataset
num_not_in_GrimAge <- length(not_in_GrimAge)

# Print the results
cat("Number of GrimAge CpG markers IN the final beta dataset:", num_common_GrimAge, "\n") # 829 CpGs
cat("Number of GrimAge CpG markers NOT in the final beta dataset:", num_not_in_GrimAge, "\n") # 216 missing = 21%


# Count how many elements in GrimAgeProbes$var start with "cg"
num_start_with_cg <- sum(grepl("^cg", GrimAgeProbes$var))

# Print the result (in total 1045 cpgs in GrimAge)
cat("Number of elements in GrimAgeProbes$var that start with 'cg':", num_start_with_cg, "\n")


# EpigeneticGProbes
common_EpigeneticG <- intersect(EpigeneticGProbes$CpG, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many EpigeneticGProbes are in the final beta dataset
num_common_EpigeneticG <- length(common_EpigeneticG)

# Find EpigeneticGProbes that are NOT in the final beta dataset
not_in_EpigeneticG <- setdiff(EpigeneticGProbes$CpG, row.names(beta.dropCpG.noob.nodupl.dropSamples))

# Count how many EpigeneticGProbes are not in the final beta dataset
num_not_in_EpigeneticG <- length(not_in_EpigeneticG)

# Print the results
cat("Number of EpigeneticG CpG markers IN the final beta dataset:", num_common_EpigeneticG, "\n") # 678199 CpGs
cat("Number of EpigeneticG CpG markers NOT in the final beta dataset:", num_not_in_EpigeneticG, "\n") # 86326 CpGs missing = 11%


#### Do Mean Imputation on Betas for NA's ####

#not necessary, cause there are no NA's in CpGs of our samples

beta.dropCpG.noob.nodupl.dropSamples.df <- as.data.frame(beta.dropCpG.noob.nodupl.dropSamples)

# Replace `your_dataframe` with the name of your dataframe
na_counts <- colSums(is.na(beta.dropCpG.noob.nodupl.dropSamples.df))

# Filter participants with 1 or more missing values
participants_with_na <- na_counts[na_counts > 0]

# Print the participants with their NA counts
print(participants_with_na) #no samples have NA


#### Cell composition Ref-free ####

#Package ‘RefFreeEWAS’ was removed from the CRAN repository.
#Formerly available versions can be obtained from the archive (https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/)

if (!requireNamespace("BiocManager", quietly = TRUE))  
  install.packages("BiocManager")  
BiocManager::install("RefFreeEWAS")  
install.packages("RefFreeEWAS")
install.packages("quadprog")
install.packages("RefFreeEWAS_2.2.tar.gz", repos = NULL, type = "source")  

library(quadprog)
library(RefFreeEWAS)
library(RefFreeEWAS)

# Set the working directory  
setwd("--/data/4_Processed Data/")  

### Saliva Reffree cell estimation----------------

# It takes some time to run the script, first try with 8 samples if the script works before running on the full sample
#beta_subset <- beta.dropCpG.noob.nodupl.dropSamples.noNA[, 1:8]

#calling the RefFreeCellMixArray function, which is used in reference-free cell mixture adjustment.
sa_cell_reffree = RefFreeCellMixArray(as.matrix(beta.dropCpG.noob.nodupl.dropSamples),Klist=5,iters=5,Yfinal=as.matrix(beta.dropCpG.noob.nodupl.dropSamples))

sa_cell_reffreeOmega = sa_cell_reffree$'5'$Omega

# So it created cell estimation for 5 cell types
head(sa_cell_reffreeOmega)
dim(sa_cell_reffreeOmega)

save(sa_cell_reffreeOmega, file = file.path("--/data/4_Processed Data/beta.dropCpG.noob.nodupl.dropSamples.sa_cell_reffreeOmega.rda"))
load("--/data/4_Processed Data/beta.dropCpG.noob.nodupl.dropSamples.sa_cell_reffreeOmega.rda")

head(sa_cell_reffreeOmega)
sa_cell_reffreeOmega.df <- as.data.frame(sa_cell_reffreeOmega)

save(sa_cell_reffreeOmega.df, file = file.path("--/data/4_Processed Data/sa_cell_reffreeOmega.df.rda"))
load("--/data/4_Processed Data/sa_cell_reffreeOmega.df.rda")

# Compute summary statistics for all columns in sa_cell_reffreeOmega.df
summary_table <- data.frame(
  CellType = colnames(sa_cell_reffreeOmega.df),
  Mean = sapply(sa_cell_reffreeOmega.df, mean, na.rm = TRUE),
  SD = sapply(sa_cell_reffreeOmega.df, sd, na.rm = TRUE),
  Min = sapply(sa_cell_reffreeOmega.df, function(x) min(x, na.rm = TRUE)),  # Ensure it captures negative values
  Max = sapply(sa_cell_reffreeOmega.df, function(x) max(x, na.rm = TRUE))
)

# Print the table
cat("Summary Table (Mean, SD, Min, Max):\n")
print(summary_table)

# Optional: Format the table nicely for presentation
library(knitr)
kable(summary_table, format = "markdown", col.names = c("Cell Type", "Mean", "SD", "Min", "Max"))


# Define a function to flag values outside 3 IQR from the median
flag_outliers <- function(x) {
  med <- median(x, na.rm = TRUE)        # Calculate the median
  iqr <- IQR(x, na.rm = TRUE)          # Calculate the interquartile range (IQR)
  lower_bound <- med - 3 * iqr         # Lower threshold
  upper_bound <- med + 3 * iqr         # Upper threshold
  x < lower_bound | x > upper_bound    # Return TRUE for outliers
}

# Create a list to store outliers for each column
outliers_by_column_cellcomp <- lapply(1:5, function(col_idx) {
  outliers <- flag_outliers(sa_cell_reffreeOmega[, col_idx])  # Check for outliers in the column
  rownames(sa_cell_reffreeOmega)[outliers]                   # Get rownames for outliers
})

# Name the list by column indices
names(outliers_by_column_cellcomp) <- paste0("Column_", 1:5)

# View outliers for each column
outliers_by_column_cellcomp

# Combine outliers across all columns into one list
all_outliers_CellComp <- unique(unlist(outliers_by_column_cellcomp))

save(outliers_by_column_cellcomp, file = file.path("--/data/4_Processed Data/outliers_by_column_cellcomp.rda"))
load("--/data/4_Processed Data/outliers_by_column_cellcomp.rda")

save(all_outliers_CellComp, file = file.path("--/data/4_Processed Data/all_outliers_CellComp.rda"))
load("--/data/4_Processed Data/all_outliers_CellComp.rda")


#### Create phenotype file including new variables ####
library(tibble)

pheno.dropSamples.df <- data.frame(pheno.dropSamples)
head(pheno.dropSamples.df)
colnames(pheno.dropSamples.df)

# Turn rowname into ID, for later merging
# Add row names as a new column while keeping them as row names
pheno.dropSamples.df$BaseID <- rownames(pheno.dropSamples.df)

# Reorder the column to make it the first column
pheno.dropSamples.df <- pheno.dropSamples.df[, c("BaseID", setdiff(names(pheno.dropSamples.df), "BaseID"))]

#Add labels to phenotype data, so its easy later to remember what they stand for
colnames(pheno.dropSamples.df)

#### Load in Basic Phenotype Data where samples are dropped ####
save(pheno.dropSamples.df, file = file.path("--/data/4_Processed Data/pheno.dropSamples.df.rda"))
load("--/data/4_Processed Data/pheno.dropSamples.df.rda")

colnames(pheno.dropSamples.df)
head(pheno.dropSamples.df)

### add cell reference to phenotypic dataset ####

colnames(pheno.dropSamples.df)

load("--/data/4_Processed Data/sa_cell_reffreeOmega.df.rda")

# Rename columns  1, 2, 3, 4, 5 in the cell composition data to celltype_1, celltype_2 etc. 
colnames(sa_cell_reffreeOmega.df) <- gsub(
  "^(\\d+)$", # Matches column names that are digits only
  "celltype_\\1", # Replaces with "celltype_" followed by the number
  colnames(sa_cell_reffreeOmega.df)
)

colnames(sa_cell_reffreeOmega.df)

# Step 1: Convert rownames of sa_cell_reffreeOmega.df to a column
sa_cell_reffreeOmega.df$BaseID <- rownames(sa_cell_reffreeOmega.df)

# Step 2: Merge cell composition 
pheno.dropSamples.df <- pheno.dropSamples.df %>%
  left_join(sa_cell_reffreeOmega.df, by = "BaseID")


#### Flag cell composition mismatch in phenotype data ####

#add outliers based on cell composition
load("--/data/4_Processed Data/all_outliers_CellComp.rda")

# Create the new variable Flag_CellComp
pheno.dropSamples.df$Flag_CellComp <- ifelse(pheno.dropSamples.df$BaseID %in% all_outliers_CellComp, 1, 0)


# View the updated dataset
#only the fifth cell type has outliers, so keep this in mind when doing sensitivity checks
# these have 31 outliers, but probably this is a very rare celltype, so don't worry about it too much
table(pheno.dropSamples.df$Flag_CellComp) 

head(pheno.dropSamples.df)
colnames(pheno.dropSamples.df)


#### Flag sex mismatch in phenotype data ####
#add flag mismatch based on sex

annotation(RGSet) #should be right annotation, loaded in earlier


manifest <- getManifest(RGSet)      
annotation <- getAnnotation(RGSet)

Meth_manifest <- data.frame(
  index = rownames(annotation),  # Probe IDs
  chr = annotation$chr           # Chromosome information
)

# Check the first few rows
head(Meth_manifest)

MSet <- preprocessRaw(RGSet)

#Get methylation
Meth <-getMeth(MSet) # M signal per probe, per sample
Unmeth <-getUnmeth(MSet) # U signal per probe, per sample



Meth <- list(
  manifest = Meth_manifest,  # The manifest you just created
  M = Meth,                # Matrix of methylated intensities
  U = Unmeth                 # Matrix of unmethylated intensities
)



# Find common probes
common_probes <- intersect(rownames(Meth$M), Meth$manifest$index)

# Check the number of common probes
length(common_probes)  # Should match the length of `Meth$manifest$index` if alignment is correct

# Subset Meth$manifest
Meth$manifest <- Meth$manifest[Meth$manifest$index %in% common_probes, ]

# Subset Meth$M and Meth$U
Meth$M <- Meth$M[common_probes, , drop = FALSE]
Meth$U <- Meth$U[common_probes, , drop = FALSE]

# Assign proper rownames to Meth$M and Meth$U
rownames(Meth$M) <- Meth$manifest$index
rownames(Meth$U) <- Meth$manifest$index

# Verify dimensions
dim(Meth$M)  # Should now match `length(Meth$manifest$index)`
dim(Meth$U)
length(Meth$manifest$index)

# Confirm that all rownames align
all(rownames(Meth$M) == Meth$manifest$index)  # Should return TRUE
all(rownames(Meth$U) == Meth$manifest$index)


# Select probes by chromosome
chrX <- Meth$manifest[Meth$manifest$chr == 'chrX', "index"]
chrY <- Meth$manifest[Meth$manifest$chr == 'chrY', "index"]
auto <- Meth$manifest[!Meth$manifest$chr %in% c("chrX", "chrY"), "index"]

# Compute total intensities
chrX <- Meth$M[chrX, , drop = FALSE] + Meth$U[chrX, , drop = FALSE]
chrY <- Meth$M[chrY, , drop = FALSE] + Meth$U[chrY, , drop = FALSE]
auto <- Meth$M[auto, , drop = FALSE] + Meth$U[auto, , drop = FALSE]

# Compute per-sample averages
chrX <- colMeans(chrX, na.rm = TRUE)
chrY <- colMeans(chrY, na.rm = TRUE)
auto <- colMeans(auto, na.rm = TRUE)

# Normalize allosomal intensities
chrX <- chrX / auto
chrY <- chrY / auto

# Store results in phenotype data
pd <- data.frame(SampleID = colnames(Meth$M))
pd$X <- chrX
pd$Y <- chrY

# View results
head(pd)

#merge two datasets
GRset <-mapToGenome(MSet)
pheno <- pData(GRset)  
# Ensure pheno is a data frame
pheno <- as.data.frame(pheno)

# Add rownames of pheno as a column for merging
pheno$SampleID <- rownames(pheno)

# Select only the desired columns for merging
pheno_subset <- pheno[, c("SampleID",  "Child.is.Female", "Code.C..child..M..mother.", "sex")]

#"Sample_Name"

# Merge the subset of pheno into pd by SampleID
pd <- merge(pd, pheno_subset, by = "SampleID", all.x = TRUE)

# View the merged dataset
head(pd)
colnames(pd)

#add sex prediction
pd$Predicted_Sex <- ifelse(abs(pd$X - pd$Y) <= 0.05, 0,  # Predict Male (0) if X and Y are close
                           ifelse(pd$X > pd$Y, 1, NA))   # Predict Female (1) if X > Y, NA otherwise




# Add a column to flag mismatches
pd$Mismatch <- pd$Predicted_Sex != pd$sex

# View mismatches
mismatches <- pd[pd$Mismatch, ]

# Count mismatches
table(pd$Mismatch)

# Find rows where Mismatch is NA
na_rows <- which(is.na(pd$Mismatch))

# View the rows with NA in Mismatch
pd_with_na <- pd[na_rows, ]

# Print the rows or save them to a file
print(pd_with_na)

#so we discovered there was a sample swop, lets merge that phenotype data here

load("--/data/4_Processed Data/pheno_Only_Swapped.rda")

colnames(pheno_Only_Swapped)

pheno_Only_Swapped_sex <- pheno_Only_Swapped[,c("BaseID",
                                                "sex")]

colnames(pheno_Only_Swapped_sex )[colnames(pheno_Only_Swapped_sex ) == "BaseID"] <- "SampleID"
colnames(pheno_Only_Swapped_sex )[colnames(pheno_Only_Swapped_sex ) == "sex"] <- "sex_correct"


pd <- pd %>%
  left_join(pheno_Only_Swapped_sex, by = "SampleID")


# Add a column to flag mismatches
pd$Mismatch2 <- pd$Predicted_Sex != pd$sex_correct

# View mismatches
mismatches2 <- pd[pd$Mismatch2, ]

# Count mismatches
table(pd$Mismatch2)

#so these are potentail sex mismatch

#ID 208250190118_R07C01 and sample ID 208290890102_R06C01


#Plot gender mismatch

library(ggplot2)

ggplot(pd, aes(x = X, y = Y, color = as.factor(Mismatch))) +
  geom_point(size = 3) +
  labs(title = "Normalized Intensities and Sex Mismatches",
       x = "Normalized Chromosome X Intensity",
       y = "Normalized Chromosome Y Intensity",
       color = "Mismatch (1=True)") +
  theme_minimal()

# flag samples that score TRUE on mismatch
# Extract SampleID where Mismatch is TRUE
mismatch_samples <- pd$SampleID[pd$Mismatch == TRUE]

# View the list of mismatched SampleIDs
print(mismatch_samples)


## these IDs are flagged ID 208250190118_R07C01 and sample ID 208290890102_R06C01

#create a flag for sex
pheno.dropSamples.df$Flag_sex <- ifelse(pheno.dropSamples.df$BaseID %in% c("208250190118_R07C01", "208290890102_R06C01"), 1, 0)

# these have 31 outliers, but probably this is a very rare celltype, so don't worry about it too much
table(pheno.dropSamples.df$Flag_sex) 


#### Load in Basic Phenotype Data where samples are dropped and cell comp is added ####
save(pheno.dropSamples.df, file = file.path("--/data/4_Processed Data/pheno.dropSamples.df_withCellComp.rda"))
load("--/data/4_Processed Data/pheno.dropSamples.df_withCellComp.rda")

pheno.dropSamples.df_CellComp_Flagg <- pheno.dropSamples.df

save(pheno.dropSamples.df_CellComp_Flagg, file = file.path("--/data/4_Processed Data/pheno.dropSamples.df_CellComp_Flagg.rda"))
load("--/data/4_Processed Data/pheno.dropSamples.df_CellComp_Flagg.rda")

colnames(pheno.dropSamples.df)

#### Sample Swap ####

# When checking the methylation age based on SkinHorvath, we saw some children are in their 20/30ies, and some mothers are around age 4.
# That cannot be, so we checked and they have potentially been swapped (see script QC_BFY_12122024_Check potential sample swaps)
# See this data file, where you see the swop within families (=Block) load("--/data/4_Processed Data/sample_swop.rda")

#so here we recreate a phenofile where we change the BaseID (=location of array) so we can later merge all the methylation dat in the correct way!

# We create two datasets, one for phenotype data (where we will change the BaseID), and one with all the epigenetic info (which we will merge again after)

# Define the columns for the first dataset
selected_columns <- c(
  "BaseID", "Sample_Name", "sampleId", "block", "Site", "Child.is.Female", 
  "sex", "Child.Sample..1.yes.", "Mother.Sample..1.yes.", 
  "Treatment.Group..blinded.", "Code.C..child..M..mother.", "Mum0Child1"
)

# Create the first dataset with selected columns, which is the pheno data
pheno_Only <- pheno.dropSamples.df[, selected_columns]
colnames(pheno_Only)

# Create the second dataset with BaseID and all remaining columns
lab_Only <- pheno.dropSamples.df[, !(colnames(pheno.dropSamples.df) %in% selected_columns) | colnames(pheno.dropSamples.df) == "BaseID"]
colnames(lab_Only)

# View the datasets
head(pheno_Only)  # this is the pheno data, e.g. sex, and location
head(lab_Only)  # this is the lab info


#in the phenodata, swop those those mother child pairs that need to be swopped
# in 11 families, mom-child samples seem to have been swapped
# List of blocks where BaseID needs to be swapped, blocks = family ID. These are the blocks in load("--/data/4_Processed Data/sample_swop.rda")
blocks_to_swap <- c(29, 56, 74, 202, 278, 368, 403, 577, 601, 691, 721)

# Iterate over each block and swap BaseID within the block
pheno_Only_Swapped <- pheno_Only %>%
  group_by(block) %>%
  mutate(
    BaseID = ifelse(block %in% blocks_to_swap, 
                    rev(BaseID),  # Reverse the BaseID order within the block
                    BaseID)
  ) %>%
  ungroup()

#Check if swap was succesful
# Filter the original and swapped datasets for the blocks to swap
original_blocks <- pheno_Only %>% filter(block %in% blocks_to_swap)
swapped_blocks <- pheno_Only_Swapped %>% filter(block %in% blocks_to_swap)

# Ensure both filtered datasets are aligned in row order
swapped_blocks <- swapped_blocks %>% arrange(block)
original_blocks <- original_blocks %>% arrange(block)

# Compare the BaseID values
comparison <- original_blocks %>%
  mutate(Swapped_BaseID = swapped_blocks$BaseID) %>%
  dplyr::select(block, BaseID, Swapped_BaseID)

# Check for differences
comparison <- comparison %>%
  mutate(Swap_Correct = BaseID != Swapped_BaseID)

# View the comparison
# All went well
comparison 

swoppedIDs <-  comparison 

#### save phenotypedata with swapped samples ####
save(swoppedIDs , file = file.path("--/data/4_Processed Data/swoppedIDs.rda"))
load("--/data/4_Processed Data/swoppedIDs.rda")


#Remerge methylation data to phenotype data
pheno_Only_Swapped  <- pheno_Only_Swapped  %>%
  mutate(BaseID = as.character(BaseID))

pheno_Only_Swapped <- pheno_Only_Swapped %>%
  left_join(lab_Only, by = "BaseID")

colnames(pheno_Only_Swapped)


save(pheno_Only_Swapped, file = file.path("--/data/4_Processed Data/pheno_Only_Swapped.rda"))
load("--/data/4_Processed Data/pheno_Only_Swapped.rda")

pheno.dropsamples.sw.df <- pheno_Only_Swapped

save(pheno.dropsamples.sw.df, file = file.path("  pheno.dropsamples.sw.df.rda"))
load("--/data/4_Processed Data/pheno_Only_Swapped.rda")


##### Additional QC specifically for EWAS ####

#### Install and load necessary packages ####
install.packages("dplyr")
install.packages("tidyr")
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("minfi", force = TRUE)
BiocManager::install("ENmix")
install.packages("locfit", dependencies = TRUE)
BiocManager::install("RnBeads")
BiocManager::install("limma", force = TRUE)
BiocManager::install("lumi", force = TRUE)
BiocManager::install("bumphunter", force = TRUE)
BiocManager::install("IlluminaHumanMethylationEPICv2manifest", force = TRUE)
install.packages("stringr")  
BiocManager::install("BiocGenerics")  
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
install.packages("remotes")  
remotes::install_github("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
install.packages("qqman")  
install.packages("ggplot2")
BiocManager::install("bacon")
if(!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("Treys925/CoMeBackV2_R_package", force = TRUE)

library(CoMeBackV2)
library(ggplot2)
library(dplyr)  # For data manipulation
library(tidyr)  # For reshaping data
library(ENmix)  # For ICC of CpGs with technical replicates
library(minfi)  # For preprocessing and analysis of Illumina methylation data
library(RnBeads)  # For comprehensive methylation data analysis and visualization
library(bumphunter)  # Dependency minfi
library(locfit)  # For local regression and smoothing
library(limma)  # For differential methylation analysis using linear models
library(lumi)  # For preprocessing and analysis of Illumina microarray data
library(IlluminaHumanMethylationEPICv2manifest)  # For annotation of EPIC v2 array probes
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38) # For annotation of EPIC v2 array probes
library(stringr)  # For string manipulation
library(BiocGenerics)  # Core Bioconductor utilities
library(remotes)  # For installing packages from remote repositories
library(qqman)  # For QQ plots and Manhattan plots in GWAS/EWAS
library(bacon) #for potential bacon adjusment for EWAS analyses
library(stats)
library(tidyverse)
library(data.table)
library(devtools)

#Before doing the EWAS, we select variable probes, probes that have an acceptable intraclass correlation and correct them for technical covariates. 
#We preregistered the following steps
#1) We select the top 50% most variable DNA Methylation sites identified based on the standard deviations per site (mothers and children separately)
#2) We will select DNA methylation sites that have an acceptable intraclass correlation coefficient (ICC) >.50 in our technical replicates (mothers and children grouped together)
#3) The selected CpG sites will be corrected for cell composition and technical covariates (plate, array) (mothers and children separately)
#4) EWAS (mothers and children separately)


### Load in Beta data for CpGs ####
# These are the betas where samples that do not meet QC are dropped, only including CpGs that meet QC standards, and betas Noob Normalized
# Note that CpG names correspond to Epic V2, so they have the _BC11 _TC11 etc annotation. So some sites have multiple CpGs annotated to same location (also refered to as duplicates)
# Note that this excludes cross reactive probes as we preregistered for the EWAS. Computation of the clocks does include cross-reactive probes!

load("beta.dropCpG.noob.dropSamples.rda")

beta_df <- as.data.frame(beta.dropCpG.noob.dropSamples)
dim(beta_df) 

#### Create beta file separately for mothers and children ####

#load in Phenodata so we know which BaseID are for kids/mums
load("data_DNAmScores_BFY_ExclTechRep_14Jan.rda")

Mothers_pheno_data <- data_DNAmScores_BFY_ExclTechRep_14Jan[data_DNAmScores_BFY_ExclTechRep_14Jan$Mum0Child1 == 0, ]
Children_pheno_data <- data_DNAmScores_BFY_ExclTechRep_14Jan[data_DNAmScores_BFY_ExclTechRep_14Jan$Mum0Child1 == 1, ]

#Create beta file only including BaseID that is in Children_pheno_data
beta_df_children <- beta_df[, colnames(beta_df) %in% Children_pheno_data$BaseID]

beta_df_mothers <- beta_df[, colnames(beta_df) %in% Mothers_pheno_data$BaseID]

#save datafiles

beta.dropCpG.noob.dropSamples.notechrep.children <- beta_df_children

save(beta.dropCpG.noob.dropSamples.notechrep.children, 
     file = "/beta.dropCpG.noob.dropSamples.notechrep.children.rda")

load("beta.dropCpG.noob.dropSamples.notechrep.children.rda")

beta.dropCpG.noob.dropSamples.notechrep.mothers <- beta_df_mothers

save(beta.dropCpG.noob.dropSamples.notechrep.mothers, 
     file = "beta.dropCpG.noob.dropSamples.notechrep.mothers.rda")

load("beta.dropCpG.noob.dropSamples.notechrep.mothers.rda")


#### Selecting top 50% most variable DNA methylation sites separately for mothers and children ####

#For Children
beta_mostvar_C <- beta.dropCpG.noob.dropSamples.notechrep.children

# Calculate Standard Deviation for Each CpG Site
beta_mostvar_C$SD <- apply(beta_mostvar_C, 1, sd)

hist(beta_mostvar_C$SD)
describe(beta_mostvar_C$SD)

CpGs_allSD_C <- data.frame(CpG_ID = rownames(beta_mostvar_C), SD = beta_mostvar_C$SD, row.names = NULL)

save(CpGs_allSD_C , 
     file = "CpGs_allSD_C .rda")

load("  CpGs_allSD_C .rda")

#Filter Top 50% Most Variable CpG Sites
threshold_C <- quantile(beta_mostvar_C$SD, 0.5)  # 50% quantile
top_50_percent_C <- beta_mostvar_C[beta_mostvar_C$SD >= threshold_C, ]

hist(top_50_percent_C$SD)
describe(top_50_percent_C$SD)

# Convert rownames (CpG sites) into a column and create a new dataframe with SD values, just for archive, and save them as csv file
CpGs_top50SD_C <- data.frame(CpG_ID = rownames(top_50_percent_C), SD = top_50_percent_C$SD, row.names = NULL)

save(CpGs_top50SD_C, 
     file = "  CpGs_top50SD_C.rda")

load("  CpGs_top50SD_C.rda")


hist(CpGs_top50SD_C$SD)
describe(CpGs_top50SD_C$SD)

#For Mothers

beta_mostvar_M <- beta.dropCpG.noob.dropSamples.notechrep.mothers

# Calculate Standard Deviation for Each CpG Site
beta_mostvar_M$SD <- apply(beta_mostvar_M, 1, sd)

hist(beta_mostvar_M$SD)
describe(beta_mostvar_M$SD)

CpGs_allSD_M <- data.frame(CpG_ID = rownames(beta_mostvar_M), SD = beta_mostvar_M$SD, row.names = NULL)

save(CpGs_allSD_M, 
     file = "  CpGs_allSD_M.rda")

load("  CpGs_allSD_M.rda")


#Filter Top 50% Most Variable CpG Sites
threshold_M <- quantile(beta_mostvar_M$SD, 0.5)  # 50% quantile
top_50_percent_M <- beta_mostvar_M[beta_mostvar_M$SD >= threshold_M, ]

hist(top_50_percent_M$SD)
describe(top_50_percent_M$SD)

# Convert rownames (CpG sites) into a column and create a new dataframe with SD values, just for archive, and save them as csv file
CpGs_top50SD_M <- data.frame(CpG_ID = rownames(top_50_percent_M), SD = top_50_percent_M$SD, row.names = NULL)

save(CpGs_top50SD_M, 
     file = "  CpGs_top50SD_M.rda")

load("  CpGs_top50SD_M.rda")

hist(CpGs_top50SD_M$SD)
describe(CpGs_top50SD_M$SD)


#### List the most 50% most variable probes  ####

load("  CpGs_top50SD_C.rda") #for Kids
ids_top50SD_C   <- CpGs_top50SD_C$CpG_ID

load("  CpGs_top50SD_M.rda") #for Mums
ids_top50SD_M   <- CpGs_top50SD_M$CpG_ID

overlap_C_M     <- intersect(ids_top50SD_C, ids_top50SD_M)
length(overlap_C_M) #410550 , 97% (410550/431674)


save(ids_top50SD_C, 
     file = "  ids_top50SD_C.rda")

save(ids_top50SD_M, 
     file = "  ids_top50SD_M.rda")



#### Listing CpGs that intraclass correlation coefficient (ICC) >.50 in our data ####

#This we do across mothers and children, as ICC is very sensitive for sample size. So we run this across both (like we pre-registered)
load("  beta.dropCpG.noob.dropSamples.rda")

beta_df <- as.data.frame(beta.dropCpG.noob.dropSamples)
dim(beta_df) #N=1528, of which N=16 technical replicates, with 863347 probes 

##read in phenotype data
#here we extract which ones are technical replicates
load("data_DNAmScores_BFY_InclTechRep_14Jan.rda")

df <- data_DNAmScores_BFY_InclTechRep_14Jan[,c("BaseID",
                                               "Sample_Name")]
head(df)

#### create dataset with just technical replicates 

# Extract unique base Sample_Names (without _repl)
df <- df %>%
  mutate(Base_Sample = gsub("_repl$", "", Sample_Name))  # Remove "_repl" to get base sample name

#Identify samples that have both original and replicate
samples_with_replicates <- df %>%
  group_by(Base_Sample) %>%
  filter(n() > 1) %>%   # Keep only those that appear more than once (original + replicate)
  ungroup()

df_filtered <- df %>%
  filter(Base_Sample %in% samples_with_replicates$Base_Sample)

TechRepl_list <- df_filtered

save(TechRepl_list , 
     file = "  TechRepl_list .rda")

#To calculate ICC, we need a beta file which just includes those samples who have a technical replicate
#And we need a file that states which sample is a technical replicate of whom


#Here we create the file stating who is a technical replicate of whom
df_filtered <- df_filtered %>%
  group_by(Base_Sample) %>%
  mutate(idx = cur_group_id()) %>%
  ungroup()

sampleid <- df_filtered[,c("BaseID",
                           "idx")]

names(sampleid)[names(sampleid) == "BaseID"] <- "id"

#Here we create the beta file only including those with a technical replicate, should have 32 obs (16 individuals with 16 technical replicates)
df_filtered$BaseID <- as.character(df_filtered$BaseID)
beta.dropCpG.noob.dropSamples_ICC <- beta.dropCpG.noob.dropSamples[, colnames(beta.dropCpG.noob.dropSamples) %in% df_filtered$BaseID]

colnames(beta.dropCpG.noob.dropSamples_ICC)
head(df_filtered)


### run ICC across CpGs ####

#we use this package: https://rdrr.io/github/xuz1/ENmix/man/repicc.html

#the Beta file needs to be a matrix
#set it as matrix 
beta.dropCpG.noob.dropSamples_ICC_matrix <- as.matrix(beta.dropCpG.noob.dropSamples_ICC)

#here we run the ICC across CpGs
icc_results <- repicc(beta.dropCpG.noob.dropSamples_ICC_matrix, sampleid)
head(icc_results)

icc_results_acrossAllCpGs <- icc_results


#Save ICC for CpGs 
save(icc_results_acrossAllCpGs , 
     file = "  icc_results_acrossAllCpGs.rda")

load("  icc_results_acrossAllCpGs.rda")

#select only CpGs that have icc of .50 or higher
icc_results_acrossAllCpGs_highICC <- icc_results_acrossAllCpGs %>%   #results in 339534 CpGs, which is 63% (of 431674)
  filter(icc >= 0.50)

#Save CpGs based on .50 or higher
save(icc_results_acrossAllCpGs_highICC , 
     file = "  icc_results_acrossAllCpGs_highICC.rda")

load("  icc_results_acrossAllCpGs_highICC.rda")

ids_highICC   <- icc_results_acrossAllCpGs_highICC$probe

save(ids_highICC , 
     file = "  ids_highICC.rda")



#### Load CpGs to be included in EWAS dataset ####

load("  ids_top50SD_C.rda") #top 50% variable probes in children

load("  ids_top50SD_M.rda") #top 50% variable probes in mothers

load("  ids_highICC.rda")  # CpGs with ICC >.50

selected_probes_C <- ids_top50SD_C[ids_top50SD_C %in% ids_highICC]

selected_probes_M <- ids_top50SD_M[ids_top50SD_M %in% ids_highICC]

save(selected_probes_C , 
     file = "  selected_probes_C.rda")

save(selected_probes_M , 
     file = "  selected_probes_M.rda")


#### Create Beta files with selected probes for mothers and children ####

#Create for children
load("  selected_probes_C.rda") #selected probes based on 50% variability and high ICC
load("  beta.dropCpG.noob.dropSamples.notechrep.children.rda") #betafile

beta.dropCpG.noob.dropSamples.notechrep.children.selectProbes <- beta.dropCpG.noob.dropSamples.notechrep.children[
  rownames(beta.dropCpG.noob.dropSamples.notechrep.children) %in% selected_probes_C,
]

save(beta.dropCpG.noob.dropSamples.notechrep.children.selectProbes, 
     file = "  beta.dropCpG.noob.dropSamples.notechrep.children.selectProbes.rda")

#create for mothers
load("  selected_probes_M.rda")
load("  beta.dropCpG.noob.dropSamples.notechrep.mothers.rda")

beta.dropCpG.noob.dropSamples.notechrep.mothers.SelectProbes <- beta.dropCpG.noob.dropSamples.notechrep.mothers[
  rownames(beta.dropCpG.noob.dropSamples.notechrep.mothers) %in% selected_probes_M,
]

save(beta.dropCpG.noob.dropSamples.notechrep.mothers.SelectProbes,
     file = "  beta.dropCpG.noob.dropSamples.notechrep.mothers.SelectProbes.rda")


### Residualize CpGs for celltype, plate and slide, separately for kids and mothers ####

### ON SELECTED PROBES

#### For Kids ####

#Phenodata (this is phenodata without tech replicates)
load("  pheno_excldpl.rda") #this is the overall phenotype data excluding technical replicates

#Beta data (these are betas selected on high variability and with ICC >.50)
load("  beta.dropCpG.noob.dropSamples.notechrep.children.selectProbes.rda") 

betas_forEwas_excldpl_kids <- beta.dropCpG.noob.dropSamples.notechrep.children.selectProbes

#selecting kids only
pheno_excldpl_kids <- pheno_excldpl %>% filter(Mum0Child1 == 1) #Make sure its just kids

#check if order is correct
identical(colnames(betas_forEwas_excldpl_kids ), pheno_excldpl_kids$BaseID) #should be TRUE

#If FALSE, run the code below
#betas_forEwas_excldpl_kids <- betas_forEwas_excldpl_kids[, pheno_excldpl_kids$BaseID]


#save pheno and beta data for kids
save(betas_forEwas_excldpl_kids, 
     file = "  betas_forEwas_excldpl_kids.rda")

save(pheno_excldpl_kids, 
     file = "  pheno_excldpl_kids.rda")

load("  betas_forEwas_excldpl_kids.rda")
load("  pheno_excldpl_kids.rda")

all(pheno_excldpl_kids$BaseID == colnames(betas_forEwas_excldpl_kids)) # should be TRUE

#Correct CpG for celltype

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_kids)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

celltype_data <- pheno_excldpl_kids[,c("BaseID",
                                       "celltype_1",
                                       "celltype_2",
                                       "celltype_3",
                                       "celltype_4",
                                       "celltype_5",
                                       "plate",
                                       "Slide")]
celltype_data <- as.data.frame(celltype_data)
rownames(celltype_data) <- celltype_data$BaseID  
celltype_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_kids[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(celltype_data), ]
  
  # Function to regress each CpG on cell types and extract residuals
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ celltype_1 + celltype_2 + 
                  celltype_3 + celltype_4 + 
                  celltype_5 + as.factor(plate) + as.factor(Slide),
                data = celltype_data)
    return(residuals(model))
  }
  
  # Apply function to batch
  residuals_batch <- apply(cpg_data_t_batch, 2, get_residuals)
  
  # Convert back so rows = CpGs, columns = participants
  residuals_list[[ceiling(i / batch_size)]] <- t(residuals_batch)
  
  # Print progress
  cat("Processed CpGs:", i, "to", batch_end, "\n")
}

# Combine all batches into a single dataframe
residuals_df <- do.call(rbind, residuals_list)

betas_forEwas_excldpl_ResCellTech_Kids <- residuals_df

# Show first few values
print(betas_forEwas_excldpl_ResCellTech_Kids[1:4, 1:4])  

save(betas_forEwas_excldpl_ResCellTech_Kids , 
     file = "  betas_forEwas_excldpl_ResCellTech_Kids.rda")

load("  betas_forEwas_excldpl_ResCellTech_Kids.rda")


#### For Mothers ####

#Phenodata (this is phenodata without tech replicates)
load("  pheno_excldpl.rda") #this is the overall phenotype data

#Beta data (these are betas selected on variability and with ICC >.50 and without technical replicates)
load("  beta.dropCpG.noob.dropSamples.notechrep.mothers.SelectProbes.rda") 

betas_forEwas_excldpl_mums <- beta.dropCpG.noob.dropSamples.notechrep.mothers.SelectProbes

#selecting mums only
pheno_excldpl_mums <- pheno_excldpl %>% filter(Mum0Child1 == 0) #Select mothers only 

#check if order is correct
identical(colnames(betas_forEwas_excldpl_mums),pheno_excldpl_mums$BaseID) #should be TRUE

#If FALSE, run the code below
#betas_forEwas_excldpl_mums <- betas_forEwas_excldpl_mums[, pheno_excldpl_mums$BaseID]

#save pheno and beta data for mums
save(betas_forEwas_excldpl_mums, 
     file = "  betas_forEwas_excldpl_mums.rda")

save(pheno_excldpl_mums, 
     file = "  pheno_excldpl_mums.rda")

load("  betas_forEwas_excldpl_mums.rda")
load("  pheno_excldpl_mums.rda")

all(pheno_excldpl_mums$BaseID == colnames(betas_forEwas_excldpl_mums)) # should be TRUE

#Correct CpG for celltype

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_mums)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

celltype_data <- pheno_excldpl_mums[,c("BaseID",
                                       "celltype_1",
                                       "celltype_2",
                                       "celltype_3",
                                       "celltype_4",
                                       "celltype_5",
                                       "plate",
                                       "Slide")]
celltype_data <- as.data.frame(celltype_data)
rownames(celltype_data) <- celltype_data$BaseID  
celltype_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_mums[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(celltype_data), ]
  
  # Function to regress each CpG on cell types and extract residuals
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ celltype_1 + celltype_2 + 
                  celltype_3 + celltype_4 + 
                  celltype_5 + as.factor(plate) + as.factor(Slide),
                data = celltype_data)
    return(residuals(model))
  }
  
  # Apply function to batch
  residuals_batch <- apply(cpg_data_t_batch, 2, get_residuals)
  
  # Convert back so rows = CpGs, columns = participants
  residuals_list[[ceiling(i / batch_size)]] <- t(residuals_batch)
  
  # Print progress
  cat("Processed CpGs:", i, "to", batch_end, "\n")
}

# Combine all batches into a single dataframe
residuals_df <- do.call(rbind, residuals_list)

betas_forEwas_excldpl_ResCellTech_mums <- residuals_df

# Show first few values
print(betas_forEwas_excldpl_ResCellTech_mums[1:4, 1:4])  # After residualizing larger sample, should be same results

save(betas_forEwas_excldpl_ResCellTech_mums , 
     file = "  betas_forEwas_excldpl_ResCellTech_mums.rda")

load("  betas_forEwas_excldpl_ResCellTech_mums.rda")

#### ON ALL PROBES ####


#### For Kids ####

load("  pheno_excldpl_kids.rda") #phenodata without tech replicates
load("  beta.dropCpG.noob.dropSamples.notechrep.children.rda") #beta data for all CpGs , so 863 347 probes

betas_forEwas_excldpl_kids_allCpGs <- beta.dropCpG.noob.dropSamples.notechrep.children

all(pheno_excldpl_kids$BaseID == colnames(betas_forEwas_excldpl_kids_allCpGs)) # should be TRUE

#If FALSE, run the code below
#betas_forEwas_excldpl_kids_allCpGs <- betas_forEwas_excldpl_kids_allCpGs[, pheno_excldpl_kids$BaseID]

#Correct CpG for celltype

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_kids_allCpGs)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

celltype_data <- pheno_excldpl_kids[,c("BaseID",
                                       "celltype_1",
                                       "celltype_2",
                                       "celltype_3",
                                       "celltype_4",
                                       "celltype_5",
                                       "plate",
                                       "Slide")]
celltype_data <- as.data.frame(celltype_data)
rownames(celltype_data) <- celltype_data$BaseID  
celltype_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_kids_allCpGs[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(celltype_data), ]
  
  # Function to regress each CpG on cell types and extract residuals
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ celltype_1 + celltype_2 + 
                  celltype_3 + celltype_4 + 
                  celltype_5 + as.factor(plate) + as.factor(Slide),
                data = celltype_data)
    return(residuals(model))
  }
  
  # Apply function to batch
  residuals_batch <- apply(cpg_data_t_batch, 2, get_residuals)
  
  # Convert back so rows = CpGs, columns = participants
  residuals_list[[ceiling(i / batch_size)]] <- t(residuals_batch)
  
  # Print progress
  cat("Processed CpGs:", i, "to", batch_end, "\n")
}

# Combine all batches into a single dataframe
residuals_df <- do.call(rbind, residuals_list)

betas_forEwas_excldpl_ResCellTech_Kids_AllCpGs <- residuals_df

# Show first few values
print(betas_forEwas_excldpl_ResCellTech_Kids_AllCpGs[1:4, 1:4])  

save(betas_forEwas_excldpl_ResCellTech_Kids_AllCpGs , 
     file = "  betas_forEwas_excldpl_ResCellTech_Kids_AllCpGs.rda")

load("  betas_forEwas_excldpl_ResCellTech_Kids_AllCpGs.rda")


#### For Mothers ####

load("  beta.dropCpG.noob.dropSamples.notechrep.mothers.rda")
load("  pheno_excldpl_mums.rda")

betas_forEwas_excldpl_mums_allCpGs <- beta.dropCpG.noob.dropSamples.notechrep.mothers

all(pheno_excldpl_mums$BaseID == colnames(betas_forEwas_excldpl_mums_allCpGs)) # should be TRUE

#If FALSE, run the code below
#betas_forEwas_excldpl_mums_allCpGs <- betas_forEwas_excldpl_mums_allCpGs[, pheno_excldpl_kids$BaseID]

#Correct CpG for celltype

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_mums_allCpGs)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

celltype_data <- pheno_excldpl_mums[,c("BaseID",
                                       "celltype_1",
                                       "celltype_2",
                                       "celltype_3",
                                       "celltype_4",
                                       "celltype_5",
                                       "plate",
                                       "Slide")]
celltype_data <- as.data.frame(celltype_data)
rownames(celltype_data) <- celltype_data$BaseID  
celltype_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_mums_allCpGs[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(celltype_data), ]
  
  # Function to regress each CpG on cell types and extract residuals
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ celltype_1 + celltype_2 + 
                  celltype_3 + celltype_4 + 
                  celltype_5 + as.factor(plate) + as.factor(Slide),
                data = celltype_data)
    return(residuals(model))
  }
  
  # Apply function to batch
  residuals_batch <- apply(cpg_data_t_batch, 2, get_residuals)
  
  # Convert back so rows = CpGs, columns = participants
  residuals_list[[ceiling(i / batch_size)]] <- t(residuals_batch)
  
  # Print progress
  cat("Processed CpGs:", i, "to", batch_end, "\n")
}

# Combine all batches into a single dataframe
residuals_df <- do.call(rbind, residuals_list)

betas_forEwas_excldpl_ResCellTech_mums_AllCpGs <- residuals_df


# Show first few values
print(betas_forEwas_excldpl_ResCellTech_mums_AllCpGs[1:4, 1:4])  # After residualizing larger sample, should be same results

save(betas_forEwas_excldpl_ResCellTech_mums_AllCpGs , 
     file = "  betas_forEwas_excldpl_ResCellTech_mums_AllCpGs.rda")

load("  betas_forEwas_excldpl_ResCellTech_mums_AllCpGs.rda")




#### P-value thresholds for EWAS ####

#In our pre-registration, we note three thresholds for p-values
# 1) Bonferroni correction α= 0.05/N CpG sites
# 2) α= 0.05/N uncorrelated CpG sites , applying Co-MethylaSonwith genomic CpG Background technique (CoMBack) to determine uncorrelated DNA methylation sites
# 3) α= 0.05/N uncorrelated CpG sites , applying Matrix Spectral DecomposiSon (MSD) to determine uncorrelated DNA methylation sites

### Calculate different P-value Thresholds for EWAS ####

install.packages("devtools")
devtools::install_github("Treys925/CoMeBackV2_R_package", force = TRUE)

library(CoMeBackV2) #Package and code below developed by Trey Smith, UMichigan
library(devtools)
library(dplyr)
library(stringr)



### Load in different Betas
#Kids/Mothers seperatly, only selected CpGs (on variability and ICC) and Residualized for cellcomp/tech and unresidualized

#load("betas_forEwas_excldpl_ResCellTech_Kids.rda")

load("betas_forEwas_excldpl_ResCellTech_mums.rda")

#Full beta files (pick which one is of interest)

#beta <- betas_forEwas_excldpl_ResCellTech_Kids

beta <- betas_forEwas_excldpl_ResCellTech_mums

#### Classic Bonferroni ####
dim(beta)

0.05 / dim(beta)[1] 


#### ComBat based on script Trey ####

#Note: script needs time to run, so you can also select 5 samples to check if script works
#beta <- betas_forEwas_excldpl[1:270000, 1:5]

#Make sure you have probes as columns (so transpose your data)
tbeta <- t(beta)

#EpicV2 data has different suffix then EpicV1 so make sure that numbers after _ are omitted and duplicate probes are removed
#Keep only one of each duplicate (based on part before underscore)
col_prefix <- str_extract(colnames(tbeta), "^[^_]+")
unique_cols <- !duplicated(col_prefix)
df_unique <- tbeta[, unique_cols]

#Remove underscore and everything after
colnames(df_unique) <- gsub("_.*", "", colnames(df_unique))

#Turn it into data frame
tbeta_nodupl <- as.data.frame(df_unique)

#Turn it into a matrix 
beta_noCTC1 <- as.matrix(tbeta_nodupl)

#### Determine uncorrelated CpGs based on Combat ####

## Construct CMRs using the defaults
## CMRs<- cmr(Mdata = beta_noCTC1, Iarray = "450K") #also possible to choose other arrays, eg 450K

## construct CMRs using fixed Spearman correlation cut-off of 0.3 and maximum probe distance of 2Kbp
#Note: it will take some time for this script to run
CMRs_Scor50_pw2K <- cmr(beta_noCTC1, Iarray = "EPICv2", cormethod = "spearman", corlo = 0.3, maxprbdst =2000)

#here you save the output
#Combat_children <- CMRs_Scor50_pw2K #this is Combat based on betas residualized for cell comb and tech cov

#save(Combat_children , 
#     file = "Combat_children.rda")

Combat_moms <- CMRs_Scor50_pw2K #this is Combat based on betas residualized for cell comb and tech cov




#Calculate the P-value
total <- 0
for(i in 1:length(Combat_children)){
  total = total + length(Combat_children[[i]])
}

print (total)
print(0.05/total)




#### Determine uncorrelated CpGs based on Matrix Spectral Decomposition	####

## Eigendecomposition
# Matrix Spectral Decomposition function
DetermineMEffLi <- function(Data,Columns){
  MeanCenteredAutoScaledData <- Data
  for(i in 1:length(Columns)){
    Entry <- Columns[i]
    Mean <- mean(Data[,Entry])
    Std <- sd(Data[,Entry])
    MeanCenteredAutoScaledData[,Entry] <- (MeanCenteredAutoScaledData[,Entry] - Mean)/Std
  }
  SVD <- svd(x=as.matrix(MeanCenteredAutoScaledData[,Columns]))
  EigenValues <- SVD$d^2 / dim(SVD$u)[1]
  M <- length(EigenValues)
  L <- M - 1
  Var <- var(EigenValues)
  MEff <- M*(1.0-(L*Var/M^2))
  IntEVals <- as.numeric(EigenValues>=1.0)
  NonIntEVals <- EigenValues - floor(EigenValues)
  MEffLi <- ceiling(sum(IntEVals) + sum(NonIntEVals))
  MEffLi
}

## Make sure CpGs are Columns!
#Note: this code takes some time to run
CpGs <- colnames(tbeta_nodupl)
MEffLiC <- DetermineMEffLi(Data=tbeta_nodupl, Columns=CpGs) 
print (MEffLiC) 

#save your output
MSD_kids <- MEffLiC
save(MSD_kids, 
     file = "MSD_kids.rda")

load("MSD_mums.rda")

#bonferonni correction 
print(0.05/MEffLiC) 

############################  End ############################################################################
