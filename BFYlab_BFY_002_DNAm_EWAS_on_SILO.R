################################################### 
########    EWAS Mothers and Children  ############ 
################################################### 

### Install Packages ####
BiocManager::install("bacon")
BiocManager::install("limma")
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")

BiocManager::install("minfi", force = TRUE)
BiocManager::install("IRanges", force = TRUE)

library(bacon)
library(limma)
library(stringr)
library(minfi)
library(qqman)

devtools::install_local("    EWAS/Rpackages/IlluminaHumanMethylationEPICv2anno.20a1.hg38-master/IlluminaHumanMethylationEPICv2anno.20a1.hg38-master", force = TRUE)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

BiocManager::install("bumphunter", force = TRUE)
library(bumphunter)

devtools::install_local("    EWAS/bumphunter-devel/bumphunter-devel")
system.file(package = "bumphunter")


### Set Working Directory
setwd("  EWAS/EWAS April 2025")

### Load in the Data from Max Planck Institute ####

##### For Mothers ##### 

#Phenotype data
load("  pheno_excldpl_mums.rda")
dim(pheno_excldpl_mums)

#CpGs with Betas for Mothers, these have been residualized for celltype, plate and slide
load("betas_forEwas_excldpl_ResCellTech_mums.rda")
dim(betas_forEwas_excldpl_ResCellTech_mums)

#### Load in baseline covariates from SILO ####
baseline_cov_all <- read_dta("  clean_data/age4_pheno_epi_merge.dta")

#select baseline covariates
baseline_cov_select <- baseline_cov_all[,c("sampleid",
                                           "treat",
                                           "magea0",
                                           "mnbiochilda0",
                                           "mgoodhealtha0",
                                           "mdepressiona0",
                                           "mhsgrada0",
                                           "msomecollegea0",
                                           "massociatesa0",
                                           "mbachelorsa0",
                                           "medunknowna0",  
                                           "mblacka0",
                                           "mracemultiplea0",
                                           "mraceothera0",
                                           "mhispanica0",
                                           "mraceunknowna0",
                                           "mcohaba0",
                                           "mmarrieda0",
                                           "mdivorceda0",
                                           "mrelateothera0", 
                                           "mrelateunknowna0",
                                           "mcigduringavgwka0",
                                           "malcduringavgwka0",
                                           "hhincomecat2a0",
                                           "hhincomecat3a0",
                                           "hhincomecat4a0",
                                           "hhincomecat5a0",
                                           "hhincomecat6a0",
                                           "hhnetworthcat2a0",
                                           "hhnetworthcat3a0",
                                           "hhnetworthcat4a0",
                                           "hhnetworthcat5a0",
                                           "hhnetworthcat6a0",
                                           "hhnadulta0",
                                           "hhbiodada0",
                                           "cfemalea0",
                                           "cweightlba0",
                                           "cgestagewksa0",
                                           "assessageinmonthsa4",
                                           "assessinterviewera4")]


#In line with earlier BFY analyses, we impute the mean for missing values on 4 covariate that have missing values for some families, 
#We also generate a dummy variable that reflects if the value was imputed in the macro covariates list


library(dplyr)

baseline_cov_select <- baseline_cov_select %>%
  mutate(across(c(mcigduringavgwka0, malcduringavgwka0, cweightlba0, cgestagewksa0),
                list(miss_ = ~ ifelse(is.na(.), mean(., na.rm = TRUE), .),
                     macro_cov_ = ~ as.integer(is.na(.))),
                .names = "{fn}{.col}"))

macro_covariates <- names(baseline_cov_select)[grepl("^macro_cov_", names(baseline_cov_select))]



#### Create Phenotype datasets forMothers ####

#We have two phenotype datasets
# pheno_excldpl_mums = One from Max Planck Institute with participant id (sampleid) and with the connector to the CpGs (BaseID, so on which part of the plate they were)
# baseline_cov_select = One from SILO including all baseline covariates and whether participants were part of treatment of control group
# We have to merge these two on sampleid

head(pheno_excldpl_mums)
head(baseline_cov_select)

#We add a P to sampleid so the sample ids are the same across datasets
pheno_excldpl_mums_select <- pheno_excldpl_mums %>%
  mutate(sampleid = paste0("P", sampleId)) %>%
  select(1:2, sampleid, everything(),-sampleId) #remove old column

class(pheno_excldpl_mums_select$sampleid)
class(baseline_cov_select$sampleid)

#select the variables to merge
pheno_excldpl_mums_select <- pheno_excldpl_mums_select[,c("BaseID",
                                                          "sampleid",
                                                          "Treatment.Group..blinded.")]

#merge the two phenofiles, with the pheno_excldpl_mums_select as the baseline
pheno_excldpl_mums_merged <- pheno_excldpl_mums_select %>%
  left_join(baseline_cov_select, by = "sampleid")

colnames(pheno_excldpl_mums_merged)
head(pheno_excldpl_mums_merged)

pheno_excldpl_mums_mergedBaseCov <- pheno_excldpl_mums_merged
describe(pheno_excldpl_mums_mergedBaseCov)

save(pheno_excldpl_mums_mergedBaseCov , 
     file = "  pheno_excldpl_mums_mergedBaseCov.rda")

load("  pheno_excldpl_mums_mergedBaseCov.rda")


#### Run EWAS ####

## for mums

#Rename Phenotype and CpG data
FinalPD <- pheno_excldpl_mums_merged
betas <-betas_forEwas_excldpl_ResCellTech_mums[, match(FinalPD$BaseID, colnames(betas_forEwas_excldpl_ResCellTech_mums)), drop = FALSE]

#we set beta file in the same order as the pheno file, so this order should be the same (otherwhise it gets scrambled)
identical(colnames(betas), FinalPD$BaseID) #should be TRUE

### Set Covars
#here we select all the baseline covariates
#those not #Continious will be set as.factor later in the script
covars <- c( "magea0",                       #Continious
             "mnbiochilda0",                 #Continious
             "mgoodhealtha0",  
             "mdepressiona0",                #Continious
             "mhsgrada0",      
             "msomecollegea0", 
             "massociatesa0",  
             "mbachelorsa0",  
             "medunknowna0",   
             "mblacka0",       
             "mracemultiplea0", 
             "mraceothera0",    
             "mhispanica0",     
             "mraceunknowna0", 
             "mcohaba0",        
             "mmarrieda0",      
             "mdivorceda0",     
             "mrelateothera0", 
             "mrelateunknowna0",
             "hhincomecat2a0",
             "hhincomecat3a0",
             "hhincomecat4a0",
             "hhincomecat5a0",
             "hhincomecat6a0",
             "hhnetworthcat2a0",
             "hhnetworthcat3a0",
             "hhnetworthcat4a0",
             "hhnetworthcat5a0",
             "hhnetworthcat6a0",
             "hhnadulta0",                  #Continious
             "hhbiodada0",
             "cfemalea0",
             "assessageinmonthsa4",         #Continious
             "miss_mcigduringavgwka0",      #Continious
             "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
             "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
             "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
             "miss_cweightlba0",            #Continious,this is the mean imputed version
             "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
             "miss_cgestagewksa0",          #Continious, this is the mean imputed version
             "macro_cov_cgestagewksa0",     #This is the dummy if mean was imputed
             "assessinterviewera4")

#Set non-continuous variables to factor
#set them all as.factor in R
FinalPD[, c(  "mgoodhealtha0",  
              "mhsgrada0",      
              "msomecollegea0", 
              "massociatesa0",  
              "mbachelorsa0",  
              "medunknowna0",   
              "mblacka0",       
              "mracemultiplea0", 
              "mraceothera0",    
              "mhispanica0",     
              "mraceunknowna0", 
              "mcohaba0",        
              "mmarrieda0",      
              "mdivorceda0",     
              "mrelateothera0", 
              "mrelateunknowna0",
              "hhincomecat2a0",
              "hhincomecat3a0",
              "hhincomecat4a0",
              "hhincomecat5a0",
              "hhincomecat6a0",
              "hhnetworthcat2a0",
              "hhnetworthcat3a0",
              "hhnetworthcat4a0",
              "hhnetworthcat5a0",
              "hhnetworthcat6a0",
              "hhbiodada0",
              "cfemalea0",
              "macro_cov_mcigduringavgwka0",
              "macro_cov_malcduringavgwka0",
              "macro_cov_cweightlba0",
              "macro_cov_cgestagewksa0",
              "assessinterviewera4")] <-
  lapply (FinalPD[, c("mgoodhealtha0",  
                      "mhsgrada0",      
                      "msomecollegea0", 
                      "massociatesa0",  
                      "mbachelorsa0",  
                      "medunknowna0",   
                      "mblacka0",       
                      "mracemultiplea0", 
                      "mraceothera0",    
                      "mhispanica0",     
                      "mraceunknowna0", 
                      "mcohaba0",        
                      "mmarrieda0",      
                      "mdivorceda0",     
                      "mrelateothera0", 
                      "mrelateunknowna0",
                      "hhincomecat2a0",
                      "hhincomecat3a0",
                      "hhincomecat4a0",
                      "hhincomecat5a0",
                      "hhincomecat6a0",
                      "hhnetworthcat2a0",
                      "hhnetworthcat3a0",
                      "hhnetworthcat4a0",
                      "hhnetworthcat5a0",
                      "hhnetworthcat6a0",
                      "hhbiodada0",
                      "cfemalea0",
                      "macro_cov_mcigduringavgwka0",
                      "macro_cov_malcduringavgwka0",
                      "macro_cov_cweightlba0",
                      "macro_cov_cgestagewksa0",
                      "assessinterviewera4")], as.factor)

#check if all went well
sapply(FinalPD, class)


#### Run EWAS for mums ####
### Set Outcome 
outcome <- "treat" #this is our outcome (treatment/control), which we set as factor in line below
FinalPD$treat <- as.factor(FinalPD$treat)

### Run EWAS 
ptmp <- FinalPD #we have no missing variables, as we imputed the means. If we had, we would use here FinalPD[complete.cases(FinalPD$x) & complete.cases(FinalPD$xxx) & complete.cases(FinalPD$xxx),] 
myregressionmodel <- paste0("~", outcome, "+", paste(covars, collapse = "+") ) %>% as.formula()
mod <- model.matrix(myregressionmodel, ptmp)
btmp <- betas[, match((ptmp$BaseID), colnames(betas))]

# Run the single site association model
out <- lmFit(btmp, mod)
out <- eBayes(out)
ss.hits <- limma::topTable(out, coef = 2, number = nrow(out))

ss.hits.mums.cov <- ss.hits

save(ss.hits.mums.cov , 
     file = "ss.hits.mums.cov.rda")

load("ss.hits.mums.cov.rda")

ss.hits.mums.cov.excel <- ss.hits.mums.cov %>%
  tibble::rownames_to_column(var = "CpGs")


write_xlsx(ss.hits.mums.cov.excel, "ss.hits.mums.cov.xlsx")

# Check if any sites are significant at first significance threshold

#See 00_DNAm_QC_general and QC EWAS for calculation p-value thresholds

#ComeBack
any(ss.hits.mums.cov$P.Value < 1.761494e-06)
num_significant <- sum(ss.hits.mums.cov$P.Value < 1.761494e-06, na.rm = TRUE)
num_significant

#MSD
any(ss.hits.mums.cov$P.Value < 5.820722e-05)
num_significant <- sum(ss.hits.mums.cov$P.Value < 5.820722e-05, na.rm = TRUE)
num_significant

#### Do Bacon Adjustment mums ####

#bacon adjustment
bc <- bacon(teststatistics = ss.hits.mums.cov$t)
PVals <- pval(bc)
ss.hits.mums.cov.bacon <- data.frame(PVals)
row.names(ss.hits.mums.cov.bacon) <- row.names(ss.hits.mums.cov)
colnames(ss.hits.mums.cov.bacon) <- "P.Value"

save(ss.hits.mums.cov.bacon , 
     file = "ss.hits.mums.cov.bacon.rda")

load("ss.hits.mums.cov.bacon.rda")

ss.hits.mums.cov.bacon.excel <-ss.hits.mums.cov.bacon %>%
  tibble::rownames_to_column(var = "CpGs")


write_xlsx(ss.hits.mums.cov.bacon.excel, "ss.hits.mums.cov.bacon.xlsx")

# Check if any sites are significant at first significance threshold
#Bonferroni
any(ss.hits.mums.cov.bacon$P.Value < 1.77072e-07)
num_significant <- sum(ss.hits.mums.cov.bacon$P.Value <  1.77072e-07, na.rm = TRUE)
num_significant

#ComeBack
any(ss.hits.mums.cov.bacon$P.Value < 1.761494e-06)
num_significant <- sum(ss.hits.mums.cov.bacon$P.Value < 1.761494e-06, na.rm = TRUE)
num_significant

#MSD
any(ss.hits.mums.cov.bacon$P.Value < 5.820722e-05)
num_significant <- sum(ss.hits.mums.cov.bacon$P.Value < 5.820722e-05, na.rm = TRUE)
num_significant



#### Create EWAS plots ####

#QQ plot - before bacon
#here you see some inflation, and lambda above 1.2 so we try bacon adjustment too

observed <- -log10(sort(ss.hits.mums.cov$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.mums.cov$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_mums.pdf") #open PDF device to save the plot

qq(ss.hits.mums.cov$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#QQ plot - after bacon adjustment

observed <- -log10(sort(ss.hits.mums.cov.bacon$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.mums.cov.bacon$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_mums_bacon.pdf") #open PDF device to save the plot

qq(ss.hits.mums.cov.bacon$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#### Volcano Plot ####

pdf("volcano_zscore_mums.pdf") #open PDF device to save the plot

ggplot(ss.hits.mums.cov, aes(x = ((t-50)/10) , y = -log10(P.Value))) +
  geom_point() +
  theme_bw() +
  labs(x = "Z-Score", y = "-log10 P Value")

dev.off() #closes the PDF device, ensuring that the file is properly saved.

pdf("volcanoplot_mums.pdf") #open PDF device to save the plot

ggplot(ss.hits.mums.cov, aes(x = logFC , y = -log10(P.Value))) +
  geom_point() +
  theme_bw() +
  labs(x = "Log Fold Change", y = "-log10 P Value", title = "Volcano Plot", subtitle = "", caption = "" )

dev.off() #closes the PDF device, ensuring that the file is properly saved.

# Manhattanplot 


#Create annotation data where CpGs are linked to location on the genome
anno2 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations
anno2 <- data.frame(anno2)
anno2$chr <- str_split(anno2$chr, "chr",n=2,simplify=T)[,2]
anno2 <- anno2[anno2$chr %in% 1:22,]
anno2$chr <- as.numeric(anno2$chr)

#anno2 includes 905565 CpGs across all EpicV2
# "ss.hits.mums.cov" is our EPIC v2 data hits for mums

forman_mums <- merge(anno2[c("chr", "pos")], ss.hits.mums.cov[c("P.Value")], by = "row.names")
forman_mums <- forman_mums[complete.cases(forman_mums),]
forman_mums$SNP <- ""


# Generate the Manhattanplot
pdf("manhattan_ewas_mums.pdf") #open PDF device to save the plot

manhattan_ewas_mums <- qqman::manhattan(
  forman_mums, 
  chr = "chr", 
  bp = "pos", 
  p = "P.Value", 
  suggestiveline = FALSE, 
  genomewideline = -log10(1.836541e-07)
)

dev.off() #closes the PDF device, ensuring that the file is properly saved.


### Residualize CpGs for baseline covariates ####

#This is for later elastic net analyses

#Phenodata (this is phenodata including baseline covariates)
load("  pheno_excldpl_mums_mergedBaseCov.rda") #this is the overall phenotype data
colnames(pheno_excldpl_mums_mergedBaseCov)
dim(pheno_excldpl_mums_mergedBaseCov)

#Beta data (these are betas selected on variability and with ICC >.50)
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_mums.rda")
dim(betas_forEwas_excldpl_ResCellTech_mums)
betas_forEwas_excldpl_ResCellTech_mums[1:5, 1:5]

#check if order is correct
identical(colnames(betas_forEwas_excldpl_ResCellTech_mums),pheno_excldpl_mums_mergedBaseCov$BaseID) #should be TRUE

all(pheno_excldpl_mums_mergedBaseCov$BaseID == colnames(betas_forEwas_excldpl_ResCellTech_mums)) # should be TRUE

#Correct CpGs for baseline covariates (they have already been residualized for celltype / plate / array)

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_ResCellTech_mums)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

baselinecov_data <- pheno_excldpl_mums_mergedBaseCov[,c("BaseID",
                                                        "magea0",                       #Continious
                                                        "mnbiochilda0",                 #Continious
                                                        "mgoodhealtha0",  
                                                        "mdepressiona0",                #Continious
                                                        "mhsgrada0",      
                                                        "msomecollegea0", 
                                                        "massociatesa0",  
                                                        "mbachelorsa0",  
                                                        "medunknowna0",   
                                                        "mblacka0",       
                                                        "mracemultiplea0", 
                                                        "mraceothera0",    
                                                        "mhispanica0",     
                                                        "mraceunknowna0", 
                                                        "mcohaba0",        
                                                        "mmarrieda0",      
                                                        "mdivorceda0",     
                                                        "mrelateothera0", 
                                                        "mrelateunknowna0",
                                                        "hhincomecat2a0",
                                                        "hhincomecat3a0",
                                                        "hhincomecat4a0",
                                                        "hhincomecat5a0",
                                                        "hhincomecat6a0",
                                                        "hhnetworthcat2a0",
                                                        "hhnetworthcat3a0",
                                                        "hhnetworthcat4a0",
                                                        "hhnetworthcat5a0",
                                                        "hhnetworthcat6a0",
                                                        "hhnadulta0",                  #Continious
                                                        "hhbiodada0",
                                                        "cfemalea0",
                                                        "assessageinmonthsa4",         #Continious
                                                        "assessinterviewera4", 
                                                        "miss_mcigduringavgwka0",      #Continious
                                                        "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
                                                        "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
                                                        "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
                                                        "miss_cweightlba0",            #Continious,this is the mean imputed version
                                                        "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
                                                        "miss_cgestagewksa0",          #Continious, this is the mean imputed version
                                                        "macro_cov_cgestagewksa0")]     #This is the dummy if mean was imputed




baselinecov_data <- as.data.frame(baselinecov_data)
rownames(baselinecov_data) <- baselinecov_data$BaseID  
baselinecov_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_ResCellTech_mums[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(baselinecov_data), ]
  
  # Function to regress each CpG on baseline covariates and extract residuals
  
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ magea0 +  mnbiochilda0 + as.factor(mgoodhealtha0) + mdepressiona0 +              
                  as.factor(mhsgrada0) +      
                  as.factor(msomecollegea0)+ 
                  as.factor(massociatesa0)+  
                  as.factor(mbachelorsa0)+  
                  as.factor(medunknowna0)+   
                  as.factor(mblacka0)+       
                  as.factor(mracemultiplea0)+ 
                  as.factor(mraceothera0)+    
                  as.factor(mhispanica0)+     
                  as.factor(mraceunknowna0)+ 
                  as.factor(mcohaba0)+        
                  as.factor(mmarrieda0)+      
                  as.factor(mdivorceda0)+     
                  as.factor(mrelateothera0)+ 
                  as.factor(mrelateunknowna0)+
                  as.factor(hhincomecat2a0)+
                  as.factor(hhincomecat3a0)+
                  as.factor(hhincomecat4a0)+
                  as.factor(hhincomecat5a0)+
                  as.factor(hhincomecat6a0)+
                  as.factor(hhnetworthcat2a0)+
                  as.factor(hhnetworthcat3a0)+
                  as.factor(hhnetworthcat4a0)+
                  as.factor(hhnetworthcat5a0)+
                  as.factor(hhnetworthcat6a0)+
                  hhnadulta0 +                                 #Continious
                  as.factor(hhbiodada0)+
                  as.factor(cfemalea0)+
                  assessageinmonthsa4+                         #Continious
                  as.factor(assessinterviewera4)+ 
                  miss_mcigduringavgwka0+                      #Continious
                  as.factor(macro_cov_mcigduringavgwka0) +     #This is the dummy if mean was imputed
                  miss_malcduringavgwka0+                      #Continious, this is the mean imputed version
                  as.factor(macro_cov_malcduringavgwka0) +     #This is the dummy if mean was imputed
                  miss_cweightlba0+                            #Continious,this is the mean imputed version
                  macro_cov_cweightlba0+                       #This is the dummy if mean was imputed
                  miss_cgestagewksa0+                          #Continious, this is the mean imputed version
                  as.factor(macro_cov_cgestagewksa0),          #This is the dummy if mean
                data = baselinecov_data)
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

betas_forEwas_excldpl_ResCellTech_ResBaseline_mums <- residuals_df

# Show first few values
print(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[1:4, 1:4]) 
dim(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums)

save(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums ,
     file ="  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_mums.rda") 


load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_mums.rda")


#### Check for outliers in top hits ####

### Load in the data

#load phenotype date with treatment variable
load("  pheno_excldpl_mums_mergedBaseCov.rda")

#load beta data residualized for baseline covariates
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_mums.rda")

#load tophits after bacon adjustment
load("  ss.hits.mums.cov.bacon.rda")
CpGhitsBacon <- ss.hits.mums.cov.bacon %>%
  arrange(P.Value)

### Create dataset with only tophits

#select 17 hits that were significant in EWAS
top_cpgs <- rownames(CpGhitsBacon)[1:17]

#create betadata just for those 17 hits
Beta_tophits <- as.data.frame(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[top_cpgs,])
head(Beta_tophits)

#transpose the data
Beta_tophits_transposed <- as.data.frame(t(Beta_tophits))

#turn rownames into collumn names
Beta_tophits_transposed$BaseID <- rownames(Beta_tophits_transposed)
rownames(Beta_tophits_transposed) <- NULL

#merge treatment variable to it
hits_treat <- Beta_tophits_transposed %>%
  left_join(pheno_excldpl_mums_mergedBaseCov %>% select(BaseID, treat), by = "BaseID")

table(hits_treat$treat) #0 is control, 1 is treatment

#### Plot methylation values ####
library(ggplot2)
colnames(hits_treat)

hits_treat$treat <- as.factor(hits_treat$treat)


ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg13412754_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg13412754_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg19137044_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg19137044_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg25576801_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg25576801_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg04188396_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg04188396_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg01144399_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg01144399_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg03786343_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg03786343_BC21 by Treatment Group") +
  theme_minimal()

#7
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg15462959_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg15462959_BC21 by Treatment Group") +
  theme_minimal()

#8
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg13412754_TC22)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg13412754_TC22 by Treatment Group") +
  theme_minimal()

#9
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg09739347_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg09739347_TC21 by Treatment Group") +
  theme_minimal()

#10
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg16717411_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg16717411_BC21 by Treatment Group") +
  theme_minimal()

#11
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg01801101_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg01801101_BC21 by Treatment Group") +
  theme_minimal()

#12
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg14413795_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg14413795_TC21 by Treatment Group") +
  theme_minimal()

#13
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg12795046_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg12795046_TC21 by Treatment Group") +
  theme_minimal()

#14
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg14641625_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg14641625_BC21 by Treatment Group") +
  theme_minimal()

#15
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg16483795_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg16483795_BC21 by Treatment Group") +
  theme_minimal()

#16
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg17495671_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg17495671_BC21 by Treatment Group") +
  theme_minimal()

#17
ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg25314817_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg25314817_BC21 by Treatment Group") +
  theme_minimal()


### now do it based on sd ####
library(ggplot2)
library(dplyr)

cpg_list <- c("cg13412754_TC21", "cg19137044_TC21", "cg25576801_TC21", "cg04188396_BC21", "cg01144399_BC21", "cg03786343_BC21", "cg15462959_BC21", "cg13412754_TC22", "cg09739347_TC21", "cg16717411_BC21", "cg01801101_BC21", "cg14413795_TC21",
              "cg12795046_TC21", "cg14641625_BC21", "cg16483795_BC21", "cg17495671_BC21", "cg25314817_BC21")

for (cpg in cpg_list) {
  
  plot_data <- hits_treat %>%
    group_by(treat) %>%
    mutate(
      mean_val = mean(.data[[cpg]], na.rm = TRUE),
      sd_val = sd(.data[[cpg]], na.rm = TRUE),
  is_outlier_sd = abs(.data[[cpg]] - mean_val) > 3 * sd_val
  ) %>%
  ungroup()
  
  p <- ggplot(plot_data, aes(x = factor(treat, levels = c(0,1), labels = c("Control", "Treatment")),
                             y = .data[[cpg]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = filter(plot_data, is_outlier_sd),
               aes(y = .data[[cpg]]), color = "black", shape = 16, size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
    labs(x = "Group", y = "Methylation Level",
         title = paste("Methylation", cpg, "by Treatment Group")) +
    theme_minimal()

print(p)
    
}



#### Run QQplot but excluding outlier CpGS ####

#check which cpgs have extreme outliers >5 sd


extreme_flags <-apply(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums, 1, function(row) {
  mu <- mean(row, na.rm =TRUE)
  sigma <- sd(row, na.rm = TRUE)
  
  if(is.na(sigma) || sigma ==0) return(FALSE)
  any(abs(row-mu) >5 *sigma, na.rm = TRUE)
})

sum(extreme_flags)

cpgs_with_extreme_outliers <- rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums)[extreme_flags]

sum(Beta_tophits %in% cpgs_with_extreme_outliers) #none of them are in the tophits


#create beta file excluding cpgs with extreme outliers
betas_nooutliers <- betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[!extreme_flags, ]


##### Run EWAS with qqplot on data without those outliers ####

load("  pheno_excldpl_mums_mergedBaseCov.rda")


#Rename Phenotype and CpG data
FinalPD <- pheno_excldpl_mums_mergedBaseCov
betas <-betas_nooutliers[, match(FinalPD$BaseID, colnames(betas_nooutliers)), drop = FALSE]

#we set beta file in the same order as the pheno file, so this order should be the same (otherwhise it gets scrambled)
identical(colnames(betas), FinalPD$BaseID) #should be TRUE

### Set Covars
#here we select all the baseline covariates
#those not #Continious will be set as.factor later in the script
covars <- c( "magea0",                       #Continious
             "mnbiochilda0",                 #Continious
             "mgoodhealtha0",  
             "mdepressiona0",                #Continious
             "mhsgrada0",      
             "msomecollegea0", 
             "massociatesa0",  
             "mbachelorsa0",  
             "medunknowna0",   
             "mblacka0",       
             "mracemultiplea0", 
             "mraceothera0",    
             "mhispanica0",     
             "mraceunknowna0", 
             "mcohaba0",        
             "mmarrieda0",      
             "mdivorceda0",     
             "mrelateothera0", 
             "mrelateunknowna0",
             "hhincomecat2a0",
             "hhincomecat3a0",
             "hhincomecat4a0",
             "hhincomecat5a0",
             "hhincomecat6a0",
             "hhnetworthcat2a0",
             "hhnetworthcat3a0",
             "hhnetworthcat4a0",
             "hhnetworthcat5a0",
             "hhnetworthcat6a0",
             "hhnadulta0",                  #Continious
             "hhbiodada0",
             "cfemalea0",
             "assessageinmonthsa4",         #Continious
             "miss_mcigduringavgwka0",      #Continious
             "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
             "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
             "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
             "miss_cweightlba0",            #Continious,this is the mean imputed version
             "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
             "miss_cgestagewksa0",          #Continious, this is the mean imputed version
             "macro_cov_cgestagewksa0",     #This is the dummy if mean was imputed
             "assessinterviewera4")

#Set non-continuous variables to factor
#set them all as.factor in R
FinalPD[, c(  "mgoodhealtha0",  
              "mhsgrada0",      
              "msomecollegea0", 
              "massociatesa0",  
              "mbachelorsa0",  
              "medunknowna0",   
              "mblacka0",       
              "mracemultiplea0", 
              "mraceothera0",    
              "mhispanica0",     
              "mraceunknowna0", 
              "mcohaba0",        
              "mmarrieda0",      
              "mdivorceda0",     
              "mrelateothera0", 
              "mrelateunknowna0",
              "hhincomecat2a0",
              "hhincomecat3a0",
              "hhincomecat4a0",
              "hhincomecat5a0",
              "hhincomecat6a0",
              "hhnetworthcat2a0",
              "hhnetworthcat3a0",
              "hhnetworthcat4a0",
              "hhnetworthcat5a0",
              "hhnetworthcat6a0",
              "hhbiodada0",
              "cfemalea0",
              "macro_cov_mcigduringavgwka0",
              "macro_cov_malcduringavgwka0",
              "macro_cov_cweightlba0",
              "macro_cov_cgestagewksa0",
              "assessinterviewera4")] <-
  lapply (FinalPD[, c("mgoodhealtha0",  
                      "mhsgrada0",      
                      "msomecollegea0", 
                      "massociatesa0",  
                      "mbachelorsa0",  
                      "medunknowna0",   
                      "mblacka0",       
                      "mracemultiplea0", 
                      "mraceothera0",    
                      "mhispanica0",     
                      "mraceunknowna0", 
                      "mcohaba0",        
                      "mmarrieda0",      
                      "mdivorceda0",     
                      "mrelateothera0", 
                      "mrelateunknowna0",
                      "hhincomecat2a0",
                      "hhincomecat3a0",
                      "hhincomecat4a0",
                      "hhincomecat5a0",
                      "hhincomecat6a0",
                      "hhnetworthcat2a0",
                      "hhnetworthcat3a0",
                      "hhnetworthcat4a0",
                      "hhnetworthcat5a0",
                      "hhnetworthcat6a0",
                      "hhbiodada0",
                      "cfemalea0",
                      "macro_cov_mcigduringavgwka0",
                      "macro_cov_malcduringavgwka0",
                      "macro_cov_cweightlba0",
                      "macro_cov_cgestagewksa0",
                      "assessinterviewera4")], as.factor)

#check if all went well
sapply(FinalPD, class)

### Set Outcome 
outcome <- "treat" #this is our outcome (treatment/control), which we set as factor in line below
FinalPD$treat <- as.factor(FinalPD$treat)

### Run EWAS 
ptmp <- FinalPD #we have no missing variables, as we imputed the means. If we had, we would use here FinalPD[complete.cases(FinalPD$x) & complete.cases(FinalPD$xxx) & complete.cases(FinalPD$xxx),] 
myregressionmodel <- paste0("~", outcome, "+", paste(covars, collapse = "+") ) %>% as.formula()
mod <- model.matrix(myregressionmodel, ptmp)
btmp <- betas[, match((ptmp$BaseID), colnames(betas))]

# Run the single site association model
out <- lmFit(btmp, mod)
out <- eBayes(out)
ss.hits <- limma::topTable(out, coef = 2, number = nrow(out))

ss.hits.mums.excOutl <- ss.hits

save(ss.hits.mums.excOutl ,
     file ="  DNAm betas/ss.hits.mums.excOutl.rda") 



#### Create EWAS plots ####

#QQ plot - before bacon
#here you see some inflation, and lambda above 1.2 so we try bacon adjustment too

observed <- -log10(sort(ss.hits.mums.excOutl$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.mums.excOutl$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_mums_noOutliers.pdf") #open PDF device to save the plot

qq(ss.hits.mums.excOutl$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.


#### For Children ####

#Phenotype data
load("pheno_excldpl_kids.rda")
dim(pheno_excldpl_kids) #n=735 kids

#CpGs with Betas for Kids, these have been residualized for celltype, plate and slide
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_Kids.rda")
dim(betas_forEwas_excldpl_ResCellTech_Kids) #N=735 kids, with 278771 CpGs 
betas_forEwas_excldpl_ResCellTech_Kids[1:3, 1:3]

#### Load in baseline covariates from SILO ####
baseline_cov_all <- read_dta("  clean_data/age4_pheno_epi_merge.dta")
dim(baseline_cov_all)#N=1000, this is in wide format with every row representing a family. 

table(baseline_cov_all$treat, baseline_cov_all$treatmentgroupblindedC) #in pheno_excldpl_kids we have the blinded version, which has a slight diff sample size then "treat variable". We will use "Treat" in analyses

#select baseline covariates
#these are the baseline covariates that are part of the BFY-preregistration
baseline_cov_select <- baseline_cov_all[,c("sampleid",
                                           "treat",
                                           "magea0",
                                           "mnbiochilda0",
                                           "mgoodhealtha0",
                                           "mdepressiona0",
                                           "mhsgrada0",
                                           "msomecollegea0",
                                           "massociatesa0",
                                           "mbachelorsa0",
                                           "medunknowna0",  
                                           "mblacka0",
                                           "mracemultiplea0",
                                           "mraceothera0",
                                           "mhispanica0",
                                           "mraceunknowna0",
                                           "mcohaba0",
                                           "mmarrieda0",
                                           "mdivorceda0",
                                           "mrelateothera0", 
                                           "mrelateunknowna0",
                                           "mcigduringavgwka0",
                                           "malcduringavgwka0",
                                           "hhincomecat2a0",
                                           "hhincomecat3a0",
                                           "hhincomecat4a0",
                                           "hhincomecat5a0",
                                           "hhincomecat6a0",
                                           "hhnetworthcat2a0",
                                           "hhnetworthcat3a0",
                                           "hhnetworthcat4a0",
                                           "hhnetworthcat5a0",
                                           "hhnetworthcat6a0",
                                           "hhnadulta0",
                                           "hhbiodada0",
                                           "cfemalea0",
                                           "cweightlba0",
                                           "cgestagewksa0",
                                           "assessageinmonthsa4",
                                           "assessinterviewera4")]


#In line with earlier BFY analyses, we impute the mean for missing values on 4 covariats that have missing values for some families (mcigduringavgwka0, malcduringavgwka0, cweightlba0, cgestagewksa0)
baseline_cov_select <- baseline_cov_select %>%
  mutate(across(c(mcigduringavgwka0, malcduringavgwka0, cweightlba0, cgestagewksa0),
                list(miss_ = ~ ifelse(is.na(.), mean(., na.rm = TRUE), .),
                     macro_cov_ = ~ as.integer(is.na(.))),
                .names = "{fn}{.col}"))

#We also generate a dummy variable that reflects if the value was imputed, and these we add to the covariates in the EWAS too
macro_covariates <- names(baseline_cov_select)[grepl("^macro_cov_", names(baseline_cov_select))]

colnames(baseline_cov_select)

#see if imputation works, by checking mean original and newly created variable (means should be very similar)
describe(baseline_cov_select$mcigduringavgwka0)
describe(baseline_cov_select$miss_mcigduringavgwka0)


#### Create Phenotype datasets for Kids ####

#We have two phenotype datasets
# pheno_excldpl_kids = One from Max Planck Institute with participant id (sampleid) and with the connector to the CpGs (BaseID, so on which part of the plate they were)
# baseline_cov_select = One from SILO including all baseline covariates and whether participants were part of treatment of control group
# We have to merge these two on sampleid

head(pheno_excldpl_kids)
head(baseline_cov_select)

#We add a P to sampleid so the sample ids are the same across datasets
pheno_excldpl_kids_select <- pheno_excldpl_kids %>%
  mutate(sampleid = paste0("P", sampleId)) %>%
  select(1:2, sampleid, everything(),-sampleId) #remove old column

class(pheno_excldpl_kids_select$sampleid)
class(baseline_cov_select$sampleid)

#select the variables to merge
pheno_excldpl_kids_select <- pheno_excldpl_kids_select[,c("BaseID",
                                                          "sampleid")]

#merge the two phenofiles, with the pheno_excldpl_kids_select as the baseline
pheno_excldpl_kids_merged <- pheno_excldpl_kids_select %>%
  left_join(baseline_cov_select, by = "sampleid")

colnames(pheno_excldpl_kids_merged)
head(pheno_excldpl_kids_merged)

describe(pheno_excldpl_kids_merged)
table(pheno_excldpl_kids_merged$treat)

pheno_excldpl_kids_mergedBaseCov <- pheno_excldpl_kids_merged

save(pheno_excldpl_kids_mergedBaseCov , 
     file = "pheno_excldpl_kids_mergedBaseCov.rda")

load("    EWAS/pheno_excldpl_kids_merged.rda")

#### Run EWAS ####

## for Kids

#Rename Phenotype and CpG data
FinalPD <- pheno_excldpl_kids_merged
betas <-betas_forEwas_excldpl_ResCellTech_Kids[, match(FinalPD$BaseID, colnames(betas_forEwas_excldpl_ResCellTech_Kids)), drop = FALSE]

#We set beta file in the same order as the pheno file, so this order should be the same (otherwhise it gets scrambled)
identical(colnames(betas), FinalPD$BaseID) #should be TRUE

### Set Covars
#here we select all the baseline covariates
#those not #Continious will be set as.factor later in the script
covars <- c( "magea0",                       #Continious
             "mnbiochilda0",                 #Continious
             "mgoodhealtha0",  
             "mdepressiona0",                #Continious
             "mhsgrada0",      
             "msomecollegea0", 
             "massociatesa0",  
             "mbachelorsa0",  
             "medunknowna0",   
             "mblacka0",       
             "mracemultiplea0", 
             "mraceothera0",    
             "mhispanica0",     
             "mraceunknowna0", 
             "mcohaba0",        
             "mmarrieda0",      
             "mdivorceda0",     
             "mrelateothera0", 
             "mrelateunknowna0",
             "hhincomecat2a0",
             "hhincomecat3a0",
             "hhincomecat4a0",
             "hhincomecat5a0",
             "hhincomecat6a0",
             "hhnetworthcat2a0",
             "hhnetworthcat3a0",
             "hhnetworthcat4a0",
             "hhnetworthcat5a0",
             "hhnetworthcat6a0",
             "hhnadulta0",                  #Continious
             "hhbiodada0",
             "cfemalea0",
             "assessageinmonthsa4",         #Continious
             "assessinterviewera4", 
             "miss_mcigduringavgwka0",      #Continious
             "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
             "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
             "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
             "miss_cweightlba0",            #Continious,this is the mean imputed version
             "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
             "miss_cgestagewksa0",          #Continious, this is the mean imputed version
             "macro_cov_cgestagewksa0")     #This is the dummy if mean was imputed


#Set non-continuous variables to factor

#set them all as.factor in R
FinalPD[, c(  "mgoodhealtha0",  
              "mhsgrada0",      
              "msomecollegea0", 
              "massociatesa0",  
              "mbachelorsa0",  
              "medunknowna0",   
              "mblacka0",       
              "mracemultiplea0", 
              "mraceothera0",    
              "mhispanica0",     
              "mraceunknowna0", 
              "mcohaba0",        
              "mmarrieda0",      
              "mdivorceda0",     
              "mrelateothera0", 
              "mrelateunknowna0",
              "hhincomecat2a0",
              "hhincomecat3a0",
              "hhincomecat4a0",
              "hhincomecat5a0",
              "hhincomecat6a0",
              "hhnetworthcat2a0",
              "hhnetworthcat3a0",
              "hhnetworthcat4a0",
              "hhnetworthcat5a0",
              "hhnetworthcat6a0",
              "hhbiodada0",
              "cfemalea0",
              "assessinterviewera4",
              "macro_cov_mcigduringavgwka0",
              "macro_cov_malcduringavgwka0",
              "macro_cov_cweightlba0",
              "macro_cov_cgestagewksa0"
)] <-
  lapply (FinalPD[, c("mgoodhealtha0",  
                      "mhsgrada0",      
                      "msomecollegea0", 
                      "massociatesa0",  
                      "mbachelorsa0",  
                      "medunknowna0",   
                      "mblacka0",       
                      "mracemultiplea0", 
                      "mraceothera0",    
                      "mhispanica0",     
                      "mraceunknowna0", 
                      "mcohaba0",        
                      "mmarrieda0",      
                      "mdivorceda0",     
                      "mrelateothera0", 
                      "mrelateunknowna0",
                      "hhincomecat2a0",
                      "hhincomecat3a0",
                      "hhincomecat4a0",
                      "hhincomecat5a0",
                      "hhincomecat6a0",
                      "hhnetworthcat2a0",
                      "hhnetworthcat3a0",
                      "hhnetworthcat4a0",
                      "hhnetworthcat5a0",
                      "hhnetworthcat6a0",
                      "hhbiodada0",
                      "cfemalea0",
                      "assessinterviewera4",
                      "macro_cov_mcigduringavgwka0",
                      "macro_cov_malcduringavgwka0",
                      "macro_cov_cweightlba0",
                      "macro_cov_cgestagewksa0"
  )], as.factor)

#check if all went well
sapply(FinalPD, class)


#### Run EWAS for KIDS ####
### Set Outcome 
outcome <- "treat" #this is our outcome (treatment/control), which we set as factor in line below
FinalPD$treat <- as.factor(FinalPD$treat)

### Run EWAS 
ptmp <- FinalPD #we have no missing variables, as we imputed the means. If we had, we would use here FinalPD[complete.cases(FinalPD$x) & complete.cases(FinalPD$xxx) & complete.cases(FinalPD$xxx),] 
myregressionmodel <- paste0("~", outcome, "+", paste(covars, collapse = "+") ) %>% as.formula()
mod <- model.matrix(myregressionmodel, ptmp)
btmp <- betas[, match((ptmp$BaseID), colnames(betas))]

# Run the single site association model
out <- lmFit(btmp, mod)
out <- eBayes(out)
ss.hits <- limma::topTable(out, coef = 2, number = nrow(out))

ss.hits.kids.cov <- ss.hits

save(ss.hits.kids.cov , 
     file = "ss.hits.kids.cov.rda")

load("ss.hits.kids.cov.rda")

ss.hits.kids.cov.excel <- ss.hits.kids.cov %>%
  tibble::rownames_to_column(var = "CpGs")


write_xlsx(ss.hits.kids.cov.excel, "ss.hits.kids.cov.xlsx")


# Check if any sites are significant at first significance threshold
#See 00_DNAm_QC_general and QC EWAS for calculation p-value thresholds

#ComeBack
any(ss.hits.kids.cov$P.Value < 1.855701e-06)
num_significant <- sum(ss.hits.kids.cov$P.Value < 1.855701e-06, na.rm = TRUE)
num_significant

#MSD
any(ss.hits.kids.cov$P.Value < 6.157635e-05)
num_significant <- sum(ss.hits.kids.cov$P.Value < 6.157635e-05, na.rm = TRUE)
num_significant

#### Do Bacon Adjustment Kids ####

#bacon adjustment
bc <- bacon(teststatistics = ss.hits.kids.cov$t)
PVals <- pval(bc)
ss.hits.kids.cov.bacon <- data.frame(PVals)
row.names(ss.hits.kids.cov.bacon) <- row.names(ss.hits.kids.cov)
colnames(ss.hits.kids.cov.bacon) <- "P.Value"

save(ss.hits.kids.cov.bacon , 
     file = "ss.hits.kids.cov.bacon.rda")

load("ss.hits.kids.cov.bacon.rda")

ss.hits.kids.cov.bacon.excel <- ss.hits.kids.cov.bacon %>%
  tibble::rownames_to_column(var = "CpGs")


write_xlsx(ss.hits.kids.cov.bacon.excel, "ss.hits.kids.cov.bacon.xlsx")

# Check if any sites are significant at first significance threshold
#Bonferroni
any(ss.hits.kids.cov.bacon$P.Value < 1.793587e-07)
num_significant <- sum(ss.hits.kids.cov.bacon$P.Value <  1.793587e-07, na.rm = TRUE)
num_significant

#ComeBack
any(ss.hits.kids.cov.bacon$P.Value < 1.855701e-06)
num_significant <- sum(ss.hits.kids.cov.bacon$P.Value < 1.855701e-06, na.rm = TRUE)
num_significant

#MSD
any(ss.hits.kids.cov.bacon$P.Value < 6.157635e-05)
num_significant <- sum(ss.hits.kids.cov.bacon$P.Value < 6.157635e-05, na.rm = TRUE)
num_significant


#### Create EWAS plots ####

#### QQ-plot ####

#QQ plot - before bacon
#here you see some inflation, and lambda above 1.2 so we try bacon adjustment too

observed <- -log10(sort(ss.hits.kids.cov$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.kids.cov$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_kids.pdf") #open PDF device to save the plot

qq(ss.hits.kids.cov$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#QQ plot - after bacon adjustment

observed <- -log10(sort(ss.hits.kids.cov.bacon$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.kids.cov.bacon$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_kids_bacon.pdf") #open PDF device to save the plot

qq(ss.hits.kids.cov.bacon$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.

#### Volcano Plot ####

pdf("volcano_zscore_kids.pdf") #open PDF device to save the plot

ggplot(ss.hits.kids.cov, aes(x = ((t-50)/10) , y = -log10(P.Value))) +
  geom_point() +
  theme_bw() +
  labs(x = "Z-Score", y = "-log10 P Value")

dev.off() #closes the PDF device, ensuring that the file is properly saved.

pdf("volcanoplot_kids.pdf") #open PDF device to save the plot

ggplot(ss.hits.kids.cov, aes(x = logFC , y = -log10(P.Value))) +
  geom_point() +
  theme_bw() +
  labs(x = "Log Fold Change", y = "-log10 P Value", title = "Volcano Plot", subtitle = "", caption = "" )

dev.off() #closes the PDF device, ensuring that the file is properly saved.


### Residualize CpGs for baseline covariates ####

#Phenodata (this is phenodata including baseline covariates)
load("  pheno_excldpl_kids_mergedBaseCov.rda") #this is the overall phenotype data
colnames(pheno_excldpl_kids_mergedBaseCov)
dim(pheno_excldpl_kids_mergedBaseCov)

#Beta data (these are betas selected on variability and with ICC >.50)
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_Kids.rda")
dim(betas_forEwas_excldpl_ResCellTech_Kids)
betas_forEwas_excldpl_ResCellTech_Kids[1:5, 1:5]

#check if order is correct
identical(colnames(betas_forEwas_excldpl_ResCellTech_Kids),pheno_excldpl_kids_mergedBaseCov$BaseID) #should be TRUE

all(pheno_excldpl_kids_mergedBaseCov$BaseID == colnames(betas_forEwas_excldpl_ResCellTech_Kids)) # should be TRUE

#Correct CpGs for baseline covariates (they have already been residualized for celltype / plate / array)

# Define batch size (e.g., process 1000 CpGs at a time), it is computational heavy so that is why we run it across 1000 CpGs (outcomes do not depend on how large these chunks are, I checked it for chunks of 10, 100, 1000 and results are the same)
batch_size <- 1000  

# Get total number of CpGs
num_cpgs <- nrow(betas_forEwas_excldpl_ResCellTech_Kids)

# Create an empty list to store results
residuals_list <- vector("list", length = ceiling(num_cpgs / batch_size))

baselinecov_data <- pheno_excldpl_kids_mergedBaseCov[,c("BaseID",
                                                        "magea0",                       #Continious
                                                        "mnbiochilda0",                 #Continious
                                                        "mgoodhealtha0",  
                                                        "mdepressiona0",                #Continious
                                                        "mhsgrada0",      
                                                        "msomecollegea0", 
                                                        "massociatesa0",  
                                                        "mbachelorsa0",  
                                                        "medunknowna0",   
                                                        "mblacka0",       
                                                        "mracemultiplea0", 
                                                        "mraceothera0",    
                                                        "mhispanica0",     
                                                        "mraceunknowna0", 
                                                        "mcohaba0",        
                                                        "mmarrieda0",      
                                                        "mdivorceda0",     
                                                        "mrelateothera0", 
                                                        "mrelateunknowna0",
                                                        "hhincomecat2a0",
                                                        "hhincomecat3a0",
                                                        "hhincomecat4a0",
                                                        "hhincomecat5a0",
                                                        "hhincomecat6a0",
                                                        "hhnetworthcat2a0",
                                                        "hhnetworthcat3a0",
                                                        "hhnetworthcat4a0",
                                                        "hhnetworthcat5a0",
                                                        "hhnetworthcat6a0",
                                                        "hhnadulta0",                  #Continious
                                                        "hhbiodada0",
                                                        "cfemalea0",
                                                        "assessageinmonthsa4",         #Continious
                                                        "assessinterviewera4", 
                                                        "miss_mcigduringavgwka0",      #Continious
                                                        "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
                                                        "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
                                                        "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
                                                        "miss_cweightlba0",            #Continious,this is the mean imputed version
                                                        "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
                                                        "miss_cgestagewksa0",          #Continious, this is the mean imputed version
                                                        "macro_cov_cgestagewksa0")]     #This is the dummy if mean was imputed




baselinecov_data <- as.data.frame(baselinecov_data)
rownames(baselinecov_data) <- baselinecov_data$BaseID  
baselinecov_data$BaseID <- NULL

# Loop over CpGs in chunks
for (i in seq(1, num_cpgs, by = batch_size)) {
  # Define CpG subset (batch)
  batch_end <- min(i + batch_size - 1, num_cpgs)  # Ensure we don't exceed total CpGs
  cpg_data_batch <- betas_forEwas_excldpl_ResCellTech_Kids[i:batch_end, ]
  
  # Transpose batch data so rows = participants, columns = CpGs
  cpg_data_t_batch <- t(cpg_data_batch)
  
  # Ensure participant order is the same
  cpg_data_t_batch <- cpg_data_t_batch[rownames(baselinecov_data), ]
  
  # Function to regress each CpG on baseline covariates and extract residuals
  
  get_residuals <- function(cpg_values) {
    model <- lm(cpg_values ~ magea0 +  mnbiochilda0 + as.factor(mgoodhealtha0) + mdepressiona0 +              
                  as.factor(mhsgrada0) +      
                  as.factor(msomecollegea0)+ 
                  as.factor(massociatesa0)+  
                  as.factor(mbachelorsa0)+  
                  as.factor(medunknowna0)+   
                  as.factor(mblacka0)+       
                  as.factor(mracemultiplea0)+ 
                  as.factor(mraceothera0)+    
                  as.factor(mhispanica0)+     
                  as.factor(mraceunknowna0)+ 
                  as.factor(mcohaba0)+        
                  as.factor(mmarrieda0)+      
                  as.factor(mdivorceda0)+     
                  as.factor(mrelateothera0)+ 
                  as.factor(mrelateunknowna0)+
                  as.factor(hhincomecat2a0)+
                  as.factor(hhincomecat3a0)+
                  as.factor(hhincomecat4a0)+
                  as.factor(hhincomecat5a0)+
                  as.factor(hhincomecat6a0)+
                  as.factor(hhnetworthcat2a0)+
                  as.factor(hhnetworthcat3a0)+
                  as.factor(hhnetworthcat4a0)+
                  as.factor(hhnetworthcat5a0)+
                  as.factor(hhnetworthcat6a0)+
                  hhnadulta0 +                                 #Continious
                  as.factor(hhbiodada0)+
                  as.factor(cfemalea0)+
                  assessageinmonthsa4+                         #Continious
                  as.factor(assessinterviewera4)+ 
                  miss_mcigduringavgwka0+                      #Continious
                  as.factor(macro_cov_mcigduringavgwka0) +     #This is the dummy if mean was imputed
                  miss_malcduringavgwka0+                      #Continious, this is the mean imputed version
                  as.factor(macro_cov_malcduringavgwka0) +     #This is the dummy if mean was imputed
                  miss_cweightlba0+                            #Continious,this is the mean imputed version
                  macro_cov_cweightlba0+                       #This is the dummy if mean was imputed
                  miss_cgestagewksa0+                          #Continious, this is the mean imputed version
                  as.factor(macro_cov_cgestagewksa0),          #This is the dummy if mean
                data = baselinecov_data)
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

betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids <- residuals_df

# Show first few values
print(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[1:4, 1:4])  

save(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids , 
     file = "  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids.rda")

#### Check for outliers in top hits ####

### Load in the data

#load phenotype date with treatment variable
load("  pheno_excldpl_kids_mergedBaseCov.rda")

#load beta data residualized for baseline covariates
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids.rda")

#load tophits after bacon adjustment
load("  ss.hits.kids.cov.bacon.rda")
CpGhitsBacon <- ss.hits.kids.cov.bacon %>%
  arrange(P.Value)

### Create dataset with only tophits

#select 15 hits that were significant in EWAS
top_cpgs <- rownames(CpGhitsBacon)[1:15]

#create betadata just for those 15 hits
Beta_tophits <- as.data.frame(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[top_cpgs,])
head(Beta_tophits)

#transpose the data
Beta_tophits_transposed <- as.data.frame(t(Beta_tophits))

#turn rownames into collumn names
Beta_tophits_transposed$BaseID <- rownames(Beta_tophits_transposed)
rownames(Beta_tophits_transposed) <- NULL

#merge treatment variable to it
hits_treat <- Beta_tophits_transposed %>%
  left_join(pheno_excldpl_kids_mergedBaseCov %>% select(BaseID, treat), by = "BaseID")

table(hits_treat$treat) #0 is control, 1 is treatment

#### Plot methylation values ####
library(ggplot2)
colnames(hits_treat)

hits_treat$treat <- as.factor(hits_treat$treat)


ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg06152865_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg06152865_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg20272170_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg20272170_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg14480531_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg14480531_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg11040181_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg11040181_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg08861115_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg08861115_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg17280975_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg17280975_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg24450643_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg24450643_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg11588197_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg11588197_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg01854039_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg01854039_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg21190595_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg21190595_TC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg01974334_BC22)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg01974334_BC22 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg26607429_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg26607429_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg21932360_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg21932360_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg15584219_BC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg15584219_BC21 by Treatment Group") +
  theme_minimal()

ggplot(hits_treat, aes(x= factor(treat, levels = c(0,1), labels = c("Control", "Treatment")), y = cg16586756_TC21)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "blue") +
  labs(x = "Group", y = "Methylation Level", title = "Methylation cg16586756_TC21 by Treatment Group") +
  theme_minimal()


### now do it based on sd ####
library(ggplot2)
library(dplyr)

cpg_list <- c("cg06152865_BC21", "cg20272170_TC21", "cg14480531_TC21", "cg11040181_BC21", "cg08861115_TC21", "cg17280975_BC21", "cg24450643_BC21", "cg11588197_BC21", "cg01854039_BC21", "cg21190595_TC21", "cg01974334_BC22", "cg26607429_BC21",
              "cg21932360_BC21", "cg15584219_BC21", "cg16586756_TC21")

for (cpg in cpg_list) {
  
  plot_data <- hits_treat %>%
    group_by(treat) %>%
    mutate(
      mean_val = mean(.data[[cpg]], na.rm = TRUE),
      sd_val = sd(.data[[cpg]], na.rm = TRUE),
      is_outlier_sd = abs(.data[[cpg]] - mean_val) > 3 * sd_val
    ) %>%
    ungroup()
  
  p <- ggplot(plot_data, aes(x = factor(treat, levels = c(0,1), labels = c("Control", "Treatment")),
                             y = .data[[cpg]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = filter(plot_data, is_outlier_sd),
               aes(y = .data[[cpg]]), color = "black", shape = 16, size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
    labs(x = "Group", y = "Methylation Level",
         title = paste("Methylation", cpg, "by Treatment Group")) +
    theme_minimal()
  
  print(p)
  
}


#### Run QQplot but excluding outlier CpGS ####

#check which cpgs have extreme outliers >5 sd


extreme_flags <-apply(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids, 1, function(row) {
  mu <- mean(row, na.rm =TRUE)
  sigma <- sd(row, na.rm = TRUE)
  
  if(is.na(sigma) || sigma ==0) return(FALSE)
  any(abs(row-mu) >5 *sigma, na.rm = TRUE)
})

sum(extreme_flags)

cpgs_with_extreme_outliers <- rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids)[extreme_flags]

sum(Beta_tophits %in% cpgs_with_extreme_outliers) #none of them are in the tophits


#create beta file excluding cpgs with extreme outliers
betas_nooutliers <- betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[!extreme_flags, ]


##### Run EWAS with qqplot on data without those outliers ####

load("  pheno_excldpl_kids_mergedBaseCov.rda")


#Rename Phenotype and CpG data
FinalPD <- pheno_excldpl_kids_mergedBaseCov
betas <-betas_nooutliers[, match(FinalPD$BaseID, colnames(betas_nooutliers)), drop = FALSE]

#we set beta file in the same order as the pheno file, so this order should be the same (otherwhise it gets scrambled)
identical(colnames(betas), FinalPD$BaseID) #should be TRUE

### Set Covars
#here we select all the baseline covariates
#those not #Continious will be set as.factor later in the script
covars <- c( "magea0",                       #Continious
             "mnbiochilda0",                 #Continious
             "mgoodhealtha0",  
             "mdepressiona0",                #Continious
             "mhsgrada0",      
             "msomecollegea0", 
             "massociatesa0",  
             "mbachelorsa0",  
             "medunknowna0",   
             "mblacka0",       
             "mracemultiplea0", 
             "mraceothera0",    
             "mhispanica0",     
             "mraceunknowna0", 
             "mcohaba0",        
             "mmarrieda0",      
             "mdivorceda0",     
             "mrelateothera0", 
             "mrelateunknowna0",
             "hhincomecat2a0",
             "hhincomecat3a0",
             "hhincomecat4a0",
             "hhincomecat5a0",
             "hhincomecat6a0",
             "hhnetworthcat2a0",
             "hhnetworthcat3a0",
             "hhnetworthcat4a0",
             "hhnetworthcat5a0",
             "hhnetworthcat6a0",
             "hhnadulta0",                  #Continious
             "hhbiodada0",
             "cfemalea0",
             "assessageinmonthsa4",         #Continious
             "miss_mcigduringavgwka0",      #Continious
             "macro_cov_mcigduringavgwka0", #This is the dummy if mean was imputed
             "miss_malcduringavgwka0",      #Continious, this is the mean imputed version
             "macro_cov_malcduringavgwka0", #This is the dummy if mean was imputed
             "miss_cweightlba0",            #Continious,this is the mean imputed version
             "macro_cov_cweightlba0",       #This is the dummy if mean was imputed
             "miss_cgestagewksa0",          #Continious, this is the mean imputed version
             "macro_cov_cgestagewksa0",     #This is the dummy if mean was imputed
             "assessinterviewera4")

#Set non-continuous variables to factor
#set them all as.factor in R
FinalPD[, c(  "mgoodhealtha0",  
              "mhsgrada0",      
              "msomecollegea0", 
              "massociatesa0",  
              "mbachelorsa0",  
              "medunknowna0",   
              "mblacka0",       
              "mracemultiplea0", 
              "mraceothera0",    
              "mhispanica0",     
              "mraceunknowna0", 
              "mcohaba0",        
              "mmarrieda0",      
              "mdivorceda0",     
              "mrelateothera0", 
              "mrelateunknowna0",
              "hhincomecat2a0",
              "hhincomecat3a0",
              "hhincomecat4a0",
              "hhincomecat5a0",
              "hhincomecat6a0",
              "hhnetworthcat2a0",
              "hhnetworthcat3a0",
              "hhnetworthcat4a0",
              "hhnetworthcat5a0",
              "hhnetworthcat6a0",
              "hhbiodada0",
              "cfemalea0",
              "macro_cov_mcigduringavgwka0",
              "macro_cov_malcduringavgwka0",
              "macro_cov_cweightlba0",
              "macro_cov_cgestagewksa0",
              "assessinterviewera4")] <-
  lapply (FinalPD[, c("mgoodhealtha0",  
                      "mhsgrada0",      
                      "msomecollegea0", 
                      "massociatesa0",  
                      "mbachelorsa0",  
                      "medunknowna0",   
                      "mblacka0",       
                      "mracemultiplea0", 
                      "mraceothera0",    
                      "mhispanica0",     
                      "mraceunknowna0", 
                      "mcohaba0",        
                      "mmarrieda0",      
                      "mdivorceda0",     
                      "mrelateothera0", 
                      "mrelateunknowna0",
                      "hhincomecat2a0",
                      "hhincomecat3a0",
                      "hhincomecat4a0",
                      "hhincomecat5a0",
                      "hhincomecat6a0",
                      "hhnetworthcat2a0",
                      "hhnetworthcat3a0",
                      "hhnetworthcat4a0",
                      "hhnetworthcat5a0",
                      "hhnetworthcat6a0",
                      "hhbiodada0",
                      "cfemalea0",
                      "macro_cov_mcigduringavgwka0",
                      "macro_cov_malcduringavgwka0",
                      "macro_cov_cweightlba0",
                      "macro_cov_cgestagewksa0",
                      "assessinterviewera4")], as.factor)

#check if all went well
sapply(FinalPD, class)

### Set Outcome 
outcome <- "treat" #this is our outcome (treatment/control), which we set as factor in line below
FinalPD$treat <- as.factor(FinalPD$treat)

### Run EWAS 
ptmp <- FinalPD #we have no missing variables, as we imputed the means. If we had, we would use here FinalPD[complete.cases(FinalPD$x) & complete.cases(FinalPD$xxx) & complete.cases(FinalPD$xxx),] 
myregressionmodel <- paste0("~", outcome, "+", paste(covars, collapse = "+") ) %>% as.formula()
mod <- model.matrix(myregressionmodel, ptmp)
btmp <- betas[, match((ptmp$BaseID), colnames(betas))]

# Run the single site association model
out <- lmFit(btmp, mod)
out <- eBayes(out)
ss.hits <- limma::topTable(out, coef = 2, number = nrow(out))

ss.hits.kids.excOutl <- ss.hits

save(ss.hits.kids.excOutl ,
     file ="  DNAm betas/ss.hits.kids.excOutl.rda") 



#### Create EWAS plots ####

#QQ plot - before bacon
#here you see some inflation, and lambda above 1.2 so we try bacon adjustment too

observed <- -log10(sort(ss.hits.kids.excOutl$P.Value, decreasing = F))
expected <- -log10(ppoints(length(ss.hits.kids.excOutl$P.Value)))
lambda <- median(observed) / median(expected)

pdf("QQplot_mums_noOutliers.pdf") #open PDF device to save the plot

qq(ss.hits.kids.excOutl$P.Value, main = sprintf("Lambda value of %.3f", lambda))

dev.off() #closes the PDF device, ensuring that the file is properly saved.


#### Sensitivity analyses for smoking ####

#####Loading the Adjusted CpGs####

#KIDS
#ss.hits.kids.cov.bacon.xlsx"

Kids_Adjusted_CpGs <- read_excel("ss.hits.kids.cov.bacon.xlsx")  

#Mothers
#  ss.hits.mums.cov.bacon.xlsx"
Moms_Adjusted_CpGs <- read_excel("ss.hits.mums.cov.bacon.xlsx")  




### Smoking CpGs (Van Dongen 2023) ####


#Loading list of CpGs from: https://doi.org/10.7554/eLife.83286
#Saved here:   Smoking/elife-Smoking_CpGs_Jenny.xlsx
#Effects of smoking on genome-wide DNA methylation profiles: A study of discordant and concordant monozygotic twin pairs

#In pairs discordant for current smoking, 13 differentially methylated CpGs were 
#found between current smoking twins and their genetically identical co-twin who never smoked.

#All 13 CpGs have been previously associated with smoking in unrelated individuals 
#and data from monozygotic pairs discordant for former smoking indicated that 
#methylation patterns are to a large extent reversible upon smoking cessation.


Smoking_CpGs <- read_excel("elife-Smoking_CpGs_Jenny.xlsx")  



# Convert Smoking_CpGs to a dataframe  
Smoking_CpGs <- as.data.frame(Smoking_CpGs)  

# Set row names of the data frame to the values in the "IlmnID" column  
rownames(Smoking_CpGs) <- Smoking_CpGs$IlmnID  

# Remove the "IlmnID" column as it's now redundant  
Smoking_CpGs$IlmnID <- NULL  



######Overlap between Kids CpGs and smoking---------

# Split the rownames of "Kids_Adjusted_CpGs" at the underscore and keep only the part before the underscore  
kids_row_names <- sapply(strsplit(rownames(Kids_Adjusted_CpGs), "_"), `[`, 1)  

# Find the common row names between "kids_row_names" and the rownames of "Smoking_CpGs"  
common_rows <- intersect(kids_row_names, rownames(Smoking_CpGs))  

# Print the number of common rows  
print(length(common_rows))  #10

#"cg21161138" "cg19089201" "cg09935388" "cg01940273" "cg21188533" "cg22132788" "cg00336149" "cg05575921" "cg01901332" "cg21566642"

common_rows_indices <- which(kids_row_names %in% common_rows)
#"cg21161138" "cg19089201" "cg09935388" "cg01940273" "cg21188533" "cg22132788" "cg00336149" "cg05575921" "cg01901332" "cg21566642"
#[1]   6067   9253  11173  36444  64841  77051  81124  85423 197741 257857



######Overlap between Moms CpGs and smoking---------

# Split the rownames of "Moms_Adjusted_CpGs" at the underscore and keep only the part before the underscore  
moms_row_names <- sapply(strsplit(rownames(Moms_Adjusted_CpGs), "_"), `[`, 1)  

# Find the common row names between "moms_row_names" and the rownames of "Smoking_CpGs"  
common_rows_mom <- intersect(moms_row_names, rownames(Smoking_CpGs))  

# Print the number of common rows  
print(length(common_rows_mom))  #10

#"cg19089201" "cg22132788" "cg05575921" "cg01901332" "cg21188533" "cg21566642" "cg21161138" "cg01940273" "cg00336149" "cg09935388"

common_rows_indices_moms <- which(moms_row_names %in% common_rows_mom)

#"cg19089201" "cg22132788" "cg05575921" "cg01901332" "cg21188533" "cg21566642" "cg21161138" "cg01940273" "cg00336149" "cg09935388"
#457   9296  10095  96403 132403 154862 161971 238386 240804 250849




### Smoking CpGs (Barcelona et al., 2019) ####


#Loading list of CpGs from: https://doi.org/10.1080/15592294.2019.1588683
#Saved here: ~/Projects/04_data_analysis/003_BY/04_EWAS/EWAS bioannotation/Smoking/Smoking_CpGs_Barcelona_2019.csv
#Novel DNA methylation sites associated with cigarette smoking among African Americans

#After controlling for age, body mass index, population structure and cell composition, 
#26 epigenome-wide significant sites (FDR q < 0.05) were identified, including the AHRR and PHF14 genes associated with atherosclerosis and lung disease, respectively.


Smoking_CpGs_2019 <- read.csv("Smoking_CpGs_Barcelona_2019.csv")  

# Set column names to the first row
colnames(Smoking_CpGs_2019) <- as.character(unlist(Smoking_CpGs_2019[1, ]))

# Remove the first row from the data (since it's now in colnames)
Smoking_CpGs_2019 <- Smoking_CpGs_2019[-1, ]

# Reset rownames (optional, to have consecutive numbers)
rownames(Smoking_CpGs_2019) <- NULL 

#Now we 26 CpGs


# Convert Smoking_CpGs to a dataframe  
Smoking_CpGs_2019 <- as.data.frame(Smoking_CpGs_2019)  

# Set row names of the data frame to the values in the "IlmnID" column  
rownames(Smoking_CpGs_2019) <- Smoking_CpGs_2019$CpG 



######Overlap between Kids CpGs and smoking---------

# Set row names of the "Kids_Adjusted_CpGs" to the values in the "CpGs" column  
rownames(Kids_Adjusted_CpGs) <- Kids_Adjusted_CpGs$CpGs 

# Split the rownames of "Kids_Adjusted_CpGs" at the underscore and keep only the part before the underscore  
kids_row_names <- sapply(strsplit(rownames(Kids_Adjusted_CpGs), "_"), `[`, 1)  

# Find the common row names between "kids_row_names" and the rownames of "Smoking_CpGs_2019"  
common_rows <- intersect(kids_row_names, rownames(Smoking_CpGs_2019))  

# Print the number of common rows  
print(length(common_rows))  #22

# [1] "cg09935388" "cg22996023" "cg12956751" "cg14051805" "cg17739917" "cg14753356" "cg01940273" "cg26703534" "cg25189904"
# [10] "cg05575921" "cg14389122" "cg27174698" "cg07439098" "cg27241845" "cg00748718" "cg16937168" "cg00073090" "cg22063959"
# [19] "cg07824483" "cg21566642" "cg19695041" "cg07741821"



common_rows_indices <- which(kids_row_names %in% common_rows)
# [1]  11173  11994  18280  25233  34076  34743  36444  39447  75416  85423 100679 106231 133639 158881 180743 216416 216570 239557
# [19] 239944 257857 259333 266460


######Overlap between Moms CpGs and smoking---------

# Set row names of the "Kids_Adjusted_CpGs" to the values in the "CpGs" column  
rownames(Moms_Adjusted_CpGs) <- Moms_Adjusted_CpGs$CpGs 

# Split the rownames of "Moms_Adjusted_CpGs" at the underscore and keep only the part before the underscore  
moms_row_names <- sapply(strsplit(rownames(Moms_Adjusted_CpGs), "_"), `[`, 1)  

# Find the common row names between "moms_row_names" and the rownames of "Smoking_CpGs_2019"  
common_rows_mom <- intersect(moms_row_names, rownames(Smoking_CpGs_2019))  

# Print the number of common rows  
print(length(common_rows_mom))  #22

# [1] "cg05575921" "cg12956751" "cg17739917" "cg14753356" "cg14389122" "cg00748718" "cg27174698" "cg25189904" "cg07741821"
# [10] "cg07439098" "cg14051805" "cg22063959" "cg27241845" "cg21566642" "cg16937168" "cg26703534" "cg00073090" "cg22996023"
# [19] "cg01940273" "cg09935388" "cg19695041" "cg07824483"


common_rows_indices_moms <- which(moms_row_names %in% common_rows_mom)

# [1]  10095  15950  53422  55821  57443  66932  84350  92311 104598 114431 132878 141070 147566 154862 158679 205342 228409 236898
# [19] 238386 250849 261203 278527


### Smoking CpGs (Dawes et al. 2021) ####


#Loading list of CpGs from: https://doi.org/10.1038/s41598-021-01088-7
#The relationship of smoking to cg05575921 methylation in blood and saliva DNA samples from several studies

#Just one CPG: cg05575921

#Rank in kids:85423
#Rank in moms:10095

#### Overlap CpGS EWAS Mothers and Kids ####

#### Overlap cpgs children and mums ####

load("/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/betas.top30000.kids.rda")
load("/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/betas.top30000.mums.rda")

betas.top30000.kids[1:5, 1:5]
betas.top30000.mums[1:5, 1:5]

#get CpG names
cpgs_kids <- rownames(betas.top30000.kids)
cpgs_mums <- rownames(betas.top30000.mums)

#find overlapping CpGs
overlappingCpGs <- intersect(cpgs_kids, cpgs_mums)


length(overlappingCpGs) #4296 / 30000


###### Correlation between tophits from EWAS ####

#### Install & load-in Packages ####
install.packages("apaTables")
install.packages("openxlsx")
install.packages("psych")
library(apaTables)
library(psych)
library(openxlsx)


#### For Mothers ####

#load beta data residualized for baseline covariates
#I picked the one residualzied for baseline covariates as this allows for cleaner correaltions (on the methylation level not confoudned by differences on baseline covariates)
load("/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_mums.rda")
betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[1:5, 1:5]


#create dataset only with top hits from the EWAS
cpgs_to_keep <- c(
  "cg13412754_TC21",
  "cg19137044_TC21",
  "cg25576801_TC21",
  "cg04188396_BC21",
  "cg01144399_BC21",
  "cg03786343_BC21",
  "cg15462959_BC21",
  "cg13412754_TC22",
  "cg09739347_TC21",
  "cg16717411_BC21",
  "cg01801101_BC21",
  "cg14413795_TC21",
  "cg12795046_TC21",
  "cg14641625_BC21",
  "cg16483795_BC21",
  "cg17495671_BC21",
  "cg25314817_BC21")

#check if all of these are in the betafile
setdiff(cpgs_to_keep, rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums))

#create subset data with betas but 
cpg_subset <- betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[cpgs_to_keep, ,drop = FALSE]
dim(cpg_subset) #17 CpGs for 777 mums

betas_tophits_mums <- cpg_subset
save(betas_tophits_mums ,
     file ="/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/betas_tophits_mums.rda") 

#Create correlation table

cpg_cor_test <- corr.test(
  t(cpg_subset),
  use = "pairwise")

out_dir <- "/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025"

apa_cor <- apa.cor.table(
  t(cpg_subset),
  table.number = 1)

apa_cor <- apa.cor.table(
  t(cpg_subset),
  table.number = 1,
  filename = file.path(out_dir, "CpGcors_tophits_Mums.rtf"))




#### For Children ####

#load beta data residualized for baseline covariates
#I picked the one residualzied for baseline covariates as this allows for cleaner correaltions (on the methylation level not confoudned by differences on baseline covariates)
load("/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids.rda")
betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[1:5, 1:5]

#create dataset only with top hits from the EWAS
cpgs_to_keep <- c(
  "cg06152865_BC21",
  "cg20272170_TC21",
  "cg14480531_TC21",
  "cg11040181_BC21",
  "cg08861115_TC21",
  "cg17280975_BC21",
  "cg24450643_BC21",
  "cg11588197_BC21",
  "cg01854039_BC21",
  "cg21190595_TC21",
  "cg01974334_BC22",
  "cg26607429_BC21",
  "cg21932360_BC21",
  "cg15584219_BC21",
  "cg16586756_TC21")

#check if all of these are in the betafile
setdiff(cpgs_to_keep, rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids))

#create subset data with betas 
cpg_subset <- betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[cpgs_to_keep, ,drop = FALSE]
dim(cpg_subset) #15 CpGs for 735 mums

betas_tophits_kids <- cpg_subset
save(betas_tophits_kids ,
     file ="/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025/betas_tophits_kids.rda") 


#Create correlation table

cpg_cor_test <- corr.test(
  t(cpg_subset),
  use = "pairwise")

out_dir <- "/data4/project/irp4/BabyFirstYear/Age4_epigenetics/Yayouk/EWAS/EWAS April 2025"

apa_cor <- apa.cor.table(
  t(cpg_subset),
  table.number = 1)

apa_cor <- apa.cor.table(
  t(cpg_subset),
  table.number = 1,
  filename = file.path(out_dir, "CpGcors_tophits_Kids.rtf"))












######################## The END ########################