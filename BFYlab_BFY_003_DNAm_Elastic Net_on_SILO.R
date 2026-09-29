################################################### 
########  ElasticNet in Mothers and Children  #####
################################################### 

# Load necessary libraries
install.packages("glmnet")
install.packages("tidyverse")
install.packages("caret")
install.packages("pROC")
library(pROC)
library(glmnet)
library(tidyverse)
library(caret)


####### For Mums ####


#### Load in Beta and Phenotype data ####


#load phenotype date with treatment variable
load("  pheno_excldpl_mums_mergedBaseCov.rda")
dim(pheno_excldpl_mums_mergedBaseCov)
pheno_excldpl_mums_mergedBaseCov[1:5, 1:5]

#load beta data residualized for cell and tech comp and baseline covariates
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_mums.rda")
dim(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums)
betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[1:5, 1:5]

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums), pheno_excldpl_mums_mergedBaseCov$BaseID) #should be True

#### Select CpGs that are both on EpicV1 and EpicV2 ####

#load in CpGs specific to Epic v1
cpg_sites_EpicV1 <- readRDS("    cpg_sites_EpicV1.rds")

cpg_sites <- cpg_sites_EpicV1

# Remove suffixes (_B11, _B12, etc.) from rownames in your beta file (cause those suffixes are specific to epicv2)
rownames_clean <- sub("_[^_]+$", "", rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums))

# Create a logical vector for matching
rows_to_keep <- rownames_clean %in% cpg_sites$x

# Subset the dataframe
betas_filtered <- betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[rows_to_keep, ]
dim(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums) #282371 CpGs
dim(betas_filtered) #22217 CpGs that are on Epic V1 and Epic V2

#select only the top 30.000 cpgs to go into the model
#first load in the top hits after bacon adjustment
ss_hits_mums_cov_bacon <- read_excel("  ss.hits.mums.cov.bacon.xlsx")

#Order hits data by p-value
ss_hits_ordered <- ss_hits_mums_cov_bacon[order(ss_hits_mums_cov_bacon$P.Value), ]

# Get the CpGs present in both datasets
common_cpgs <- intersect(rownames(betas_filtered), ss_hits_ordered$CpGs)

# Keep only the hits that are in the beta file
hits_in_betas <- ss_hits_ordered[ss_hits_ordered$CpGs %in% common_cpgs, ]

#Select the top 30,000 CpGs from those that match
ss_top30000_hits <- hits_in_betas[1:30000, ]

# Subset the beta file to these top 30,000 CpGs
betas_top30000 <- betas_filtered[rownames(betas_filtered) %in% ss_top30000_hits$CpGs, ]
dim(betas_top30000)

#check if a couple of tophits from EWAS are in these top 30.000
"cg13412754_TC21" %in% rownames(betas_top30000)
"cg04188396_BC21" %in% rownames(betas_top30000) #all True, meaning top EWAS hits on epicv1/v2 are in the betafile
"cg25576801_TC21" %in% rownames(betas_top30000)
"cg13412754_TC22" %in% rownames(betas_top30000)
"cg16717411_BC21" %in% rownames(betas_top30000)
"cg01801101_BC21" %in% rownames(betas_top30000)
"cg17495671_BC21" %in% rownames(betas_top30000)

#while epicV2 specific hits should not be in there
"cg19137044_TC21" %in% rownames(betas_top30000) #these yield False, cause indeed these are specific to epic v2 and we dont want them in our beta data
"cg14413795_TC21" %in% rownames(betas_top30000)

#### Perform stratified Random Sampling ####

#make sure its set as a factor
table(pheno_excldpl_mums_mergedBaseCov$treat) #0= control, 1=treatment
pheno_excldpl_mums_mergedBaseCov$treat <- as.factor(pheno_excldpl_mums_mergedBaseCov$treat)

pheno = pheno_excldpl_mums_mergedBaseCov
beta = betas_top30000
dim(beta)


# Create stratified partition
train_index <- createDataPartition(pheno$treat, p = 0.75, list = FALSE) #change this to actual treatment variable 

# Subset the data into training and tune sets
pheno_train <- pheno[train_index, ]
pheno_tune <- pheno[-train_index, ]

# Check distribution in the original dataset
prop.table(table(pheno$treat))

# Check distribution of treatment in the training set
prop.table(table(pheno_train$treat))

# Check distribution of treatment in the tune set
prop.table(table(pheno_tune$treat))

#create beta files for the training and tune data

# Subset the beta matrix for training samples
beta_train <- beta[, train_index]
dim(beta_train)
dim(beta)

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(beta_train), pheno_train$BaseID) #should be True

# Subset the beta matrix for tune samples
beta_tune <- beta[, -train_index]

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(beta_tune), pheno_tune$BaseID)

dim(beta_tune)
dim(beta_train)

##### Running Elastic Net ####

#### Assing 10-fold cross-tune ####
foldcount <- 10
set.seed(123)
foldAssignments <- sample(rep(1:foldcount, length.out = ncol(beta_train)))

#### Run Elastic Net model ####

cv.fit <- cv.glmnet(
  x = t(beta_train),
  y = pheno_train$treat,
  foldid = foldAssignments,
  family = "binomial",
  type.measure = "auc",   # or "class"
  alpha = 0.5,
  nfolds = foldcount,
  parallel = FALSE,       # set to TRUE if using doMC
  trace.it = TRUE,
  nlambda = 100
)

head(cv.fit)

#### Extract Selected CpGs of your MPS ####

# Coefficients at lambda.min

pred.vals <- predict(cv.fit, newx = t(beta_tune), s = cv.fit$lambda.min)
colnames(pred.vals)[1] = "DNAm_Cash_score"

##merge with pheno file
pred.vals <- pred.vals %>% as.data.frame() %>% tibble::rownames_to_column("BaseID")

#merge
pheno_tune_withDNAm <- merge(pred.vals, pheno_tune, by= "BaseID")

colnames(pheno_tune_withDNAm)


#Run Regression
logit_model <- glm(treat ~ DNAm_Cash_score, data = pheno_tune_withDNAm, family = binomial())
summary(logit_model)
confint(logit_model)

#Run Regression (with standardized DNAm score)
logit_model <- glm(treat ~ scale(DNAm_Cash_score), data = pheno_tune_withDNAm, family = binomial())
summary(logit_model)
confint(logit_model)

#Check Mean differences
library(ggplot2)


#pheno_tune_withDNAm <- pheno.tune.withDNAm.mums  #load this data in via link below

ggplot(pheno_tune_withDNAm, aes(x = factor(treat), y = DNAm_Cash_score, fill = factor(treat))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Control", "1" = "Treatment")) +
  labs(x = "Treatment Group", y = "Methylation Profile Score (MPS)", fill = "Group") +
  theme_minimal()

t_test_result <- t.test(DNAm_Cash_score ~ factor(treat), data = pheno_tune_withDNAm, var.equal = TRUE)
print(t_test_result)

table(pheno_tune_withDNAm$treat)


# Plot ROC
roc_obj <- roc(pheno_tune_withDNAm$treat, pheno_tune_withDNAm$`DNAm_Cash_score`)
plot(roc_obj, col = "blue", main = paste0("AUC = ", round(auc(roc_obj), 3)))

## get weights of probes ###

# Extract the coefficients of the model at lambda.min
coefficients <- coef(cv.fit, s = "lambda.min")

# Keep only non-zero coefficients (i.e., the selected CpGs)
selected_probes <- coefficients[coefficients[,1] != 0,]

# Convert the selected probes into a data frame
enet_weights <- as.data.frame(selected_probes) %>%
  rownames_to_column("probe")

# View the first few rows of the selected probes
head(enet_weights)

#remove _TC11 suffixes as they are specific to epicv2
enet_weights$probe <- sub("_[^_]+$", "", enet_weights$probe)

#check how many of these selected probes were in the top hits of EWAS
#make sure there are no EpicV2 attachments
load("  ss.hits.mums.cov.bacon.rda")

enet_weights$probe_clean <- sub("_[^_]+$", "", enet_weights$probe)
ss.hits.mums.cov.bacon$probe_clean <- sub("_[^_]+$", "", rownames(ss.hits.mums.cov.bacon))

#Take top 17 hits EWAS
#make sure they are ordered
ss.hits.mums.cov.bacon <- ss.hits.mums.cov.bacon[order(ss.hits.mums.cov.bacon$P.Value),]

top17_CpGs <- head(ss.hits.mums.cov.bacon, 17)

overlap_EWAS_ML <- sum(top17_CpGs$probe_clean %in% enet_weights$probe_clean)
print(overlap_EWAS_ML) # N of top CpGs in ML


##### Data to save so it can be shared with others for replication ####
cv.fit.mums <- cv.fit
enet.weights.mums <- enet_weights
betas.top30000.mums <- as.data.frame(betas_top30000)
pheno.tune.withDNAm.mums <- pheno_tune_withDNAm

save(cv.fit.mums, file = "  cv.fit.mums.rda")
save(enet.weights.mums, file= "  enet.weights.mums.rda")
save(betas.top30000.mums, file = "  betas.top30000.mums.rda")
save(pheno.tune.withDNAm.mums, file = "  pheno.tune.withDNAm.mums.rda")

load("  enet.weights.mums.rda")
load("  cv.fit.mums.rda")

load("  pheno.tune.withDNAm.mums.rda")
colnames(pheno.tune.withDNAm.mums)

#Save participants who were part of tuning (in Manuscript testing set! is more accurate label!) so these can be used again as testing set at age 6
tuning_set_mums_elastic_net <- pheno.tune.withDNAm.mums["BaseID"]
save(tuning_set_mums_elastic_net, file = "  tuning_set_mums_elastic_net")

#### Sensitivity Analyses ####

#ch.19.748535F , ch.22.48274842R are two probes in the enet model, but they look different from other probe names
#check if they are in the betadata

targets <- c("ch.19.748535F", "ch.22.48274842R")
cols <- colnames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums)
matches <- sapply(targets, function(t) any(grepl(t, cols)))
names(matches) <- targets
print(matches)



targets %in% colnames(betas_forEwas_excldpl_ResCellTech_ResBaseline_mums)

betas_forEwas_excldpl_ResCellTech_ResBaseline_mums[1:5, 1:5]



### Clean weights ####

#create updated enet.weights 

enet.weights.mums.July <- enet.weights.mums[!enet.weights.mums$probe %in% c("ch.19.748535F", "ch.22.48274842R"), ]

#check duplicates
enet_cpgs <- enet.weights.mums.July$probe_clean
meanbeta_cpgs <- mean_beta_df_unique$CpG  

enet_cpgs <- enet_cpgs[enet_cpgs != "(Intercept)"]
missing_cpgs <- enet_cpgs[!enet_cpgs %in% meanbeta_cpgs]

print(missing_cpgs)

any(duplicated(enet.weights.mums.July$probe))
unique(enet.weights.mums.July$probe[duplicated(enet.weights.mums.July$probe)])



####### For Kids ####

#### Load in Beta and Phenotype data ####


#load phenotype date with treatment variable
load("  pheno_excldpl_kids_mergedBaseCov.rda")
dim(pheno_excldpl_kids_mergedBaseCov)
pheno_excldpl_kids_mergedBaseCov[1:5, 1:5]

#load beta data residualized for baseline covariates
load("  DNAm betas/betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids.rda")
dim(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids)
betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[1:5, 1:5]

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids), pheno_excldpl_kids_mergedBaseCov$BaseID) #should be True

#### Select CpGs that are both on EpicV1 and EpicV2 ####

#load in CpGs specific to Epic v1
cpg_sites_EpicV1 <- readRDS("    cpg_sites_EpicV1.rds")
cpg_sites <- cpg_sites_EpicV1

# Remove suffixes (_B11, _B12, etc.) from rownames in your beta file (cause those suffixes are specific to epicv2)
rownames_clean <- sub("_[^_]+$", "", rownames(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids))

# Create a logical vector for matching
rows_to_keep <- rownames_clean %in% cpg_sites$x

# Subset the dataframe
betas_filtered <- betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids[rows_to_keep, ]
dim(betas_forEwas_excldpl_ResCellTech_ResBaseline_Kids) #278771 CpGs
dim(betas_filtered) #218902 CpGs that are on Epic V1 and Epic V2

#select only the top 30.000 cpgs to go into the model
#first load in the top hits after bacon adjustment
ss_hits_kids_cov_bacon <- read_excel("  ss.hits.kids.cov.bacon.xlsx")

#Order hits data by p-value
ss_hits_ordered <- ss_hits_kids_cov_bacon[order(ss_hits_kids_cov_bacon$P.Value), ]

# Get the CpGs present in both datasets
common_cpgs <- intersect(rownames(betas_filtered), ss_hits_ordered$CpGs)

# Keep only the hits that are in the beta file
hits_in_betas <- ss_hits_ordered[ss_hits_ordered$CpGs %in% common_cpgs, ]

#Select the top 30,000 CpGs from those that match
ss_top30000_hits <- hits_in_betas[1:30000, ]

# Subset the beta file to these top 30,000 CpGs
betas_top30000 <- betas_filtered[rownames(betas_filtered) %in% ss_top30000_hits$CpGs, ]
dim(betas_top30000)

#check if a couple of tophits from EWAS are in these top 30.000
"cg06152865_BC21" %in% rownames(betas_top30000)
"cg14480531_TC21" %in% rownames(betas_top30000) #all True, meaning top EWAS hits on epicv1/v2 are in the betafile
"cg08861115_TC21" %in% rownames(betas_top30000)
"cg17280975_BC21" %in% rownames(betas_top30000)
"cg11588197_BC21" %in% rownames(betas_top30000)
"cg06470558_BC21" %in% rownames(betas_top30000)


#while epicV2 specific hits should not be in there (cause we only select CpGs present at epicV1 and epicV2)
"cg24450643_BC21" %in% rownames(betas_top30000) #these yield False, cause indeed these are specific to epic v2 and we dont want them in our beta data
"cg01854039_BC21" %in% rownames(betas_top30000)
"cg01854039_BC21" %in% rownames(betas_top30000)




#### Perform stratified Random Sampling ####

#make sure its set as a factor
table(pheno_excldpl_kids_mergedBaseCov$treat) #0= control, 1=treatment
pheno_excldpl_kids_mergedBaseCov$treat <- as.factor(pheno_excldpl_kids_mergedBaseCov$treat)

#set your pheno and beta data
pheno = pheno_excldpl_kids_mergedBaseCov
beta = betas_top30000
dim(beta)

# Create stratified partition
train_index <- createDataPartition(pheno$treat, p = 0.75, list = FALSE) #change this to actual treatment variable 

# Subset the data into training and tune sets
pheno_train <- pheno[train_index, ]
pheno_tune <- pheno[-train_index, ]

# Check distribution in the original dataset
prop.table(table(pheno$treat))

# Check distribution of treatment in the training set
prop.table(table(pheno_train$treat))

# Check distribution of treatment in the tune set
prop.table(table(pheno_tune$treat))

#create beta files for the training and tune data

# Subset the beta matrix for training samples
beta_train <- beta[, train_index]
dim(beta_train)
dim(beta)

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(beta_train), pheno_train$BaseID) #should be True

# Subset the beta matrix for tune samples
beta_tune <- beta[, -train_index]

# Check if the column names of beta_train and sample_id in pheno_train are identical and in the same order
identical(colnames(beta_tune), pheno_tune$BaseID)


##### Running Elastic Net ####

#### Assing 10-fold cross-tune ####
foldcount <- 10
set.seed(123)
foldAssignments <- sample(rep(1:foldcount, length.out = ncol(beta_train)))

#### Run Elastic Net model ####

cv.fit <- cv.glmnet(
  x = t(beta_train),
  y = pheno_train$treat,
  foldid = foldAssignments,
  family = "binomial",
  type.measure = "auc",   # or "class"
  alpha = 0.5,
  nfolds = foldcount,
  parallel = FALSE,       # set to TRUE if using doMC
  trace.it = TRUE,
  nlambda = 100
)

head(cv.fit)

#### Extract Selected CpGs of your MPS ####

# Coefficients at lambda.min

pred.vals <- predict(cv.fit, newx = t(beta_tune), s = cv.fit$lambda.min)
colnames(pred.vals)[1] = "DNAm_Cash_score"

##merge with pheno file
pred.vals <- pred.vals %>% as.data.frame() %>% tibble::rownames_to_column("BaseID")

#merge
pheno_tune_withDNAm <- merge(pred.vals, pheno_tune, by= "BaseID")

colnames(pheno_tune_withDNAm)


#Run Regression
logit_model <- glm(treat ~ DNAm_Cash_score, data = pheno_tune_withDNAm, family = binomial())
summary(logit_model)
confint(logit_model)

#Run Regression , standardized DNAm cash score
logit_model <- glm(treat ~ scale(DNAm_Cash_score), data = pheno_tune_withDNAm, family = binomial())
summary(logit_model)
confint(logit_model)

#Check Mean differences
library(ggplot2)

ggplot(pheno_tune_withDNAm, aes(x = factor(treat), y = DNAm_Cash_score, fill = factor(treat))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Control", "1" = "Treatment")) +
  labs(x = "Treatment Group", y = "Methylation Profile Score (MPS)", fill = "Group") +
  theme_minimal()

t_test_result <- t.test(DNAm_Cash_score ~ factor(treat), data = pheno_tune_withDNAm, var.equal = TRUE)
print(t_test_result)

table(pheno_tune_withDNAm$treat)
pheno_tune_withDNAm$treat <- as.factor(pheno_tune_withDNAm$treat)

# Plot ROC
roc_obj <- roc(pheno_tune_withDNAm$treat, pheno_tune_withDNAm$`DNAm_Cash_score`)
plot(roc_obj, col = "blue", main = paste0("AUC = ", round(auc(roc_obj), 3)))

## get weights of probes ###

# Extract the coefficients of the model at lambda.min
coefficients <- coef(cv.fit, s = "lambda.min")

# Keep only non-zero coefficients (i.e., the selected CpGs)
selected_probes <- coefficients[coefficients[,1] != 0,]

# Convert the selected probes into a data frame
enet_weights <- as.data.frame(selected_probes) %>%
  rownames_to_column("probe")

# View the first few rows of the selected probes
head(enet_weights)

#remove _TC11 suffixes as they are specific to epicv2
enet_weights$probe <- sub("_[^_]+$", "", enet_weights$probe)

#check how many of these selected probes were in the top hits of EWAS
#make sure there are no EpicV2 attachments
enet_weights$probe_clean <- sub("_[^_]+$", "", enet_weights$probe)
ss.hits.kids.cov.bacon$probe_clean <- sub("_[^_]+$", "", rownames(ss.hits.kids.cov.bacon))

#Take top 15 hits EWAS
#make sure they are ordered
ss.hits.kids.cov.bacon <- ss.hits.kids.cov.bacon[order(ss.hits.kids.cov.bacon$P.Value),]

top15_CpGs <- head(ss.hits.kids.cov.bacon, 15)

overlap_EWAS_ML <- sum(top15_CpGs$probe_clean %in% enet_weights$probe_clean)
print(overlap_EWAS_ML) # N of top CpGs in ML

##### Data to save so it can be shared with others for replication ####
cv.fit.kids <- cv.fit
enet.weights.kids <- enet_weights
betas.top30000.kids <- as.data.frame(betas_top30000)

save(cv.fit.kids, file = "  cv.fit.kids.rda")
save(enet.weights.kids, file = "  enet.weights.kids.rda")
save(betas.top30000.kids, file = "  betas.top30000.kids.rda")
save(betas.top30000.kids, file = "  FromSILOtoMPI/November/betas.top30000.kids.rda")
save(pheno_tune_withDNAm, file = "  pheno_tune_withDNAm.rda")


load("  cv.fit.kids.rda")
load("  enet.weights.kids.rda")
load("  betas.top30000.kids.rda")
load("  pheno_tune_withDNAm.rda")

#save tuning set (in manuscript referred to as testing set! that's more correct) so we can use these participants in follow-up again as testing set at age 6
tuning_set_kids_elastic_net <- pheno_tune_withDNAm["BaseID"]
save(tuning_set_kids_elastic_net, file = "  tuning_set_kids_elastic_net")

### CpG overlap EpiCash and Epigenetic aging probes ###################


#load in CpGs that are part of epigenetic aging probes. 

DunedinProbes <- read.csv("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/DunedinPACE_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
GrimAgeProbes  <- read.csv("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/GrimAgev1_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
PhenoAgeProbes  <- read.csv("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/PhenoAge_probeslist.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(GrimAgeProbes)[colnames(GrimAgeProbes) == "var"] <- "CpGmarker"

#load in CpGs included in EpiCash Moms
DNAm_CashGift_Comp_Mums_Excel <- read_excel("MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/DNAm.CashGift.Comp.Mums.Excel.xlsx")

#Select CpGs from EpiCash Moms
cash_cpgs <- DNAm_CashGift_Comp_Mums_Excel$probe             
cash_cpgs <- unique(na.omit(cash_cpgs))
cash_cpgs <- cash_cpgs[cash_cpgs != "(Intercept)"]

#Check overlap with epigenetic aging probes

dunedin <- unique(na.omit(DunedinProbes$CpGmarker))
grim    <- unique(na.omit(GrimAgeProbes$CpGmarker))
pheno   <- unique(na.omit(PhenoAgeProbes$CpGmarker))

ol_dunedin <- intersect(dunedin, cash_cpgs)
ol_grim    <- intersect(grim,    cash_cpgs)
ol_pheno   <- intersect(pheno,   cash_cpgs)

c(Dunedin = length(ol_dunedin),  
  GrimAge = length(ol_grim),
  PhenoAge = length(ol_pheno))

#see which CpGs overlap
ol_dunedin # 1
ol_grim    # 2
ol_pheno   # 1


#load in CpGs included in EpiCash Moms
DNAm_CashGift_Comp_Kids_Excel <- read_excel("MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/DNAm.CashGift.Comp.Kids.Excel.xlsx")

#Select CpGs from EpiCash Moms
cash_cpgs <- DNAm_CashGift_Comp_Kids_Excel$probe             
cash_cpgs <- unique(na.omit(cash_cpgs))
cash_cpgs <- cash_cpgs[cash_cpgs != "(Intercept)"]

#Check overlap with epigenetic aging probes

dunedin <- unique(na.omit(DunedinProbes$CpGmarker))
grim    <- unique(na.omit(GrimAgeProbes$CpGmarker))
pheno   <- unique(na.omit(PhenoAgeProbes$CpGmarker))

ol_dunedin <- intersect(dunedin, cash_cpgs)
ol_grim    <- intersect(grim,    cash_cpgs)
ol_pheno   <- intersect(pheno,   cash_cpgs)

c(Dunedin = length(ol_dunedin),  
  GrimAge = length(ol_grim),
  PhenoAge = length(ol_pheno))

#see which CpGs overlap
ol_dunedin # 0
ol_grim    # 0
ol_pheno   # 0



#### EWAS + Elastic net Bio Annotation and Enrichment Analysis (GO & KEGG)

# Install required packages (run once if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # Installs Bioconductor manager
}
BiocManager::install("missMethyl")  # For gene ontology and pathway enrichment analysis
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")  # EPICv2 array annotation data
install.packages("readxl")          # For reading Excel files
BiocManager::install("org.Hs.eg.db") # For mapping gene IDs (Entrez, SYMBOL)
install.packages("clusterProfiler")  # For additional enrichment and pathway analysis tools

# Load libraries
library(missMethyl)                         # GO/KEGG enrichment for methylation data
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)  # EPICv2 CpG annotation
library(readxl)                             # Read Excel sheet input data
library(org.Hs.eg.db)                       # Map gene identifiers
library(clusterProfiler)                    # More advanced pathway/GO analysis (optional)





# ------------------------------- EWAS hits annotation & enrichmnet -------------------------------------

### Load list of CpGs ####

#####Unadjusted CpGs (moms + Kids)####

#KIDS

Kids_Unadjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/EWAS_Hits/ss.hits.kids.cov.xlsx")  

#Mothers
Moms_Unadjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/EWAS_Hits/ss.hits.mums.cov.xlsx")  


#####Adjusted CpGs (moms + Kids)####

#KIDS

Kids_Adjusted_CpGs <- read_excel("s~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/EWAS_Hits/s.hits.kids.cov.bacon.xlsx")  

#Mothers

Moms_Adjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/EWAS_Hits/ss.hits.mums.cov.bacon.xlsx")  


# Load EWAS results

### Analyses for kids ####
#####Unadjusted significant CpGs #####

##MSD (148)

MSD <- 6.157635e-05 

Kids_MSD_survived <- subset(Kids_Unadjusted_CpGs, P.Value <= MSD)  

write_xlsx(Kids_MSD_survived, "Kids_MSD_survived_unadjusted.xlsx")



#####Bacon adjusted significant CpGs #####

##MSD (15)

MSD <- 6.157635e-05 

Kids_MSD_survived_adjusted <- subset(Kids_Adjusted_CpGs, P.Value <= MSD)  

write_xlsx(Kids_MSD_survived_adjusted, "Kids_MSD_survived_adjusted.xlsx")



# Kids_Adjusted_CpGs: all tested CpGs (background, n = 278,771)
# Kids_MSD_survived: significant CpGs, 148 top hits (from unadjusted p-value list)


# Prepare input data
sigCpGs <- Kids_MSD_survived$CpGs       # Significant CpGs (148 hits)
allCpGs <- Kids_Adjusted_CpGs$CpGs      # All tested CpGs (278,771 total)

#GO MSD 148 unadjusted hits -------------------------
# GO Enrichment Analysis
#-------------------------
go_results_KIDS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.

# Explanation:
# The bias plot displayed by gometh shows the relationship between the number of CpGs per gene 
# and the probability that a gene is called significant. If the plot slopes upward, it means 
# genes with more CpGs are more likely to be called significant by chance (probe-number bias).
# gometh corrects for this bias using the Wallenius non-central hypergeometric test, ensuring
# that enrichment results are not inflated by technical design of the array.

# Extract significant GO terms (FDR < 0.05)
sig_go_KIDS <- go_results_KIDS[go_results_KIDS$FDR < 0.05, ]
# Result: No GO terms were significant after multiple testing correction (FDR < 0.05).
# Most P.DE (raw p-values) and all FDRs were 1, indicating no statistical evidence for 
# enrichment of any GO term among your significant CpGs. This means your CpGs did not 
# cluster in any known biological processes more than expected by chance.

#KEGG MSD 148 unadjusted hits -------------------------
# KEGG Pathway Analysis
#-------------------------
kegg_results_KIDS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_KIDS <- kegg_results_KIDS[kegg_results_KIDS$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).
# As with GO, most P.DE values and all FDRs were 1, indicating no evidence for pathway-level 
# enrichment among your significant CpGs.

#-------------------------
# Interpretation and Reporting
#-------------------------
# - The absence of significant GO terms or KEGG pathways (all FDR = 1) means that your 
#   differentially methylated CpGs did not show functional clustering in any known biological 
#   processes or pathways, as defined by these databases.
# - This can happen if the input list is small, the CpGs are dispersed across many genes, 
#   or there is genuinely no pathway-level signal in your data.
# - The results are valid and reflect the current data and thresholds. It is standard to 
#   report null results in enrichment analyses.
# - If you wish to increase power, consider relaxing your significance threshold to include 
#   more CpGs, but always interpret such findings with caution.

write.csv(go_results_KIDS, "GO_enrichment_results_Kids_148.csv", row.names = FALSE)
write.csv(kegg_results_KIDS, "KEGG_enrichment_resultsKids_148.csv", row.names = FALSE)


# Prepare input data
sigCpGs <- Kids_MSD_survived_adjusted$CpGs       # Significant CpGs (15 hits)
allCpGs <- Kids_Adjusted_CpGs$CpGs      # All tested CpGs (278,771 total)

#GO MSD 15 adjusted hits -------------------------
# GO Enrichment Analysis
#-------------------------
go_results_KIDS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_KIDS <- go_results_KIDS[go_results_KIDS$FDR < 0.05, ] #0


#KEGG MSD 15 adjusted hits -------------------------
# KEGG Pathway Analysis
#-------------------------
kegg_results_KIDS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_KIDS <- kegg_results_KIDS[kegg_results_KIDS$FDR < 0.05, ] #0



##### Manhattan plot with adjusted p.values in KIDS #####
library(stringr)
library(qqman)
# Create annotation data where CpGs are linked to location on the genome  
anno2 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations  
anno2 <- data.frame(anno2) 

# Clean up chromosome column  
anno2$chr <- str_split(anno2$chr, "chr", n = 2, simplify = TRUE)[, 2]  
# Keep all chromosomes (1-22, X, Y) 
anno2 <- anno2[anno2$chr %in% c(1:22, "X", "Y"), ]  # Include X and Y chromosomes  

# Create a mapping function for chromosomes  
map_chr_to_num <- function(chr) {  
  if (chr %in% as.character(1:22)) {  
    return(as.numeric(chr))  
  } else if (chr == "X") {  
    return(23)  
  } else if (chr == "Y") {  
    return(24)  
  } else {  
    return(NA)  # For any other chromosomes (like MT)  
  }  
}  


# Apply the mapping to convert chromosome identifiers to numeric  
anno2$chr <- sapply(anno2$chr, map_chr_to_num)  

#anno2$chr <- as.numeric(anno2$chr)  # Convert to numeric 

#anno2 includes 930075 CpGs across all EpicV2
# "ss.hits.kids" is our EPIC v2 data hits for kids

# check the data before merging with the annotation file, same row names (CpG ids)
Kids_Adjusted_CpGs <- as.data.frame(Kids_Adjusted_CpGs)
row.names(Kids_Adjusted_CpGs) <- Kids_Adjusted_CpGs$CpGs 


# Merge with EPIC v2 data hits for kids  
# Kids_Unadjusted_CpGs contains P.Value  
forman_kids_Adjusted <- merge(anno2[c("chr", "pos")], Kids_Adjusted_CpGs[c("P.Value")], by = "row.names", all = TRUE)  
forman_kids_Adjusted <- forman_kids_Adjusted[complete.cases(forman_kids_Adjusted), ]  # Remove rows with NA values, 277041 CpGs
forman_kids_Adjusted$SNP <- ""  # Create an empty column for SNPs 


save(forman_kids_Adjusted , file = "forman_kids_adjusted_XY.rda")



chrom_colors <- c("darkblue", "lightblue4")
gw_thresh <- -log10(1.855701e-06)
suggestive_thresh <- -log10(6.157635e-05)

manhattan(
  forman_kids_Adjusted,
  chr = "chr",
  bp = "pos",
  p = "P.Value",
  col = chrom_colors,                     # Black and gray alternation by chromosome
  chrlabs = c(1:22, "X", "Y"),            # Label all chromosomes
  suggestiveline = FALSE,                 # Do not draw suggestive line automatically
  genomewideline = FALSE,                 # Do not draw genomewide line automatically
  main = "",
  cex.main = 2
)

# Genomewide threshold: dashed (long dashed, _ _ _)
abline(h = gw_thresh, col = "black", lty = 2, lwd = 3)             # Dashed line
# Suggestive threshold: solid (full line)
abline(h = suggestive_thresh, col = "black", lty = 1, lwd = 3)     # Solid line
# Optional: add light gray grid for easier reading
abline(h = seq(0, max(-log10(forman_kids_Adjusted$P.Value), na.rm = TRUE), by = 1),
       col = "lightgray", lty = "dotted")











### Analyses for Moms ####
#####Unadjusted significant CpGs #####

##MSD (61)

MSD <- 5.820722e-05

Moms_MSD_survived <- subset(Moms_Unadjusted_CpGs, P.Value <= MSD)  


#####Bacon adjusted significant CpGs #####

##MSD (17)

MSD <- 5.820722e-05

Moms_MSD_survived_adjusted <- subset(Moms_Adjusted_CpGs, P.Value <= MSD)  



#GO MSD 61 unadjusted hits -------------------------
# GO Enrichment Analysis

# Prepare input data
sigCpGs <- Moms_MSD_survived$CpGs       # Significant CpGs (61 hits)
allCpGs <- Moms_Adjusted_CpGs$CpGs      # All tested CpGs (282,371 total)

go_results_MOMS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot here:not enough data to observe a bias (a flat blue line)
)
# The analysis tested 21,762 GO terms for enrichment.

# Extract bin ranges from the plot data
plot_data <- attr(go_results_MOMS, "bias.plot")
high_bins <- plot_data[plot_data$Proportion > 0.004, ]  # Adjust threshold
print(high_bins) #0


# Extract significant GO terms (FDR < 0.05)
sig_go_MOMS <- go_results_MOMS[go_results_MOMS$FDR < 0.05, ]


#KEGG MSD 61 unadjusted hits -------------------------
# KEGG Pathway Analysis
#-------------------------
kegg_results_MOMS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_MOMS <- kegg_results_MOMS[kegg_results_MOMS$FDR < 0.05, ]


write.csv(go_results_KIDS, "GO_enrichment_results_Moms_61.csv", row.names = FALSE)
write.csv(kegg_results_KIDS, "KEGG_enrichment_results_Moms_61.csv", row.names = FALSE)



#GO MSD 17 adjusted hits -------------------------
# GO Enrichment Analysis

# Load EWAS results
# Moms_Adjusted_CpGs: all tested CpGs (background, n = 282,371)
# Moms_MSD_survived_adjusted: significant CpGs, 17 top hits (from unadjusted p-value list)

# Prepare input data
sigCpGs <- Moms_MSD_survived_adjusted$CpGs       # Significant CpGs (17 hits)
allCpGs <- Moms_Adjusted_CpGs$CpGs      # All tested CpGs (282,371 total)

go_results_MOMS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot here:not enough data to observe a bias (a flat blue line)
)
# The analysis tested 21,751 GO terms for enrichment.

# Extract bin ranges from the plot data
plot_data <- attr(go_results_MOMS, "bias.plot")
high_bins <- plot_data[plot_data$Proportion > 0.004, ]  # Adjust threshold
print(high_bins) #0


# Extract significant GO terms (FDR < 0.05)
sig_go_MOMS <- go_results_MOMS[go_results_MOMS$FDR < 0.05, ]

#KEGG MSD 17 adjusted hits -------------------------
# KEGG Pathway Analysis
#-------------------------
kegg_results_MOMS <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 366 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_MOMS <- kegg_results_MOMS[kegg_results_MOMS$FDR < 0.05, ]


write.csv(go_results_KIDS, "GO_enrichment_results_Moms_17.csv", row.names = FALSE)
write.csv(kegg_results_KIDS, "KEGG_enrichment_results_Moms_17.csv", row.names = FALSE)



##### Manhattan plot with adjusted p.values in MOMS #####

# Create annotation data where CpGs are linked to location on the genome  
anno2 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations  
anno2 <- data.frame(anno2) 

# Clean up chromosome column  
anno2$chr <- str_split(anno2$chr, "chr", n = 2, simplify = TRUE)[, 2]  
# Keep all chromosomes (1-22, X, Y) 
anno2 <- anno2[anno2$chr %in% c(1:22, "X"), ]  # Include X and Y chromosomes  

# Create a mapping function for chromosomes  
map_chr_to_num <- function(chr) {  
  if (chr %in% as.character(1:22)) {  
    return(as.numeric(chr))  
  } else if (chr == "X") {  
    return(23)  
  } else {  
    return(NA)  # For any other chromosomes (like MT)  
  }  
}  


# Apply the mapping to convert chromosome identifiers to numeric  
anno2$chr <- sapply(anno2$chr, map_chr_to_num)  

#anno2 includes 905565 CpGs across all EpicV2
# "Moms_Adjusted_CpGs" is our EPIC v2 data hits for kids

# check the data before merging with the annotation file, same row names (CpG ids)
Moms_Adjusted_CpGs <- as.data.frame(Moms_Adjusted_CpGs)
row.names(Moms_Adjusted_CpGs) <- Moms_Adjusted_CpGs$CpGs 


# Merge with EPIC v2 data hits for kids  
# Moms_Adjusted_CpGs contains P.Value  
forman_moms_adjusted <- merge(anno2[c("chr", "pos")], Moms_Adjusted_CpGs[c("P.Value")], by = "row.names", all = TRUE)  
forman_moms_adjusted <- forman_moms_adjusted[complete.cases(forman_moms_adjusted), ]  # Remove rows with NA values, 267856 CpGs
forman_moms_adjusted$SNP <- ""  # Create an empty column for SNPs 


save(forman_moms_adjusted , file = "forman_moms_adjusted.rda")


chrom_colors <- c("darkblue", "lightblue4")
gw_thresh <- -log10(1.761494e-06)
suggestive_thresh <- -log10(5.820722e-05)

manhattan(
  forman_moms_adjusted,
  chr = "chr",
  bp = "pos",
  p = "P.Value",
  col = chrom_colors,                     # Black and gray alternation by chromosome
  chrlabs = c(1:22, "X"),            # Label all chromosomes
  suggestiveline = FALSE,                 # Do not draw suggestive line automatically
  genomewideline = FALSE,                 # Do not draw genomewide line automatically
  main = "",
  cex.main = 2
)

# Genomewide threshold: dashed (long dashed, _ _ _)
abline(h = gw_thresh, col = "black", lty = 2, lwd = 3)             # Dashed line
# Suggestive threshold: solid (full line)
abline(h = suggestive_thresh, col = "black", lty = 1, lwd = 3)     # Solid line
# Optional: add light gray grid for easier reading
abline(h = seq(0, max(-log10(forman_kids_Adjusted$P.Value), na.rm = TRUE), by = 1),
       col = "lightgray", lty = "dotted")




#Enrichment Analysis for Elastic net model CpGs using the whole QCed CpGs------------------
## 1.Using EPICv2 annotation ------

# EPICv2 annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


## Prepare input data for Moms --------

#Load the Elastic net CpGs
#Moms CpGs (749 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_mums_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_mums_EpicV2labels$MPS_mums_EpicV2labels       # 749 CpGs 


#Load all the CpGs
#all.cpg ->	All CpGs tested in the initial EWAS (after QC/filtering, before elastic net)
mums_Adjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/GO_KEGG_Enrichmnet/ss.hits.mums.cov.bacon.xlsx") # 
allCpGs <- mums_Adjusted_CpGs$CpGs      # All tested CpGs (282,371 total)


### GO Enrichment Analysis Moms--------

go_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Moms <- go_results_Moms[go_results_Moms$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Moms------

kegg_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Moms <- kegg_results_Moms[kegg_results_Moms$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Moms, "GO_enrichment_results_Moms_Elasticnet_Probes_Epicv2.csv", row.names = FALSE)
write.csv(kegg_results_Moms, "KEGG_enrichment_results_Moms_Elasticnet_Probes_Epicv2.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Moms[order(go_results_Moms$P.DE), ][1:20, ]
library(dplyr)
top20 <- go_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Moms[order(kegg_results_Moms$P.DE), ][1:20, ]
top20 <- kegg_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)



## Prepare input data for Kids --------

#Load the Elastic net CpGs
#Kids CpGs (721 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_kids_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_kids_EpicV2labels$cpg_withEpicV2       # 721 CpGs 


#Load all the CpGs
#all.cpg ->	All CpGs tested in the initial EWAS (after QC/filtering, before elastic net)
kids_Adjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/GO_KEGG_Enrichmnet/ss.hits.kids.cov.bacon.xlsx") # 
allCpGs <- kids_Adjusted_CpGs$CpGs      # All tested CpGs (278,771 total)


### GO Enrichment Analysis Kids--------

go_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Kids <- go_results_Kids[go_results_Kids$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Kids------

kegg_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Kids <- kegg_results_Kids[kegg_results_Kids$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Kids, "GO_enrichment_results_Kids_Elasticnet_Probes_Epicv2.csv", row.names = FALSE)
write.csv(kegg_results_Kids, "KEGG_enrichment_results_Kids_Elasticnet_Probes_Epicv2.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Kids[order(go_results_Kids$P.DE), ][1:20, ]
top20 <- go_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Kids[order(kegg_results_Kids$P.DE), ][1:20, ]
top20 <- kegg_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)


## 2.Using EPICv1 annotation ------

#This time we try with the epicv1 cpg ids
# EPICv1 annotation
anno1 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)



## Prepare input data for Moms --------

#Load the Elastic net CpGs
#Moms CpGs (749 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_mums_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_mums_EpicV2labels$probe_clean       # 749 CpGs, use the epicv1 matching cpg ids


#Load all the CpGs
#all.cpg ->	All CpGs tested in the initial EWAS (after QC/filtering, before elastic net)
mums_Adjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/GO_KEGG_Enrichmnet/ss.hits.mums.cov.bacon.xlsx") # 
#edit cpg ids to match the epicv1
mums_Adjusted_CpGs$CpGs <- sub("_.*", "", mums_Adjusted_CpGs$CpGs)

allCpGs <- mums_Adjusted_CpGs$CpGs      # All tested CpGs (282,371 total)


### GO Enrichment Analysis Moms--------

go_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC",  # Use EPIC (v1) annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Moms <- go_results_Moms[go_results_Moms$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Moms------

kegg_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Moms <- kegg_results_Moms[kegg_results_Moms$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05)

write.csv(go_results_Moms, "GO_enrichment_results_Moms_Elasticnet_Probes_Epicv1.csv", row.names = FALSE)
write.csv(kegg_results_Moms, "KEGG_enrichment_results_Moms_Elasticnet_Probes_Epicv1.csv", row.names = FALSE)


#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Moms[order(go_results_Moms$P.DE), ][1:20, ]
top20 <- go_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Moms[order(kegg_results_Moms$P.DE), ][1:20, ]
top20 <- kegg_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)








## Prepare input data for Kids --------

#Load the Elastic net CpGs
#Kids CpGs (721 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_kids_EpicV2labels.rda")
#(MPS_kids_EpicV2labels)
#edit cpg ids to match the epicv1
MPS_kids_EpicV2labels$CpGs <- sub("_.*", "", MPS_kids_EpicV2labels$cpg_withEpicV2) 

sigCpGs <- MPS_kids_EpicV2labels$CpGs       # 721 CpGs 


#Load all the CpGs
#all.cpg ->	All CpGs tested in the initial EWAS (after QC/filtering, before elastic net)
kids_Adjusted_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/GO_KEGG_Enrichmnet/ss.hits.kids.cov.bacon.xlsx")

#edit cpg ids to match the epicv1
kids_Adjusted_CpGs$CpGs <- sub("_.*", "", kids_Adjusted_CpGs$CpGs) 
allCpGs <- kids_Adjusted_CpGs$CpGs      # All tested CpGs (278,771 total)


### GO Enrichment Analysis Kids--------

go_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC",  # Use EPICv1 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Kids <- go_results_Kids[go_results_Kids$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Kids------

kegg_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Kids <- kegg_results_Kids[kegg_results_Kids$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Kids, "GO_enrichment_results_Kids_Elasticnet_Probes_Epicv1.csv", row.names = FALSE)
write.csv(kegg_results_Kids, "KEGG_enrichment_results_Kids_Elasticnet_Probes_Epicv1.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Kids[order(go_results_Kids$P.DE), ][1:20, ]
top20 <- go_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Kids[order(kegg_results_Kids$P.DE), ][1:20, ]
top20 <- kegg_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)


#Enrichment Analysis for Elastic net model CpGs using the top 30000 EWAS CpGs------------------
##1.Using EPICv2 annotation ------

# EPICv2 annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


## Prepare input data for Moms --------

#Load the Elastic net CpGs
#Moms CpGs (749 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_mums_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_mums_EpicV2labels$MPS_mums_EpicV2labels       # 749 CpGs 


#Load all the CpGs
#all.cpg ->	top 30000 EWAS CpGs
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/betas.top30000.mums/betas.top30000.mums.rda") # 

moms_top_cpgs <- data.frame(CpGs = rownames(betas.top30000.mums))

allCpGs <- moms_top_cpgs$CpGs      # All tested CpGs (30000 total)


### GO Enrichment Analysis Moms--------

go_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Moms <- go_results_Moms[go_results_Moms$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Moms------

kegg_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Moms <- kegg_results_Moms[kegg_results_Moms$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Moms, "GO_enrichment_results_Moms_EWAS_Elasticnet_Epicv2.csv", row.names = FALSE)
write.csv(kegg_results_Moms, "KEGG_enrichment_results_Moms_EWAS_Elasticnet_Epicv2.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Moms[order(go_results_Moms$P.DE), ][1:20, ]
library(dplyr)
top20 <- go_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Moms[order(kegg_results_Moms$P.DE), ][1:20, ]
top20 <- kegg_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)



## Prepare input data for Kids --------

#Load the Elastic net CpGs
#Kids CpGs (721 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_kids_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_kids_EpicV2labels$cpg_withEpicV2       # 721 CpGs 


#Load all the CpGs

load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/betas.top30000.mums/betas.top30000.kids.rda") # 

kids_top_cpgs <- data.frame(CpGs = rownames(betas.top30000.kids))

allCpGs <- kids_top_cpgs$CpGs      # All tested CpGs (30000 total)



### GO Enrichment Analysis Kids--------

go_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC_V2",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Kids <- go_results_Kids[go_results_Kids$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Kids------

kegg_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC_V2",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Kids <- kegg_results_Kids[kegg_results_Kids$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Kids, "GO_enrichment_results_Kids_EWAS_Elasticnet_Epicv2.csv", row.names = FALSE)
write.csv(kegg_results_Kids, "KEGG_enrichment_results_Kids_EWAS_Elasticnet_Epicv2.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Kids[order(go_results_Kids$P.DE), ][1:20, ]
top20 <- go_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Kids[order(kegg_results_Kids$P.DE), ][1:20, ]
top20 <- kegg_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)


##2.Using EPICv1 annotation ------

# EPICv2 annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


## Prepare input data for Moms --------

#Load the Elastic net CpGs
#Moms CpGs (749 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_mums_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
sigCpGs <- MPS_mums_EpicV2labels$probe_clean       # 749 CpGs 


#Load all the CpGs
#all.cpg ->	top 30000 EWAS CpGs
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/betas.top30000.mums/betas.top30000.mums.rda") # 

moms_top_cpgs <- data.frame(CpGs = rownames(betas.top30000.mums))
moms_top_cpgs$CpGs <- sub("_.*", "", moms_top_cpgs$CpGs)


allCpGs <- moms_top_cpgs$CpGs      # All tested CpGs (30000 total)


### GO Enrichment Analysis Moms--------

go_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC",  # Use EPICv2 annotation for correct mapping
  plot.bias = TRUE         # Show probe-number bias plot
)



# Extract significant GO terms (FDR < 0.05)
sig_go_Moms <- go_results_Moms[go_results_Moms$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Moms------

kegg_results_Moms <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Moms <- kegg_results_Moms[kegg_results_Moms$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Moms, "GO_enrichment_results_Moms_EWAS_Elasticnet_Epicv1.csv", row.names = FALSE)
write.csv(kegg_results_Moms, "KEGG_enrichment_results_Moms_EWAS_Elasticnet_Epicv1.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Moms[order(go_results_Moms$P.DE), ][1:20, ]
library(dplyr)
top20 <- go_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Moms[order(kegg_results_Moms$P.DE), ][1:20, ]
top20 <- kegg_results_Moms %>%
  arrange(P.DE) %>%
  slice(1:20)



## Prepare input data for Kids --------

#Load the Elastic net CpGs
#Kids CpGs (721 CpGs)
load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Elastic Net Computation kids and moms/MPS_kids_EpicV2labels.rda")
#(MPS_mums_EpicV2labels)
MPS_kids_EpicV2labels$CpGs <- sub("_.*", "", MPS_kids_EpicV2labels$cpg_withEpicV2)

sigCpGs <- MPS_kids_EpicV2labels$CpGs       # 721 CpGs 


#Load all the CpGs

load("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/betas.top30000.mums/betas.top30000.kids.rda") # 

kids_top_cpgs <- data.frame(CpGs = rownames(betas.top30000.kids))

kids_top_cpgs$CpGs <- sub("_.*", "", kids_top_cpgs$CpGs)


allCpGs <- kids_top_cpgs$CpGs      # All tested CpGs (30000 total)



### GO Enrichment Analysis Kids--------

go_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "GO",
  array.type = "EPIC",  # 
  plot.bias = TRUE         # Show probe-number bias plot
)
# The analysis tested 21,751 GO terms for enrichment.


# Extract significant GO terms (FDR < 0.05)
sig_go_Kids <- go_results_Kids[go_results_Kids$FDR < 0.05, ]
#0 obs.


### KEGG Pathway Analysis Kids------

kegg_results_Kids <- gometh(
  sig.cpg = sigCpGs,
  all.cpg = allCpGs,
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = TRUE
)
# The analysis tested 368 KEGG pathways for enrichment.

# Extract significant KEGG pathways (FDR < 0.05)
sig_kegg_Kids <- kegg_results_Kids[kegg_results_Kids$FDR < 0.05, ]
# Result: No KEGG pathways were significant after multiple testing correction (FDR < 0.05).

write.csv(go_results_Kids, "GO_enrichment_results_Kids_EWAS_Elasticnet_Epicv1.csv", row.names = FALSE)
write.csv(kegg_results_Kids, "KEGG_enrichment_results_Kids_EWAS_Elasticnet_Epicv1.csv", row.names = FALSE)

#top Gene Ontology (GO) terms with the smallest (most significant) raw p-values
top20 <- go_results_Kids[order(go_results_Kids$P.DE), ][1:20, ]
top20 <- go_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)

#the top pathways with the smallest raw p-values from testing your CpGs for overrepresentation in KEGG pathways
top20 <- kegg_results_Kids[order(kegg_results_Kids$P.DE), ][1:20, ]
top20 <- kegg_results_Kids %>%
  arrange(P.DE) %>%
  slice(1:20)




### Smoking CpGs (Van Dongen 2023) ####


#Loading list of CpGs from: https://doi.org/10.7554/eLife.83286
#Effects of smoking on genome-wide DNA methylation profiles: A study of discordant and concordant monozygotic twin pairs

#In pairs discordant for current smoking, 13 differentially methylated CpGs were 
#found between current smoking twins and their genetically identical co-twin who never smoked.

#All 13 CpGs have been previously associated with smoking in unrelated individuals 
#and data from monozygotic pairs discordant for former smoking indicated that 
#methylation patterns are to a large extent reversible upon smoking cessation.


Smoking_CpGs <- read_excel("~/MPIB-SRT/1001-BFY/private/data analysis/02_Elastic Net/Bioannotation/EWAS_Hits/elife-Smoking_CpGs_Jenny.xlsx")  



# Convert Smoking_CpGs to a dataframe  
Smoking_CpGs <- as.data.frame(Smoking_CpGs)  

# Set row names of the data frame to the values in the "IlmnID" column  
rownames(Smoking_CpGs) <- Smoking_CpGs$IlmnID  

# Remove the "IlmnID" column as it's now redundant  
Smoking_CpGs$IlmnID <- NULL  



#####Overlap between Kids CpGs and smoking---------

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



#####Overlap between Moms CpGs and smoking---------

# Split the rownames of "Moms_Adjusted_CpGs" at the underscore and keep only the part before the underscore  
moms_row_names <- sapply(strsplit(Moms_Adjusted_CpGs$CpGs, "_"), `[`, 1)

# Find the common row names between "moms_row_names" and the rownames of "Smoking_CpGs"  
common_rows_mom <- intersect(moms_row_names, rownames(Smoking_CpGs))  

# Print the number of common rows  
print(length(common_rows_mom))  #10

#"cg19089201" "cg22132788" "cg05575921" "cg01901332" "cg21188533" "cg21566642" "cg21161138" "cg01940273" "cg00336149" "cg09935388"

common_rows_indices_moms <- which(moms_row_names %in% common_rows_mom)

#"cg19089201" "cg22132788" "cg05575921" "cg01901332" "cg21188533" "cg21566642" "cg21161138" "cg01940273" "cg00336149" "cg09935388"
#457   9296  10095  96403 132403 154862 161971 238386 240804 250849



################ The END ################


