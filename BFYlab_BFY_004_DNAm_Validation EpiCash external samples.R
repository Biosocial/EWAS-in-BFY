#############################################################
#### Validation EpiCash in external samples.             ####
#############################################################

#Packages
library(psych)
library(dplyr)
library(readr) 

#for regressions
library(lme4)
library(lmerTest)


######################### TTP ##############################

#### 1. Calculate BFY_MPS in TTP ####

#read in elastic net results (computed on SILO)
load("  Elasticnet_Validation_TTP/DNAm.CashGift.Comp.Kids.rda")

enet.weights.kids <- DNAm.CashGift.Comp.Kids

#change column names so they make sense
names(enet.weights.kids)[2] <- "weights"

#read in beta data of TTP that we want to compute MPS-BFY
#these are betas from Texas Twin Project (TTP)
load("~/MPIB-SRT/998-TTP/private/data/002_TTP_processed_data/1_Methylation data_Abby_20240915/beta_noobcombat.qc.240912.rda")

combat.beta[1:4, 1:4]

# Assume:
# enet.weights.kids is a data frame with columns: 'probe' and 'weights'
# combat.beta is a matrix or data frame with CpGs as rownames and samples as columns

# 1. Extract the intercept (if present)
intercept <- 0
if ("(Intercept)" %in% enet.weights.kids$probe) {
  intercept <- enet.weights.kids$weights[enet.weights.kids$probe == "(Intercept)"]
}
# Example: Intercept = -1.468742

# 2. Remove the intercept row to keep only CpGs
enet.weights.filtered <- enet.weights.kids[enet.weights.kids$probe != "(Intercept)", ]

# 3. Get the CpGs present in both the weights and the beta matrix
model_cpgs <- enet.weights.filtered$probe
available_cpgs <- rownames(combat.beta)
common_cpgs <- intersect(model_cpgs, available_cpgs)

# 4. Report any missing CpGs
missing_cpgs <- setdiff(model_cpgs, available_cpgs)
cat("Missing CpGs:", length(missing_cpgs), "out of", length(model_cpgs),
    "(", round(100 * length(missing_cpgs) / length(model_cpgs), 2), "%)\n")

#Missing CpGs: 69 out of 721 ( 9.57 %) , which is less than 20% so we can compute the score

# 5. Subset and reorder combat.beta and weights to match only common CpGs
beta.sub <- combat.beta[common_cpgs, ]
weights.sub <- enet.weights.filtered$weights[match(common_cpgs, enet.weights.filtered$probe)]

# 5b. Coverage of model weight retained (useful to judge impact of missing CpGs)
total_weight_abs <- sum(abs(enet.weights.filtered$weights), na.rm = TRUE)
kept_weight_abs  <- sum(abs(weights.sub), na.rm = TRUE)
coverage_abs_pct <- 100 * kept_weight_abs / total_weight_abs

total_weight_signed <- sum(enet.weights.filtered$weights, na.rm = TRUE)
kept_weight_signed  <- sum(weights.sub, na.rm = TRUE)
coverage_signed_pct <- 100 * kept_weight_signed / total_weight_signed

cat(sprintf("Retained weight (|w|): %.2f%%\n", coverage_abs_pct))

# 6. Ensure beta.sub is a matrix (just in case it's a data.frame)
beta.sub <- as.matrix(beta.sub)

# 7. Compute the DNAm score per person (sample)
DNAm_BFY <- as.numeric(crossprod(weights.sub, beta.sub)) + intercept  # 1 x n_samples

# 8. Create the final data frame with sample IDs and scores
DNAm_score_df <- data.frame(
  sample_id = colnames(beta.sub),
  DNAm_BFY = DNAm_BFY,
  row.names = NULL
)

# 9. Preview
head(DNAm_score_df)
hist(DNAm_score_df$DNAm_BFY)
describe(DNAm_score_df$DNAm_BFY)


# vars    n   mean   sd median trimmed  mad    min   max range skew kurtosis  se
# X1    1 1836 -23.61 4.22 -23.52  -23.66 4.24 -35.33 -9.92 25.41 0.13     -0.2 0.1

DNAm_BFY_kids_TTP <- DNAm_score_df
DNAm_BFY_kids_TTP_July <- DNAm_BFY_kids_TTP 

# 10. Save data
save(DNAm_BFY_kids_TTP_July, file = "  Elasticnet_Validation_TTP/DNAm_BFY_kids_TTP_July.rda")

load("  Elasticnet_Validation_TTP/DNAm_BFY_kids_TTP_July.rda")
DNAm_BFY_kids_TTP <- DNAm_BFY_kids_TTP_July
describe(DNAm_BFY_kids_TTP$DNAm_BFY)
hist(DNAm_BFY_kids_TTP$DNAm_BFY)

#### 2.Merge data ####

#With DNAm BFY
load("~/MPIB-SRT/1001-BFY/private/data/4_Processed Data/DNAm_BFY_kids_TTP.rda")
colnames(DNAm_BFY_kids_TTP)
DNAm_BFY_kids_TTP[1:2, 1:2]

#With SES variables
load("~/MPIB-SRT/998-TTP/private/data/002_TTP_processed_data/SES for BFY/ses_ttp_250617.rda")
colnames(ses_ttp)
ses_ttp[1:5, 1:5]

#with age and sex
load('~/MPIB-SRT/998-TTP/private/data/002_TTP_processed_data/0_Epigenetic clocks_Abby_20240915/ttp_epigenetic_240915.rda')
colnames(ttp_epigenetic)

#check correlations across celltype
library(psych)
library(apaTables)

# select the variables
vars <- c("Epi", "IC", "Fib", "cellf5_1", "cellf5_2", "cellf5_3", "cellf5_4", "cellf5_5")

df <- ttp_epigenetic[, vars]

# produce APA style correlation table
apa.cor.table(df, filename = "cor_table.doc")  # saves a Word doc


ttp_sex <- ttp_epigenetic[,c("Basename",
                                "sex")]

#change sample_id to to Basename for merging
colnames(DNAm_BFY_kids_TTP)[colnames(DNAm_BFY_kids_TTP) == "sample_id"] <- "Basename"


ses_ttp_DNAm <- left_join(ses_ttp, DNAm_BFY_kids_TTP, by = "Basename")
colnames(ses_ttp_DNAm)

ses_ttp_DNAm <- left_join(ses_ttp_DNAm, ttp_sex, by = "Basename")
colnames(ses_ttp_DNAm)


# Create DNAm residualized for tech covariates and Epithelial cells (Epi)
cor.test(ses_ttp_DNAm$Epi, ses_ttp_DNAm$IC) #this is -1 so we only pick one to avoid multicollinearity

#for BFY
reg = lm (DNAm_BFY ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + Epi
            , data=ses_ttp_DNAm)

ses_ttp_DNAm$DNAm_BFY_res_Epi = residuals(reg) 
ses_ttp_DNAm$DNAm_BFY_res_Epi_Z = as.numeric(scale(ses_ttp_DNAm$DNAm_BFY_res_Epi))

hist(ses_ttp_DNAm$DNAm_BFY_res_Epi_Z)
describe(ses_ttp_DNAm$DNAm_BFY_res_Epi_Z)
#vars    n mean sd median trimmed  mad   min  max range skew kurtosis   se
#X1    1 1836    0  1  -0.06   -0.02 0.96 -3.05 4.08  7.13 0.24     0.31 0.02

#for DunedinPACE
reg = lm (DunedinPACE ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + Epi,
             data=ses_ttp_DNAm)

ses_ttp_DNAm$DunedinPACE_res_Epi = residuals(reg) 
ses_ttp_DNAm$DunedinPACE_res_Epi_Z = as.numeric(scale(ses_ttp_DNAm$DunedinPACE_res_Epi))

#for Pheno
reg = lm (PCPhenoAge ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + oragene_age +Epi
            , data=ses_ttp_DNAm)

ses_ttp_DNAm$PCPhenoAge_res_Accel_Epi = residuals(reg) 
ses_ttp_DNAm$PCPhenoAge_res_Accel_Epi_Z = as.numeric(scale(ses_ttp_DNAm$PCPhenoAge_res_Accel_Epi))

#For GrimAge
reg = lm (PCGrimAge ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + oragene_age + Epi
          , data=ses_ttp_DNAm)

ses_ttp_DNAm$PCGrimAge_res_Accel_Epi = residuals(reg) 
ses_ttp_DNAm$PCGrimAge_res_Accel_Epi_Z = as.numeric(scale(ses_ttp_DNAm$PCGrimAge_res_Accel_Epi))

#create DNAm residualized for tech covariates but not cell type

#for DNAm_BFY
reg = lm (DNAm_BFY ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide)
            , data=ses_ttp_DNAm)

ses_ttp_DNAm$DNAm_BFY_res = residuals(reg) 
ses_ttp_DNAm$DNAm_BFY_res_Z = as.numeric(scale(ses_ttp_DNAm$DNAm_BFY_res))

#for DunedinPace
reg = lm (DunedinPACE ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide)
          , data=ses_ttp_DNAm)

ses_ttp_DNAm$DunedinPACE_res = residuals(reg) 
ses_ttp_DNAm$DunedinPACE_res_Z = as.numeric(scale(ses_ttp_DNAm$DunedinPACE_res))

#for PhenoAge
reg = lm (PCPhenoAge ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + oragene_age
          , data=ses_ttp_DNAm)

ses_ttp_DNAm$PCPhenoAge_res_Accel = residuals(reg) 
ses_ttp_DNAm$PCPhenoAge_res_Z_Accel = as.numeric(scale(ses_ttp_DNAm$PCPhenoAge_res_Accel))

#for GrimAge
reg = lm (PCGrimAge ~ as.factor(methyl_batch) + as.factor(array) + as.factor(slide) + oragene_age
          , data=ses_ttp_DNAm)

ses_ttp_DNAm$PCGrimAge_res_Accel = residuals(reg) 
ses_ttp_DNAm$PCGrimAge_res_Z_Accel = as.numeric(scale(ses_ttp_DNAm$PCGrimAge_res_Accel))


save(ses_ttp_DNAm, file = "  Elasticnet_Validation_TTP/ses_ttp_DNAm.rda")
load("  Elasticnet_Validation_TTP/ses_ttp_DNAm.rda")


#Add latest pheno data 
load("  Elasticnet_Validation_TTP/batch3_dnam_pheno.rda")

colnames(batch3_dnam_pheno) # n= 1838

# Match "batch3_dnam_pheno" with "ses_ttp_DNAm" by the Basename

# Finding matching and different Basenames  
matching <- batch3_dnam_pheno %>%  
  inner_join(ses_ttp_DNAm, by = "Basename")  #1836

# Merging the dataframes to add specified columns: BMI, pubertal development score(pds), poverty (ipr), externalizing and internalizing  
ses_ttp_DNAm <- ses_ttp_DNAm %>%  
  left_join(batch3_dnam_pheno %>% select(Basename, event,Sample_Name, bmi, cbcl_int, cbcl_ext, pds, pds_cat),   
            by = "Basename")  





save(ses_ttp_DNAm, file = "  Elasticnet_Validation_TTP/ses_ttp_DNAm.rda")
load("  Elasticnet_Validation_TTP/ses_ttp_DNAm.rda")


#### Add Sample Name so we can check technical replicates

#save dataset for ICC
ses_ttp_DNAm_icc <- ses_ttp_DNAm %>%
  left_join(batch3_dnam_pheno %>% select(Basename, Sample_Name), by = "Basename")

colnames(ses_ttp_DNAm_icc)

#get only technical replicates

library(dplyr)
library(stringr)

ses_ttp_DNAm_icc_pairs <- ses_ttp_DNAm_icc %>%
  mutate(base_id = str_remove(Sample_Name, "_R$")) %>%   # strip "_R" if present
  group_by(base_id) %>%
  filter(n() > 1) %>%   # keep only groups with >1 (i.e., both replicate + original exist)
  ungroup()

ses_ttp_DNAm_icc_pairs <- ses_ttp_DNAm_icc_pairs[,c("Sample_Name",
                                                      "DNAm_BFY",
                                                      "DNAm_BFY_res",
                                                      "DNAm_BFY_res_Z",
                                                     "DNAm_BFY_res_Epi")]


save(ses_ttp_DNAm_icc_pairs, file = "  Elasticnet_Validation_TTP/ses_ttp_DNAm_icc_pairs.rda")
load("  Elasticnet_Validation_TTP/ses_ttp_DNAm_icc_pairs.rda")
                                                      



#save dataset where technical replicates are removed
save(ses_ttp_DNAm_noTechRep, file = "  Elasticnet_Validation_TTP/ses_ttp_DNAm_noTechRep.rda")

ses_ttp_DNAm_noTechRep <- ses_ttp_DNAm_icc %>%
  filter(!grepl("_R", Sample_Name))



#### 3. Association Tests ####

# Using mixed-effects regression (lmer) -> for accounting for non-independence due to family structure and study batch effects
# Accounting for covariates like age, sex, and random effects for family and study
# Testing sex*age interaction


load("  Elasticnet_Validation_TTP/ses_ttp_DNAm_noTechRep.rda")

ses_ttp_DNAm <- ses_ttp_DNAm_noTechRep
colnames(ses_ttp_DNAm)

#### Demographics ####
#Sex and Age

#age
model <- lmer(DNAm_BFY_res_Epi_Z ~  scale(oragene_age) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 

#sex
model <- lmer(DNAm_BFY_res_Epi_Z ~  as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 
table(ses_ttp_DNAm$sex) #0=female, 1=male
tapply(ses_ttp_DNAm$DNAm_BFY_res_Epi_Z, ses_ttp_DNAm$sex, mean, na.rm = TRUE)


#### Socioeconomic context ####

#####Parental SES 
#Regression DNAm_BFY residualized for celltype and tech cov
model <- lmer(DNAm_BFY_res_Epi_Z ~  par_ses_z + oragene_age + as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 

#plus age*sex interaction
model <- lmer(DNAm_BFY_res_Epi_Z ~  par_ses_z + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(DNAm_BFY_res_Epi_Z ~  par_ses_z + oragene_age + as.factor(sex) + oragene_age*par_ses_z + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(DNAm_BFY_res_Epi_Z ~  par_ses_z + oragene_age + as.factor(sex) + par_ses_z*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

#####Poverty 

#Poverty Category: ses_ttp_DNAm$ipr_category (with categories: "in poverty", "near poverty", "above poverty")
table(ses_ttp_DNAm$ipr_category)

model <- lmer(DNAm_BFY_res_Epi_Z ~ as.factor(ipr_category) + oragene_age + as.factor(sex) +
                (1|family_id) + (1|study_id), data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 

#plus age*sex interaction
model <- lmer(DNAm_BFY_res_Epi_Z ~ as.factor(ipr_category) + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) +
                (1|family_id) + (1|study_id), data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(DNAm_BFY_res_Epi_Z ~ as.factor(ipr_category) + oragene_age + as.factor(sex) + as.factor(ipr_category)*as.factor(sex) +
                (1|family_id) + (1|study_id), data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(DNAm_BFY_res_Epi_Z ~ as.factor(ipr_category) + oragene_age + as.factor(sex) + as.factor(ipr_category)*oragene_age +
                (1|family_id) + (1|study_id), data = ses_ttp_DNAm, REML=F)
summary(model)

#### Mental Health ####

##Internalizing  behavior 
model <- lmer( DNAm_BFY_res_Epi_Z ~ scale(cbcl_int)   + oragene_age + as.factor(sex) + 
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)
confint(model, level = 0.95, method = "profile") 

#with age and sex interaction
model <- lmer(  scale(cbcl_int) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) +
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)

model <- lmer(  scale(cbcl_int) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + oragene_age*DNAm_BFY_res_Epi_Z +
                  (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

model <- lmer(  scale(cbcl_int) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + as.factor(sex)*DNAm_BFY_res_Epi_Z +
                  (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

##Externalizing  behavior
model <- lmer(scale(cbcl_ext) ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) + 
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)
confint(model, level = 0.95, method = "profile") 

#with sex and age intraction
model <- lmer(scale(cbcl_ext) ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) +  oragene_age*as.factor(sex)+
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

model <- lmer(scale(cbcl_ext) ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) +  oragene_age*DNAm_BFY_res_Epi_Z+
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

model <- lmer(scale(cbcl_ext) ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) +  as.factor(sex)*DNAm_BFY_res_Epi_Z+
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)

#### Physical Health ####

###Pubertal development
describe(ses_ttp_DNAm$pds)

#pubertal timing (not just pubertal status), residualize pds for age within each sex (separately for girls and boys)
ses_ttp_DNAm$sex <- as.factor(ses_ttp_DNAm$sex)

ses_ttp_DNAm$pds_res_age <- NA

# Loop over the levels of sex
for(s in levels(ses_ttp_DNAm$sex)){
  # Get indices for this sex
  idx <- which(ses_ttp_DNAm$sex == s)
  # Subset without NAs in crucial variables
  temp_data <- ses_ttp_DNAm[idx, c("pds", "oragene_age")]
  complete_cases <- complete.cases(temp_data)
  # Only run regression on complete cases for this sex
  fit <- lm(pds ~ oragene_age, data = temp_data[complete_cases, ])
  # The position in the full data
  idx_full <- idx[complete_cases]
  # Get residuals for those in the regression (complete cases)
  ses_ttp_DNAm$pds_res_age[idx_full] <- residuals(fit)
}


model_timing <- lmer(scale(pds_res_age) ~ DNAm_BFY_res_Epi_Z+ as.factor(sex) + 
                       (1|family_id) + (1|study_id), data = ses_ttp_DNAm, REML=F)

summary(model_timing)
confint(model_timing, level = 0.95, method = "profile") 

#####BMI 
describe(ses_ttp_DNAm$bmi)

model <- lmer(scale(bmi) ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) + 
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)
confint(model, level = 0.95, method = "profile") 

#with age and sex interaction
model <- lmer(scale(bmi) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) +
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

model <- lmer(scale(bmi) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*as.factor(sex) +
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

model <- lmer(scale(bmi) ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*oragene_age +
                (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)

summary(model)


#### Childhood Disease ####
healthDat <- read.csv("~/MPIB-SRT/998-TTP/private/data/002_TTP_processed_data/SES for BFY/TTP_Dep_Alc_Health.csv")
colnames(healthDat)


ses_ttp_DNAm_health <- left_join(ses_ttp_DNAm, healthDat, by = c("study_id", "event"))
colnames(ses_ttp_DNAm_health)

table(ses_ttp_DNAm_health$Disease_Twins)

model <- glmer(
  Disease_Twins ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) + (1|study_id) +
    (1 | family_id),
  data = ses_ttp_DNAm_health,
  family = binomial(link = "logit")
)

summary(model)
exp(0.21903) #odds ratio
exp(confint(model, parm = "DNAm_BFY_res_Epi_Z", method = "Wald"))


#with age*sex interaction

model <- glmer(
  Disease_Twins ~ DNAm_BFY_res_Z + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) +
    (1 | family_id),
  data = ses_ttp_DNAm_health,
  family = binomial(link = "logit")
)

model <- glmer(
  Disease_Twins ~ DNAm_BFY_res_Z + oragene_age + as.factor(sex) + oragene_age*DNAm_BFY_res_Z +
    (1 | family_id),
  data = ses_ttp_DNAm_health,
  family = binomial(link = "logit")
)

model <- glmer(
  Disease_Twins ~ DNAm_BFY_res_Z + oragene_age + as.factor(sex) + as.factor(sex)*DNAm_BFY_res_Z +
    (1 | family_id),
  data = ses_ttp_DNAm_health,
  family = binomial(link = "logit")
)

summary(model)


#### DNA methylation measurs of biological aging #####
colnames(ses_ttp_DNAm)

dat <- ses_ttp_DNAm[, c(
  "Epi",
  "DNAm_BFY_res",
  "DunedinPACE_res",
  "PCGrimAge_res_Accel",
  "PCPhenoAge_res_Accel",
  "DNAm_BFY_res_Epi",
  "DunedinPACE_res_Epi",
  "PCGrimAge_res_Accel_Epi",
  "PCPhenoAge_res_Accel_Epi")]

# correlations and p-values
cor_mat <- cor(dat, use = "pairwise.complete.obs")
p_mat <- psych::corr.test(dat)$p

# function to format APA-style (2 decimals + sig stars)
format_cor <- function(r, p) {
  stars <- ifelse(p < .001, "***",
                  ifelse(p < .01, "**",
                         ifelse(p < .05, "*", "")))
  paste0(sprintf("%.2f", r), stars)
}

apa_cor <- matrix("", nrow = ncol(dat), ncol = ncol(dat))
for (i in 2:ncol(dat)) {
  for (j in 1:(i-1)) {
    apa_cor[i, j] <- format_cor(cor_mat[i,j], p_mat[i,j])
  }
}
diag(apa_cor) <- "—"   # APA leaves diagonal as dashes

colnames(apa_cor) <- rownames(apa_cor) <- colnames(dat)
apa_cor

apaTables::apa.cor.table(dat, filename = "  Elasticnet_Validation_TTP/cor_table.doc")


#DunedinPace
model <- lmer(DunedinPACE_res_Epi_Z ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 

#with age*sex interaction
model <- lmer( DunedinPACE_res_Epi_Z ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer( DunedinPACE_res_Epi_Z ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer( DunedinPACE_res_Epi_Z ~ DNAm_BFY_res_Epi_Z  + oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*oragene_age + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)


#GrimAge Accel
model <- lmer(scale(PCGrimAge_res_Accel_Epi) ~ scale(DNAm_BFY_res_Epi)+ oragene_age + as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 


#with sex*age interaction
model <- lmer(PCGrimAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z+ oragene_age + as.factor(sex) + oragene_age*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(PCGrimAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z+ oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(PCGrimAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z+ oragene_age + as.factor(sex) + DNAm_BFY_res_Epi_Z*oragene_age + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)


#PhenoAge Accel
model <- lmer(PCPhenoAge_res_Accel_Epi_Z  ~ DNAm_BFY_res_Epi_Z + oragene_age + as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 

#with age*sex interaction
model <- lmer(PCPhenoAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z    + oragene_age + as.factor(sex) + oragene_age*as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(PCPhenoAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z    + oragene_age + as.factor(sex) + oragene_age*DNAm_BFY_res_Epi_Z + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

model <- lmer(PCPhenoAge_res_Accel_Epi_Z ~ DNAm_BFY_res_Epi_Z    + oragene_age + as.factor(sex) + as.factor(sex)*DNAm_BFY_res_Epi_Z + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm, REML=F)
summary(model)

### Add cognition

CogDat <- read_csv("MPIB-SRT/998-TTP/private/data/002_TTP_processed_data/Cognition_Education/20230329_Data/3-Fraemke_2023 Data.csv")
colnames(CogDat)


## choose relevant variables ###
vars <- c("study_id","event",
          
          #cognition
          "fsiq4","fsiq2","vci_comp","pri_comp")
          
CogDat_sub <- subset(CogDat, 
             select = c(vars))

## Merge with data

ses_ttp_DNAm_final <- left_join(ses_ttp_DNAm, CogDat_sub, by = c("study_id", "event"))

model <- lmer(scale(fsiq4) ~ DNAm_BFY_res_Epi_Z    + oragene_age + as.factor(sex) + (1|family_id) + (1|study_id) , data = ses_ttp_DNAm_final, REML=F)
summary(model)
confint(model, level = 0.95, method = "profile") 



################# SOEP #####################

##### Load in Data ####

#load in betas from SOEP
load("  /DNAm_datasets/noob_normalized_betas.rda")

#Load in info to compute MPS-mums in SOEP
load("  /Elastic Net Computation kids and moms/DNAm.CashGift.Comp.Mums.rda")
enet.weights.mums <- DNAm.CashGift.Comp.Mums


# Check:
# enet.weights.mums is a data frame with columns: 'probe' and 'weights'
colnames(enet.weights.mums)
colnames(enet.weights.mums)[colnames(enet.weights.mums) == "Weights"] <- "weights"

# beta.n is a matrix or data frame with CpGs as rownames and samples as columns
beta.n[1:5, 1:5]

# 1. Extract the intercept 
intercept <- 0
if ("(Intercept)" %in% enet.weights.mums$probe) {
  intercept <- enet.weights.mums$weights[enet.weights.mums$probe == "(Intercept)"]
}
# Intercept = -1.604

# 2. Remove the intercept row to keep only CpGs
enet.weights.filtered <- enet.weights.mums[enet.weights.mums$probe != "(Intercept)", ]

# 3. Get the CpGs present in both the weights and the beta matrix
model_cpgs <- enet.weights.filtered$probe
available_cpgs <- rownames(beta.n)
common_cpgs <- intersect(model_cpgs, available_cpgs)

# 4. Report any missing CpGs (info only)
missing_cpgs <- setdiff(model_cpgs, available_cpgs)
cat("Missing CpGs:", length(missing_cpgs), "out of", length(model_cpgs), 
    "(", round(100 * length(missing_cpgs) / length(model_cpgs), 2), "%)\n")
#Missing CpGs: 1 out of 749 ( 0.13 %)

# 5. Subset and reorder beta.n and weights to the common CpGs only (no imputation)
beta.sub <- beta.n[common_cpgs, , drop = FALSE]
weights.sub <- enet.weights.filtered$weights[match(common_cpgs, enet.weights.filtered$probe)]

# 6. Ensure beta.sub is a matrix
beta.sub <- as.matrix(beta.sub)

# 7. Compute the DNAm score per person (sample)
DNAm_BFY <- as.numeric(crossprod(weights.sub, beta.sub)) + intercept  # 1 x n_samples

# 8. Create the final data frame with sample IDs and scores
DNAm_score_df <- data.frame(
  sample_id = colnames(beta.sub),
  DNAm_BFY = DNAm_BFY,
  row.names = NULL
)

# 9. Preview
head(DNAm_score_df)
hist(DNAm_score_df$DNAm_BFY)
describe(DNAm_score_df$DNAm_BFY)

DNAm_BFY_mums_SOEP <- DNAm_score_df


save(DNAm_BFY_mums_SOEP, file = "  DNAm_BFY_mums_SOEP.rda")
load("  DNAm_BFY_mums_SOEP.rda")

##### Merge with phenotypes ####
PhenoDat <- read_csv("FullDat_Var2.csv")
colnames(PhenoDat)

#create sample ID so you can merge the data
library(dplyr)

PhenoDat <- PhenoDat %>%
  mutate(sample_id = paste(Slide, Array, sep = "_"))

colnames(PhenoDat)
colnames(DNAm_BFY_mums_SOEP)

soep_DNAm_bfy <- left_join(PhenoDat, DNAm_BFY_mums_SOEP, by = "sample_id")
describe(soep_DNAm_bfy$DNAm_BFY)

#only select rows with no NA on DNAm_BFY
soep_DNAm_bfy <- soep_DNAm_bfy[!is.na(soep_DNAm_bfy$DNAm_BFY), ]



### Residualize for technical covariates ####
#we residualize for technical covariates only in DNAm bfy , others have already been residualized
#no need for cell composition cause those with off cell composition have been excluded

reg = lm (DNAm_BFY ~ as.factor(Array) + as.factor(Slide), data=soep_DNAm_bfy)

soep_DNAm_bfy$DNAm_BFY_res = residuals(reg) 
soep_DNAm_bfy$DNAm_BFY_res_Z = as.numeric(scale(soep_DNAm_bfy$DNAm_BFY_res))

save(soep_DNAm_bfy, file = "  soep_DNAm_bfy.rda")
load("  soep_DNAm_bfy.rda")

##### Run associations with variables of interest ####

### Demographics
#Age
model = lm(DNAm_BFY_res_Z ~ Age_std , data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

#Sex
model = lm(DNAm_BFY_res_Z ~ female , data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

### SES
model = lm(DNAm_BFY_res_Z ~ Zeduinc + Age_std + female, data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

#with age*sex interaction
model = lm(DNAm_BFY_res_Z ~ Zeduinc + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model) 

model = lm(DNAm_BFY_res_Z ~ Zeduinc + Age_std + female + Age_std*Zeduinc, data = soep_DNAm_bfy)
summary(model) 

model = lm(DNAm_BFY_res_Z ~ Zeduinc + Age_std + female + female*Zeduinc, data = soep_DNAm_bfy)
summary(model) 

#### Add parental SES
### Add SES parents
bioparen <- read_dta("bioparen.dta")

#vsbil – LevelOfEducationFather
#msbil – LevelOfEducationMother

table(bioparen$bioyear, bioparen$vsbil)

EaPar <- bioparen[,c("pid",
                     "vsbil",
                     "msbil")]

#recode so that high scores mean higher EA

EaPar <- EaPar %>%
  mutate(
    vsbil_rank = case_when(
      vsbil %in% c(-1, -2, -3, -4, -5, -6, 0, 5) ~ NA_real_,
      vsbil %in% c(6, 7) ~ 0,
      vsbil == 1 ~ 1,
      vsbil == 2 ~ 2,
      vsbil == 3 ~ 3,
      vsbil == 4 ~ 4,
      TRUE ~ NA_real_
    ),
    msbil_rank = case_when(
      msbil %in% c(-1, -2, -3, -4, -5, -6, 0, 5) ~ NA_real_,
      msbil %in% c(6, 7) ~ 0,
      msbil == 1 ~ 1,
      msbil == 2 ~ 2,
      msbil == 3 ~ 3,
      msbil == 4 ~ 4,
      TRUE ~ NA_real_
    ),
    parentalEA = rowMeans(across(c(vsbil_rank, msbil_rank)), na.rm = TRUE)
  )

soep_DNAm_bfy <- left_join(soep_DNAm_bfy, EaPar, by = "pid")

model = lm(DNAm_BFY_res_Z ~ scale(parentalEA) + Age_std + female, data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 


#### Health

#recreate variables like in paper Qiao
soep_DNAm_bfy$SRH_D <- ifelse(soep_DNAm_bfy$SRH_ >= 2, 1, 0)
soep_DNAm_bfy$CDB_D <- ifelse(soep_DNAm_bfy$CDB_ >= 2, 1, 0)
soep_DNAm_bfy$ChronicDis_D <- ifelse(soep_DNAm_bfy$ChronicDis >= 2, 1, 0)
soep_DNAm_bfy$LDA_day_D <- ifelse(soep_DNAm_bfy$LDA_day <= 2, 1, 0)

#self-reported health
#in table = self-reported disease
model = lm(scale(Unhealthy) ~ DNAm_BFY_res_Z + Age_std + female, data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

#for those with binary outcomes (health)

#multimorbidity
model_logit <- glm(
  ChronicDis_D ~ DNAm_BFY_res_Z + Age_std + female,
  data = soep_DNAm_bfy,
  family = binomial(link = "logit")
)
summary(model_logit)
confint(model_logit, level = 0.95, method = "profile") 
table(soep_DNAm_bfy$ChronicDis_D)
exp(coef(model_logit))
exp(confint(model_logit))


#in table = functional limitations daily living
model_logit <- glm(
  LDA_day_D ~ DNAm_BFY_res_Z + Age_std + female,
  data = soep_DNAm_bfy,
  family = binomial(link = "logit")
)
summary(model_logit)
confint(model_logit, level = 0.95, method = "profile") 
table(soep_DNAm_bfy$LDA_day_D)
exp(coef(model_logit))
exp(confint(model_logit))

#in table = BMI
model = lm(BMI_agesex_std ~ DNAm_BFY_res_Z   + Age_std + female, data = soep_DNAm_bfy)
summary(model)
confint(model, level = 0.95, method = "profile") 

#with age*sex interactions
model = lm(DNAm_BFY_res_Z ~ scale(Unhealthy) + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model) 

model = lm(DNAm_BFY_res_Z ~ scale(Unhealthy) + Age_std + female + Age_std*scale(Unhealthy), data = soep_DNAm_bfy)
summary(model)

model = lm(DNAm_BFY_res_Z ~ scale(Unhealthy) + Age_std + female + female*scale(Unhealthy), data = soep_DNAm_bfy)
summary(model)

model = lm(DNAm_BFY_res_Z ~ BMI_agesex_std + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model)

model = lm(DNAm_BFY_res_Z ~ BMI_agesex_std + Age_std + female + Age_std*BMI_agesex_std, data = soep_DNAm_bfy)
summary(model)

model = lm(DNAm_BFY_res_Z ~ BMI_agesex_std + Age_std + female + BMI_agesex_std*female, data = soep_DNAm_bfy)
summary(model)

model_logit <- glm(
  ChronicDis_D ~ DNAm_BFY_res_Z + Age_std + female + Age_std*female,
  data = soep_DNAm_bfy,
  family = binomial(link = "logit")
)
summary(model_logit)

model_logit <- glm(
  LDA_day_D ~ DNAm_BFY_res_Z + Age_std + female + Age_std*female,
  data = soep_DNAm_bfy,
  family = binomial(link = "logit")
)
summary(model_logit)


##### DNA methylation measures of biological aging 

#DunedinPACE
model = lm(scale(DunedinPACE_res) ~ DNAm_BFY_res_Z + Age_std + female, data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

#with age*sex interaction
model = lm( scale(DunedinPACE_res) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model) 

model = lm( scale(DunedinPACE_res) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*DNAm_BFY_res_Z, data = soep_DNAm_bfy)
summary(model) 

model = lm( scale(DunedinPACE_res) ~ DNAm_BFY_res_Z + Age_std + female + DNAm_BFY_res_Z*female, data = soep_DNAm_bfy)
summary(model) 


#PhenoAge
model = lm(scale(PCPhenoAge_res_age)~ DNAm_BFY_res_Z + Age_std + female, data = soep_DNAm_bfy)
summary(model) 

#with age*sex interaction
model = lm(scale(PCPhenoAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model) 

model = lm(scale(PCPhenoAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*DNAm_BFY_res_Z, data = soep_DNAm_bfy)
summary(model) 

model = lm(scale(PCPhenoAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + DNAm_BFY_res_Z*female, data = soep_DNAm_bfy)
summary(model) 


#GrimAge
model = lm( DNAm_BFY_res_Z ~ scale(PCGrimAge_res_age) + Age_std + female, data = soep_DNAm_bfy)
summary(model) 
confint(model, level = 0.95, method = "profile") 

#with age*sex interaction
model = lm(  scale(PCGrimAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*female, data = soep_DNAm_bfy)
summary(model) 

model = lm(  scale(PCGrimAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + Age_std*DNAm_BFY_res_Z, data = soep_DNAm_bfy)
summary(model) 

model = lm(  scale(PCGrimAge_res_age) ~ DNAm_BFY_res_Z + Age_std + female + DNAm_BFY_res_Z*female, data = soep_DNAm_bfy)
summary(model) 


##### Depression ####

#what about the depression item
#see https://paneldata.org/soep-is/datasets/inno/sim1201

inno <- read_dta("MPIB-SRT/997-SOEP/private/data/001_SOEP_rawdata/002_SOEP_Pheno/001_SOEP-IS/soep-is.2020_stata_en/inno.dta")


table(inno$sim1201, inno$syear)
table(inno$sim1202)

library(dplyr)
library(haven)

inno_dep <- inno %>%
  filter(
    syear == 2018,                  # keep only 2018
    !sim1201 %in% c(-1, -5),        # remove -1, -5 from sim1201
    !sim1204 %in% c(-1, -5),
    !sim1203 %in% c(-1, -5),
    !sim1202 %in% c(-1, -5)         # remove -1, -5 from sim1202
  ) %>%
  select(pid, sim1201, sim1202, sim1203, sim1204)     # keep only pid + variables of interest

soep_DNAm_bfy <- left_join(soep_DNAm_bfy, inno_dep, by = "pid")

sum(!is.na(soep_DNAm_bfy$DNAm_BFY_res) & 
      !is.na(soep_DNAm_bfy$sim1202))

soep_DNAm_bfy <- soep_DNAm_bfy %>%
  mutate(DepAnx = rowMeans(select(., sim1201:sim1204), na.rm = TRUE))


colnames(soep_DNAm_bfy)

#with depression/anxiety
model = lm(scale(DepAnx) ~ DNAm_BFY_res_Z + Age_std + female, data = soep_DNAm_bfy)
summary(model)
confint(model, level = 0.95, method = "profile")
describe(soep_DNAm_bfy$DepAnx)

soep_DNAm_bfy_final <- soep_DNAm_bfy

#loading in final dataset
save(soep_DNAm_bfy_final, file = "  soep_DNAm_bfy_final.rda")
load("  soep_DNAm_bfy_final.rda")

#if you use these do  soep_DNAm_bfy <- soep_DNAm_bfy_final


######### FFCW #########

##################################################
#### Goal: cross-check Elastic net model kids ####
##################################################

#read in elastic net results (computed on SILO)
load("C:/Users/treysm/Downloads/enet.weights.kids.rda")

#change column names so they make sense
names(enet.weights.kids)[2] <- "weights"

#read in beta data of sample you want to compute MPS-BFY
#these are betas from FFCW
#combat.beta <- readRDS("O:/PNG-DNAsecure/Jonah/FFCWData/450k/beta.rds")
combat.beta <- readRDS("O:/PNG-DNAsecure/Jonah/FFCWData/epic/beta.rds")
#combat.beta[1:4, 1:4]

#### create MPS in FFCW ####

# Assume:
# enet.weights.kids is a data frame with columns: 'probe' and 'weights'
# combat.beta is a matrix or data frame with CpGs as rownames and samples as columns

# 1. Extract the intercept (if present)
intercept <- 0
if ("(Intercept)" %in% enet.weights.kids$probe) {
  intercept <- enet.weights.kids$weights[enet.weights.kids$probe == "(Intercept)"]
}

# 2. Remove the intercept row to keep only CpGs
enet.weights.filtered <- enet.weights.kids[enet.weights.kids$probe != "(Intercept)", ]

# 3. Get the CpGs present in both the weights and the beta matrix
model_cpgs <- enet.weights.filtered$probe
available_cpgs <- rownames(combat.beta)
common_cpgs <- intersect(model_cpgs, available_cpgs)

# 4. Report any missing CpGs
missing_cpgs <- setdiff(model_cpgs, available_cpgs)
cat("Missing CpGs:", length(missing_cpgs), "out of", length(model_cpgs), 
    "(", round(100 * length(missing_cpgs) / length(model_cpgs), 2), "%)\n")

# Missing CpGs: 68 out of 721 ( 9.43 %)

# 5. Subset and reorder combat.beta and weights to match
beta.sub <- combat.beta[common_cpgs, ]
weights.sub <- enet.weights.filtered$weights[match(common_cpgs, enet.weights.filtered$probe)]

# 6. Ensure beta.sub is a matrix (just in case it's a data.frame)
beta.sub <- as.matrix(beta.sub)

# 7. Compute the DNAm score per person (sample)
DNAm_BFY <- as.numeric(crossprod(weights.sub, beta.sub)) + intercept  # 1 x n_samples

# 8. Create the final data frame with sample IDs and scores
DNAm_score_df <- data.frame(
  sample_id = colnames(beta.sub),
  DNAm_BFY = DNAm_BFY,
  row.names = NULL
)

# 9. Preview
head(DNAm_score_df)
hist(DNAm_score_df$DNAm_BFY)
summary(DNAm_score_df$DNAm_BFY)



# 10. Safe data
#save(DNAm_score_df, file = "BBFYscore_ff450k.rda")
save(DNAm_score_df, file = "BBFYscore_ffEPIC.rda")
load("C:/Users/treysm/Documents/BBFYscore_ffEPIC.rda")


##### Merge data ####



library(data.table)

pd <- fread("O:/PNG-DNAsecure/Jonah/FFCWData/epic/pd.csv")
dat <- merge(pd, DNAm_score_df, by.x = "methid", by.y = "sample_id")
load("C:/Users/treysm/Downloads/PGS and full.rdata")
pov <- PGS_and_full[, c("idnum", "cf1povca", "cf1inpov", "cm1inpov","cf2povca", "cf2povcob", "cf3povco", "cf3povcob")]
pov$idnum <- as.numeric(pov$idnum)
ses_ttp_DNAm<- merge(dat, pov, by = "idnum")

clocks <- fread("O:/PNG-DNAsecure/Methylation/FFCW/Clocks/allclocks_bothstudies_withcells.csv")
ses_ttp_DNAm<- merge(ses_ttp_DNAm, clocks, by.x = "methid", by.y = "MethID")

age9 <- subset(ses_ttp_DNAm, group == "k5" & replicate == "")
age15 <- subset(ses_ttp_DNAm, group == "k6" & replicate == "")

age9$cf1inpov <-  ifelse(age9$cf1inpov < 0 , NA, age9$cf1inpov)
age9$cf1povca <-  ifelse(age9$cf1povca < 0 , NA, age9$cf1povca)
age9$cm1inpov <-  ifelse(age9$cm1inpov < 0 , NA, age9$cm1inpov)
age9$cf2povca <-  ifelse(age9$cf2povca < 0 , NA, age9$cf2povca)
age9$cf2povcob <-  ifelse(age9$cf2povcob < 0 , NA, age9$cf2povcob)
age9$cf3povco <-  ifelse(age9$cf3povco < 0 , NA, age9$cf3povco)
age9$cf3povcob <-  ifelse(age9$cf3povcob < 0 , NA, age9$cf3povcob)

age15$cf1inpov <-  ifelse(age15$cf1inpov < 0 , NA, age15$cf1inpov)
age15$cf1povca <-  ifelse(age15$cf1povca < 0 , NA, age15$cf1povca)
age15$cm1inpov <-  ifelse(age15$cm1inpov < 0 , NA, age15$cm1inpov)
age15$cf2povca <-  ifelse(age15$cf2povca < 0 , NA, age15$cf2povca)
age15$cf2povcob <-  ifelse(age15$cf2povcob < 0 , NA, age15$cf2povcob)
age15$cf3povco <-  ifelse(age15$cf3povco < 0 , NA, age15$cf3povco)
age15$cf3povcob <-  ifelse(age15$cf3povcob < 0 , NA, age15$cf3povcob)

summary(DNAm_BFY)

##Poverty at birth (categorical)
model <- lm(DNAm_BFY ~  as.factor(cf1povca), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  as.factor(cf1povca), data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  as.factor(cf1povca) + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  as.factor(cf1povca) + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at birth (father)
model <- lm(DNAm_BFY ~  cf1inpov, data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf1inpov, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  cf1inpov + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf1inpov + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at birth (mother)
model <- lm(DNAm_BFY ~  cm1inpov, data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cm1inpov, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  cm1inpov + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cm1inpov + Age + as.factor(sex), data = age15)
summary(model)

##Poverty at Year 1 (categorical)
model <- lm(DNAm_BFY ~  as.factor(cf2povca), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  as.factor(cf2povca), data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  as.factor(cf2povca) + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  as.factor(cf2povca) + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at Year 1 (mother)
model <- lm(DNAm_BFY ~  cf2povcob, data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf2povcob, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  cf2povcob + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf2povcob + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at Year 3 (mother)
model <- lm(DNAm_BFY ~  cf3povcob, data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf3povcob, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  cf3povcob + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf3povcob + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty threshold at Year 3 (mother)
model <- lm(DNAm_BFY ~  cf3povco, data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf3povco, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY ~  cf3povco + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY ~  cf3povco + Age + as.factor(sex), data = age15)
summary(model)


#Clocks
model <- lm(DNAm_BFY ~  PhenoAge, data = age9)
summary(model)

model <- lm(DNAm_BFY ~  DNAmGrimAge , data = age9)
summary(model)

model <- lm(DNAm_BFY ~  PoAm45 , data = age9)
summary(model)

model <- lm(DNAm_BFY ~  PhenoAge, data = age15)
summary(model)

model <- lm(DNAm_BFY ~  DNAmGrimAge , data = age15)
summary(model)

model <- lm(DNAm_BFY ~  PoAm45 , data = age15)
summary(model)

#With Age and Sex

model <- lm(DNAm_BFY ~  PhenoAge + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY ~  DNAmGrimAge + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY ~  PoAm45 + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY ~  PhenoAge + Age + as.factor(sex) , data = age15)
summary(model)

model <- lm(DNAm_BFY ~  DNAmGrimAge + Age + as.factor(sex) , data = age15)
summary(model)

model <- lm(DNAm_BFY ~  PoAm45 + Age + as.factor(sex) , data = age15)
summary(model)


library(stringr)

#create DNAm residualized for tech cov
age9$slide <- str_split_i(age9$methid, "_", i = 1)
age9$array <- str_split_i(age9$methid, "_", i = 2)
age15$slide <- str_split_i(age15$methid, "_", i = 1)
age15$array <- str_split_i(age15$methid, "_", i = 2)

age9 <- lapply(age9, unlist)
age15 <- lapply(age15, unlist)

rg = lm (DNAm_BFY ~ as.factor(plate) + as.factor(array) + as.factor(slide) 
         , data=age9)

age9$DNAm_BFY_res_tech = residuals(rg) 
age9$DNAm_BFY_res_tech_Z = as.numeric(scale(age9$DNAm_BFY_res_tech))

hist(age9$DNAm_BFY_res_cell_Z)
summary(age9$DNAm_BFY_res_cell_Z)

rg1 = lm (DNAm_BFY ~ as.factor(plate) + as.factor(array) + as.factor(slide) 
          , data=age15)

age15$DNAm_BFY_res_tech = residuals(rg1) 
age15$DNAm_BFY_res_tech_Z = as.numeric(scale(age15$DNAm_BFY_res_tech))

#create DNAm residualized for tech cov and cell-type
# age9$slide <- str_split_i(age9$methid, "_", i = 1)
# age9$array <- str_split_i(age9$methid, "_", i = 2)
# age15$slide <- str_split_i(age15$methid, "_", i = 1)
# age15$array <- str_split_i(age15$methid, "_", i = 2)

age9 <- lapply(age9, unlist)
age15 <- lapply(age15, unlist)

reg = lm (DNAm_BFY ~ as.factor(plate) + as.factor(array) + as.factor(slide) +
            IC + Epi, data=age9)

age9$DNAm_BFY_res_cell = residuals(reg) 
age9$DNAm_BFY_res_cell_Z = as.numeric(scale(age9$DNAm_BFY_res_cell))

age9 <- lapply(age9, unlist)
age15 <- lapply(age15, unlist)

hist(age9$DNAm_BFY_res_cell_Z)
summary(age9$DNAm_BFY_res_cell_Z)

reg1 = lm (DNAm_BFY ~ as.factor(plate) + as.factor(array) + as.factor(slide) +
             IC + Epi, data=age15)

age15$DNAm_BFY_res_cell = residuals(reg1) 
age15$DNAm_BFY_res_cell_Z = as.numeric(scale(age15$DNAm_BFY_res_cell))

age9 <- lapply(age9, unlist)
age15 <- lapply(age15, unlist)

hist(age15$DNAm_BFY_res_cell_Z)
summary(age15$DNAm_BFY_res_cell_Z)
summary(age9$DNAm_BFY_res_tech_Z)
summary(age15$DNAm_BFY_res_tech_Z)
#run regressions

library(lme4)
install.packages("lmerTest")
library(lmerTest)

#SES
#Regression DNAm_BFY residualized for celltype and tech cov

##Poverty at birth (categorical)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf1povca), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf1povca), data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf1povca) + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf1povca) + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at birth (father)
model <- lm(DNAm_BFY_res_cell_Z ~  cf1inpov, data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf1inpov, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  cf1inpov + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf1inpov + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at birth (mother)
model <- lm(DNAm_BFY_res_cell_Z ~  cm1inpov, data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cm1inpov, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  cm1inpov + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cm1inpov + Age + as.factor(sex), data = age15)
summary(model)

##Poverty at Year 1 (categorical)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf2povca), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf2povca), data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf2povca) + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  as.factor(cf2povca) + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at Year 1 (mother)
model <- lm(DNAm_BFY_res_cell_Z ~  cf2povcob, data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf2povcob, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  cf2povcob + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf2povcob + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty ratio at Year 3 (mother)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povcob, data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povcob, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povcob + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povcob + Age + as.factor(sex), data = age15)
summary(model)

##Income to Poverty threshold at Year 3 (mother)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povco, data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povco, data = age15)
summary(model)
#With Age and Sex
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povco + Age + as.factor(sex), data = age9)
summary(model)
model <- lm(DNAm_BFY_res_cell_Z ~  cf3povco + Age + as.factor(sex), data = age15)
summary(model)

#Clocks
model <- lm(DNAm_BFY_res_cell_Z ~  PhenoAge, data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  DNAmGrimAge , data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PoAm45 , data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PhenoAge, data = age15)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  DNAmGrimAge , data = age15)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PoAm45 , data = age15)
summary(model)

#With Age and Sex

model <- lm(DNAm_BFY_res_cell_Z ~  PhenoAge + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  DNAmGrimAge + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PoAm45 + Age + as.factor(sex) , data = age9)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PhenoAge + Age + as.factor(sex) , data = age15)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  DNAmGrimAge + Age + as.factor(sex) , data = age15)
summary(model)

model <- lm(DNAm_BFY_res_cell_Z ~  PoAm45 + Age + as.factor(sex) , data = age15)
summary(model)

#### HRS #####
#see STATA file


####### PLOTS of association across independent cohorts ####

library(tidyverse)

# ---- Paste your data block exactly as-is between the triple quotes ----
txt <- "
Variable	    beta	  se	p_value	ci_low	ci_high	Cohort
Childhood socioeconomic disadvantage	-0,04	0,02	0,021	-0,08	-0,01	HRS
Childhood low educational attainment	-0,01	0,03	0,819	-0,07	0,06	SOEP
Adulthood socioeconomic disadvantage	-0,09	0,04	0,000	-0,13	-0,05	HRS
Adulthood socioeconomic disadvantage	-0,10	0,04	0,006	-0,16	-0,03	SOEP
Adulthood poverty	-0,21	0,06	0,000	-0,32	-0,09	HRS
Depressive symptoms	0,00	0,02	0,966	-0,04	0,05	HRS
Depressive symptoms	-0,14	0,09	0,112	-0,32	0,03	SOEP
Cognitive dysfunction	-0,07	0,02	0,001	-0,11	-0,03	HRS
Body mass index	-0,05	0,02	0,057	-0,09	0,00	HRS
Body mass index	-0,04	0,04	0,298	-0,12	0,04	SOEP
Multimorbidity	-0,04	0,02	0,038	-0,09	0,00	HRS
Functional limitations daily living	-0,02	0,02	0,334	-0,07	0,02	HRS
Functional limitations instrumental	-0,05	0,02	0,053	-0,09	0,00	HRS
Physiological dysregulation	-0,09	0,03	0,001	-0,14	-0,04	HRS
Self-reported disease	-0,07	0,04	0,078	-0,14	0,01	SOEP
DunedinPACE	-0,07	0,02	0,002	-0,12	-0,02	HRS
DunedinPACE	-0,36	0,03	0,000	-0,42	-0,30	SOEP
GrimAge Acceleration	-0,10	0,02	0,000	-0,14	-0,06	HRS
GrimAge Acceleration	-0,13	0,03	0,000	-0,19	-0,07	SOEP
PhenoAge Acceleration	0,06	0,02	0,019	0,01	0,11	HRS
PhenoAge Acceleration	-0,13	0,03	0,000	-0,19	-0,06	SOEP
Age	-0,01	0,00	0,000	-0,02	-0,01	HRS
Age	-0,14	0,03	0,000	-0,20	-0,08	SOEP
"

# ---- Read as tab-delimited with decimal commas ----
df <- read_delim(
  file = I(txt),
  delim = "\t",
  locale = locale(decimal_mark = ","),
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("HRS", "SOEP"))
  )

# Order variables nicely (keep your input order, top-to-bottom in plot)
var_order <- df %>%
  distinct(Variable) %>%
  pull(Variable) %>%
  rev()

df <- df %>%
  mutate(Variable = factor(Variable, levels = var_order))

# ---- Forest-plot style ----
ggplot(df, aes(x = beta, y = Variable, color = Cohort)) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 4.5      # increase this for bigger dots
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.6),
    height = 0.25,
    linewidth = 1.1   # increase this for thicker lines
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(min(df$ci_low), 0.20)) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.title  = element_text(size = 18, colour = "black"),
    legend.text  = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16, colour = "black")
  ) +
  labs(
    x = "Standardized regression estimate and 95% CI",
    y = "Variable",
    color = "Cohort"
  )


#### For Children ####

library(tidyverse)

# ---- Paste your data block exactly as-is between the triple quotes ----
txt <- "
Variable	Standardized Coefficients	SE	p-value	ci_low	ci_high	Cohort
Socioeconomic disadvantage	-0,024	0,027	0,383	-0,077	0,030	TTP
Socioeconomic disadvantage	0,024	0,037	0,521	-0,049	0,096	FFCW 9
Socioeconomic disadvantage	0,064	0,033	0,052	0,000	0,129	FFCW 15
Poverty	-0,047	0,137	0,732	-0,215	0,122	TTP
Poverty	-0,043	0,060	0,475	-0,161	0,075	FFCW 9
Poverty	-0,074	0,057	0,196	-0,186	0,037	FFCW 15
Internalizing	0,011	0,023	0,633	-0,050	0,060	TTP
Internalizing	0,031	0,033	0,343	-0,033	0,095	FFCW 9
Internalizing	-0,001	0,032	0,987	-0,063	0,062	FFCW 15
Externalizing	-0,032	0,026	0,227	-0,083	0,020	TTP
Externalizing	-0,005	0,033	0,872	-0,069	0,058	FFCW 9
Externalizing	0,000	0,030	0,999	-0,059	0,059	FFCW 15
Cognition	-0,014	0,024	0,554	-0,060	0,032	TTP
Cognition	-0,012	0,031	0,712	-0,073	0,050	FFCW 9
Pubertal Development	-0,041	0,028	0,135	-0,073	0,094	TTP
Pubertal Development	0,022	0,030	0,456	-0,036	0,080	FFCW 9
Pubertal Development	0,003	0,031	0,063	-0,057	0,063	FFCW 15
Body mass index	0,013	0,021	0,550	-0,029	0,054	TTP
Body mass index	0,017	0,031	0,591	-0,044	0,078	FFCW 9
Body mass index	-0,004	0,033	0,906	-0,068	0,060	FFCW 15
Childhood disease	-0,044	0,032	0,160	-0,106	0,018	FFCW 9
Childhood disease	-0,012	0,032	0,709	-0,075	0,051	FFCW 15
DunedinPACE	-0,005	0,025	0,850	-0,054	0,044	TTP
DunedinPACE	0,007	0,031	0,827	-0,054	0,068	FFCW 9
DunedinPACE	0,059	0,031	0,056	-0,001	0,120	FFCW 15
GrimAge Acceleration	-0,051	0,022	0,019	-0,094	-0,008	TTP
GrimAge Acceleration	-0,078	0,031	0,012	-0,139	-0,017	FFCW 9
GrimAge Acceleration	-0,036	0,032	0,266	-0,098	0,027	FFCW 15
PhenoAge Acceleration	0,109	0,022	0,000	0,065	0,152	TTP
PhenoAge Acceleration	0,066	0,031	0,033	0,005	0,126	FFCW 9
PhenoAge Acceleration	0,149	0,031	0,000	0,088	0,210	FFCW 15
Age	0,194	0,026	0,000	0,142	0,245	TTP
Age	0,040	0,0295	0,173	-0,018	0,098	FFCW 9
Age	0,001	0,0295	0,736	-0,057	0,059	FFCW 15
"

# ---- Read as tab-delimited with decimal commas ----
df <- read_delim(
  file = I(txt),
  delim = "\t",
  locale = locale(decimal_mark = ","),
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
  # rename the coefficient column to 'beta' to reuse plotting logic
  rename(
    beta    = `Standardized Coefficients`,
    se      = SE,
    p_value = `p-value`
  ) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("TTP", "FFCW 9", "FFCW 15"))
  )

# ---- Keep variables in input order (top-to-bottom), but flipped for plotting ----
var_order <- df %>%
  distinct(Variable) %>%
  pull(Variable) %>%
  rev()

df <- df %>%
  mutate(Variable = factor(Variable, levels = var_order))

# ---- Forest-plot style ----
ggplot(df, aes(x = beta, y = Variable, color = Cohort)) +
  geom_point(
    position = position_dodge(width = 0.7),
    size = 4.5
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.7),
    height = 0.25,
    linewidth = 1.1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "TTP"     = "lightgreen",  # dark green
      "FFCW 9"  = "purple",  # medium green
      "FFCW 15" = "blue"   # light green
    )
  ) + 
  scale_x_continuous(limits = c(-0.35, 0.30)) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.title  = element_text(size = 18, colour = "black"),
    legend.text  = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16, colour = "black")
  ) +
  labs(
    x = "Standardized regression estimate and 95% CI",
    y = "Variable",
    color = "Cohort"
  )


################################ The END ################################










