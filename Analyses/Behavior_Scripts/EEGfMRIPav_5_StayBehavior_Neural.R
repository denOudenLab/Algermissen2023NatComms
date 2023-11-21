#' EEGfMRIPav_5_StayBehavior_Neural
#' 
#' Perform mixed-effects linear and logistic regression testing how 
#' trial-by-trial BOLD amplitude and EEG power predict stay behavior
#' as reported behavioral results section and in Supplementary Material S15.
#' Interactive script, run step-by-step, NOT in one Go.
#' 
#' Requires that .csv files have been recoded and concatenated using 
#' EEGfMRIPav_1_Read_Initial.R
#' 
#' Calls EEGfMRIPav_2_Preprocess_Automated.R
#' Calls EEGfMRIPav_2_Cue_Position.R
#' Calls EEGfMRIPav_2_Stay.R
#' Calls EEGfMRIPav_0_CreateBarPlots.R
#'
#' Mind to adjust rootDir.
#'
#' INPUTS:
#' none
#'
#' OUTPUTS:
#' no outputs, fits (and saves) models.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2023.

# ========================================================================================================== #
#### Load packages: ####

library(ggplot2) # Version '3.3.4'
library(car) # Version '3.0.8'
library(lattice) # Version '0.20-44'
library(ltm) # Version '1.1-1'
library(pastecs) # Version '1.3-18'
library(plyr) # Version ‘1.8.6’
library(stringr) # Version '1.4.0'

# Analyses:
library(lme4) # Version ‘1.1.26’ # for lmer and glmer
library(afex) # Version ‘0.28.1’ # for mixed
library(emmeans) # Version ‘1.5.5.1’ # for lsmeans and glht
library(effects) # Version ‘4.2.0’ # for effects plots

options(scipen = 20)
options(contrasts = c("contr.sum", "contr.poly"))
rm(list = ls())

# ========================================================================================================== #
#### 1) Load full data set: ####

## Set directories:
rootDir   <- "/project/3017042.02/" # adjust root directory to local folder structure.
scriptDir <- paste0(rootDir, "Analyses/Behavior_Scripts/")
behavDir <- paste0(rootDir, "Log/Behavior/")
modDir <- paste0(behavDir, "Behavior_Models/")
plotDir <- paste0(behavDir, "Plots_R/")
fmriDir <- paste0(rootDir, "Log/fMRI/")

## Load functions:
source(paste0(scriptDir, "EEGfMRIPav_1_Read_Initial.R")) # Load variables and bring to R data frame format
source(paste0(scriptDir, "EEGfMRIPav_2_Preprocess_Automated.R")) # Recode variables
source(paste0(scriptDir, "EEGfMRIPav_2_Cue_Position.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_2_Stay.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_0_CreateBarPlots.R")) # # Retrieve cue position

## Apply functions:
mydata <- EEGfMRIPav_1_Read_Initial() # Initial reading, converting to factors
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.
mydata <- EEGfMRIPav_2_Cue_Position(mydata) # Compute cue position
# mydata$cue_pos
mydata <- EEGfMRIPav_2_Stay(mydata) # Compute stay behavior

## Inspect:
table(mydata$stay_n)
sum(is.na(mydata$stay_n)) # 576, i.e. 16 per subject, i.e. last trial for each cue
mean(mydata[is.na(mydata$stay_n), "cue_pos"] == 40) # 1, i.e. stay always missing for last cue repetition

# ========================================================================================================== #
#### 2) Load outcome-locked BOLD trial-by-trial data: ####

## Get number of subjects:
nSub <- length(unique(mydata$PPN_n))

## Loop and retrieve all ROIs
varNamesAll <- c("GLM1ACCConjMan", "GLM1PCCConj", "GLM1vmPFCConjMan", "GLM1StriatumConj")
for (iROI in 1:length(varNamesAll)){
  
  varName <- varNamesAll[iROI]
  
  ## Set directory for respective ROI:
  fmriDir <- paste0(rootDir, "Log/EEG/OutcomeLockedResults/TAfT_Betas/TAfT_BOLD_", varName, "/")
  
  ## Load:
  mydata[,varName] <- NA # initialize
  
  ## Loop over subject files:
  for (iSub in 1:nSub){ # iSub <- 1
    cat(paste0("Start subject ",iSub, "\n"))
    fileName <- paste0(fmriDir, "TAfT_BOLD_", varName, "_sub", str_pad(iSub, 3, side = "left", pad = 0), ".csv")
    tmp <- as.matrix(read.csv(fileName, header = F))
    if (length(tmp) != 640){stop(paste0("Subject ", iSub, ": Input not of length 640 trials"))}
    subIdx <- which(mydata$PPN_n == iSub) 
    mydata[subIdx, varName] <- tmp
  }
}

## Check:
names(mydata)

# ========================================================================================================== #
#### 3a) Mixed-effects model with trial-by-trial BOLD as IV: ####

## Select subjects:
invalidSubs <- c(15, 25) # failed registration
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)
length(unique(modData$PPN_n))

## Select focal predictor:
varName <- "GLM1ACCConjMan"
varName <- "GLM1PCCConj"
varName <- "GLM1vmPFCConjMan"
varName <- "GLM1StriatumConj"

## Assign and standardize:
modData$fMRIPred <- modData[, varName]
modData$fMRIPred_z <- as.numeric(scale(modData$fMRIPred))

# ----------------------------------------------------------------------------------- #
### Run model:

# Formula:
formula <- "stay_n ~ fMRIPred_z + (fMRIPred_z|PPN_f)"

## Fit model:
mod <- glmer(formula = formula, modData, family = binomial(), 
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

# Plot:
plot(effect("fMRIPred_z", mod))

## p-value with LRTs:
mod_LRT <- mixed(mod, modData, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

## ACC:
#             Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)  0.96135    0.10914   8.809 <0.0000000000000002 ***
# fMRIPred_z  -0.01883    0.01588  -1.186               0.236 
# --> not significant, if anything negative association

## PCC:
#             Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)  0.96167    0.10916   8.810 <0.0000000000000002 ***
# fMRIPred_z  -0.03550    0.01608  -2.207              0.0273 *  
# --> just significant, negative association

## vmPFC:
#             Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)  0.96268    0.10929   8.808 < 0.0000000000000002 ***
# fMRIPred_z  -0.07642    0.01742  -4.387            0.0000115 ***
# --> clear negative association

## Striatum:
#             Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)  0.96556    0.10988   8.788 < 0.0000000000000002 ***
# fMRIPred_z   0.06695    0.02433   2.752              0.00593 ** 
# --> positive association

# ------------------------------------------------------------- #
### LRTs:
## ACC:
#            Df  Chisq Chi Df Pr(>Chisq)
# fMRIPred_z  4 1.2939      1     0.2553

## PCC:
#            Df Chisq Chi Df Pr(>Chisq)  
# fMRIPred_z  4 3.421      1    0.06437 .

## vmPFC:
#            Df  Chisq Chi Df Pr(>Chisq)    
# fMRIPred_z  4 15.559      1 0.00007996 ***

## Striatum:
#            Df  Chisq Chi Df Pr(>Chisq)   
# fMRIPred_z  4 9.0512      1   0.002625 **

# ========================================================================================================== #
#### 3b) Plot p(Stay) ~ trial-by-trial BOLD amplitude: ####

## Select subjects:
invalidSubs <- c(15, 25) # failed registration
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)
length(unique(modData$PPN_n))

names(modData)

allVar <- c("GLM1StriatumConj", "GLM1vmPFCConjMan", "GLM1PCCConj", "GLM1ACCConjMan")
allCol <- c("#FFC000", "#0070C0", "03BEBE", "#FF0000")

## Loop:
for (iVar in 1:length(allVar)){ # iVar <- 1
  for (iPerc in c(3, 4, 5)){ # iPerc <- 3
    nPerc <- iPerc; selVar <- allVar[iVar]; percVar <- paste0(selVar, "_", nPerc, "Perc");
    modData <- create_percentiles(modData, nPerc = nPerc, selVar)
    custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                    xLab = "BOLD (percentiles)", yLab = "p(Stay)", selCol = allCol[iVar],
                    isPoint = F, yLim = c(0.65, 0.75), savePNG = T, saveEPS = F)
    # With dots:
    custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                    xLab = "BOLD (percentiles)", yLab = "p(Stay)", selCol = allCol[iVar],
                    isPoint = T, yLim = c(0.25, 0.95), savePNG = T, saveEPS = F)
      }
}

# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
#### 4) Load outcome-locked EEG trial-by-trial data: ####

nSub <- length(unique(mydata$PPN_n))

## Loop and retrieve all EEG variables:
allVarNames <- c("Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta",
                 "Preferred_CzFCzFz_beta", "Action_CzFCzFz_lowalpha")
for (iROI in 1:length(allVarNames)){
  
  varName <- allVarNames[iROI]
  
  ## Set appropriate EEG directory:
  eegDir <- paste0(rootDir, "Log/EEG/OutcomeLockedResults/TF_singleTrial/", varName, "/")
  
  ## Load:
  mydata[,varName] <- NA # initialize
  for (iSub in 1:nSub){ # iSub <- 1
    cat(paste0("Start subject ", iSub, "\n"))
    fileName <- paste0(eegDir, "TF_singleTrial_", varName, "_sub", str_pad(iSub, 3, side = "left", pad = 0), ".csv")
    tmp <- as.matrix(read.csv(fileName, header = F))
    if (length(tmp) != 640){stop(paste0("Subject ", iSub, ": Input not of length 640 trials"))}
    subIdx <- which(mydata$PPN_n == iSub) 
    mydata[subIdx, varName] <- tmp
  }
}

### Select data:
invalidSubs <- c(11, 12, 23, 30) ## People excluded after EEG
# invalidSubs <- c(11, 12, 15, 23, 25, 30) ## People excluded after BOLD and EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

### Log transform:
allVarNames <- c("Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta",
                 "Preferred_CzFCzFz_beta", "Action_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_beta", "ActionLong_CzFCzFz_broad")
for (iROI in 1:length(allVarNames)){
  varName <- allVarNames[iROI]
  newVarName <- paste0(varName, "_Log")
  modData[, newVarName] <- 10*log(modData[, varName] - min(modData[, varName]) + 0.1)
}

#Check:
names(modData)

# ========================================================================================================== #
### Prepare EEG variable:

## Select subjects for EEG: 
invalidSubs <- c(11, 12, 23, 30) # subjects excluded after EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)
length(unique(modData$PPN_n))

## Select variable:
varName <- "Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta"
# varName <- "Preferred_CzFCzFz_beta"
# varName <- "Action_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_beta"
# varName <- "ActionLong_CzFCzFz_broad"

modData$EEGPred <- modData[,varName] # put onto EEGPred variable
## Inspect:
densityplot(modData$EEGPred)
## Re-standardize:
modData$EEGPred_z <- as.numeric(scale(modData$EEGPred))
## Log-transform (conversion to decibel):
modData$EEGPredLog <- 10*log(modData$EEGPred - min(modData$EEGPred) + 0.1)
mean(is.na(modData$EEGPredLog))
densityplot(modData$EEGPredLog)
## Re-standardize:
modData$EEGPredLog_z <- as.numeric(scale(modData$EEGPredLog))
densityplot(modData$EEGPredLog_z)

# ========================================================================================================== #
#### 5a) Mixed-effects models with trial-by-trial EEG power as IV: ####

## Selects data:
invalidSubs <- c(11, 12, 23, 30) ## People excluded after EEG
# invalidSubs <- c(11, 12, 15, 23, 25, 26, 30) ## People excluded after BOLD and EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)
length(unique(modData$PPN_f))

## Select variable:
varName <- "Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta"
# varName <- "Preferred_CzFCzFz_beta"
# varName <- "Action_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_beta"
# varName <- "ActionLong_CzFCzFz_broad"

modData$EEGPred <- modData[, varName] # put onto EEGPred variable

## Inspect:
densityplot(modData$EEGPred)
## Re-standardize:
modData$EEGPred_z <- as.numeric(scale(modData$EEGPred))
## Log-transform (conversion to decibel):
modData$EEGPredLog <- 10*log(modData$EEGPred - min(modData$EEGPred) + 0.1)
mean(is.na(modData$EEGPredLog))
densityplot(modData$EEGPredLog)
## Re-standardize:
modData$EEGPredLog_z <- as.numeric(scale(modData$EEGPredLog))
densityplot(modData$EEGPredLog_z)

# ----------------------------------- #

## Formula:
formula <- "stay_n ~ EEGPredLog_z + (EEGPredLog_z|PPN_f)"

## Fit model:
mod <- glmer(formula = formula, modData, family = binomial(), 
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Plot:
plot(effect("EEGPredLog_z", mod))

## p-values via LRT:
mod_LRT <- mixed(mod, modData, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

### With lme4:
## Low alpha:
#             Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)  0.97123    0.11135   8.722 <0.0000000000000002 ***
# EEGPred_z    0.09179    0.06695   1.371                0.17 
# --> non-significant positive association

## Theta:
#              Estimate Std. Error z value            Pr(>|z|)    
# (Intercept)   1.06608    0.10760   9.908 <0.0000000000000002 ***
# EEGPredLog_z -0.09898    0.04719  -2.097               0.036 *  
# --> negative association?

## Beta:
#             Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)  1.06092    0.10705   9.910 < 0.0000000000000002 ***
# EEGPred_z    0.14527    0.04057   3.581             0.000343 ***
# --> significant positive association

# --------------------------------- #
### With LRT:
## Low alpha:
#           Df  Chisq Chi Df Pr(>Chisq)
# EEGPred_z  4 1.2771      1     0.2584
## Theta:
#              Df  Chisq Chi Df Pr(>Chisq)  
# EEGPredLog_z  4 4.1787      1    0.04094 *
## Beta:
#           Df  Chisq Chi Df Pr(>Chisq)    
# EEGPred_z  4 11.886      1  0.0005657 ***

# ========================================================================================================== #
#### 5b) Plot p(Stay) ~ EEG TF Power: ####

allVar <- c("Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta_Log",
            "Preferred_CzFCzFz_beta_Log", "Action_CzFCzFz_lowalpha_Log", "ActionLong_CzFCzFz_lowalpha_Log", "ActionLong_CzFCzFz_beta_Log", "ActionLong_CzFCzFz_broad")
allCol <- c("#0070C0", "#FFC000", "#FF0000", "#FF0000", "#FF0000", "#FF0000")

## Loop:
for (iVar in 1:length(allVar)){
  for (iPerc in c(3, 4, 5)){
    nPerc <- iPerc; selVar <- allVar[iVar]; percVar <- paste0(selVar, "_", nPerc, "Perc");
    modData <- create_percentiles(modData, nPerc = nPerc, selVar)
    custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                    xLab = "TF power (percentiles)", yLab = "p(Stay)", selCol = allCol[iVar],
                    isPoint = F, yLim = c(0.65, 0.75), savePNG = T, saveEPS = F)
    # With dots:
    custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                    xLab = "TF power (percentiles)", yLab = "p(Stay)", selCol = allCol[iVar],
                    isPoint = T, yLim = c(0.25, 0.95), savePNG = T, saveEPS = F)
  }
}

# END OF FILE.
