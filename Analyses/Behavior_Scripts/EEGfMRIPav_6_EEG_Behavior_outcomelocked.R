#' EEGfMRIPav_6_EEG_Behavior_outcomelocked
#' 
#' Fit mixed-effects models testing whether trial-by-trial EEG power changes over time 
#' and whether effects of behavioral variables interact with time
#' as reported in Supplementary Material S14.
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

## Analyses:
library(lme4) # Version ‘1.1.26’ # for lmer and glmer
library(afex) # Version ‘0.28.1’ # for mixed
library(emmeans) # Version ‘1.5.5.1’ # for lsmeans and glht
library(effects) # Version ‘4.2.0’ # for effects plots

options(scipen = 20)
options(contrasts = c("contr.sum", "contr.poly"))
rm(list = ls())

# ========================================================================================================== #
#### 1) Load full data set ####

## Set directories:
rootDir   <- "/project/3017042.02/"
scriptDir <- paste0(rootDir, "Analyses/Behavior_Scripts/")
modDir <- paste0(rootDir, "Log/Behavior/Behavior_Models/")
fmriDir <- paste0(rootDir, "Log/fMRI/")

## Load functions:
setwd(scriptDir)
source(paste0(scriptDir, "EEGfMRIPav_1_Read_Initial.R")) # Load variables and bring to R data frame format
source(paste0(scriptDir, "EEGfMRIPav_2_Preprocess_Automated.R")) # Recode variables
source(paste0(scriptDir, "EEGfMRIPav_2_Cue_Position.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_2_Stay.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_0_CreateBarPlots.R")) # # Retrieve cue position

## Apply functions:
# mydata <- EEGfMRIPav_1_Read_Initial() # Initial reading, converting to factors 
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.
mydata <- EEGfMRIPav_2_Cue_Position(mydata) # Compute cue position
mydata <- EEGfMRIPav_2_Stay(mydata) # Compute stay behavior

table(mydata$stay_n)
sum(is.na(mydata$stay_n)) # 576, i.e. 16 per subject, i.e. last trial for each cue
mean(mydata[is.na(mydata$stay_n), "cue_pos"] == 40) # 1, i.e. stay always missing for last cue repetition

# setwd(modDir) # where to store models

# ========================================================================================================== #
#### 2) Load outcome-locked EEG trial-by-trial data ####

nSub <- length(unique(mydata$PPN_n))
# nSub <- 36

## Name:
allVarNames <- c("Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta",
                 "Preferred_CzFCzFz_beta", "Action_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_beta", "ActionLong_CzFCzFz_broad")

for (varName in allVarNames){
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

# ========================================================================================================== #
#### 3) Prepare EEG variable and data set: ####

## Select subjects for EEG: 
invalidSubs <- c(11, 12, 23, 30) # subjects excluded in EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

## Select variable:
# varName <- "Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta"
# varName <- "Preferred_CzFCzFz_beta"

varName <- "Action_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_beta"
# varName <- "ActionLong_CzFCzFz_broad"

modData$EEGPred <- modData[, varName] # put onto EEGPred variable

## Inspect:
densityplot(modData$EEGPred)
plot(density(modData$EEGPred))

## Re-standardize:
modData$EEGPred_z <- as.numeric(scale(modData$EEGPred))

## Log-transform (conversion to decibel):
# negative numbers produces NaNs; very small numbers produce -Inf 
# --> first make all positive by subtracting mean; add 0.1 to avoid very small numbers?
modData$EEGPredLog <- 10*log(modData$EEGPred - min(modData$EEGPred) + 0.1)
mean(is.na(modData$EEGPredLog))
densityplot(modData$EEGPredLog)

## Re-standardize:
modData$EEGPredLog_z <- as.numeric(scale(modData$EEGPredLog))
densityplot(modData$EEGPredLog_z)

# ========================================================================================================== #
#### 4a) Mixed-effects model: EEG as function of action and block: #####

varName <- "Action_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_beta"
# varName <- "ActionLong_CzFCzFz_broad"

## Select variable:
modData$EEGPred <- modData[, varName] # put onto EEGPred variable
modData$EEGPred_z <- as.numeric(scale(modData$EEGPred))
modData$EEGPredLog <- 10 * log(modData$EEGPred - min(modData$EEGPred) + 0.1) # log-transform non-negative numbers
modData$EEGPredLog_z <- as.numeric(scale(modData$EEGPredLog)) # re-standardize

## Standardize block:
modData$block_z <- as.numeric(scale(modData$block_n))

formula <- "EEGPredLog_z ~ is_go_cor_f*block_z + (is_go_cor_f*block_z|PPN_f)"
mod <- lmer(formula = formula, modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## p-value with LRTs:
mod_LRT <- mixed(mod, modData, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
anova(mod_LRT)

## Action_CzFCzFz_lowalpha:
#                        Estimate Std. Error         df t value Pr(>|t|)  
# (Intercept)           -0.007955   0.048645  31.033705  -0.164   0.8712  
# is_go_cor_f1           0.034635   0.014592  31.425687   2.374   0.0239 *
# block_z                0.051753   0.019329  31.128919   2.677   0.0117 *
# is_go_cor_f1:block_z   0.003000   0.007520 118.289204   0.399   0.6907  
#                     Df  Chisq Chi Df Pr(>Chisq)   
# is_go_cor_f         14 5.3495      1   0.020728 * 
# block_z             14 6.6452      1   0.009943 **
# is_go_cor_f:block_z 14 0.1562      1   0.692699   

## ActionLong_CzFCzFz_lowalpha:
#                       Estimate Std. Error        df t value Pr(>|t|)  
# (Intercept)          -0.014487   0.044680 31.037203  -0.324   0.7479  
# is_go_cor_f1          0.068181   0.029744 31.129156   2.292   0.0288 *
# block_z               0.071962   0.029095 31.059237   2.473   0.0191 *
# is_go_cor_f1:block_z  0.009634   0.008825 50.366885   1.092   0.2802
#                     Df  Chisq Chi Df Pr(>Chisq)  
# is_go_cor_f         14 5.0091      1    0.02521 *
# block_z             14 5.7574      1    0.01642 *
# is_go_cor_f:block_z 14 1.1837      1    0.27661  

## ActionLong_CzFCzFz_beta:
#                        Estimate Std. Error         df t value Pr(>|t|)  
# (Intercept)           -0.018010   0.046078  30.999247  -0.391   0.6986  
# is_go_cor_f1           0.082701   0.031811  31.117960   2.600   0.0141 *
# block_z               -0.041771   0.020568  30.998929  -2.031   0.0509 .
#                     Df  Chisq Chi Df Pr(>Chisq)  
# is_go_cor_f         14 6.3014      1    0.01206 *
# block_z             14 4.0071      1    0.04531 *
# is_go_cor_f:block_z 14 0.0295      1    0.86366  

# ========================================================================================================== #
#### 4b) Action x block: custom_barplots2 ####

## Select subjects for EEG: 
invalidSubs <- c(11, 12, 23, 30) # subjects excluded in EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

# varName <- "StrongWeakBias_CzFCzFz_lowalpha"
# varName <- "Action_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_lowalpha"
# varName <- "ActionLong_CzFCzFz_beta"
# varName <- "ActionLong_CzFCzFz_broad"

plotDir <- "/project/3017042.02/Log/OutcomeLockedPaperPlots/"

selVarNames <- c("Action_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_beta")
for (iVar in 1:length(selVarNames)){
  selVar <- selVarNames[iVar]
  custom_barplot2(modData, xVar = "block_n", yVar = selVar, zVar = "is_go_cor_f", subVar = "PPN_f",
                xLab = "Block number", yLab = "TF power (dB)", zLab = "Performed action", selCol = c("red", "blue"),
                  isPoint = F, savePNG = T, saveEPS = F)
  # With dots:
  custom_barplot2(modData, xVar = "block_n", yVar = selVar, zVar = "is_go_cor_f", subVar = "PPN_f",
                  xLab = "Block number", yLab = "TF power (dB)", zLab = "Performed action", selCol = c("red", "blue"),
                  isPoint = T, savePNG = T, saveEPS = F)
}

# ========================================================================================================== #
#### 5a) Mixed-effects model: EEG as function of action and trial number: #####

# modData$trial_nr_n # from 1-320
modData$trial_nr_all_n <- (modData$session_n - 1) * 320 + modData$trial_nr_n # from 1-640
table(modData$trial_nr_all_n)  
modData$trial_nr_all_z <- as.numeric(scale(modData$trial_nr_all_n))

formula <- "EEGPredLog_z ~ is_go_cor_f*trial_nr_all_n + (is_go_cor_f*trial_nr_all_n|PPN_f)"
mod <- lmer(formula = formula, modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Action_CzFCzFz_lowalpha:
#                                Estimate  Std. Error          df t value Pr(>|t|)  
# (Intercept)                 -0.10589148  0.08248820 27.09565154  -1.284    0.210  
# is_go_cor_f1                 0.03071568  0.02035976 18.47531225   1.509    0.148  
# trial_nr_all_n               0.00030464  0.00012224 23.85639581   2.492    0.020 *
# is_go_cor_f1:trial_nr_all_n  0.00001344  0.00004252 63.93889570   0.316    0.753  

plot(effect("is_go_cor_f", mod))
plot(effect("trial_nr_all_n", mod))
plot(effect("is_go_cor_f:trial_nr_all_n", mod))

# ========================================================================================================== #
#### 5b) Plot with percentiles: #####

nPerc <- 5
varName <- "trial_nr_all_n"
newVarName <- paste0(varName, "_", nPerc, "Perc") # new variable name
modData$selVar <- modData[, varName]
modData[,newVarName] <- with(modData, cut(selVar, 
                                          breaks = quantile(selVar, probs = seq(0, 1, by = 1/nPerc), na.rm=TRUE), 
                                          include.lowest = TRUE))
modData[, newVarName] <- as.numeric(modData[, newVarName]) # into numeric

## Plot:
plotData <- modData
names(modData)
ggplot(plotData,aes(trial_nr_all_n_5Perc, StrongWeakBias_CzFCzFz_lowalpha, fill = is_go_cor_f)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.6) + 
  labs(x = "Trial number (percentiles)", fill = "Executed action") + 
  ggtitle("EEG signal as function of trial number and performed action") + 
  scale_fill_manual(values = c("blue", "red")) +
  theme_classic() + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20), 
        title = element_text(size = 15)
        )

# END OF FILE.
