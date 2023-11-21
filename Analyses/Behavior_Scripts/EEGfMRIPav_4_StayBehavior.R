#' EEGfMRIPav_4_StayBehavior()
#' 
#' Perform mixed-effects linear and logistic regression testing effects of 
#' outcome valence, outcome salience, and performed action on stay behavior.
#' as reported in behavioral results section.
#' Interactive script, run step-by-step, NOT in one Go.
#' 
#' Requires that .csv files have been recoded and concatenated using 
#' EEGfMRIPav_1_Read_Initial.R
#' 
#' Calls EEGfMRIPav_2_Preprocess_Automated.R
#' Calls EEGfMRIPav_2_Cue_Position.R
#' Calls EEGfMRIPav_2_Stay.R
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

# ================================================================================================================================================ #
#### Load packages: ####

require(plyr) # Version ‘1.8.6’ # for ddply
library(ggplot2) # Version '3.3.4' # for ggplot
library(stringr) # Version '1.4.0'

# Analyses:
library(lme4) # Version ‘1.1.26’ # for lmer and glmer
library(afex) # Version ‘0.28.1’ # for mixed
library(emmeans) # Version ‘1.5.5.1’ # for lsmeans and glht
library(effects) # Version ‘4.2.0’ # for effects plots

options(scipen = 20)
options(contrasts = c("contr.sum", "contr.poly"))
rm(list = ls())

# ================================================================================================================================================ #
#### 1) Load full data set ####

## Set directories:
rootDir   <- "/project/3017042.02/" # adjust root directory to local folder structure.
scriptDir <- paste0(rootDir, "Analyses/Behavior_Scripts")
modDir <- paste0(rootDir, "Log/Behavior/Behavior_Models/")

## Load functions:
source(paste0(scriptDir, "EEGfMRIPav_1_Read_Initial.R")) # Load variables and bring to R data frame format
source(paste0(scriptDir, "EEGfMRIPav_2_Preprocess_Automated.R")) # Recode variables
source(paste0(scriptDir, "EEGfMRIPav_2_Cue_Position.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_2_Stay.R")) # # Retrieve cue position

## Apply functions:
fullData <- EEGfMRIPav_1_Read_Initial() # Initial reading, converting to factors 
fullData <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.
fullData <- EEGfMRIPav_2_Cue_Position(fullData) # Compute cue position
fullData <- EEGfMRIPav_2_Stay(fullData) # Compute stay behavior

# ================================================================================================================================================ #
#### 2) Select data ####

## All subjects:
modData <- fullData

## Exclude certain subjects:
# modData <- subset(fullData, !(PPN_n %in% c(11, 12, 15, 23, 25, 26, 30))) # outcome-locked (N = 29)

# ================================================================================================================================================ #
#### 3) Mixed-effects model across all trials: ####

## Fit model:
mod <- glmer(stay_n ~ is_go_cor_f * outcome_salience_f * outcome_valence_f + (is_go_cor_f * outcome_salience_f * outcome_valence_f|PPN_f),
             modData, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod, paste0(modDir, "mod_stay_action_salience_valence.rds"))
mod <- readRDS(paste0(modDir, "mod_stay_action_salience_valence.rds"))
summary(mod)
#                                                      Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)                                          0.991923   0.120679   8.220 < 0.0000000000000002 ***
# is_go_cor_f1                                        -0.005394   0.056750  -0.095               0.9243    
# outcome_salience_f1                                  0.057602   0.025035   2.301               0.0214 *  
# outcome_valence_f1                                  -0.503843   0.053152  -9.479 < 0.0000000000000002 ***
# is_go_cor_f1:outcome_salience_f1                    -0.011031   0.031718  -0.348               0.7280    
# is_go_cor_f1:outcome_valence_f1                     -0.061850   0.033686  -1.836               0.0663 .  
# outcome_salience_f1:outcome_valence_f1               0.450365   0.063998   7.037     0.00000000000196 ***
# is_go_cor_f1:outcome_salience_f1:outcome_valence_f1  0.247688   0.047697   5.193     0.00000020694970 ***

## With LRT:
mod_LRT <- mixed(mod, modData, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT, paste0(modDir, "LRT_mod_stay_action_salience_valence.rds"))
mod_LRT <- readRDS(paste0(modDir, "LRT_mod_stay_action_salience_valence.rds"))
anova(mod_LRT)
#                                                  Df   Chisq Chi Df       Pr(>Chisq)    
# is_go_cor_f                                      43  0.0090      1          0.92453    
# outcome_salience_f                               43  5.1482      1          0.02327 *  
# outcome_valence_f                                43 45.5948      1 0.00000000001454 ***
# is_go_cor_f:outcome_salience_f                   43  0.1200      1          0.72905    
# is_go_cor_f:outcome_valence_f                    43  3.2367      1          0.07201 .  
# outcome_salience_f:outcome_valence_f             43 30.9545      1 0.00000002641475 ***
# is_go_cor_f:outcome_salience_f:outcome_valence_f 43 19.7324      1 0.00000890800787 ***

## Plots:
plot(effect("is_go_cor_f", mod)) # no effect
plot(effect("outcome_salience_f", mod)) # higher for neutral than salient (only sligthly)
plot(effect("outcome_valence_f", mod)) # higher for positive than negative
plot(effect("is_go_cor_f:outcome_valence_f", mod)) # stronger valence effect for Go (marginally significant)
plot(effect("is_go_cor_f:outcome_salience_f:outcome_valence_f", mod)) # stronger valence effect for Go, but only if salient

# ================================================================================================================================================ #
#### 3a) Go trials only: ####

## Select subset of data:
modData_Go <- subset(modData, is_go_cor_f == "Go")

## Fit model:
mod <- glmer(stay_n ~ outcome_salience_f * outcome_valence_f + (outcome_salience_f * outcome_valence_f|PPN_f),
             modData_Go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod, paste0(modDir, "mod_stay_action_salience_valence_onlyGo.rds"))
mod <- readRDS(paste0(modDir, "mod_stay_action_salience_valence_onlyGo.rds"))
summary(mod)
#                                        Estimate Std. Error z value            Pr(>|z|)
# (Intercept)                             0.98538    0.11857   8.311 <0.0000000000000002
# outcome_salience_f1                     0.04638    0.03632   1.277               0.202
# outcome_valence_f1                     -0.56614    0.06323  -8.954 <0.0000000000000002
# outcome_salience_f1:outcome_valence_f1  0.70017    0.06789  10.313 <0.0000000000000002
plot(effect("outcome_salience_f:outcome_valence_f", mod)) # valence effect only for Go

# LRTs:
mod_LRT <- mixed(mod, modData_Go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                                      Df   Chisq Chi Df        Pr(>Chisq)    
# outcome_salience_f                   13  1.5632      1            0.2112    
# outcome_valence_f                    13 43.2832      1 0.000000000047363 ***
# outcome_salience_f:outcome_valence_f 13 49.8219      1 0.000000000001684 ***

# ================================================================================================================================================ #
#### 3b) NoGo trials only: ####

## Select subset of data:
modData_NoGo <- subset(modData, is_go_cor_f == "NoGo")

## Fit model:
mod <- glmer(stay_n ~ outcome_salience_f * outcome_valence_f + (outcome_salience_f * outcome_valence_f|PPN_f),
             modData_NoGo, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod, paste0(modDir, "mod_stay_action_salience_valence_onlyNoGo.rds"))
mod <- readRDS(paste0(modDir, "mod_stay_action_salience_valence_onlyNoGo.rds"))
summary(mod)
#                                        Estimate Std. Error z value         Pr(>|z|)    
# (Intercept)                             0.98249    0.14847   6.617 0.00000000003654 ***
# outcome_salience_f1                     0.06665    0.04209   1.584           0.1133    
# outcome_valence_f1                     -0.44388    0.06250  -7.102 0.00000000000123 ***
# outcome_salience_f1:outcome_valence_f1  0.19444    0.09305   2.090           0.0366 * 

# LRTs:
mod_LRT <- mixed(mod, modData_NoGo, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                                      Df   Chisq Chi Df    Pr(>Chisq)    
# outcome_salience_f                   13  2.4136      1       0.12029    
# outcome_valence_f                    13 31.5002      1 0.00000001994 ***
# outcome_salience_f:outcome_valence_f 13  4.0145      1       0.04511 * 

# ================================================================================================================================================ #
#### 3c) Salient trials only: ####

## Select subset of data:
modData_Salient <- subset(modData, outcome_salience_f == "Salient")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f * is_go_cor_f + (outcome_valence_f * is_go_cor_f|PPN_f),
             modData_Salient, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod, paste0(modDir, "mod_stay_action_salience_valence_onlySalient.rds"))
mod <- readRDS(paste0(modDir, "mod_stay_action_salience_valence_onlySalient.rds"))
summary(mod)
#                                  Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)                      0.935629   0.117828   7.941  0.00000000000000201 ***
# outcome_valence_f1              -0.957350   0.098929  -9.677 < 0.0000000000000002 ***
# is_go_cor_f1                     0.004229   0.067929   0.062                 0.95    
# outcome_valence_f1:is_go_cor_f1 -0.308215   0.063622  -4.844  0.00000126946981369 ***
plot(effect("outcome_valence_f:is_go_cor_f", mod)) # stronger valence effect for Go than NoGo

# LRTs:
mod_LRT <- mixed(mod, modData_Salient, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                               Df   Chisq Chi Df        Pr(>Chisq)    
# outcome_valence_f             13 46.3554      1 0.000000000009863 ***
# is_go_cor_f                   13  0.0039      1            0.9501    
# outcome_valence_f:is_go_cor_f 13 17.7975      1 0.000024570573798 ***

# ================================================================================================================================================ #
#### 3d) Neutral trials only: ####

## Select subset of data:
modData_Neutral <- subset(modData, outcome_salience_f == "Neutral")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f * is_go_cor_f + (outcome_valence_f * is_go_cor_f|PPN_f),
             modData_Neutral, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod, paste0(modDir, "mod_stay_action_salience_valence_onlyNeutral.rds"))
mod <- readRDS(paste0(modDir, "mod_stay_action_salience_valence_onlyNeutral.rds"))
summary(mod)
#                                 Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)                      1.05213    0.12842   8.193 0.000000000000000255 ***
# outcome_valence_f1              -0.05597    0.06116  -0.915             0.360148    
# is_go_cor_f1                    -0.01913    0.05975  -0.320             0.748793    
# outcome_valence_f1:is_go_cor_f1  0.18774    0.04933   3.806             0.000141 ***
plot(effect("outcome_valence_f:is_go_cor_f", mod)) # Negative > Positive for Go, Positive > Negative for NoGo

# LRTs:
mod_LRT <- mixed(mod, modData_Neutral, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                               Df   Chisq Chi Df Pr(>Chisq)    
# outcome_valence_f             13  0.8295      1  0.3624217    
# is_go_cor_f                   13  0.1019      1  0.7495232    
# outcome_valence_f:is_go_cor_f 13 12.3201      1  0.0004481 ***

# ================================================================================================================================================ #
#### 3e) Salient Go trials only: ####

## Select subset of data:
modData_GoSalient <- subset(modData, is_go_cor_f == "Go" & outcome_salience_f == "Salient")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f + (outcome_valence_f|PPN_f),
             modData_GoSalient, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
# Fixed effects:
#                    Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)          0.9403     0.1217   7.729   0.0000000000000108 ***
# outcome_valence_f1  -1.2764     0.1152 -11.079 < 0.0000000000000002 ***
plot(effect("outcome_valence_f", mod)) # Positive > Negative

# LRTs:
mod_LRT <- mixed(mod, modData_GoSalient, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                   Df  Chisq Chi Df         Pr(>Chisq)    
# outcome_valence_f  4 53.932      1 0.0000000000002075 ***

# ------------------------------------------- #
#### 3f) Salient NoGo trials only: ####

## Select subset of data:
modData_NoGoSalient <- subset(modData, is_go_cor_f == "NoGo" & outcome_salience_f == "Salient")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f + (outcome_valence_f|PPN_f),
             modData_NoGoSalient, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
# Fixed effects:
#                    Estimate Std. Error z value      Pr(>|z|)    
# (Intercept)          0.9105     0.1526   5.966 0.00000000242 ***
# outcome_valence_f1  -0.6373     0.1274  -5.004 0.00000056163 ***
plot(effect("outcome_valence_f", mod)) # Positive > Negative

# LRTs:
mod_LRT <- mixed(mod, modData_NoGoSalient, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                   Df  Chisq Chi Df Pr(>Chisq)    
# outcome_valence_f  4 18.228      1  0.0000196 ***

# ------------------------------------------- #
#### 3g) Neutral Go trials only: ####

## Select subset of data:
modData_GoNeutral <- subset(modData, is_go_cor_f == "Go" & outcome_salience_f == "Neutral")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f + (outcome_valence_f|PPN_f),
             modData_GoNeutral, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
#                    Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)         1.03031    0.12742   8.086 0.000000000000000616 ***
# outcome_valence_f1  0.13370    0.06825   1.959               0.0501 . 
plot(effect("outcome_valence_f", mod)) # Negative > Positive

# LRTs:
mod_LRT <- mixed(mod, modData_GoNeutral, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                   Df  Chisq Chi Df Pr(>Chisq)  
# outcome_valence_f  4 3.5978      1    0.05786 .

# ------------------------------------------- #
#### 3h) Neutral NoGo trials only: ####

## Select subset of data:
modData_NoGoNeutral <- subset(modData, is_go_cor_f == "NoGo" & outcome_salience_f == "Neutral")

## Fit model:
mod <- glmer(stay_n ~ outcome_valence_f + (outcome_valence_f|PPN_f),
             modData_NoGoNeutral, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
#                    Estimate Std. Error z value         Pr(>|z|)    
# (Intercept)         1.05745    0.15425   6.855 0.00000000000711 ***
# outcome_valence_f1 -0.24414    0.08715  -2.801          0.00509 ** 
plot(effect("outcome_valence_f", mod)) # Positive > Negative

# LRTs:
mod_LRT <- mixed(mod, modData_NoGoNeutral, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                   Df  Chisq Chi Df Pr(>Chisq)   
# outcome_valence_f  4 7.2095      1   0.007252 **

# END OF FILE.