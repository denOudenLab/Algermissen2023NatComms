#'% FigureS22.R
#'
#' Plots for Supplementary Figure S22.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2023.
#'
#' we are here:
#' /project/3017042.02/Analyses/FiguresOutcomeLocked

# ========================================================================================================== #
#### Load packages: ####

library(ggplot2) # Version '3.3.4'
library(car) # Version '3.0.8'
library(lattice) # Version '0.20-44'
library(ltm) # Version '1.1-1'
library(pastecs) # Version '1.3-18'
library(plyr) # Version ‘1.8.6’
library(stringr) # Version '1.4.0'

# rm(list = ls())

# ========================================================================================================== #
#### 1) Load full data set: ####

## Set directories:
rootDir   <- "/project/3017042.02/" # adjust root directory to local folder structure.
scriptDir <- paste0(rootDir, "Analyses/Behavior_Scripts/")
behavDir <- paste0(rootDir, "Log/Behavior/")
fmriDir <- paste0(rootDir, "Log/fMRI/")
plotDir <- paste0(rootDir, "Log/OutcomeLockedPaperPlots/")

## Load functions:
source(paste0(scriptDir, "EEGfMRIPav_2_Preprocess_Automated.R")) # Recode variables
source(paste0(scriptDir, "EEGfMRIPav_2_Cue_Position.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_2_Stay.R")) # # Retrieve cue position
source(paste0(scriptDir, "EEGfMRIPav_0_CreateBarPlots.R")) # # Retrieve cue position

## Apply functions:
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.
mydata <- EEGfMRIPav_2_Cue_Position(mydata) # Compute cue position
mydata <- EEGfMRIPav_2_Stay(mydata) # Compute stay behavior

# ========================================================================================================== #
#### 2) Load outcome-locked BOLD trial-by-trial data: ####

## Get number of subjects:
nSub <- length(unique(mydata$PPN_n))

## Loop and retrieve all ROIs
roiVec <- c("GLM1ACCConjMan", "GLM1vmPFCConjMan", "GLM1StriatumConj")
for (iROI in 1:length(roiVec)){
  
  varName <- roiVec[iROI]
  
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

## Select subjects:
invalidSubs <- c(15, 25) # failed registration
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)
length(unique(modData$PPN_n))

# ========================================================================================================== #
#### 3) Plot p(Stay) ~ trial-by-trial BOLD amplitude: ####

## Set percentiles:
nPerc <- 5

## Variables to loop over:
roiVec <- c("GLM1ACCConjMan", "GLM1vmPFCConjMan", "GLM1StriatumConj")
colVec <- c("#FF0000", "#0070C0", "#FFC000")
letterVec <- c("A", "B", "C")

## Loop:
for (iROI in 1:length(roiVec)){ # iROI <- 1
  
  ## Create percentiles:
  selVar <- roiVec[iROI]; percVar <- paste0(selVar, "_", nPerc, "Perc");
  modData <- create_percentiles(modData, nPerc = nPerc, selVar)
  
  # Create plot with dots:
  p <- custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                  xLab = "BOLD (percentiles)", yLab = "p(Stay)", selCol = colVec[iROI],
                  isPoint = T, yLim = c(0.25, 0.95), savePNG = F, saveEPS = F)
  
  ## Save plot:
  png(paste0(plotDir, "FigS22", letterVec[iROI], ".png"), width = 480, height = 480)
  print(p)
  dev.off()
  
  ## Aggregate data:
  modData$percentile_n <- modData[, percVar]
  aggrData <- ddply(modData, .(PPN_f, percentile_n), function(x){
    pStay_n <- mean(x$stay_n, na.rm = T)
    return(data.frame(pStay_n))
    dev.off()})
  
  ## Save source data file:
  write.csv(aggrData, paste0(plotDir, "FigS22", letterVec[iROI], ".csv"), row.names = F)
  
}

# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
# ========================================================================================================== #
#### 4) Load outcome-locked EEG trial-by-trial data: ####

nSub <- length(unique(mydata$PPN_n))

## Set all EEG variables:
roiVec <- c("Action_CzFCzFz_lowalpha",
                 "Preferred_AF3AF4AF7AF8F1F2F3F4F5F6F7F8FC1FC2FC3FC4FC5FC6FCzFp1Fp2FpzFz_deltatheta",
                 "Preferred_CzFCzFz_beta")

### Loop over ROIs and retrieve all EEG variables:
for (iROI in 1:length(roiVec)){
  
  varName <- roiVec[iROI]
  
  ## Set appropriate EEG directory:
  eegDir <- paste0(rootDir, "Log/EEG/OutcomeLockedResults/TF_singleTrial/", varName, "/")
  
  ## Load:
  mydata[, varName] <- NA # initialize
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
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

### Log-transform EEG variables:
for (iROI in 1:length(roiVec)){
  varName <- roiVec[iROI]
  newVarName <- paste0(varName, "_Log")
  modData[, newVarName] <- 10*log(modData[, varName] - min(modData[, varName]) + 0.1)
}

# ========================================================================================================== #
#### 5) Plot p(Stay) ~ EEG TF Power: ####

## Set percentile:
nPerc <- 5

## Update ROI names:
roiLogVec <- paste0(roiVec, "_Log")
colVec <- c("#FF0000", "#0070C0", "#FFC000")
letterVec <- c("D", "E", "F")

## Loop:
for (iROI in 1:length(roiLogVec)){
  
  ## Select data:
  selVar <- roiLogVec[iROI]; percVar <- paste0(selVar, "_", nPerc, "Perc");
  modData <- create_percentiles(modData, nPerc = nPerc, selVar)
  
  # Create plot with dots:
  p <- custom_barplot1(modData, xVar = percVar, yVar = "stay_n", subVar = "PPN_f",
                  xLab = "TF power (percentiles)", yLab = "p(Stay)", selCol = colVec[iROI],
                  isPoint = T, yLim = c(0.25, 0.95), savePNG = F, saveEPS = F)
  ## Save plot:
  png(paste0(plotDir, "FigS22", letterVec[iROI], ".png"), width = 480, height = 480)
  print(p)
  dev.off()
  
  ## Aggregate data:
  modData$percentile_n <- modData[, percVar]
  aggrData <- ddply(modData, .(PPN_f, percentile_n), function(x){
    pStay_n <- mean(x$stay_n, na.rm = T)
    return(data.frame(pStay_n))
    dev.off()})
  
  ## Save source data file:
  write.csv(aggrData, paste0(plotDir, "FigS22", letterVec[iROI], ".csv"), row.names = F)
}

# END OF FILE.