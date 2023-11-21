#'% FigureS21ABC.R
#'
#' Plots for Supplementary Figure S21, panels A, B, and C.
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
#### 1) Load full data set ####

## Set directories:
rootDir   <- "/project/3017042.02/"
scriptDir <- paste0(rootDir, "Analyses/Behavior_Scripts/")
fmriDir <- paste0(rootDir, "Log/fMRI/")
plotDir <- paste0(rootDir, "Log/OutcomeLockedPaperPlots/")

## Load functions:
setwd(scriptDir)
source(paste0(scriptDir, "EEGfMRIPav_1_Read_Initial.R")) # Load variables and bring to R data frame format
source(paste0(scriptDir, "EEGfMRIPav_2_Preprocess_Automated.R")) # Recode variables
source(paste0(scriptDir, "EEGfMRIPav_0_CreateBarPlots.R")) # # Retrieve cue position

## Apply functions to load behavior:
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.

# ========================================================================================================== #
#### 2) Load outcome-locked EEG trial-by-trial data: ####

## Number subjects:
nSub <- length(unique(mydata$PPN_n))

## ROI names:
roiVec <- c("Action_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_lowalpha", "ActionLong_CzFCzFz_beta")

for (varName in roiVec){
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

# ========================================================================================================== #
#### 3) Create plots: ####

## Select subjects for EEG: 
invalidSubs <- c(11, 12, 23, 30) ## People excluded after EEG
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

## Names for plotting:
letterVec <- c("A", "B", "C")

## Loop over ROIs:
for (iROI in 1:length(roiVec)){
  selVar <- roiVec[iROI]
  
  ## Create plot:
  p <- custom_barplot2(modData, xVar = "block_n", yVar = selVar, zVar = "is_go_cor_f", subVar = "PPN_f",
                  xLab = "Block number", yLab = "TF power (dB)", zLab = "Performed action", selCol = c("red", "blue"),
                  isPoint = F, savePNG = F, saveEPS = F)
  ## Save:
  png(paste0(plotDir, "FigS21", letterVec[iROI], ".png"), width = 480, height = 480)
  print(p)
  dev.off()
  
  ## Aggregate data:
  aggrData <- ddply(modData, .(PPN_f, block_n, is_go_cor_f), function(x){
    EEG <- mean(x[, selVar], na.rm = T)
    return(data.frame(EEG))
    dev.off()})
  
  ## Save source data file:
  write.csv(aggrData, paste0(plotDir, "FigS21", letterVec[iROI], ".csv"), row.names = F)
  
}

# END OF FILE.