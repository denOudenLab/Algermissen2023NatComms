#' preprocess_PEs_outcomeLocked.R
#' 
#' Load .txt file output from featquery, extract relevant summary statistic,
#' rearrange conditions, save for all subjects in one .csv file. 
#' 
#' Requires that Run_ROI/run_ROI_PEs_2ndlevel.sh has been run and .txt files 
#' are available under e.g. sub-001/GLM1/featquery_GLM2ACCAction
#'
#' Mind to adjust rootDir.
#'
#' INPUTS:
#' none
#'
#' OUTPUTS:
#' no outputs, saves .txt file to Log/fMRI/fMRI_ROIs/ROIPlots.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2023.

# ========================================================================================================================================
#### General settings: ####

## GLM to extract from:
GLMID <- "2" # GLM where PEs were extracted from

## Parameter estimates and their names:
peVec <- c(1, 3, 2, 4, 5, 6, 7, 8) # use for selecting and reordering PEs
nPE <- length(peVec) # number of PEs to extract

# Valid subjects:
invalidSubs <- c(15, 25) # subjects to exclude (no PEs available)
validSubs <- which(!(1:36 %in% invalidSubs))
nSub <- length(validSubs)

# ======================================================================================================================================== #
#### Set directories: ####

rootDir <- "/project/3017042.02/"
plotDir <- paste0(rootDir, "Log/fMRI/fMRI_ROIs/ROIPlots/")

# Directory dependent on GLMID:
targetDir <- paste0(rootDir, "Log/fMRI/fMRI_ROIs/GLM", GLMID, "_ROIData/")
if (!dir.exists(targetDir)) {dir.create(targetDir)}

# ======================================================================================================================================== #
#### Create function create_wide() to create wide data frame: ####

create_wide <- function(GLMID, ROIID){

  cat(paste0("Load PEs from GLM ", GLMID, ", ROI ", ROIID, "\n"))
  
  # 1) Initialize empty data frame:
  ROI_wide <- data.frame(PPN_n = numeric(nSub),
                         PE1 = numeric(nSub), 
                         PE2 = numeric(nSub), 
                         PE3 = numeric(nSub), 
                         PE4 = numeric(nSub), 
                         PE5 = numeric(nSub), 
                         PE6 = numeric(nSub), 
                         PE7 = numeric(nSub), 
                         PE8 = numeric(nSub), 
                         stringsAsFactors = FALSE) 
  
  # 2) Create subject ID:
  ROI_wide$PPN_n <- validSubs # skip subjects 15 and 25 because of unreliable registrations
  ROI_wide$PPN_f <- as.factor(ROI_wide$PPN_n) # as factor
  ROI_wide$group <- 1 # fake x-axis for plots
  
  # 3) Loop through subjects and load:
  iCount <- 0 # initialize count # fill subjects consecutively, skip subject IDs not contained in validSubs
  for (iSub in unique(ROI_wide$PPN_n)){ # iSub <- 1 # loop over valid subjects
    iCount <- iCount + 1;

    # Set directory to respective subject/GLM/ROI:
    dataDir <- sprintf("%s/fMRI/sub-%03d/GLM%s/GLM%s_sub%03d.feat/",rootDir, iSub, GLMID, GLMID, iSub)
    fullFileName <- paste0(dataDir, "featquery_GLM", GLMID, ROIID, "/report.txt")
    
    # Columns are: name, number of voxels, min, 10%, mean, median, max, stddev, 3 columns vox(fMRI space), 3 columns mm(standard space) 
    report <- read.delim(fullFileName, sep = " ", header = F) # load
    names(report) <- c("X", "PEname", "nVoxels", 
                       "min", "10%", "mean", "median", "90%", "max", "stddev", 
                       "vox1", "vox2", "vox3", 
                       "mm1", "mm2", "mm3")
    
    # Extract mean (or different summary statistic):
    for (iPE in 1:nPE){
      ROI_wide[iCount, paste0("PE", iPE)] <- report[iPE, "mean"]
    }
  }
  return(ROI_wide)
}

# ======================================================================================================================================== #
#### Loop over ROIs to use, execute create_wide(): ####

# GLMID <- 2 # set above

for (ROIID in c("ACCAction", "vmPFCValence", "StriatumValence")){
  
  # Load and extract:
  ROI_wide <- create_wide(GLMID, ROIID)
  
  # Write as .txt file for Matlab: 
  output <- t(as.matrix(ROI_wide[, c("PE1", "PE2", "PE3", "PE4", "PE5", "PE6", "PE7", "PE8")]))
  outFileName <- paste0(targetDir, "GLM", GLMID, ROIID, ".txt")
  write.table(output, outFileName, row.names = F, col.names = F)
  
}

cat("Finished :-)\n")

# END OF FILE.