# EEGfMRIPav_1_Read_Initial

EEGfMRIPav_1_Read_Initial <- function(){
  #' EEGfMRIPav_1_Read_Initial()
  #' 
  #' Reads to .csv files obtained by executing EEGfMRIPav_extract_rawdata.m,
  #' give names to columns, recode required action and response, convert to factors, 
  #' save as EEGfMRIPav_all.csv
  #'
  #' Mind to adjust rootDir and targetDir.
  #'
  #' INPUTS:
  #' none
  #'
  #' OUTPUTS:
  #' @return saves all data of all subjects under EEGfMRIPav_all.csv in targetDir
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2023.

  # ============================================================================================== #  
  #### 1) Load and combine all data frames ####

  ## Set directories:
  rootDir <- "/project/3017042.02/" # root directory--needs to be adapted to users' folder structure
  inputDir <- paste0(rootDir, "Log/Behavior/Data_beh_csv")
  targetDir <- paste0(rootDir, "Log/Behavior/Behavior_concatenated/")

  ## Load all data files and concatenate them 
  inputData <- do.call("rbind", lapply(list.files(inputDir, full = TRUE), read.csv, header = F))
  # should be 36*2*320 = 23040 entries of 10 variables

  # ============================================================================================== #  
  #### 2) Name variables; ####
  
  cat("Set variable names\n")
  names(inputData) <- c("PPN_n", "session_n", "trial_nr_n", "stim_n", "req_action_n", "RT", "ACC_n", "response_n", "outcome_n", "is_go_n")

  # ============================================================================================== #  
  #### 3) Convert subject/ session/ trial variables to factors; ####
  
  cat("Participant, session, and trial number from numeric to factor\n")
  inputData$PPN_f <- factor(inputData$PPN_n)
  inputData$session_f <- factor(inputData$session_n)
  inputData$trial_nr_f <- factor(inputData$trial_nr_n)
  inputData$trial_nr_all_n <- (inputData$session_n -1 ) * 320 + inputData$trial_nr_n # trial number across sessions
  inputData$PPN_Ses_Trial_f <- factor(paste(inputData$PPN_f, inputData$session_f, inputData$trial_nr_n, sep = "_"))
  
  cat("Create block variable\n")
  # if (!(inputData[iRow, "trial_nr_n"]) %in% c(110,220,320,430,540,640)){
    
  inputData$block_n <- ifelse(inputData$session_n == 1 & inputData$trial_nr_n %in% 1:110, 1,
                           ifelse(inputData$session_n == 1 & inputData$trial_nr_n %in% 111:220, 2,
                                  ifelse(inputData$session_n == 1 & inputData$trial_nr_n %in% 221:320, 3,
                                         ifelse(inputData$session_n == 2 & inputData$trial_nr_n %in% 1:110, 4,
                                                ifelse(inputData$session_n == 2 &inputData$trial_nr_n %in% 111:220, 5,
                                                       ifelse(inputData$session_n == 2 & inputData$trial_nr_n %in% 221:320, 6,
                                                              NA))))))
  table(inputData$block_n)
  inputData$block_f <- factor(inputData$block_n)
  
  # ============================================================================================== #  
  #### 4) Convert stimulus features to factors: ####
  
  cat("Numeric stimulus variables to factors\n")
  inputData$stim_win_f <- factor(ifelse(inputData$stim_n %in% c(1, 2, 5, 6), "Win", "Avoid"))
  inputData$stim_win_n <- ifelse(inputData$stim_n %in% c(1, 2, 5, 6), 1, 0)
  inputData$stim_go_f <- factor(ifelse(inputData$stim_n %in% 1:4, "Go", "NoGo"))
  inputData$stim_go_n <- ifelse(inputData$stim_n %in% 1:4, 1, 0)
  inputData$stim_cong_f <- factor(ifelse(inputData$stim_win_n == inputData$stim_go_n, "congruent", "incongruent"))
  inputData$stim_cong_n <- ifelse(inputData$stim_win_n == inputData$stim_go_n, 1, 0)

  ## Ordinal factors:
  cat("Standard factors to ordered factors\n")
  inputData$stim_win_o <- ordered(inputData$stim_win_f, levels = c("Win", "Avoid"))
  inputData$stim_go_o <- ordered(inputData$stim_go_f, levels = c("Go", "NoGo"))
  
  # ============================================================================================== #  
  #### 5) Recode required actions (look up in pavParams.m: 101 = left, 97 = right): ####
  
  cat("Required action from numeric to factor\n")
  inputData$req_action_f <- factor(ifelse(inputData$req_action_n == 101, "Left",
                                       ifelse(inputData$req_action_n == 97, "Right", "NoGo")))
  
  # ============================================================================================== #  
  #### 6) Recode recorded response, convert into factor: ####

  cat("Uninstructed keys\n")
  inputData$is_uninstructed_key <- ifelse(inputData$response_n %in% c(0, 97, 101), 0, 1)
  
  cat("Recode instructed keys to respective correct key of same response side\n")
  inputData$response_n <- ifelse(inputData$response_n %in% c(101:104, 69:72), 101,
                              ifelse(inputData$response_n %in% c(97:100, 65:68), 97,
                                     ifelse(inputData$response_n == 0, 0,
                                            0)))

  cat("Performed response from numeric to factor\n")
  inputData$response_f <- factor(ifelse(inputData$response_n == 101, "Left",
                                     ifelse(inputData$response_n == 97, "Right",
                                            ifelse(inputData$response_n == 0, "NoGo", 
                                                   NA))))
  
  # ============================================================================================== #  
  #### 7) Convert accuracy to factor: ####
  
  cat("Accuracy from numeric to factor\n")
  inputData$ACC_f <- factor(ifelse(inputData$ACC_n == 1, "Correct", "Incorrect"))
  
  # ============================================================================================== #  
  #### 8) Convert outcome to factor: ####
  
  cat("Outcome from numeric to factor\n")
  inputData$outcome_f <- factor(ifelse(inputData$outcome_n == 1, "Reward",
                                    ifelse(inputData$outcome_n == 0, "Neutral",
                                           "Punishment")))
  
  # ============================================================================================== #  
  #### 9) Reorder rows: ####
  
  cat("Reorder rows based on participant number, session number, trial number\n")
  outputData <- inputData[order(inputData$PPN_n, inputData$session_n, inputData$trial_nr_n), ]
  
  # ============================================================================================== #  
  #### 10) Save: ####
  
  ## Save concatenated data:
  cat("Save data as EEGfMRIPav_all.csv in ", targetDir, "\n")
  write.csv(outputData, paste0(targetDir, "EEGfMRIPav_all.csv"), row.names = F)
  
}
# END OF FILE.
