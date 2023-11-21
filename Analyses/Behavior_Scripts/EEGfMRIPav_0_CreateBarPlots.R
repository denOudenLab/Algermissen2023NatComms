# EEGfMRIPav_0_CreateBarPlots.R

# Various functions for aggregation and bar plotting.

# =============================================================================================== #
####  Create percentiles ####

create_percentiles <- function(data, varName, nPerc = 5, perSub = F){
  #' Create new variable of percentiles based on numeric variable in data frame.
  #' @param data      data frame with trial-level data.
  #' @param varName   string, variable for which to create percentiles.
  #' @param nPerc     scalar integer, number of percentiles to create.
  #' @param perSub    Boolean, create percentiles for each subject (T) or across all subjects (F), default F.
  #' @return data     same data frame with all character variables turned into proper factors.
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2023.

  ## Set subject variable:
  subVar <- "PPN_f" # name of subject identifier variable
  nSub <- length(unique(data[, subVar])) # number of subjects
  
  ## Copy over selected variable:
  data$selVar <- data[, varName] # selected variable
  
  if (perSub == F){ # percentiles across subjects
    
    # a) Across entire data set:
    newVarName <- paste0(varName, "_", nPerc, "Perc") # create new variable name
    data[, newVarName] <- NA # initialize new variable
    
    data[, newVarName] <- with(data, cut(selVar, 
                                        breaks = quantile(selVar, probs = seq(0, 1, by = 1/nPerc), na.rm = TRUE), 
                                        include.lowest = TRUE))
    
  } else { # percentiles for each subject separately
    
    ## b) For each subject separately:
    newVarName <- paste0(varName, "_", nPerc, "PercSub") # create new variable name
    data[, newVarName] <- NA # initialize new variable
    
    # Loop over subjects:
    for (iSub in 1:nSub){ # iSub <- 1
      subIdx <- which(data[, subVar] == iSub) # rows of data from this subject
      subData <- data[subIdx, ] # select data from this subject
      subQuantiles <- quantile(data$selVar, probs = seq(0, 1, by = 1/nPerc), na.rm = TRUE) # quantile cutoffs
      data[subIdx, newVarName] <- with(subData, cut(selVar, breaks = subQuantiles, 
                                                   include.lowest = TRUE))
    }
  }
  
  # Turn into numeric:
  data[, newVarName] <- as.numeric(data[, newVarName]) # into numeric
  cat(paste0("Take variable", varName, ", add variable ", newVarName, "\n"))
  
  # Give bin sizes:
  cat("Bin sizes:\n")
  print(table(data[, newVarName])) # bin sizes
  # Output:
  cat("Finished :-)\n")
  return(data)
}

# ================================================================================================================================================ #
#### 2) Aggregate per condition per subject, plot (1 IV on x-axis): ####

custom_barplot1 <- function(data, xVar = NULL, yVar = NULL, subVar = "PPN_f", 
                            xLab = "BOLD signal (a.u.)", yLab = "p(Stay)", selCol = "red",
                            isPoint = T, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make bar plot with error bars and individual-subject data points. 1 independent variable (x-axis).
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param subVar string, name of variable containing subject identifier.
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param selCol vector of strings (HEX colors), colors for bars (default: "red").
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param savePNG Boolean, save as .png file (default: false).
  #' @param saveEPS Boolean, save as .eps file (default: false).
  #' @return creates (and saves) plot.
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2023.

  ## Load required packages:
  require(plyr) # Version ‘1.8.6’ # for ddply
  require(Rmisc) # Version '1.5' # for summarySEwithin
  require(ggthemes) # Version '4.2.0' # for custom theme
  
  ## Fixed plotting settings:
  LWD <- 1.5
  FTS <- 20
  dodgeVal <- 0.6
  colAlpha <- 1
  nRound <- 3
  
  ## Create variables under standardized names:
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$subject <- data[, subVar]
  
  # ------------------------------------------------------------------------- #
  ### 1) Aggregate data:
  
  aggrData <- ddply(data, .(subject, x), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
    dev.off()})
  cat(paste0("Min = ", round(min(aggrData$y), nRound), "; Max = ", round(max(aggrData$y), nRound)), "\n")
  
  ## Add jittered x-axis variable for points:
  aggrData$xj <- as.numeric(aggrData$x) # to numeric
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = .09) # jitter
  aggrData$xj <- aggrData$xj + aggrData$j # add jitter
  
  # ------------------------------------------------------------------------- #
  ### 2) Aggregate across subjects with Rmisc:
  
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x"))
  
  
  # ------------------------------------------------------------------------- #
  ### 3) ggplot:
  
  ## Figure name for saving:
  plotName <- paste0(yVar, "_", xVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(plotDir, plotName, ".png"), width = 480, height = 480)}
  
  ## Start plot:
  p <- ggplot(summary_d,aes(x, y))
  
  ## Bars of means:
  p <- stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.6,
                    lwd = LWD, fill = selCol, color = "black") 
  
  ## Error bars:
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = 0.1,
                         lwd = LWD, color = "black", alpha = 1)
  
  ## Individual data points:
  if (isPoint){
    p <- p + geom_point(data = aggrData, aes(x = xj), shape = 1, size = 2, stroke = 1, # size = 0.6, 
                        color = "black", alpha = colAlpha)
  }
  
  ## Further settings:
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim = yLim)}
  p <- p + labs(x = xLab, y = yLab) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30), 
          title = element_text(size = 30),
          axis.line = element_line(colour = 'black', size = LWD)
          ) # fixed font sizes
  
  ## Print ans save:
  print(p)
  if(savePNG | saveEPS){dev.off()}
  return(p)
}

# ================================================================================================================================================ #
#### 3) Aggregate per condition per subject, plot (2 IVs, 1 on x-axis, on as color/ adjacent bars): ####

custom_barplot2 <- function(data, xVar = NULL, yVar = NULL, zVar = NULL, subVar="PPN_f", 
                           xLab = "BOLD signal (a.u.)", yLab = "p(Stay)", zLab = "Action", selCol = c("blue", "red"),
                           isPoint = T, yLim = NULL, savePNG = F, saveEPS = F){
  #' Make bar plot with error bars and individual-subject data points. 2 independent variables (x-axis and colour).
  #' @param data data frame, trial-by-trial data.
  #' @param xVar string, name of variable that goes on x-axis. If numeric, it will be converted to an (ordered) factor.
  #' @param yVar string, name of variable that goes on y-axis. Needs to be numeric.
  #' @param zVar string, name of variable that determines bar coloring. Needs to be a factor.
  #' @param subVar string, name of variable containing subject identifier.
  #' @param xLab string, label for x-axis (default: "x").
  #' @param yLab string, label for y-axis (default: "y").
  #' @param zLab string, label for color legend (default: "z").
  #' @param selCol vector of strings (HEX colors), colors for input levels of zVar (default: c("blue", "red")).
  #' @param yLim vector of two numbers, y-axis (default: automatically determined by ggplot).
  #' @param isPoint Boolean, plot individual data points per condition as small points (default: TRUE).
  #' @param savePNG Boolean, save as .png file (default: false).
  #' @param saveEPS Boolean, save as .eps file (default: false).
  #' @return creates (and saves) plot.
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2023.
  
  ## Load packages:
  require(plyr) # Version ‘1.8.6’ # for ddply
  require(Rmisc) # Version '1.5' # for summarySEwithin
  library(ggplot2) # Version '3.3.4' # for ggplot
  require(ggthemes) # Version '4.2.0' # for custom theme
  
  ## Fixed plotting settings:
  LWD <- 1.5
  FTS <- 20
  dodgeVal <- 0.6
  colAlpha <- 1
  
  ## Create variables under standardized names:
  data$x <- data[, xVar]
  data$y <- data[, yVar]
  data$z <- data[, zVar]
  data$subject <- data[, subVar]
  
  # ------------------------------------------------------------------------- #
  ### 1) Aggregate data:
  
  aggrData <- ddply(data, .(subject, x, z), function(x){
    y <- mean(x$y, na.rm = T)
    return(data.frame(y))
  dev.off()})
  
  ## Add jittered x-axis for points:
  aggrData$xj <- as.numeric(aggrData$x) # to numeric
  aggrData$xj <- aggrData$xj + (as.numeric(aggrData$z) - 1.5) * 2 * 0.2 # convert to [1 2], to [-0.5,0.5], * 2 so [-1 1], scale
  aggrData$j <- jitter(rep(0, nrow(aggrData)), amount = .09) # jitter
  aggrData$xj <- aggrData$xj + aggrData$j # add jitter
  
  # ------------------------------------------------------------------------- #
  ### 2) Aggregate across subjects with Rmisc:
  
  summary_d <- summarySEwithin(aggrData, measurevar = "y", idvar = "subject", na.rm = T,
                               withinvars = c("x", "z"))
  
  # ------------------------------------------------------------------------- #
  ### 3) ggplot:
  
  ## Figure name for saving:
  plotName <- paste0(yVar, "_", xVar, "_", zVar)
  if (isPoint){plotName <- paste0(plotName, "_points")} 
  
  ## Saving:
  if (saveEPS){cat("Save as eps\n"); setEPS(); postscript(paste0(plotDir, plotName, ".eps"), width = 480, height = 480)}
  if (savePNG){cat("Save as png\n"); png(paste0(plotDir, plotName, ".png"), width = 480, height = 480)}
  
  ## Start plot:
  p <- ggplot(summary_d,aes(x, y, fill = z))
  
  ## Bars of means:
  p <- p + stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.6,
                        lwd = LWD, color = "black") 
  ## Error bars:
  p <- p + geom_errorbar(data = summary_d,
                         aes(x = x, y = y, ymin = y - se, ymax = y + se),
                         position = position_dodge(width = dodgeVal), width = 0.3,
                         lwd = LWD, color = "black", alpha = 1)
  
  ## Individual data points:
  if (isPoint){
    p <- p + geom_point(data = aggrData, aes(x = xj), shape = 0, size = 0.2, 
                        alpha = colAlpha)
  }
  
  ## Further settings:
  if(!(is.null(yLim))){p <- p + coord_cartesian(ylim = yLim)}
  
  p <- p + labs(x=xLab, y = yLab, fill = zLab) +
    scale_fill_manual(values = selCol) + 
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30), 
          title = element_text(size = 30),
          legend.title = element_blank(), 
          legend.position = "none")
  
  ## Print and save:
  print(p)
  if(savePNG | saveEPS){dev.off()}
  return(p)
}

# END OF FILE.
