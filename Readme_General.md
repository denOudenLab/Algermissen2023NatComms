# README General

## Content ##
This collection contains the following files:
- **Analyses**: Scripts for all behavioral, EEG, and fMRI analyses required to reproduce the reported results
- **Log**: behavioral, EEG, and fMRI data required to reproduce the reported results

**Raw data** for this project are available under https://doi.org/10.34973/pezs-pw62. In line with requirements of the Ethics Committee and the Radboud University security officer, potentially identifying data (such as imaging data) can only be shared to identifiable researchers. Hence, researchers requesting access to the data have to register and accept a data user agreement; access will then automatically be granted via a click-through procedure (without involvement of authors or data stewards).

For more details, see the modality-specific READMEs.

## Root directory ##
Note that for most files, the **root directory** of this project needs to be adjusted. The folders in *rootDir/Analyses/EEG_Scripts/OutcomeLockedAnalyses* have one central file where the root directory can be set centrally (*XX_set_rootDir.m*). In all other files, check whether the root directory must be adjusted at the beginning.

## Figures ##
Under rootDir/Analyses/FiguresOutcomeLocked, see individual analysis scripts for reproducing figures in the manuscript and supplementary files.

## Source data files ##
Under rootDir/Log/SourceDataFiles, .xlsx spreadsheet files with the data underlying each figure in the manuscript and supplementary files are provided.

END OF FILE.
