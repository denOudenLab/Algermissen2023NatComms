#!/bin/bash

# FigureS11.sh

# Plots for Supplementary Figure S12.
# Just executes ./extract_ROI_HOA_GLM1_TAfT.sh.
# Similar script as for Supplementary Figure S11.
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2023.

# we are here:
# cd /project/3017042.02/Analyses/FiguresOutcomeLocked

# Inventory cortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-cort-maxprob-thr25-2mm.xml
# Inventory subcortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-sub-maxprob-thr25-1mm.xml

# Set root directory:
rootDir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure

./${rootDir}/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_TAfT.sh

# END OF FILE.
