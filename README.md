# matlab-codes-IOP-2026
# MATLAB codes for "A Frequency-Optimized Optogenetic Study of Network-Level Potentiation in Cortical Cultures on Microelectrode Arrays"

This repository contains the MATLAB codes used in the paper:

*"A Frequency-Optimized Optogenetic Study of Network-Level Potentiation in Cortical Cultures on Microelectrode Arrays"*  
Authors: Matteo Dominici, Ilya Auslender, Clara Zaccaria, Yasaman Heydari, and Lorenzo Pavesi  
Journal: IOP Publishing, Year

## Contents
- **test_stimulus_analysis.m** (script) : runs the full MEA analysis pipeline (spike detection, Tr/Tb spike-rate + 8×8 diagnostic maps, and sigmoid fitting of selected channels)
                                          across multiple HDF5 files and saves the final fitResults structure
- **boxplots_Tr_Tb.m** (script) : post-processing and aggregation of sigmoid fit results
- **analysis_LTP.m** (script) : performs spike-based analysis of MEA HDF5 recordings from LTP protocols, computing PSTH-based efficacy changes, classifying electrodes into potentiated
                                and non-potentiated groups, and saving results for group-level analysis
- **LTP_average.m** (script) : aggregates ΔA (efficacy) measurements from multiple LTP experiments, computes population averages and variability, and plots the mean temporal evolution for blue and red electrodes
- **AnalogDataExtr.m** (function) : extracts the time vector, analogue signals (µV), and channel identifiers from an HDF5 recording file over a specified sample range
- **DetectSpikes_PTSD_best.m** (function) : detects spikes in filtered MEA signals using adaptive peak-to-peak thresholding
- **AnalysisFile.m** (function) : performs band-pass filtering, spike detection, and stimulation onset extraction from HDF5 electrophysiological recording files
- **ComputePSTH.m** (function) : calculates a PSTH from spike times aligned to stimulation events using user-defined binning and analysis windows
- **fitSigmoidForSelectedChannels** (function) : performs sigmoid fitting of raw and smoothed spike-rate responses for selected MEA channels, returning fit parameters in physical units and goodness-of-fit metrics





## Requirements
- MATLAB R20xx

## Notes
The codes are provided for reproducibility of the results presented in the paper.
