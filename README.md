# matlab-codes-IOP-2026
# MATLAB codes for "A Frequency-Optimized Optogenetic Study of Network-Level Potentiation in Cortical Cultures on Microelectrode Arrays"

This repository contains the MATLAB codes used in the paper:

"A Frequency-Optimized Optogenetic Study of Network-Level Potentiation in Cortical Cultures on Microelectrode Arrays"  
Authors: Matteo Dominici, Ilya Auslender, Clara Zaccaria, Yasaman Heydari, and Lorenzo Pavesi  
Journal: IOP Publishing, Year

## Contents
- script1.m : description
- LTP_average.m (script) : aggregates ΔA (efficacy) measurements from multiple LTP experiments, computes population averages and variability, and plots the mean temporal evolution for blue and red electrodes
- boxplots_Tr_Tb.m (script) : post-processing and aggregation of sigmoid fit results from single-electrode LTP experiments
- AnalogDataExtr.m (function) : extracts the time vector, analogue signals (µV), and channel identifiers from an HDF5 recording file over a specified sample range


## Requirements
- MATLAB R20xx

## Notes
The codes are provided for reproducibility of the results presented in the paper.
