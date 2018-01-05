# Probabilistic-RfD
Data and code for probabilistic RfD calculations conducted by Chiu et al. (2018) 

## Code
HDMI_functions.R: R functions to do probabilistic RfD calculations
HDMI_calculations.R: R script to for implementing automated probabilistic RfD workflow

## Inputs
RfD_InputData.csv: Input CSV file of curated RfD database, used by HDMI_calculations.R
RfD_InputLabels.csv: CSV file describing columns in RfD_InputData.csv

## Outputs of HDMI_calculations.R
RfDEffectTypes.csv: Output CSV file with frequency of effect types in RfD_InputData.csv
RfDPODs.csv: Output CSV file with frequency of different types of points of departure (PODs) in RfD_InputData.csv, separated by source
RfDUFcomps.csv: Output CSV file with frequency of composite uncertainty factors (UFs) in RfD_InputData.csv, separated by source
RfDEndpoints.csv: Output CSV file with endpoint-specific RfDs and corresponding new composite UFs
RfD_HDMI_results.csv: Output CSV file with calculated probabilistic RfDs
RfD_HDMI_mc_results.csv: Output CSV file with calculated probabilistic RfDs, including values calculated by Monte Carlo simulation
RfD_OutputLabels.csv: CSV file describing the columns of RfD_HDMI_results.csv and RfD_HDMI_mc_results.csv
