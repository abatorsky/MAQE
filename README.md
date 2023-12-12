# MAQE
This repository contains code to implement the simulations described in the publication "Integrating randomized and observational studies to estimate optimal dynamic
treatment regimes"

This folder contains R scripts and a .csv file for running the simulation experiments described in the manuscript.
The scripts and functions can be adapted for using real data or a different subset of simulation parameters. 

param.csv contains the simulation parameters for the different experiments, one row per experiment
Functions.R contains the functions used in RunMAQE.R and must be sourced to run that script
GenerateTestSet.R generates the 20,000-participant test set data using the simulation parameters described in param.csv
	It is sourced in RunMAQE.R before running the loop to run the simulations
RunMAQE.R is the script that runs the simulations and saves the results (PCC, value) for each simulation experiment
	Working directory needs to be specified by user at the top of the code
	Input: param.csv
	Output: CC-.csv file with %correctly classified for each estimator for each simulation run
		values-.csv file with values for each estimator for each simulation run
		parameterlist-.csv file that lists simulation parameters to make sure the parameters are those intended for that run

Options:

variable simnum can be changed to any number of simulations desired

variable selection
lists of vectors mu1var and mu0var can be changed to include different subsets of variables to estimate regression parameters to calculate potential outcomes; must include intercept
lists of vectors Hk can be changed to include different subsets of variables to estimate regression parameters in Q-function; must not include intercept
	MEk must be the same list of variables as Hk but must also include intercept


To use real data instead of running the simulations on simulated datasets, the user can utilize the function MAQE in file Functions.R
which takes in arguments as a list, where the first vector of the list is for stage 1 and the second vector is for stage 2 inputs
mu1var: variables for modeling potential outcomes where treatment=1
mu0var: variables for modeling potential outcomes where treatment=0
txind: treatment indicator variable (i.e., a1 or a2 in the simulated data)
Hk: variables for modeling the Q-function; do not include intercept
Yk: outcome variable (i.e., y1 or y2 in the simulated data)
Ak: treatment variable (i.e., a1 or a2 in the simulated data)
pi1k: probability of receiving treatment A=1 in trial
pi0k: probability of receiving treatment A=0 in trial
rctind: trial indicator variable (i.e., in simulated data if smart=1 that participant was in the trial, else if smart=0 they were in the OS)
MEk: must be the same as Hk for each stage but MUST include intercept also
RCTy: outcome variable in the trial (i.e., y1 or y2 in the simulated trial data)
RCTa: treatment variable in the trial (i.e., a1 or a2 in the simulated trial data)
group: indicates which data to use for potential outcome models; possible values: {RCT, OS, RCTOS}

This function can also be edited to change the definition of q, such that the regression parameters (DTR) estimated using the MAQE can be weighted more or less heavily by trial or OS data
