# HepBModel

A detailed description of the model structure and parameters can be found in the published paper: https://pubmed.ncbi.nlm.nih.gov/31419332/
In addition, Ben Cowie's thesis also contains details of the original study and data used – some of which is still in use within the current model: https://minerva-access.unimelb.edu.au/items/372f724d-24ac-5152-8a62-fadec1366b94

Contained within the file "Inputs" are all the necessary inputs to execute the model. Provided below are the required Matlab functions to perform a single simulation of the model

i.	SingleRun_HBV2021.m -  this is the main file that calls the model files and simulates a single simulation of the HBV model. When running this make sure all the files listed below are in the same folder.
ii.	HBV2021_struct_Nat.m – this creates a structure containing all the relevant files to run a single simulation or Latin hypercube sampling experiments. It gathers the files containing parameters, the model equations, the initial state of the model etc. 
iii. Nat_params.m - This is where parameters are defined and their default values are set
iv.	Nat_init_state.m – defines the initial model state
v.	Nat_stats.m - statistics/output for the model to report
vi.	HBV2021_Eqs_report.m - This file describes the ordinary differential equations (ODEs) to be solved
vii.	Single_simulation.m, solve.m and solve_continuously.m - are functions required to solve the ODEs.
viii.	ABS_MigrationEstimates_2004onwards.m – this file collates required data to determine the number of people migrating into the population by age and disease state for the years 2004 to 2050. 
ix.	ABS_MigrationEstimates_1991to2003.m - this file collates required data to determine the number of people migrating into the population by age and disease state for the years 1991 to 2003.
x.	MigrationEstimates_1951to1990.m - this file collates required data to determine the number of people migrating into the population by age and disease state for the years 1951 to 1990.
