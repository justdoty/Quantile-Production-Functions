# Estimating Quantile Production Functions: A Control Function Approach
by Justin Doty and Suyong Song (2021)
## Abstract
We propose a new approach to estimating production functions when elasticities are heterogeneous across the conditional distribution of output. We provide a quantile regression framework which controls for unobserved productivity as an extension of the control function approaches of [Ackerberg *et al.* (2015)](https://doi.org/10.3982/ECTA13408) and [Levinsohn and Petrin (2003)](https://doi.org/10.1111/1467-937X.00246). Production function parameters are estimated in a simple two-stage approach which relies a location-shift assumption on unobserved productivity. We show that this method allows us to capture heterogeneity in output elasticities which may not be found in the conditional mean estimates. We provide small-sample evidence in a Monte Carlo study to show that this approach is robust compared to other production function estimators. The method is applied to firm and plant-level manufacturing data from the U.S., Chile, and Colombia. The findings confirm that the proposed method captures unobserved heterogeneity in output elasticities.

## Software Implementation
All code is written in R language and evaluated in batches on an HPC system.

## Replication Files
The following folders contain replication files for the tables and figures in the paper.
1. [Environments](/Environments): A folder containing all R environments used in the analysis
2. [Figures](/Figures): A folder containing all figures in the paper
	- [US Tables and Figures](/Figures/US/Tables_and_Figures_US.R): File used to produce all tables and figures from the US application
	- [Chile Tables and Figures](/Figures/Chile/Tables_and_Figures_CHL.R): File used to produce all tables and figures from the Chile application
	- [Colombia Tables and Figures](/Figures/Colombia/Tables_and_Figures_COL.R): File used to produce all tables and figures from the Colombia application
3. [Functions](/Functions): A folder containing all the replication files
	- [Auxiliary Files](/Functions/Aux_Fun.R): Contains auxiliary files used in estimation procedures
	- [US Data Cleaning](/Functions/Compustat_Cleaning.R): File detailing methodology for cleaning data from Compustat
	- [Chile Data Cleaning](/Functions/ENIA_Cleaning.R): File detailing methodology for cleaning data from Chile
	- [Colombia Data Cleaning](/Functions/COL_Cleaning.R): File detailing methodology for cleaning data from Colombia
	- [Monte Carlo ACF](/Functions/Monte_Carlo_ACF.R): Monte Carlo code for ACF data generating processes
	- [Monte Carlo LP](/Functions/Monte_Carlo_LP.R): Monte Carlo code for LP data generating processes
	- [DS ACF Estimator](/Functions/QACF_Boot.R): Estimation file when ACF is used as an initial estimate of productivity
	- [DS LP Estimator](/Functions/QLP_Boot.R): Estimation file when LP is used as an initial estimate of productivity for gross-output production function
	- The remaining files found in this folder are used to send the code to the HPC system

