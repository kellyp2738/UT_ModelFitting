# UT_ModelFitting

This repository contains the scripts used for modeling transmission of *Ehrlichia chaffeensis* in tick and wildlife populations.

###MCMC_multivariate_logscale_metropolis.py

**Dependency** VeggieDeath_multivariate_logscale.py

This script imports the SIR model script (below) for *E. chaffeensis* transmission and uses a Bayesian Markov Chain Monte Carlo procedure to estimate three key parameters of interest: preference of ticks for host animals, probability of *E. chaffeensis* transmission from deer to ticks, and probability of *E. chaffeensis* transmission from ticks to deer.

Parameter proposals are generated from a multivariate lognormal proposal distribution.

For usage information, see

	python MCMC_multivariate_logscale_metropolis.py --help

###VeggieDeath_multivariate_logscale.py

**Do not run this script directly.** It is set up to be imported as a module into other helper-scripts.

This script contains the SIR model itself, as well as functions for plotting single runs and calculating equilibrium disease prevalence.

The model includes the following compartments:

1. On-Host Ticks
	- susceptible larvae
	- susceptible nymphs
	- susceptible adults
	- infected nymphs
	- infected adults
	
2. Off-Host ('Veggie') Ticks
	- susceptible questing larvae
	- susceptible replete larvae/questing nymphs
	- susceptible replete nymphs/questing adults
	- susceptible replete adults
	- infected replete larvae/questing nymphs
	- infected replete nymphs/questing adults
	- infected replete adults
3. Deer
	- susceptible
	- infected
	- recovered

Previous versions of the MCMC script supported other model modules. The existence of this script as a separate file from the MCMC script is a legacy from those previous versions.


###toy_model_prevalence.py

**Dependency**
VeggieDeath_multivariate_logscale.py

This script outputs as CSV files the expected disease prevalence values from  a set of hard-coded preference and transmission parameter values (no parameter inference).

	python toy_model_prevalence.py
	
###RevisedDissertationPlots.r

**Note 1: this is not a stand-alone script. It should be run interactively in R or R Studio.**

**Note 2: currently the associated script *ticks.r* and the empirical disease prevalence data are not online -- these components are necessary for plotting empirical data**. 

The data (model outputs) are large files not available online. However, the python scripts can be used to generate new model fitting data that can be visualized with this file. Changes to directory paths for data will be necessary in order for this script to run.

Outputs from full model fitting runs are parsed and visualized in the following ways:

- Data from multiple MCMC chains are reorganized and checked for evidence of convergence using the Gelman-Rubin convergence statistic.
- Marginal posterior distributions for estimated parameters are visualized.
- Marginal posterior distributions of disease prevalence are visualized in conjunction with empirial prevalence distributions to check that the model outputs are consistent with empirical data. **see Note 2 above**
- Scatterplots showing the correlations between estimated parameter values are  generated.

###plot_toy_model.r

**Note: run interactively.**

This script plots the output of toy model runs for different preference values and relative abundance scenarios.