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