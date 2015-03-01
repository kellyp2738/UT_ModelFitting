#!/usr/bin/python

import VeggieDeath_multivariate_logscale as VD

#run_pars=[0.2, 0.06, 0.26]
initial_pops=[10,10,10,10,10,10,10,10,10,10,190,5,5,1000,100,100,100,100,100,100]
racc_pop=1000
prefs=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]

for p in prefs:
	run_pars=[p, 0.06, 0.26]
	r_all=VD.SIRprevalences(run_pars, initial_pops, 1000, 0, 1, racc_pop)
	r=r_all[0]
	model_output=r[-1:][0] # take the last entry from the model time series
	#print(model_output)
	prevs=VD.calculate_prev(model_output)
	print(prevs)                    