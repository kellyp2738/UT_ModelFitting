#!/usr/bin/python

import numpy as np
import math
import VeggieDeath_multivariate_logscale as VD

initial_pops=[10,10,10,10,10,10,10,10,10,10,490,5,5,1000,100,100,100,100,100,100]
run_pars=[0.5, 0.06, 0.26]
r_all=VD.SIRprevalences(run_pars, initial_pops, 1000, 0, 1, 500)
VD.mkPlot(r_all[0])
		