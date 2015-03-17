#!/usr/bin/python

import numpy as np
import math
import VeggieDeath_multivariate_logscale as VD

def calculate_prev(model_output): # take the population sizes and turn them into prevalences. dependent on the ordering of the compartments in the model output
    #print(model_output)
    nymphal_prev=round(((model_output[3]+model_output[8]+model_output[17])/(model_output[1]+model_output[3]+model_output[6]+model_output[8]+model_output[14]+model_output[17])), 5)
    adult_prev=round(((model_output[4]+model_output[9]+model_output[18]+model_output[19])/(model_output[2]+model_output[4]+model_output[7]+model_output[9]+model_output[15]+model_output[16]+model_output[18]+model_output[19])), 5) 
    deer_prev=round(((model_output[11]/(model_output[10]+model_output[11]+model_output[12]))), 5)
    deer_AB_prev=round((model_output[11]+model_output[12])/(model_output[10]+model_output[11]+model_output[12]), 5)
                
    all_prevs=[nymphal_prev, adult_prev, deer_prev, deer_AB_prev]
    #print(adult_prev, round(adult_prev, 2))
    
    return all_prevs

#run_pars=[0.2, 0.06, 0.26]
racc_pops=[900, 700, 500]
deer_pops=[90, 290, 490]

#prefs=[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
prefs=np.arange(0.01,1,0.01)

header=['prefs', 'Scenario_1', 'Scenario_2', 'Scenario_3']
save_tick=np.array([[0.0]*len(header)]*len(prefs))
save_deer=np.array([[0.0]*len(header)]*len(prefs))
save_deerAB=np.array([[0.0]*len(header)]*len(prefs))
for j in range(len(prefs)):
	p=prefs[j]
	save_tick[j,0]=p
	save_deer[j,0]=p
	save_deerAB[j,0]=p
	for i in range(3):
		initial_pops=[10,10,10,10,10,10,10,10,10,10,deer_pops[i],5,5,1000,100,100,100,100,100,100]
		run_pars=[p, 0.06, 0.26]
		r_all=VD.SIRprevalences(run_pars, initial_pops, 1000, 0, 1, racc_pops[i])
		#VD.mkPlot(r_all[0])
		r=r_all[0]
		model_output=r[-1:][0] # take the last entry from the model time series
		prevs=calculate_prev(model_output)
		save_tick[j,i+1]=prevs[1]
		save_deer[j,i+1]=prevs[2]
		save_deerAB[j,i+1]=prevs[3]
                 
#save_final_tick=np.vstack((header, save_tick))
#save_final_deer=np.vstack((header, save_deer))
#save_final_deerAB=np.vstack((header, save_deerAB))
#np.savetxt('toy_model_scenarios_tick_out.csv', save_final_tick, delimiter=',', fmt='%20s')
#np.savetxt('toy_model_scenarios_deer_out.csv', save_final_deer, delimiter=',', fmt='%20s')
np.savetxt('toy_model_scenarios_deerAB_out.csv', save_final_deerAB, delimiter=',', fmt='%20s')