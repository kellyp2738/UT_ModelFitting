#!/usr/bin/python

#from __future__ import division #because otherwise 1/2 = 0. Thanks, python!
import pdb
import copy
import random as rand
#rand.seed(0) # set the seed for the random number generator 
import numpy as np
#np.random.seed(0) #set the seed for the random number generator (different seed from above)
from scipy import integrate
from scipy import stats
import scipy
import math
import pylab
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.figure as fig
from decimal import *
import decimal
import time
import csv
import os
import sys, getopt
#import pandas as pd

# ----------------------------------------------------------------------------------------------
# -- Read in parameters file
# ----------------------------------------------------------------------------------------------   
   
def acquire_params(filename):    
    params=[]
    #filename=raw_input("Please enter the name of the parameter file. Do not use quotes. Ensure that file to load is in the current directory and follows the format expected by this program! ")
    with open(filename, 'rU') as csvfile:
        param_object=csv.reader(csvfile, delimiter=",")
        for row in param_object:
            params.append(row)
    if len(params[0])>31:
        print 'PARAMETER ERROR: This version of the model does not include a second reservoir host. Please check the parameter configuration.'
        return
    for i in range(len(params[0])): #for each column, convert the values in rows 1-4 to floats (row 0 remains a str)
        if len(params) > 2:
            params[1][i]=float(params[1][i]) #starting value
            params[2][i]=float(params[2][i]) #prior mean
            params[3][i]=float(params[3][i]) #prior SD
        else:
            params[1][i]=float(params[1][i]) #starting value
            
    return params


# ----------------------------------------------------------------------------------------------
# -- Run SIRS model
# ----------------------------------------------------------------------------------------------

def SIRprevalences(run_params, R, finalTime, counter, num_prefs, racc_pop): 

    '''
    run_params = vector of input parameters (either fixed or estimated)
    R = vector or array of initial conditions
    finalTime = last time step at which to evaluate the system of equations
    counter = tracks how many times the function has been called recurisvely in an attempt to reach equilibrium. counter needs to be zero on the first call, it is subsequently updated by the function to keep track of the number of recursive calls

    This function contains an embedded function which is the model itself.
    This function returns the end point prevalences -- the infection prevalence where ever the model stops evalutating.
    The model (the embedded function) is hard-coded to run for 2000 time steps. Changing this requires a source code edit.
    The complete model output (time series) can be retrieved by editing this function to also return the complete output.'''

    # -- set counter to track number of recurions

    if counter == 0:
       x = R
       #print 'starting populations', x
    else:
       x=R[-1:] #get the last value in the numpy array
       x.tolist() #convert it from numpy array to regular list b/c odeint can't handle numpy arrays
       x=x[0] #the converted list is the wrong dimension -- it's [[x]], and it needs to be [x]

    #print 'Evaluating Prevalence'
    #pdb.set_trace()
    #print 'run pars in SIR model', run_params

    # -- determine which preference scheme is in use
    
    phiD=run_params[0]
    phiR=run_params[1]
    #phiD=math.exp(run_params[0])   # for a single pref, specify only phiDTL
    #phiR=1-math.exp(run_params[0])
    #print phiD
    
    # -- define demographic and epidemiologic constants
    # -- raccoon epi parameters have been removed
    
    rhoTD = run_params[2]
    rhoDT = run_params[3]
    #rhoTD = math.exp(run_params[1])
    #rhoDT = math.exp(run_params[2])
    #print phiD, rhoTD, rhoDT
    
    birthT=10.75# rate at which adults lay eggs and those eggs hatch
    
    lambdaD = 0.23 # transition from infected to recovered, est. =1 from data but that's not a good ode param, so set to 0.9 instead
    nuD = 0.05 # recovered to susceptible
    kappaD = 0 #recovery w/o immunity
    
    feedL = 0.9 # rate at which L feed and drop off
    feedN = 0.9 # rate at which N feed and drop off
    feedA = 0.5 # rate at which A feed and drop off
    
    densityR = racc_pop # separate param
    densityD = 20 # par persists in the params file, but doesn't show up anywhere in the equations
    densityDT = 100
    densityRT = 100
    
    deathL_on = 0.62
    deathN_on = 0.69
    deathA_on = 0.21
    
    findL=1
    findN=1
    findA=1
    
    enviro_capacity = 200000 #e_cap # separate parameter
    
    deathL_off = 0.02
    deathN_off = 0.01
    deathA_off = 0.001
    death_A_engorged = 0.08
      
    # ----------------------------------------------------------------------------------------------
    # -- define the system of differential equations
    # ----------------------------------------------------------------------------------------------
      
    def derv(x,t):
        y=np.zeros(20);
        deer_pop=x[10]+x[11]+x[12]
        racc_pop_size=densityR
        deer_cap=densityDT*deer_pop
        racc_cap=densityRT*racc_pop_size
        host_pop=deer_pop+racc_pop_size
        
        # global on-host tick K
        tick_on_deer_pop=x[0]+x[1]+x[2]+x[3]+x[4]
        tick_on_racc_pop=x[5]+x[6]+x[7]+x[8]+x[9]
        
        Questing_tick_pop = x[13] + x[14] + x[15] + x[16] + x[17] + x[18] + x[19]
        
        ## ON DEER
        #susceptible ticks
        y[0]=x[13]*findL*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[0]*feedL - x[0]*deathL_on # feeding larvae
        y[1]=x[14]*findN*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[1]*feedN - x[1]*deathN_on # feeding nymphs
        y[2]=x[15]*findA*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[2]*feedA - x[2]*deathA_on # feeding adults 
        #infected ticks
        y[3]=x[17]*findN*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[3]*feedN - x[3]*deathN_on  # feeding inf N
        y[4]=x[18]*findA*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[4]*feedA - x[4]*deathA_on  # feeding inf A
        
        ## ON RACCOONS
        #susceptible
        y[5]=x[13]*findL*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[5]*feedL - x[5]*deathL_on # feeding larvae
        y[6]=x[14]*findN*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[6]*feedN - x[6]*deathN_on # feeding nymphs 
        y[7]=x[15]*findA*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[7]*feedA - x[7]*deathA_on # feeding adults 
        #infected ticks (inf picked up elsewhere)
        y[8]=x[17]*findN*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[8]*feedN - x[8]*deathN_on  # feeding inf N  
        y[9]=x[18]*findA*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[9]*feedA - x[9]*deathA_on  # feeding inf A  
        
        #host 1: SIRS model
        # susceptible:
        y[10]=math.log(1-rhoDT)*x[10]*(x[3]+x[4])+nuD*x[12]+kappaD*x[11]  #susceptible; host density is constant
        # infected:
        y[11]=-math.log(1-rhoDT)*x[10]*(x[3]+x[4])-lambdaD*x[11]-kappaD*x[11] # infected
        # recovered:
        y[12]=lambdaD*x[11]-nuD*x[12] #recovered
        
        ## OFF HOSTS
        #susceptible ticks
        y[13] = (x[16] + x[19])*birthT*((enviro_capacity-Questing_tick_pop)/enviro_capacity) - x[13]*findL*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[13]*findL*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[13]*deathL_off
        y[14] = (x[5]*feedL + x[0]*feedL) - x[0]*feedL*(rhoTD*x[11]/deer_pop) - x[14]*findN*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[14]*findN*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[14]*deathN_off  #engorged larvae + questing nymphs
        y[15] = (x[6]*feedN + x[1]*feedN) - x[1]*feedN*(rhoTD*x[11]/deer_pop) - x[15]*findA*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[15]*findA*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[15]*deathA_off  #engorged nymphs + questing adults
        y[16] = (x[7]*feedA + x[2]*feedA) - x[2]*feedA*(rhoTD*x[11]/deer_pop) - x[16]*death_A_engorged # engorged A death rate ...deathA_off #engorged adults ready to reproduce 
        #infected ticks
        y[17] = x[0]*feedL*(rhoTD*x[11]/deer_pop) - x[17]*findN*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[17]*findN*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[17]*deathN_off # infected engorged larvae and questing nymphs
        y[18] = (x[8]*feedN + x[3]*feedN + x[1]*feedN*(rhoTD*x[11]/deer_pop)) - x[18]*findA*phiD*(deer_pop/host_pop)*((deer_cap-tick_on_deer_pop)/deer_cap) - x[18]*findA*phiR*(racc_pop_size/host_pop)*((racc_cap-tick_on_racc_pop)/racc_cap) - x[18]*deathA_off # infected engorged nymphs and questing adults
        y[19] = (x[9]*feedA + x[4]*feedA + x[2]*feedA*(rhoTD*x[11]/deer_pop)) - x[19]*death_A_engorged # engorged A death rate ...deathA_off # infected engorged adults ready to reproduce. 

        return y
        

    # ----------------------------------------------------------------------------------------------
    # -- integrate
    # ----------------------------------------------------------------------------------------------

    time=np.arange(0,finalTime,0.01)

    r=scipy.integrate.odeint(derv, x, time)

    # ----------------------------------------------------------------------------------------------
    # -- check equilibrium
    # ----------------------------------------------------------------------------------------------

    # -- start a counter
    c = counter
    
    # -- get the last values from the SIRS time series
    end_values=r[-1:][0]
    #print(end_values)
    
    # -- calculate and store the mean population size for each compartment over the last 100 iterations
    means=[] #each value here will correspond to a compartment of the model
    for column in range(len(r[1])):
        ave_last100 = np.mean(r[-100:,column]) #slice out the last 100 values in each column and take their mean
        means.append(ave_last100) #add it to an array
        
    # -- how different is the last value from the time series compared to the mean of the last 100 iterations?
    mean_diffs=means-end_values
    abs_diffs=abs(mean_diffs)
    max_diff=np.amax(abs_diffs)
    
    # -- if the difference is big, run the SIRS model again, starting where it left off
    if max_diff > 0.1:
        if c<=50: # cap number of recursions at 20
            c+=1
            full_out=np.row_stack((R,r))
            
            return SIRprevalences(run_params, full_out, 100, c, num_prefs, racc_pop) #if we're not in equilibrium, run the model another 100 time steps starting with the last values produced by the model
        
        else:
            r = 'failed to reach equilibrium'
            return r

    else: # if we reached equilibrium
        full_out=np.row_stack((R, r))
        if c == 0:
            time_to_eq = finalTime
        else:
            time_to_eq = 100 + (c*100) # 100 for initial run 
        return full_out, time_to_eq

# ----------------------------------------------------------------------------------------------
# -- calculate prevalence for ticks and deer
# ----------------------------------------------------------------------------------------------
    
def calculate_prev(model_output): # take the population sizes and turn them into prevalences. dependent on the ordering of the compartments in the model output
    #print(model_output)
    nymphal_prev=((model_output[3]+model_output[8]+model_output[17])/(model_output[1]+model_output[3]+model_output[6]+model_output[8]+model_output[14]+model_output[17]))
    adult_prev=((model_output[4]+model_output[9]+model_output[18]+model_output[19])/(model_output[2]+model_output[4]+model_output[7]+model_output[9]+model_output[15]+model_output[16]+model_output[18]+model_output[19])) 
    deer_prev=((model_output[11]/(model_output[10]+model_output[11]+model_output[12])))
    deer_AB_prev=(model_output[11]+model_output[12])/(model_output[10]+model_output[11]+model_output[12])
                
    all_prevs=[nymphal_prev, adult_prev, deer_prev, deer_AB_prev]
    
    return all_prevs

# ----------------------------------------------------------------------------------------------
# -- plot the resulting epi curves for single runs
# ----------------------------------------------------------------------------------------------

def mkPlot(ts):

    print 'Plotting'

    time = np.arange(0, len(ts[:,1]))
    
    #plt.figure(figsize=(10,10), dpi=300)
    
    plt.subplot(2,2,1)
    plt.plot(time, ts[:,0], "r", time, ts[:,1], "g", time, ts[:,2], "b", time, ts[:,3], "k", time, ts[:,4], "c")
    #plt.plot(time, (ts[:,0]+ts[:,1]+ts[:,2]+ts[:,3]+ts[:,4]))
    #plt.legend(("S(L)", "S(N)", "S(A)", "I(N)", "I(A)"), loc=0)
    plt.ylabel("Total Number of Ticks on Deer")
    plt.xlabel("Time")
    plt.xticks(rotation=45)
    #title("SIR Model")
    
    plt.subplot(2,2,2)
    #plt.plot(time, ts[:,5], "r", time, ts[:,6], "g", time, ts[:,7], "b", time, ts[:,8], "k", time, ts[:,9], "c")
    plt.plot(time, (ts[:,5]+ts[:,6]+ts[:,7]+ts[:,8]+ts[:,9]))
    #plt.legend(("S(L)", "S(N)", "S(A)", "I(N)", "I(A)"), loc=0)
    plt.ylabel("Total Number of Ticks on Raccoons")
    plt.xlabel("Time")
    plt.xticks(rotation=45)
    #title("SIR Model")
    
    plt.subplot(2,2,3)
    plt.plot(time, ts[:,13], "r", time, ts[:,14], "g", time, ts[:,15], "b", time, ts[:,16], "k", time, ts[:,17], "c", time, ts[:,18], "m", time, ts[:,19], "y")
    #plt.plot(time, (ts[:,14]+ts[:,15]+ts[:,16]+ts[:,17]+ts[:,18]+ts[:,19]+ts[:,20]))
    #plt.legend(("S(L) off", "S(L/N) off", "S(N/A) off", "S(A) off", "I(L/N) off", "I(N/A) off", "I(A) off"), loc=0)
    plt.ylabel("Number of Ticks off Hosts")
    plt.xlabel("Time")
    plt.xticks(rotation=45)
    
    plt.subplot(2,2,4)
    plt.plot(time, ts[:,10], "r", time, ts[:,11], "g", time, ts[:,12], "b")
    #plt.ylim([0,30])
    #plt.legend(("Deer", "Raccoons"), loc=0)
    #plt.legend(("SDeer", "IDeer", "RDeer"), loc=0)
    plt.ylabel("Number of Hosts")
    plt.xlabel("Time")
    plt.xticks(rotation=45)
    
    #fig.set_tight_layout(True)
    plt.show()
    
    #plt.savefig('VD_odeint_200D_200R_Find2.png')
