## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Revised Plotting Scripts for Model Fitting Dissertation Chapter
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

## Where are the data?
##  ~/Desktop/Multivar_Metropolis_200D_200R
##  ~/Desktop/Multivar_Metropolis_400D_200R
##  ~/Desktop/Multivar_Metropolis_200D_400R
##  ~/Desktop has some of the plots
##  ~/Dropbox/ModelFitting/FutureProof/Multivariate_Metropolis_Plots_Copy has copies of the plots

source("/Users/kellypierce/Dropbox/ModelFitting/R Scripts/MCMC_Source.r")
library(scales)
library(coda)
library(magicaxis)

## -----------------------------------------------------------------------------------
## 1. Reorganize data: there are multiple chains for each relative abundance scenario.
##    Convergence diagnostics require the parameter chains to be grouped together for the different replicates
##    That is, we need dataframes that are pref1, pref2, ... instead of pref1, rhoTD1, rhoDT1, ...
## -----------------------------------------------------------------------------------

# this for-loop takes the multiple chains stored in folders 'Multivar_Metropolis...' and resorts the data
# parsed data is written to file
setwd('~/Desktop')
#dirs<-grep('Multivar_Metropolis', list.files(getwd()), value=TRUE)
dirs<-grep('MV_', list.files(getwd()), value=TRUE)
for(dir in dirs){
  #setwd('~/Desktop/BlockUpdate_ChangeKE')
  #setwd('~/Desktop/FinalModelMCMCRuns/SymTrans')
  setwd(file.path('~/Desktop/', dir))
  #chain.dirs<-grep('Run', list.files(getwd()), value=TRUE)
  chain.dirs<-grep('Multivar', list.files(getwd()), value=TRUE)
  #chain.dirs<-grep('rhoTD_chain', list.files(getwd()), value=TRUE)
  pref.chains<-c()
  rhoDT.chains<-c()
  rhoTD.chains<-c()
  pref.chains.burn<-c()
  rhoDT.chains.burn<-c()
  rhoTD.chains.burn<-c()
  loglik.burn<-c()
  aPrev.burn<-c()
  dPrev.burn<-c()
  dABprev.burn<-c()
  for(chain in chain.dirs){ #a bunch of chain 1 files are messed up for some reason...
    #setwd(file.path('~/Desktop/FinalModelMCMCRuns/SymTrans', dir))
    setwd(file.path('~/Desktop/', dir))
    print(chain)
    print(getwd())
    setwd(file.path(getwd(), chain))
    dir.contents<-list.files(getwd())
    print(getwd())
    file<-dir.contents[1]
    data<-read.csv(file)
    
    #sort data without removing burnin (for plotting chain comparisons)
    #thin.file<-thin.chain(data)
    est.file<-estimated.pars(data, 7)
    pref.chains<-cbind(pref.chains, est.file[,1])
    rhoTD.chains<-cbind(rhoTD.chains, est.file[,2])
    rhoDT.chains<-cbind(rhoDT.chains, est.file[,3])
    
    #sort data with 20k iteration 'burnin' removed (for convergence diagnostics)
    #thin.file.burn<-thin.chain(data[20000:100000,])
    est.file.burn<-data[20000:120000,]
    pref.chains.burn<-cbind(pref.chains.burn, est.file.burn[,1])
    rhoTD.chains.burn<-cbind(rhoTD.chains.burn, est.file.burn[,2])
    rhoDT.chains.burn<-cbind(rhoDT.chains.burn, est.file.burn[,3])
    loglik.burn<-cbind(loglik.burn, est.file.burn[,28])
    aPrev.burn<-c(aPrev.burn, est.file.burn[,25])
    dPrev.burn<-c(dPrev.burn, est.file.burn[,26])
    dABprev.burn<-c(dABprev.burn, est.file.burn[,27])
  }
  #write the sorted data to file
  setwd(file.path('~/Desktop/', dir))
  write.csv(pref.chains, file='preference_multi_chains.csv')
  write.csv(rhoTD.chains, file='rhoTD_multi_chains.csv')
  write.csv(rhoDT.chains, file='rhoDT_multi_chains.csv')
  write.csv(pref.chains.burn, file='preference_burned_multi_chains.csv')
  write.csv(rhoTD.chains.burn, file='rhoTD_burned_multi_chains.csv')
  write.csv(rhoDT.chains.burn, file='rhoDT_burned_multi_chains.csv')
  write.csv(loglik.burn, file='loglik_burned_multi_chains.csv')
  write.csv(aPrev.burn, file='aPrev_burned_multi_chains.csv')
  write.csv(dPrev.burn, file='dPrev_burned_multi_chains.csv')
  write.csv(dABprev.burn, file='dABprev_burned_multi_chains.csv')
}

## -----------------------------------------------------------------------------------
## 2. Plot chains and perform convergence diagnostics on the re-organized data
## -----------------------------------------------------------------------------------

# read in replicate chains, coerce to MCMC objects, and perform the gelman-rubin convergence test
#pdf(file='~/Desktop/Multivariate_Metropolis_MultpleChains.pdf', height=5, width=7)
for(dir in dirs){
  # change to appropriate directory
  setwd('~/Desktop/')
  setwd(file.path(getwd(), dir))
  print(dir)
  
  # read data, select appropriate columns, exponentiate, and set column names
  pref<-read.csv("preference_multi_chains.csv")
  pref<-exp(pref[,2:length(pref[1,])])
  names(pref)<-rep('pref', length(pref[1,]))
  
  rhoDT<-read.csv("rhoDT_multi_chains.csv")
  rhoDT<-exp(rhoDT[,2:length(rhoDT[1,])])
  names(rhoDT)<-rep('rhoDT', length(rhoDT[1,]))
  
  rhoTD<-read.csv("rhoTD_multi_chains.csv")
  rhoTD<-exp(rhoTD[,2:length(rhoTD[1,])])
  names(rhoTD)<-rep('rhoTD', length(rhoTD[1,]))
  
  # plot the chains
  line.colors<-c('violetred3', 'salmon4', 'blue1', 'paleturquoise3', 'yellow3')
  par(mar=c(4,4,4,2))
  
  ## preference ##
  png(file=file.path('~/Desktop/',paste(dir,'Chains.png')), height=30, width=30, units='cm', res=300)
  par(mfrow=c(3,1), mar=c(5,6,2,2))
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  plot(1:length(pref[,1]), pref[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
       ylab='', las=1)
  mtext(expression(phi['D']), side=2, line=4, cex=1.5)
  for(i in 2:length(pref[1,])){
    lines(1:length(pref[,i]), pref[,i], col=line.colors[i])
  }
  abline(v=10000, lty=2, lwd=3)
  #A=hist(pref[,1], plot=FALSE, breaks=100)
  #plot(y=A$mids, x=A$counts, type='l', lwd=2, col=alpha(line.colors[i]), 
  #     xlim=c(0,max(A$counts)), ylim=c(0,1), xlab='Counts', yaxt='n', ylab="")
  #for(i in 2:length(pref[1,])){
  #  B=hist(pref[,i], plot=FALSE, breaks=100)
  #  lines(y=B$mids, x=B$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  
  ## rhoTD ##
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  plot(1:length(rhoTD[,1]), rhoTD[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
       xlab='Iteration',ylab='')
  mtext(expression(rho['TD']), side=2, line=4, cex=1.5)
  for(i in 2:length(rhoTD[1,])){
    lines(1:length(rhoTD[,i]), rhoTD[,i], col=line.colors[i])
  }
  abline(v=10000, lty=2, lwd=3)
  #C=hist(rhoTD[,1], plot=FALSE, breaks=25)
  #plot(y=C$mids, x=C$counts, type='l', lwd=2, col=alpha(line.colors[1]), 
  #     xlim=c(0, max(C$counts)), ylim=c(0,1), yaxt='n', xlab='Counts')
  #for(i in 2:length(rhoTD[1,])){
  #  D=hist(rhoTD[,i], plot=FALSE, breaks=25)
  #  lines(y=D$mids, x=D$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  
  ## rhoDT ##
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  plot(1:length(rhoDT[,1]), rhoDT[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
       xlab='Iteration',ylab='')
  mtext(expression(rho['DT']), side=2, line=4, cex=1.5)
  for(i in 2:length(rhoDT[1,])){
    lines(1:length(rhoDT[,i]), rhoDT[,i], col=line.colors[i])
  }
  abline(v=10000, lty=2, lwd=3)
  #E=hist(rhoDT[,1], plot=FALSE, breaks=30)
  #plot(y=E$mids, x=E$counts, type='l', lwd=2, col=alpha(line.colors[1]), 
  #     xlim=c(0, max(E$counts)), ylim=c(0,1), yaxt='n', xlab='Counts')
  #for(i in 2:length(rhoDT[1,])){
  #  G=hist(rhoTD[,i], plot=FALSE, breaks=30)
  #  lines(y=G$mids, x=G$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  dev.off()
  # convert chains into mcmc objects and create mcmc lists
  rhoDT.mcmc<-list()
  rhoTD.mcmc<-list()
  pref.mcmc<-list()
  for(i in 1:length(pref[1,])){
    pref.mcmc[[i]]<-as.mcmc(pref[20000:120000,i])
    rhoDT.mcmc[[i]]<-as.mcmc(rhoDT[20000:120000,i])
    rhoTD.mcmc[[i]]<-as.mcmc(rhoTD[20000:120000,i])
  }
  pref.mcmc<-mcmc.list(pref.mcmc)
  rhoDT.mcmc<-mcmc.list(rhoDT.mcmc)
  rhoTD.mcmc<-mcmc.list(rhoTD.mcmc)
  
  # run the gelman-rubin diagnostic  
  print(gelman.diag(pref.mcmc))
  print(gelman.diag(rhoDT.mcmc))
  print(gelman.diag(rhoTD.mcmc))
  gelman.plot(pref.mcmc, main=paste('preference', dir))
  gelman.plot(rhoDT.mcmc, main=paste('rhoDT', dir))
  gelman.plot(rhoTD.mcmc, main=paste('rhoTD', dir))
}


## -----------------------------------------------------------------------------------
## 3. Look at marginal posterior probabilities for the estimated parameters
## -----------------------------------------------------------------------------------

## Preference Posterior

#pref2.2<-read.csv('~/Desktop/Multivar_Metropolis_200D_200R/preference_burned_multi_chains.csv')
#pref2.4<-read.csv('~/Desktop/Multivar_Metropolis_200D_400R/preference_burned_multi_chains.csv')
#pref4.2<-read.csv('~/Desktop/Multivar_Metropolis_400D_200R/preference_burned_multi_chains.csv')

pref2.2<-read.csv('~/Desktop/MV_200_200_fixed/preference_burned_multi_chains.csv')
pref2.4<-read.csv('~/Desktop/MV_200_400_fixed/preference_burned_multi_chains.csv')
pref4.2<-read.csv('~/Desktop/MV_400_200_fixed/preference_burned_multi_chains.csv')

pref2.2.hist<-hist(exp(pref2.2[,2]), breaks=70, plot=FALSE)
pref2.4.hist<-hist(exp(pref2.4[,2]), breaks=70, plot=FALSE)
pref4.2.hist<-hist(exp(pref4.2[,2]), breaks=70, plot=FALSE)

pref2.2.95<-quantile(exp(pref2.2[,2]), probs=c(0.025, 0.975))
pref2.4.95<-quantile(exp(pref2.4[,2]), probs=c(0.025, 0.975))
pref4.2.95<-quantile(exp(pref4.2[,2]), probs=c(0.025, 0.975))

pref2.2.95.hist.data<-cbind(pref2.2.hist$mids, pref2.2.hist$density)
pref2.2.95.hist<-subset(pref2.2.95.hist.data, pref2.2.95.hist.data[,1] > pref2.2.95[1] & pref2.2.95.hist.data[,1] < pref2.2.95[2])
pref2.4.95.hist.data<-cbind(pref2.4.hist$mids, pref2.4.hist$density)
pref2.4.95.hist<-subset(pref2.4.95.hist.data, pref2.4.95.hist.data[,1] > pref2.4.95[1] & pref2.4.95.hist.data[,1] < pref2.4.95[2])
pref4.2.95.hist.data<-cbind(pref4.2.hist$mids, pref4.2.hist$density)
pref4.2.95.hist<-subset(pref4.2.95.hist.data, pref4.2.95.hist.data[,1] > pref4.2.95[1] & pref4.2.95.hist.data[,1] < pref2.2.95[2])



## rhoTD posterior

rhoTD2.2<-read.csv('~/Desktop/MV_200_200_fixed/rhoTD_burned_multi_chains.csv')
rhoTD2.4<-read.csv('~/Desktop/MV_200_400_fixed/rhoTD_burned_multi_chains.csv')
rhoTD4.2<-read.csv('~/Desktop/MV_400_200_fixed/rhoTD_burned_multi_chains.csv')
rhoTD2.2.hist<-hist(exp(rhoTD2.2[,2]), breaks=25, plot=FALSE)
rhoTD2.4.hist<-hist(exp(rhoTD2.4[,2]), breaks=25, plot=FALSE)
rhoTD4.2.hist<-hist(exp(rhoTD4.2[,2]), breaks=25, plot=FALSE)

rhoTD2.2.95<-quantile(exp(rhoTD2.2[,2]), probs=c(0.025, 0.975))
rhoTD2.4.95<-quantile(exp(rhoTD2.4[,2]), probs=c(0.025, 0.975))
rhoTD4.2.95<-quantile(exp(rhoTD4.2[,2]), probs=c(0.025, 0.975))

rhoTD2.2.95.hist.data<-cbind(rhoTD2.2.hist$mids, rhoTD2.2.hist$density)
rhoTD2.2.95.hist<-subset(rhoTD2.2.95.hist.data, rhoTD2.2.95.hist.data[,1] > rhoTD2.2.95[1] & rhoTD2.2.95.hist.data[,1] < rhoTD2.2.95[2])
rhoTD2.4.95.hist.data<-cbind(rhoTD2.4.hist$mids, rhoTD2.4.hist$density)
rhoTD2.4.95.hist<-subset(rhoTD2.4.95.hist.data, rhoTD2.4.95.hist.data[,1] > rhoTD2.4.95[1] & rhoTD2.4.95.hist.data[,1] < rhoTD2.4.95[2])
rhoTD4.2.95.hist.data<-cbind(rhoTD4.2.hist$mids, rhoTD4.2.hist$density)
rhoTD4.2.95.hist<-subset(rhoTD4.2.95.hist.data, rhoTD4.2.95.hist.data[,1] > rhoTD4.2.95[1] & rhoTD4.2.95.hist.data[,1] < rhoTD2.2.95[2])


## rhoDT posterior

rhoDT2.2<-read.csv('~/Desktop/MV_200_200_fixed/rhoDT_burned_multi_chains.csv')
rhoDT2.4<-read.csv('~/Desktop/MV_200_400_fixed/rhoDT_burned_multi_chains.csv')
rhoDT4.2<-read.csv('~/Desktop/MV_400_200_fixed/rhoDT_burned_multi_chains.csv')
rhoDT2.2.hist<-hist(exp(rhoDT2.2[,2]), breaks=30, plot=FALSE)
rhoDT2.4.hist<-hist(exp(rhoDT2.4[,2]), breaks=30, plot=FALSE)
rhoDT4.2.hist<-hist(exp(rhoDT4.2[,2]), breaks=30, plot=FALSE)

rhoDT2.2.95<-quantile(exp(rhoDT2.2[,2]), probs=c(0.025, 0.975))
rhoDT2.4.95<-quantile(exp(rhoDT2.4[,2]), probs=c(0.025, 0.975))
rhoDT4.2.95<-quantile(exp(rhoDT4.2[,2]), probs=c(0.025, 0.975))

rhoDT2.2.95.hist.data<-cbind(rhoDT2.2.hist$mids, rhoDT2.2.hist$density)
rhoDT2.2.95.hist<-subset(rhoDT2.2.95.hist.data, rhoDT2.2.95.hist.data[,1] > rhoDT2.2.95[1] & rhoDT2.2.95.hist.data[,1] < rhoDT2.2.95[2])
rhoDT2.4.95.hist.data<-cbind(rhoDT2.4.hist$mids, rhoDT2.4.hist$density)
rhoDT2.4.95.hist<-subset(rhoDT2.4.95.hist.data, rhoDT2.4.95.hist.data[,1] > rhoDT2.4.95[1] & rhoDT2.4.95.hist.data[,1] < rhoDT2.4.95[2])
rhoDT4.2.95.hist.data<-cbind(rhoDT4.2.hist$mids, rhoDT4.2.hist$density)
rhoDT4.2.95.hist<-subset(rhoDT4.2.95.hist.data, rhoDT4.2.95.hist.data[,1] > rhoDT4.2.95[1] & rhoDT4.2.95.hist.data[,1] < rhoDT2.2.95[2])


## Marginal Posterior Plots - Overlapping curves ##

png(file='~/Desktop/Multivariate_Metropolis_MarginalPosteriors_RoughDraft.png', height=10, width=30, units='cm', res=300)
par(mfrow=c(1,3))
plot(pref2.2.hist$mids, pref2.2.hist$counts, type='l', 
     ylim=c(0,max(c(pref2.2.hist$counts, pref2.4.hist$counts, pref4.2.hist$counts))),
     xlab='Posterior Preference', ylab='', las=1)
lines(pref2.4.hist$mids, pref2.4.hist$counts, type='l')
lines(pref4.2.hist$mids, pref4.2.hist$counts, type='l')
polygon(pref2.2.hist$mids, pref2.2.hist$counts, col=alpha('blue', 0.25))
polygon(pref2.4.hist$mids, pref2.4.hist$counts, col=alpha('yellow', 0.25))
polygon(x=c(0,pref4.2.hist$mids), y=c(0,pref4.2.hist$counts), col=alpha('red', 0.25))
legend(x='topright', legend=c('1:2', '1:1', '2:1'), 
       fill=c(alpha('yellow', 0.25), alpha('blue', 0.25), alpha('red', 0.25)),
       bty='n')

plot(rhoTD2.2.hist$mids, rhoTD2.2.hist$counts, type='l', 
     ylim=c(0,max(c(rhoTD2.2.hist$counts, rhoTD2.4.hist$counts, rhoTD4.2.hist$counts))),
     xlab='Posterior rhoTD', ylab='', las=1)
lines(rhoTD2.4.hist$mids, rhoTD2.4.hist$counts, type='l')
lines(rhoTD4.2.hist$mids, rhoTD4.2.hist$counts, type='l')
polygon(x=c(0,rhoTD2.2.hist$mids, max(rhoTD2.2.hist$mids)), y=c(0,rhoTD2.2.hist$counts,0), col=alpha('blue', 0.25))
polygon(x=c(0,rhoTD2.4.hist$mids, max(rhoTD2.4.hist$mids)), y=c(0,rhoTD2.4.hist$counts,0), col=alpha('yellow', 0.25))
polygon(x=c(0,rhoTD4.2.hist$mids, max(rhoTD4.2.hist$mids)), y=c(0,rhoTD4.2.hist$counts,0), col=alpha('red', 0.25))
legend(x='topright', legend=c('1:2', '1:1', '2:1'), 
       fill=c(alpha('yellow', 0.25), alpha('blue', 0.25), alpha('red', 0.25)),
       bty='n')

plot(rhoDT2.2.hist$mids, rhoDT2.2.hist$counts, type='l', 
     ylim=c(0,max(c(rhoDT2.2.hist$counts, rhoDT2.4.hist$counts, rhoDT4.2.hist$counts))),
     xlab='Posterior rhoDT', ylab='', las=1)
lines(rhoDT2.4.hist$mids, rhoDT2.4.hist$counts, type='l')
lines(rhoDT4.2.hist$mids, rhoDT4.2.hist$counts, type='l')
polygon(x=c(0,rhoDT2.2.hist$mids, max(rhoDT2.2.hist$mids)), y=c(0,rhoDT2.2.hist$counts,0), col=alpha('blue', 0.25))
polygon(x=c(0,rhoDT2.4.hist$mids, max(rhoDT2.4.hist$mids)), y=c(0,rhoDT2.4.hist$counts,0), col=alpha('yellow', 0.25))
polygon(x=c(0,rhoDT4.2.hist$mids, max(rhoDT4.2.hist$mids)), y=c(0,rhoDT4.2.hist$counts,0), col=alpha('red', 0.25))
legend(x='topright', legend=c('1:2', '1:1', '2:1'), 
       fill=c(alpha('yellow', 0.25), alpha('blue', 0.25), alpha('red', 0.25)),
       bty='n')

dev.off()

## Marginal posterior plots -- 3x3 matrix w/95%CI shading

png(file='~/Desktop/MV_DensityDT_Fixed_MarginalPosteriors_3x3Matrix.png', height=30, width=30, units='cm', res=300)
par(mfcol=c(3,3), mar=c(5,4,2,2))

# preference plots
plot(pref2.4.hist$mids, pref2.4.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
mtext(text='1:2', cex=1.5, side=2, line=1)
axis(side=1, cex.axis=2)
polygon(c(min(pref2.4.hist$mids), pref2.4.hist$mids, max(pref2.4.hist$mids)), 
        c(0,pref2.4.hist$density,0), col=alpha('darkolivegreen4', 0.5))
polygon(c(min(pref2.4.95.hist[,1]), pref2.4.95.hist[,1], max(pref2.4.95.hist[,1])), 
        c(0,pref2.4.95.hist[,2],0), col=alpha('darkolivegreen4', 0.5))
abline(v=0.5, lty=2, lwd=3, col='red')

plot(pref2.2.hist$mids, pref2.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
mtext(text='1:1', cex=1.5, side=2, line=1)
axis(side=1, cex.axis=2)
polygon(c(min(pref2.2.hist$mids), pref2.2.hist$mids, max(pref2.2.hist$mids)), 
        c(0,pref2.2.hist$density,0), col=alpha('darkolivegreen4', 0.5))
polygon(c(min(pref2.2.95.hist[,1]), pref2.2.95.hist[,1], max(pref2.2.95.hist[,1])), 
        c(0,pref2.2.95.hist[,2],0), col=alpha('darkolivegreen4', 0.5))
abline(v=0.5, lty=2, lwd=3, col='red')

plot(pref4.2.hist$mids, pref4.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
mtext(text='2:1', cex=1.5, side=2, line=1)
axis(side=1, cex.axis=2)
polygon(c(min(pref4.2.hist$mids), pref4.2.hist$mids, max(pref4.2.hist$mids)), 
        c(0,pref4.2.hist$density,0), col=alpha('darkolivegreen4', 0.5))
polygon(c(min(pref4.2.95.hist[,1]), pref4.2.95.hist[,1], max(pref4.2.95.hist[,1])), 
        c(0,pref4.2.95.hist[,2],0), col=alpha('darkolivegreen4', 0.5))
abline(v=0.5, lty=2, lwd=3, col='red')
mtext(text=expression(phi['D']), cex=1.5, line=4, side=1)

# rhoDT plots

plot(rhoDT2.4.hist$mids, rhoDT2.4.hist$density, type='l', axes=FALSE, main='', xlim=c(0,0.4),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoDT2.4.hist$mids), rhoDT2.4.hist$mids, max(rhoDT2.4.hist$mids)), 
        c(0,rhoDT2.4.hist$density,0), col=alpha('darkorange2', 0.5))
polygon(c(min(rhoDT2.4.95.hist[,1]), rhoDT2.4.95.hist[,1], max(rhoDT2.4.95.hist[,1])), 
        c(0,rhoDT2.4.95.hist[,2],0), col=alpha('darkorange2', 0.5))
abline(v=0.26, lty=2, lwd=3, col='black')

plot(rhoDT2.2.hist$mids, rhoDT2.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,0.4),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoDT2.2.hist$mids), rhoDT2.2.hist$mids, max(rhoDT2.2.hist$mids)), 
        c(0,rhoDT2.2.hist$density,0), col=alpha('darkorange2', 0.5))
polygon(c(min(rhoDT2.2.95.hist[,1]), rhoDT2.2.95.hist[,1], max(rhoDT2.2.95.hist[,1])), 
        c(0,rhoDT2.2.95.hist[,2],0), col=alpha('darkorange2', 0.5))
abline(v=0.26, lty=2, lwd=3, col='black')

plot(rhoDT4.2.hist$mids, rhoDT4.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,0.4),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoDT4.2.hist$mids), rhoDT4.2.hist$mids, max(rhoDT4.2.hist$mids)), 
        c(0,rhoDT4.2.hist$density,0), col=alpha('darkorange2', 0.5))
polygon(c(min(rhoDT4.2.95.hist[,1]), rhoDT4.2.95.hist[,1], max(rhoDT4.2.95.hist[,1])), 
        c(0,rhoDT4.2.95.hist[,2],0), col=alpha('darkorange2', 0.5))
abline(v=0.26, lty=2, lwd=3, col='black')
mtext(text=expression(rho['DT']), cex=1.5, line=4, side=1)

# rhoTD plots

plot(rhoTD2.4.hist$mids, rhoTD2.4.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoTD2.4.hist$mids), rhoTD2.4.hist$mids, max(rhoTD2.4.hist$mids)), 
        c(0,rhoTD2.4.hist$density,0), col=alpha('orchid4', 0.5))
polygon(c(min(rhoTD2.4.95.hist[,1]), rhoTD2.4.95.hist[,1], max(rhoTD2.4.95.hist[,1])), 
        c(0,rhoTD2.4.95.hist[,2],0), col=alpha('orchid4', 0.5))
abline(v=0.02, lty=2, lwd=3, col='black')

plot(rhoTD2.2.hist$mids, rhoTD2.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoTD2.2.hist$mids), rhoTD2.2.hist$mids, max(rhoTD2.2.hist$mids)), 
        c(0,rhoTD2.2.hist$density,0), col=alpha('orchid4', 0.5))
polygon(c(min(rhoTD2.2.95.hist[,1]), rhoTD2.2.95.hist[,1], max(rhoTD2.2.95.hist[,1])), 
        c(0,rhoTD2.2.95.hist[,2],0), col=alpha('orchid4', 0.5))
abline(v=0.06, lty=2, lwd=3, col='black')

plot(rhoTD4.2.hist$mids, rhoTD4.2.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
axis(side=1, cex.axis=2)
polygon(c(min(rhoTD4.2.hist$mids), rhoTD4.2.hist$mids, max(rhoTD4.2.hist$mids)), 
        c(0,rhoTD4.2.hist$density,0), col=alpha('orchid4', 0.5))
polygon(c(min(rhoTD4.2.95.hist[,1]), rhoTD4.2.95.hist[,1], max(rhoTD4.2.95.hist[,1])), 
        c(0,rhoTD4.2.95.hist[,2],0), col=alpha('orchid4', 0.5))
abline(v=0.06, lty=2, lwd=3, col='black')
mtext(text=expression(rho['TD']), cex=1.5, line=4, side=1)

dev.off()

## split chains and look at mean 1st half vs. mean 2nd half
pref2.2.means1<-colMeans(pref2.2[1:50000,2:6])
pref2.2.means2<-colMeans(pref2.2[50001:100000,2:6])
pref2.2.means1-pref2.2.means2
pref2.4.means1<-colMeans(pref2.4[1:50000,2:6])
pref2.4.means2<-colMeans(pref2.4[50001:100000,2:6])
pref2.4.means1-pref2.4.means2
pref4.2.means1<-colMeans(pref4.2[1:50000,2:6])
pref4.2.means2<-colMeans(pref4.2[50001:100000,2:6])
pref4.2.means1-pref4.2.means2

rhoTD2.2.means1<-colMeans(rhoTD2.2[1:50000,2:6])
rhoTD2.2.means2<-colMeans(rhoTD2.2[50001:100000,2:6])
rhoTD2.2.means1-rhoTD2.2.means2
rhoTD2.4.means1<-colMeans(rhoTD2.4[1:50000,2:6])
rhoTD2.4.means2<-colMeans(rhoTD2.4[50001:100000,2:6])
rhoTD2.4.means1-rhoTD2.4.means2
rhoTD4.2.means1<-colMeans(rhoTD4.2[1:50000,2:6])
rhoTD4.2.means2<-colMeans(rhoTD4.2[50001:100000,2:6])
rhoTD4.2.means1-rhoTD4.2.means2

rhoDT2.2.means1<-colMeans(rhoDT2.2[1:50000,2:6])
rhoDT2.2.means2<-colMeans(rhoDT2.2[50001:100000,2:6])
rhoDT2.2.means1-rhoDT2.2.means2
rhoDT2.4.means1<-colMeans(rhoDT2.4[1:50000,2:6])
rhoDT2.4.means2<-colMeans(rhoDT2.4[50001:100000,2:6])
rhoDT2.4.means1-rhoDT2.4.means2
rhoDT4.2.means1<-colMeans(rhoDT4.2[1:50000,2:6])
rhoDT4.2.means2<-colMeans(rhoDT4.2[50001:100000,2:6])
rhoDT4.2.means1-rhoDT4.2.means2

## -----------------------------------------------------------------------------------
## 4. Look at parameter correlations in the different relative abundance scenarios by
##    making a series of pairs plots (1 master plot with 9 panels) EXPONENTIATED
## -----------------------------------------------------------------------------------

#these are all "chain1", which is the same as prefX.X[,2], rhoTDX.X[,2], rhoDTX.X[,2]
#the full chains have the posteriors
chain2.2<-read.csv('~/Desktop/MV_200_200_fixed/Multivar_Metropolis_200D_200R_1/1 MCMC Accepted Iteration Log 2014-10-21 22-38 486.csv')
chain2.4<-read.csv('~/Desktop/MV_200_400_fixed/Multivar_Metropolis_200D_400R_1/1 MCMC Accepted Iteration Log 2014-10-21 22-38 293.csv')
chain4.2<-read.csv('~/Desktop/MV_400_200_fixed/Multivar_Metropolis_400D_200R_1/1 MCMC Accepted Iteration Log 2014-10-21 22-38 887.csv')
chain2.2.burn<-chain2.2[20000:length(chain2.2[,1]),]
chain2.4.burn<-chain2.4[20000:length(chain2.4[,1]),]
chain4.2.burn<-chain4.2[20000:length(chain4.2[,1]),]
post2.4colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(chain2.4.burn$posterior)))))
post2.2colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(chain2.2.burn$posterior)))))
post4.2colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(chain4.2.burn$posterior)))))

png(file='~/Desktop/DensityDT_Fixed_PairsPlotLegend.png', height=10, width=5, units='cm', res=300)
plot.new()
legend_image <- as.raster(matrix(post2.2colors, ncol=1))
rasterImage(legend_image, 0, 0, 1,1)
mtext(text = c('poor'), side=2, line=0, at=c(0), las=1, cex=1.5)
mtext(text = c('good'), side=2, line=0, at=c(1), las=1, cex=1.5)
mtext(text = c('Model Fit'), side=2, line=2, cex=1.5)
dev.off()

png(file='~/Desktop/DensityDT_Fixed_Multivariate_Metropolis_ExpPairs_Colored.png', height=30, width=30, units='cm', res=300)
par(mfrow=c(3,3), mar=c(5,7,2,2))

# row 1 (1:2 rel. abundance)
plot(exp(pref2.4[,2]), exp(rhoTD2.4[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['TD']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.4colors[abs(round(chain2.4.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
mtext(side=2, text='1:2', cex=1.5, line=5)
#title('A', adj=0)
plot(exp(pref2.4[,2]), exp(rhoDT2.4[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.4colors[abs(round(chain2.4.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('B', adj=0)
plot(exp(rhoTD2.4[,2]), exp(rhoDT2.4[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(rho['TD']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.4colors[abs(round(chain2.4.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('C', adj=0)
# row 2 (1:1 rel. abundance)
plot(exp(pref2.2[,2]), exp(rhoTD2.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['TD']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.2colors[abs(round(chain2.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
mtext(side=2, text='1:1', cex=1.5, line=5)
#title('D', adj=0)
plot(exp(pref2.2[,2]), exp(rhoDT2.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.2colors[abs(round(chain2.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('E', adj=0)
plot(exp(rhoTD2.2[,2]), exp(rhoDT2.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(rho['TD']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post2.2colors[abs(round(chain2.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('F', adj=0)
# row 3 (2:1 rel. abundance)
plot(exp(pref4.2[,2]), exp(rhoTD4.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['TD']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post4.2colors[abs(round(chain4.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
mtext(side=2, text='2:1', cex=1.5, line=5)
#title('G', adj=0)
plot(exp(pref4.2[,2]), exp(rhoDT4.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post4.2colors[abs(round(chain4.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('H', adj=0)
plot(exp(rhoTD4.2[,2]), exp(rhoDT4.2[,2]), xlim=c(0,1), ylim=c(0,1), axes=FALSE,
     xlab=expression(rho['TD']), ylab=expression(rho['DT']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post4.2colors[abs(round(chain4.2.burn$posterior))])
magaxis(side=c(1,2), las=1, cex.axis=1.5)
#title('I', adj=0)
dev.off()

## -----------------------------------------------------------------------------------
## 5. Meta-analysis empirical prevalence posteriors compared to model output prev
## -----------------------------------------------------------------------------------

# read in the data generated in the ticks.R script
empirical<-read.csv('~/Dropbox/ModelFitting/FutureProof/TickPrevGLMM_MetaAnalysis_NoRacc_Results_LogitTransform.csv')
ilogit = function(x) 1/{1+exp(-x)}

ap2.2<-read.csv('~/Desktop/Multivar_Metropolis_200D_200R/aPrev_burned_multi_chains.csv')
ap2.4<-read.csv('~/Desktop/Multivar_Metropolis_200D_400R/aPrev_burned_multi_chains.csv')
ap4.2<-read.csv('~/Desktop/Multivar_Metropolis_400D_200R/aPrev_burned_multi_chains.csv')
ap2.2.hist<-hist(ap2.2[,2], breaks=25, plot=FALSE)
ap2.4.hist<-hist(ap2.4[,2], breaks=25, plot=FALSE)
ap4.2.hist<-hist(ap4.2[,2], breaks=25, plot=FALSE)
empirical.ap.hist<-hist(ilogit(empirical[,2]), breaks=25, plot=FALSE)
ymax<-max(c(ap2.2.hist$density, ap2.4.hist$density, ap4.2.hist$density, empirical.ap.hist$density))

dp2.2<-read.csv('~/Desktop/Multivar_Metropolis_200D_200R/dPrev_burned_multi_chains.csv')
dp2.4<-read.csv('~/Desktop/Multivar_Metropolis_200D_400R/dPrev_burned_multi_chains.csv')
dp4.2<-read.csv('~/Desktop/Multivar_Metropolis_400D_200R/dPrev_burned_multi_chains.csv')
dp2.2.hist<-hist(dp2.2[,2], breaks=25, plot=FALSE)
dp2.4.hist<-hist(dp2.4[,2], breaks=25, plot=FALSE)
dp4.2.hist<-hist(dp4.2[,2], breaks=25, plot=FALSE)
empirical.dp.hist<-hist(ilogit(empirical[,3]), breaks=25, plot=FALSE)
ymax.d<-max(c(dp2.2.hist$density, dp2.4.hist$density, dp4.2.hist$density, empirical.dp.hist$density))

dabp2.2<-read.csv('~/Desktop/Multivar_Metropolis_200D_200R/dABprev_burned_multi_chains.csv')
dabp2.4<-read.csv('~/Desktop/Multivar_Metropolis_200D_400R/dABprev_burned_multi_chains.csv')
dabp4.2<-read.csv('~/Desktop/Multivar_Metropolis_400D_200R/dABprev_burned_multi_chains.csv')
dabp2.2.hist<-hist(dabp2.2[,2], breaks=25, plot=FALSE)
dabp2.4.hist<-hist(dabp2.4[,2], breaks=25 plot=FALSE)
dabp4.2.hist<-hist(dabp4.2[,2], breaks=25 plot=FALSE)

empirical.dabp.hist<-hist(ilogit(empirical[,4]), breaks=25, plot=FALSE)
ymax.dab<-max(c(dabp2.2.hist$density, dabp2.4.hist$density, dabp4.2.hist$density, empirical.dabp.hist$density))

png(file='~/Desktop/Multivar_Metropolis_Prevalence_Comparison.png', height=10, width=30, units='cm', res=300)
par(mfrow=c(1,3), mar=c(8,2,4,2))
plot(empirical.ap.hist$mids, empirical.ap.hist$density, ylim=c(0,ymax), type='l', 
     axes=FALSE, xlab='', ylab='', las=1, xlim=c(0,0.15))
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nAdult Ticks', side=1, line=6, cex=1.5)
polygon(empirical.ap.hist$mids, empirical.ap.hist$density, col=alpha('violetred3',0.25))
polygon(ap2.2.hist$mids, ap2.2.hist$density, col=alpha('violetred3',0.25))
polygon(ap2.4.hist$mids, ap2.4.hist$density, col=alpha('violetred3',0.25))
polygon(ap4.2.hist$mids, ap4.2.hist$density, col=alpha('violetred3',0.25))
title('A', adj=0, cex.main=3)

plot(c(min(empirical.dp.hist$mids), empirical.dp.hist$mids, max(empirical.dp.hist$mids)), 
     c(0,empirical.dp.hist$density,0), ylim=c(0,ymax.d), type='l', xlim=c(0,1),
     axes=FALSE, xlab='', ylab='', las=1)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(min(empirical.dp.hist$mids), empirical.dp.hist$mids, max(empirical.dp.hist$mids)),
        c(0,empirical.dp.hist$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp2.2.hist$mids),dp2.2.hist$mids, max(dp2.2.hist$mids)), 
        c(0,dp2.2.hist$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp2.4.hist$mids),dp2.4.hist$mids,max(dp2.4.hist$mids)), 
        c(0,dp2.4.hist$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp4.2.hist$mids),dp4.2.hist$mids,max(dp4.2.hist$mids)),
        c(0,dp4.2.hist$density,0), col=alpha('blue1',0.25))
title('B', adj=0, cex.main=3)

plot(c(min(empirical.dabp.hist$mids), empirical.dabp.hist$mids, max(empirical.dabp.hist$mids)), 
     c(0,empirical.dabp.hist$density,0), ylim=c(0,ymax.dab), type='l', xlim=c(0,1),
     axes=FALSE, xlab='', ylab='', las=1)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Seroprevalence \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(min(empirical.dabp.hist$mids), empirical.dabp.hist$mids, max(empirical.dabp.hist$mids)),
        c(0,empirical.dabp.hist$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp2.2.hist$mids),dabp2.2.hist$mids, max(dabp2.2.hist$mids)), 
        c(0,dabp2.2.hist$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp2.4.hist$mids),dabp2.4.hist$mids,max(dabp2.4.hist$mids)), 
        c(0,dabp2.4.hist$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp4.2.hist$mids),dabp4.2.hist$mids,max(dabp4.2.hist$mids)),
        c(0,dabp4.2.hist$density,0), col=alpha('yellow3',0.25))
title('C', adj=0, cex.main=3)

dev.off()

## -----------------------------------------------------------------------------------
## 6. Meta-analysis distributions w/original empirical data
## -----------------------------------------------------------------------------------

ticks = read.csv("~/Dropbox/ModelFitting/FutureProof/CompletePrevsFixed_2.csv", header=FALSE)
names(ticks) = c("species", "ncases", "ntrials")
ticks = ticks[order(ticks$species),]

tick.prev<-cbind(ticks, (ticks$ncases/ticks$ntrials))
names(tick.prev)<-c('species', 'ncases', 'ntrials', 'prev')

ilogit = function(x) 1/{1+exp(-x)}
line.colors<-c('violetred3', 'salmon4', 'blue1', 'paleturquoise3', 'yellow3')

empirical<-read.csv('~/Dropbox/ModelFitting/FutureProof/TickPrevGLMM_MetaAnalysis_NoRacc_Results_LogitTransform.csv')
empirical.ap.hist<-hist(ilogit(empirical[,2]), breaks=40, plot=FALSE)
empirical.dp.hist<-hist(ilogit(empirical[,3]), breaks=30, plot=FALSE)
empirical.dabp.hist<-hist(ilogit(empirical[,4]), breaks=28, plot=FALSE)

ap.95<-quantile(ilogit(empirical[,2]), probs=c(0.025, 0.975))
ap.hist.data<-cbind(empirical.ap.hist$mids, empirical.ap.hist$density)
ap.hist.data.95<-subset(ap.hist.data, ap.hist.data[,1] > ap.95[1] & ap.hist.data[,1] < ap.95[2])

dp.95<-quantile(ilogit(empirical[,3]), probs=c(0.025, 0.975))
dp.hist.data<-cbind(empirical.dp.hist$mids, empirical.dp.hist$density)
dp.hist.data.95<-subset(dp.hist.data, dp.hist.data[,1] > dp.95[1] & dp.hist.data[,1] < dp.95[2])

dabp.95<-quantile(ilogit(empirical[,4]), probs=c(0.025, 0.975))
dabp.hist.data<-cbind(empirical.dabp.hist$mids, empirical.dabp.hist$density)
dabp.hist.data.95<-subset(dabp.hist.data, dabp.hist.data[,1] > dabp.95[1] & dabp.hist.data[,1] < dabp.95[2])

png(file='~/Desktop/EmpiricalPrev_Posterior_wRug.png', height=10, width=30, res=300, units='cm')
par(mfrow=c(1,3), mar=c(8,2,4,2))
plot(empirical.ap.hist$mids, empirical.ap.hist$density, type='l', 
     xlim=c(0, max(subset(tick.prev, species=='adult')$prev)+.05),
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2,
     axes=FALSE)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nAdult Ticks', side=1, line=6, cex=1.5)
polygon(c(min(empirical.ap.hist$mids),empirical.ap.hist$mids,max(empirical.ap.hist$mids)), 
        c(0,empirical.ap.hist$density,0), col=alpha('violetred3', 0.5))
polygon(c(min(ap.hist.data.95[,1]),ap.hist.data.95[,1],max(ap.hist.data.95[,1])), 
        c(0,ap.hist.data.95[,2],0), col=alpha('violetred3', 0.5))
rug(subset(tick.prev, species=='adult')$prev, ticksize=0.03, lwd=2.5, col='red')
title('A', adj=0, cex.main=3)

plot(empirical.dp.hist$mids, empirical.dp.hist$density, type='l', axes=FALSE, xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(0,empirical.dp.hist$mids), c(0,empirical.dp.hist$density), col=alpha('blue1', 0.5))
polygon(c(min(dp.hist.data.95[,1]),dp.hist.data.95[,1],max(dp.hist.data.95[,1])), 
        c(0,dp.hist.data.95[,2],0), col=alpha('blue1', 0.5))
rug(subset(tick.prev, species=='deer')$prev, ticksize=0.03, lwd=2.5, col='red')
title('B', adj=0, cex.main=3)

plot(empirical.dabp.hist$mids, empirical.dabp.hist$density, type='l', xlim=c(0,1),axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Seroprevalence \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(min(empirical.dabp.hist$mids), empirical.dabp.hist$mids, max(empirical.dabp.hist$mids)), 
        c(0,empirical.dabp.hist$density,0), col=alpha('yellow3', 0.5))
polygon(c(min(dabp.hist.data.95[,1]),dabp.hist.data.95[,1],max(dabp.hist.data.95[,1])), 
        c(0,dabp.hist.data.95[,2],0), col=alpha('yellow3', 0.5))
rug(subset(tick.prev, species=='deerAB')$prev, ticksize=0.03, lwd=2.5, col='red')
title('C', adj=0, cex.main=3)

dev.off()

## -------------------------------------------------------------------------------------
## 7. Tick Burdens
## -------------------------------------------------------------------------------------

setwd('~/Desktop/')
data.2.2<-read.csv('Multivar_Metropolis_200D_200R/Run1_200D_200R/1 MCMC Accepted Iteration Log 2014-10-11 17-46 839.csv')
data.2.2b<-data.2.2[20000:100001,]
deer.burden<-rowSums(data.2.2b[,4:8])
hist(deer.burden/200)
alt.burden<-rowSums(data.2.2b[,9:13])
hist(alt.burden/200)
all.ticks<-rowSums(data.2.2b[,4:23])
hist(all.ticks)
on.host<-rowSums(cbind(deer.burden, alt.burden))
hist(on.host/400)
on.veg<-rowSums(data.2.2b[,14:23])
hist(on.veg)

data.2.4<-read.csv('Multivar_Metropolis_200D_400R/Run1_200D_400R/1 MCMC Accepted Iteration Log 2014-10-11 17-46 643.csv')
data.2.4b<-data.2.4[20000:100001,]
deer.burden24<-rowSums(data.2.4b[,4:8])
alt.burden24<-rowSums(data.2.4b[,9:13])
all.ticks24<-rowSums(data.2.4b[,4:23])
on.host24<-rowSums(cbind(deer.burden24, alt.burden24))
on.veg24<-rowSums(data.2.4b[,14:23])
hist(deer.burden24/200)
hist(alt.burden24/400)
hist(on.veg24)
hist(on.host24/600)

data.4.2<-read.csv('Multivar_Metropolis_400D_200R/Run1_400D_200R/1 MCMC Accepted Iteration Log 2014-10-11 17-46 389.csv')
data.4.2b<-data.4.2[20000:100001,]
deer.burden42<-rowSums(data.4.2b[,4:8])
alt.burden42<-rowSums(data.4.2b[,9:13])

par(mfrow=c(3,2))
hist(deer.burden24/200, main='1:2', xlim=c(0,35))
hist(alt.burden24/400, xlim=c(0,35))
hist(deer.burden/200, main='1:1', xlim=c(0,35), breaks=25)
hist(alt.burden/200, xlim=c(0,35), breaks=25)
hist(deer.burden42/400, main='2:1', xlim=c(0,35))
hist(alt.burden42/200, xlim=c(0,35))
