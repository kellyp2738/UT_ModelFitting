## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## Revised Plotting Scripts for Model Fitting Dissertation Chapter
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

## Where are the data?
##  ON THE BACKUP HARD DRIVE FOR THE WORK COMPUTER
##  ~/Desktop/Multivar_Metropolis_200D_200R
##  ~/Desktop/Multivar_Metropolis_400D_200R
##  ~/Desktop/Multivar_Metropolis_200D_400R
##  ~/Desktop has some of the plots
##  ~/Dropbox/ModelFitting/FutureProof/Multivariate_Metropolis_Plots_Copy has copies of the plots
##  ON THE NEW LAPTOP
##  ~/Desktop/MV_UnequalPopSizes/MV_200_600_fixed
##  ~/Desktop/MV_UnequalPopSizes/MV_600_200_fixed
##  ~/Desktop/MV_UnequalPopSizes/MV_200_800_fixed
##  ~/Desktop/MV_UnequalPopSizes/MV_800_200_fixed

## note: home directory on new laptop is 'kelly' not 'kellypierce':
source("/Users/kelly/Dropbox/2014&Older/ModelFitting/R Scripts/MCMC_Source.r")
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
#setwd('~/Desktop/MV_UnequalPopSizes/')
setwd('~/Desktop/MV_ConstantHostPop/')
#dirs<-grep('Multivar_Metropolis', list.files(getwd()), value=TRUE)
dirs<-grep('MV_', list.files(getwd()), value=TRUE)
for(dir in dirs[4:5]){ 
  #setwd('~/Desktop/BlockUpdate_ChangeKE')
  #setwd('~/Desktop/FinalModelMCMCRuns/SymTrans')
  #setwd(file.path('~/Desktop/MV_UnequalPopSizes', dir))
  setwd(file.path('~/Desktop/MV_ConstantHostPop/', dir))
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
  deer.burden<-c()
  alt.burden<-c()
  for(chain in chain.dirs){ #a bunch of chain 1 files are messed up for some reason...
    #setwd(file.path('~/Desktop/FinalModelMCMCRuns/SymTrans', dir))
    #setwd(file.path('~/Desktop/MV_UnequalPopSizes', dir))
    setwd(file.path('~/Desktop/MV_ConstantHostPop/', dir))
    print(chain)
    print(getwd())
    setwd(file.path(getwd(), chain))
    dir.contents<-list.files(getwd())
    print(getwd())
    file<-dir.contents[1] # the proper MCMC chain is the first file in each directory
    data<-read.csv(file)
    
    #sort data without removing burnin (for plotting chain comparisons)
    #thin.file<-thin.chain(data)
    est.file<-estimated.pars(data, 7)
    pref.chains<-cbind(pref.chains, est.file[,1])
    rhoTD.chains<-cbind(rhoTD.chains, est.file[,2])
    rhoDT.chains<-cbind(rhoDT.chains, est.file[,3])
    
    #sort data with 20k iteration 'burnin' removed (for convergence diagnostics)
    #thin.file.burn<-thin.chain(data[20000:100000,])
    est.file.burn<-data[100000:200000,]
    pref.chains.burn<-cbind(pref.chains.burn, est.file.burn[,1])
    rhoTD.chains.burn<-cbind(rhoTD.chains.burn, est.file.burn[,2])
    rhoDT.chains.burn<-cbind(rhoDT.chains.burn, est.file.burn[,3])
    loglik.burn<-cbind(loglik.burn, est.file.burn[,28])
    aPrev.burn<-cbind(aPrev.burn, est.file.burn[,25])
    dPrev.burn<-cbind(dPrev.burn, est.file.burn[,26])
    dABprev.burn<-cbind(dABprev.burn, est.file.burn[,27])
    deer.burden<-cbind(deer.burden, rowSums(est.file.burn[,4:8]))
    alt.burden<-cbind(alt.burden, rowSums(est.file.burn[,9:13]))
  }
  #write the sorted data to file
  #setwd(file.path('~/Desktop/MV_UnequalPopSizes', dir))
  setwd(file.path('~/Desktop/MV_ConstantHostPop/', dir))
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
  write.csv(deer.burden, file='deer_burden_chains.csv')
  write.csv(alt.burden, file='alt_burden_chains.csv')
}

## -----------------------------------------------------------------------------------
## 2. Plot chains and perform convergence diagnostics on the re-organized data
## -----------------------------------------------------------------------------------

# read in replicate chains, coerce to MCMC objects, and perform the gelman-rubin convergence test
#pdf(file='~/Desktop/Multivariate_Metropolis_MultpleChains.pdf', height=5, width=7)
for(dir in dirs[4:5]){
  # change to appropriate directory
  #setwd('~/Desktop/MV_UnequalPopSizes/')
  setwd(file.path('~/Desktop/MV_ConstantHostPop/', dir))
  #setwd(file.path(getwd(), dir))
  print(dir)
  
  # read data, select appropriate columns (first column is iteration number, which isn't necessary), exponentiate, and set column names
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
  # how many chains?
  nc<-length(pref)
  png(file=file.path(getwd(),paste(dir,'Pref_Chains.png')), height=10*nc, width=60, units='cm', res=300)
  par(mfrow=c(nc,1), mar=c(5,6,2,2))
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  for(i in seq(1,nc)){
    plot(1:length(pref[,i]), pref[,i], type='l', col=line.colors[i], 
         main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
         ylab='', las=1)
    mtext(expression(phi['D']), side=2, line=4, cex=1.5)
    abline(v=100000, lty=2, lwd=3)
    abline(v=50000, lty=2, lwd=3)
  }
  dev.off()
  ## plot chains on top of each other
  #png(file=file.path('~/Desktop/',paste(dir,'Chains.png')), height=30, width=30, units='cm', res=300)
  #par(mfrow=c(3,1), mar=c(5,6,2,2))
  #plot(1:length(pref[,1]), pref[,1], type='l', col=line.colors[1], 
  #     main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
  #     ylab='', las=1)
  #mtext(expression(phi['D']), side=2, line=4, cex=1.5)
  #for(i in 2:length(pref[1,])){
  #  lines(1:length(pref[,i]), pref[,i], col=line.colors[i])
  #}
  #abline(v=10000, lty=2, lwd=3)
  #A=hist(pref[,1], plot=FALSE, breaks=100)
  #plot(y=A$mids, x=A$counts, type='l', lwd=2, col=alpha(line.colors[i]), 
  #     xlim=c(0,max(A$counts)), ylim=c(0,1), xlab='Counts', yaxt='n', ylab="")
  #for(i in 2:length(pref[1,])){
  #  B=hist(pref[,i], plot=FALSE, breaks=100)
  #  lines(y=B$mids, x=B$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  
  ## rhoTD ##
  nt<-length(rhoTD)
  png(file=file.path(getwd(),paste(dir,'RhoTD_Chains.png')), height=10*nt, width=60, units='cm', res=300)
  par(mfrow=c(nt,1), mar=c(5,6,2,2))
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  for(i in seq(1,nt)){
    plot(1:length(rhoTD[,i]), rhoTD[,i], type='l', col=line.colors[i], 
         main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
         ylab='', las=1)
    mtext(expression(rho['TD']), side=2, line=4, cex=1.5)
    abline(v=100000, lty=2, lwd=3)
    abline(v=50000, lty=2, lwd=3)
  }
  dev.off()
  # plot chains on top of each other (requires uncommenting lines 141-142 that call the png() function)
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  #plot(1:length(rhoTD[,1]), rhoTD[,1], type='l', col=line.colors[1], 
  #     main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
  #     xlab='Iteration',ylab='')
  #mtext(expression(rho['TD']), side=2, line=4, cex=1.5)
  #for(i in 2:length(rhoTD[1,])){
  #  lines(1:length(rhoTD[,i]), rhoTD[,i], col=line.colors[i])
  #}
  #abline(v=10000, lty=2, lwd=3)
  #C=hist(rhoTD[,1], plot=FALSE, breaks=25)
  #plot(y=C$mids, x=C$counts, type='l', lwd=2, col=alpha(line.colors[1]), 
  #     xlim=c(0, max(C$counts)), ylim=c(0,1), yaxt='n', xlab='Counts')
  #for(i in 2:length(rhoTD[1,])){
  #  D=hist(rhoTD[,i], plot=FALSE, breaks=25)
  #  lines(y=D$mids, x=D$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  
  ## rhoDT ##
  nd<-length(rhoDT)
  png(file=file.path(getwd(),paste(dir,'RhoDT_Chains.png')), height=10*nd, width=60, units='cm', res=300)
  par(mfrow=c(nd,1), mar=c(5,6,2,2))
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  for(i in seq(1,nd)){
    plot(1:length(rhoDT[,i]), rhoDT[,i], type='l', col=line.colors[i], 
         main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
         ylab='', las=1)
    mtext(expression(rho['DT']), side=2, line=4, cex=1.5)
    abline(v=100000, lty=2, lwd=3)
    abline(v=50000, lty=2, lwd=3)
  }
  dev.off()
  # plot chains on top of each other (requires uncommenting lines 141-142 that call the png() function)
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  #plot(1:length(rhoDT[,1]), rhoDT[,1], type='l', col=line.colors[1], 
  #     main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
  #     xlab='Iteration',ylab='')
  #mtext(expression(rho['DT']), side=2, line=4, cex=1.5)
  #for(i in 2:length(rhoDT[1,])){
  #  lines(1:length(rhoDT[,i]), rhoDT[,i], col=line.colors[i])
  #}
  #abline(v=10000, lty=2, lwd=3)
  #E=hist(rhoDT[,1], plot=FALSE, breaks=30)
  #plot(y=E$mids, x=E$counts, type='l', lwd=2, col=alpha(line.colors[1]), 
  #     xlim=c(0, max(E$counts)), ylim=c(0,1), yaxt='n', xlab='Counts')
  #for(i in 2:length(rhoDT[1,])){
  #  G=hist(rhoTD[,i], plot=FALSE, breaks=30)
  #  lines(y=G$mids, x=G$counts, lwd=2, col=alpha(line.colors[i], 0.9))   
  #}
  #dev.off()
  
  # convert chains into mcmc objects and create mcmc lists
  rhoDT.mcmc<-list()
  rhoTD.mcmc<-list()
  pref.mcmc<-list()
  for(i in 1:length(pref[1,])){
    pref.mcmc[[i]]<-as.mcmc(pref[100000:200000,i])
    rhoDT.mcmc[[i]]<-as.mcmc(rhoDT[100000:200000,i])
    rhoTD.mcmc[[i]]<-as.mcmc(rhoTD[100000:200000,i])
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

## Decent evidence for convergence in 100/900, 300/700, and 500/500
## Will only consider scenarios where deer << alternative hosts (maybe universially true, esp. w/respect to counts of animals if not biomass)

## overlay the chains (instead of printing each chain separately)
for(dir in dirs[1:3]){
  # change to appropriate directory
  setwd('~/Desktop/MV_ConstantHostPop/')
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
  png(file=file.path('~/Desktop/MV_ConstantHostPop/',paste(dir,'Overlaid_Chains.png')), height=30, width=30, units='cm', res=300)
  par(mfrow=c(3,1), mar=c(5,6,2,2))
  #layout(matrix(c(1,2),1,2), widths=c(5,2))
  plot(1:length(pref[,1]), pref[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), xlab='Iteration', cex.lab=1.5, cex.axis=1.5,
       ylab='', las=1)
  mtext(expression(phi['D']), side=2, line=4, cex=1.5)
  for(i in 2:length(pref[1,])){
    lines(1:length(pref[,i]), pref[,i], col=line.colors[i])
  }
  abline(v=50000, lty=2, lwd=3)
  abline(v=100000, lty=2, lwd=3)

  plot(1:length(rhoTD[,1]), rhoTD[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
       xlab='Iteration',ylab='')
  mtext(expression(rho['TD']), side=2, line=4, cex=1.5)
  for(i in 2:length(rhoTD[1,])){
    lines(1:length(rhoTD[,i]), rhoTD[,i], col=line.colors[i])
  }
  abline(v=50000, lty=2, lwd=3)
  abline(v=100000, lty=2, lwd=3)

  plot(1:length(rhoDT[,1]), rhoDT[,1], type='l', col=line.colors[1], 
       main='', ylim=c(0,1), las=1, cex.lab=1.5, cex.axis=1.5,
       xlab='Iteration',ylab='')
  mtext(expression(rho['DT']), side=2, line=4, cex=1.5)
  for(i in 2:length(rhoDT[1,])){
    lines(1:length(rhoDT[,i]), rhoDT[,i], col=line.colors[i])
  }
  abline(v=50000, lty=2, lwd=3)
  abline(v=100000, lty=2, lwd=3)

  dev.off()
}

## -----------------------------------------------------------------------------------
## 3. Look at marginal posterior probabilities for the estimated parameters
## -----------------------------------------------------------------------------------

## Preference Posterior (see DissertationPlots.r for the code from the dissertation; those lines have been deleted from this code)

## Read in the relevant data
pref1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/preference_burned_multi_chains.csv')
pref3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/preference_burned_multi_chains.csv')
pref5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/preference_burned_multi_chains.csv')
pref7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/preference_burned_multi_chains.csv')
pref9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/preference_burned_multi_chains.csv')

rhoDT1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/rhoDT_burned_multi_chains.csv')
rhoDT3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/rhoDT_burned_multi_chains.csv')
rhoDT5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/rhoDT_burned_multi_chains.csv')
rhoDT7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/rhoDT_burned_multi_chains.csv')
rhoDT9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/rhoDT_burned_multi_chains.csv')

rhoTD1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/rhoTD_burned_multi_chains.csv')
rhoTD3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/rhoTD_burned_multi_chains.csv')
rhoTD5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/rhoTD_burned_multi_chains.csv')
rhoTD7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/rhoTD_burned_multi_chains.csv')
rhoTD9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/rhoTD_burned_multi_chains.csv')

## --- 100/900 scenario --- 

## plots to show convergence/overlap of posteriors for different chains 

png('~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/PosteriorOverlap-May2015.png', height=30, width=30, units='cm', res=300)
par(mfrow=c(3,3), mar=c(5,4,2,2))

#png('~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/PosteriorOverlap-100_900.png', height=10, width=30, units='cm', res=300)
#par(mfrow=c(1,3), mar=c(5,4,2,2))

#PREFERENCE

post.plot<-function(param.data, param.name, y.label, color, add.line=NULL, add.line.color=NULL){
  hist.list<-list()
  for(chain in 2:6){
    hist.list[[chain-1]]<-hist(exp(param.data[,chain]), plot=FALSE)
  }
  ymax=0
  for(hist in 1:5){
    if(max(hist.list[[hist]]$counts)>ymax){
      ymax=max(hist.list[[hist]]$counts)
    }
  }
  plot(hist.list[[1]]$mids, hist.list[[1]]$counts, type='l', xlim=c(0,1), ylim=c(0, ymax),
       xlab=param.name, axes=FALSE, cex.lab=1.5, ylab='')
  for(hist in 1:5){
    polygon(x=c(min(hist.list[[hist]]$mids), hist.list[[hist]]$mids, max(hist.list[[hist]]$mids)),
            y=c(0, hist.list[[hist]]$counts, 0), col=alpha(color, 0.25))
  }
  axis(side=1, cex.axis=1.5)
  mtext(text=y.label, side=2, line=0.5, cex=1.5)
  if(!is.null(add.line)){
    abline(v=add.line, col=add.line.color, lty=2, lwd=3)
  }
  #return(hist.list)
}

par(mfrow=c(3,5))
blues<-brewer.pal(5, 'Blues') #original color from dissertation: 'darkolive4'
post.plot(pref1.9, expression(phi['D']), 'Scenario 1', blues[1], add.line=0.5, add.line.color='red')
post.plot(pref3.7, expression(phi['D']), 'Scenario 2', blues[2], add.line=0.5, add.line.color='red')
post.plot(pref5.5, expression(phi['D']), 'Scenario 3', blues[3], add.line=0.5, add.line.color='red')
post.plot(pref7.3, expression(phi['D']), 'Scenario 4', blues[4], add.line=0.5, add.line.color='red')
post.plot(pref9.1, expression(phi['D']), 'Scenario 5', blues[5], add.line=0.5, add.line.color='red')

#par(mfrow=c(1,5))
reds<-brewer.pal(5, 'Reds') #original color from dissertation: 'orchid4'
post.plot(rhoTD1.9, expression(rho['D'%->%'T']), 'Scenario 1', reds[1], add.line=0.06, add.line.color='black')
post.plot(rhoTD3.7, expression(rho['D'%->%'T']), 'Scenario 2', reds[2], add.line=0.06, add.line.color='black')
post.plot(rhoTD5.5, expression(rho['D'%->%'T']), 'Scenario 3', reds[3], add.line=0.06, add.line.color='black')
post.plot(rhoTD7.3, expression(rho['D'%->%'T']), 'Scenario 4', reds[4], add.line=0.06, add.line.color='black')
post.plot(rhoTD9.1, expression(rho['D'%->%'T']), 'Scenario 5', reds[5], add.line=0.06, add.line.color='black')

greens<-brewer.pal(5, 'Greens') #original color from dissertation: 'darkorange2'
post.plot(rhoDT1.9, expression(rho['T'%->%'D']), 'Scenario 1', greens[1], add.line=0.26, add.line.color='black')
post.plot(rhoDT3.7, expression(rho['T'%->%'D']), 'Scenario 2', greens[2], add.line=0.26, add.line.color='black')
post.plot(rhoDT5.5, expression(rho['T'%->%'D']), 'Scenario 3', greens[3], add.line=0.26, add.line.color='black')
post.plot(rhoDT7.3, expression(rho['T'%->%'D']), 'Scenario 4', greens[4], add.line=0.26, add.line.color='black')
post.plot(rhoDT9.1, expression(rho['T'%->%'D']), 'Scenario 5', greens[5], add.line=0.26, add.line.color='black')

## example code for how you would find the 95% CI and shade in that section of the distribution

pref2.2.95<-quantile(exp(pref2.2[,2]), probs=c(0.025, 0.975))
pref2.2.95.hist.data<-cbind(pref2.2.hist$mids, pref2.2.hist$density)
pref2.2.95.hist<-subset(pref2.2.95.hist.data, pref2.2.95.hist.data[,1] > pref2.2.95[1] & pref2.2.95.hist.data[,1] < pref2.2.95[2])
plot(pref2.4.hist$mids, pref2.4.hist$density, type='l', axes=FALSE, main='', xlim=c(0,1),
     xlab='', ylab='', las=1, cex.lab=2)
mtext(text='1:2', cex=1.5, side=2, line=1)
axis(side=1, cex.axis=2)
polygon(c(min(pref2.4.hist$mids), pref2.4.hist$mids, max(pref2.4.hist$mids)), 
        c(0,pref2.4.hist$density,0), col=alpha('darkolivegreen4', 0.5))
polygon(c(min(pref2.4.95.hist[,1]), pref2.4.95.hist[,1], max(pref2.4.95.hist[,1])), 
        c(0,pref2.4.95.hist[,2],0), col=alpha('darkolivegreen4', 0.5))
abline(v=0.5, lty=2, lwd=3, col='red')

# each param on a single graph, all scenarios
get.hist.coords<-function(data.vector){
  mids<-hist(data.vector, breaks=60, plot=FALSE)$mids
  dens<-hist(data.vector, breaks=60, plot=FALSE)$density
  coords<-cbind(mids, dens)
  return(coords)
}

pref.list<-data.frame(cbind(exp(pref1.9[,2]), exp(pref3.7[,2]), exp(pref5.5[,2]), exp(pref7.3[,2])))
rhoTD.list<-data.frame(cbind(exp(rhoTD1.9[,2]), exp(rhoTD3.7[,2]), exp(rhoTD5.5[,2]), exp(rhoTD7.3[,2])))
rhoDT.list<-data.frame(cbind(exp(rhoDT1.9[,2]), exp(rhoDT3.7[,2]), exp(rhoDT5.5[,2]), exp(rhoDT7.3[,2])))

scenario.colors<-c('#a1dab4', '#41b6c4', '#225ea8')

png(file="~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/parameter_transitions_June2015.png", 
    height=30, width=20, units='cm', res=300)
par(mar=c(5,3,3,2), mfrow=c(3,1))
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,5), axes=FALSE, ylab="", xlab="")
mtext(side=1, text=expression(phi['D']), line=3, cex=1.5)
for(i in 1:3){
  i.coords<-get.hist.coords(pref.list[,i])
  #lines(i.coords[,1], i.coords[,2])
  extra.points<-seq(0,1,0.001)
  predict.posterior<-predict(loess(i.coords[,2]~i.coords[,1], span=0.25), 
                    newdata=extra.points)
  pos.posterior<-(predict.posterior + abs(predict.posterior))/2
  p1<-cbind(extra.points, pos.posterior)
  p2<-p1[complete.cases(p1),]
  #print(head(p2))
  p3<-rbind(c(abs(min(p2[,1])),0), abs(p2), c(max(p2[,1]), 0))
  polygon(p3, col=alpha(scenario.colors[i], 0.5))
}
title('A', adj=0, cex.main=3)
legend(x='topright', bty="n", legend=c('Scenario 1 (100:900)', 
                             'Scenario 2 (300:700)', 
                             'Scenario 3 (500:500)'), 
       fill=alpha(scenario.colors, 0.25), cex=1.5)
axis(side=1, cex.axis=1.5)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,4), axes=FALSE, ylab="", xlab="")
mtext(side=1, text=expression(rho['D'%->%'T']), line=3, cex=1.5)
for(i in 1:3){
  i.coords<-get.hist.coords(rhoTD.list[,i])
  #lines(i.coords[,1], i.coords[,2])
  extra.points<-seq(0,1,0.001)
  predict.posterior<-predict(loess(i.coords[,2]~i.coords[,1], span=0.25), 
                             newdata=extra.points)
  pos.posterior<-(predict.posterior + abs(predict.posterior))/2
  p1<-cbind(extra.points, pos.posterior)
  p2<-p1[complete.cases(p1),]
  #print(head(p2))
  p3<-rbind(c(abs(min(p2[,1])),0), abs(p2), c(max(p2[,1]), 0))
  polygon(p3, col=alpha(scenario.colors[i], 0.5))
}
title('B', adj=0, cex.main=3)
abline(v=0.06, lty=2, lwd=3)
axis(side=1, cex.axis=1.5)
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,170), axes=FALSE, ylab="", xlab="")
mtext(side=1, text=expression(rho['T'%->%'D']), line=3, cex=1.5)
for(i in 1:3){
  i.coords<-get.hist.coords(rhoDT.list[,i])
  #extra.points<-seq(0,1,0.001)
  predict.posterior<-predict(loess(i.coords[,2]~i.coords[,1], span=0.25), 
                             newdata=extra.points)
  pos.posterior<-(predict.posterior + abs(predict.posterior))/2
  p1<-cbind(extra.points, pos.posterior)
  p2<-p1[complete.cases(p1),]
  #print(head(p2))
  p3<-rbind(c(abs(min(p2[,1])),0), abs(p2), c(max(p2[,1]), 0))
  polygon(p3, col=alpha(scenario.colors[i], 0.5))
}
title('C', adj=0, cex.main=3)
abline(v=0.26, lty=2, lwd=3)
axis(side=1, cex.axis=1.5)
dev.off()

## -----------------------------------------------------------------------------------
## 4. Look at parameter correlations in the different relative abundance scenarios by
##    making a series of pairs plots (1 master plot with 9 panels) EXPONENTIATED
## -----------------------------------------------------------------------------------

# plot characteristics:
# - chains thinned to 10% (every 10th row is plotted)
# - all chains are plotted in a given scatterplot (not a representative chain as in the original dissertation)
# - color interpolation is based on maximum of all of the chains, not just a single chain
# - log likelihoods are identical to posteriors because uniform priors are used in the model fitting (likelihood = posterior probability)

ll1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/loglik_burned_multi_chains.csv')
ll3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/loglik_burned_multi_chains.csv')
ll5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/loglik_burned_multi_chains.csv')
ll7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/loglik_burned_multi_chains.csv')
ll9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/loglik_burned_multi_chains.csv')

png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Multivariate_Metropolis_ExpPairs_grayscale_Legend.png', height=10, width=5, units='cm', res=300)
plot.new()
post.colors<-colorRampPalette(c('gray20', 'gray90'))(length(1:max(abs(round(ll1.9[,2:length(ll1.9[1,])])))))
legend_image <- as.raster(matrix(post.colors, ncol=1))
rasterImage(legend_image, 0, 0, 1,1)
mtext(text = c('poor'), side=2, line=0, at=c(0), las=1, cex=1.5)
mtext(text = c('good'), side=2, line=0, at=c(1), las=1, cex=1.5)
mtext(text = c('Model Fit'), side=2, line=2, cex=1.5)
dev.off()

#png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Multivariate_Metropolis_ExpPairs_Colored_May2015.png', height=30, width=30, units='cm', res=300)
png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Multivariate_Metropolis_ExpPairs_Colored_June2015.png', height=30, width=30, units='cm', res=300)
#par(mfrow=c(3,3), mar=c(5,7,2,2))

plot.param.corr<-function(param.data.1, param.data.2, 
                          param.1.name, param.2.name,
                          lik.data, thin.by, y.label, 
                          incl.axis1=TRUE, incl.axis2=TRUE, alpha.lvl=1){
  #par(mar=c(6,8,2,2))
  thin.by.idx<-seq(1, 100001, thin.by)
  # color ramp was originally gray20 to yellow
  post.colors<-colorRampPalette(c('gray20', 'gray90'))(length(1:max(abs(round(lik.data[,2:length(lik.data)])))))
  plot(param.data.1[thin.by.idx,2], param.data.2[thin.by.idx,2], xlim=c(log(10^-2.5), log(1)), 
       ylim=c(log(10^-2.5), log(1)), axes=FALSE, pch=16,
       xlab='', ylab=param.2.name, las=1, cex.axis=1.5, 
       cex.lab=1.5, col=alpha(post.colors[abs(round(lik.data[thin.by.idx,2]))], alpha.lvl))
  mtext(side=1, text=param.1.name, line=2)
  for(chain in 3:6){
    points(param.data.1[thin.by.idx,chain], param.data.2[thin.by.idx,chain], pch=16,
           col=alpha(post.colors[abs(round(lik.data[thin.by.idx,chain]))], alpha.lvl))
  }
  if(incl.axis1==TRUE){
    magaxis(side=c(1), las=1, cex.axis=1.5, unlog=TRUE)
  }
  if(incl.axis2==TRUE){
    magaxis(side=c(2), las=1, cex.axis=1.5, unlog=TRUE)
  }
  mtext(side=2, text=y.label, cex=1.3, line=4.2)
}

par(mfcol=c(3,3))
par(mar=c(3,6,0.5,0.5))
plot.param.corr(pref1.9, rhoDT1.9, '', expression(rho['T'%->%'D']), 
                ll1.9, 10, 'Scenario 1', incl.axis1=FALSE, incl.axis2=TRUE)
plot.param.corr(pref3.7, rhoDT3.7, '', expression(rho['T'%->%'D']), 
                ll3.7, 10, 'Scenario 2', incl.axis1=FALSE)
plot.param.corr(pref5.5, rhoDT5.5, expression(phi['D']), expression(rho['T'%->%'D']), 
                ll5.5, 10, 'Scenario 3')
#plot.param.corr(pref7.3, rhoDT7.3, expression(phi['D']), expression(rho['T'%->%'D']), ll7.3, 10, 'Scenario 4')
#plot.param.corr(pref9.1, rhoDT9.1, expression(phi['D']), expression(rho['T'%->%'D']), ll9.1, 10, 'Scenario 5')

plot.param.corr(pref1.9, rhoTD1.9, '', expression(rho['D'%->%'T']), 
                ll1.9, 10, '', incl.axis1=FALSE)
plot.param.corr(pref3.7, rhoTD3.7, '', expression(rho['D'%->%'T']), 
                ll3.7, 10, '', incl.axis1=FALSE)
plot.param.corr(pref5.5, rhoTD5.5, expression(phi['D']), expression(rho['D'%->%'T']), 
                ll5.5, 10, '')
#plot.param.corr(pref7.3, rhoTD7.3, expression(phi['D']), expression(rho['D'%->%'T']), ll7.3, 10, 'Scenario 4')
#plot.param.corr(pref9.1, rhoTD9.1, expression(phi['D']), expression(rho['D'%->%'T']), ll9.1, 10, 'Scenario 5')

plot.param.corr(rhoTD1.9, rhoDT1.9, '', expression(rho['T'%->%'D']), 
                ll1.9, 10, '', incl.axis1=FALSE)
plot.param.corr(rhoTD3.7, rhoDT3.7, '', expression(rho['T'%->%'D']), 
                ll3.7, 10, '', incl.axis1=FALSE)
plot.param.corr(rhoTD5.5, rhoDT5.5, expression(rho['D'%->%'T']), expression(rho['T'%->%'D']), 
                ll5.5, 10, '')
#plot.param.corr(rhoTD7.3, rhoDT7.3, expression(rho['D'%->%'T']), expression(rho['T'%->%'D']), ll7.3, 10, 'Scenario 4')
#plot.param.corr(rhoTD9.1, rhoDT9.1, expression(rho['D'%->%'T']), expression(rho['T'%->%'D']), ll9.1, 10, 'Scenario 5')

dev.off()

thin.by.idx<-seq(1, 100001, 10)
cor(pref1.9[thin.by.idx,2], rhoDT1.9[thin.by.idx,2])
cor(pref1.9[thin.by.idx,2], rhoTD1.9[thin.by.idx,2])
cor(rhoTD1.9[thin.by.idx,2], rhoDT1.9[thin.by.idx,2])

cor(pref3.7[thin.by.idx,2], rhoDT3.7[thin.by.idx,2])
cor(pref3.7[thin.by.idx,2], rhoTD3.7[thin.by.idx,2])
cor(rhoTD3.7[thin.by.idx,2], rhoDT3.7[thin.by.idx,2])

cor(pref5.5[thin.by.idx,2], rhoDT5.5[thin.by.idx,2])
cor(pref5.5[thin.by.idx,2], rhoTD5.5[thin.by.idx,2])
cor(rhoTD5.5[thin.by.idx,2], rhoDT5.5[thin.by.idx,2])


## -----------------------------------------------------------------------------------
## 5. Meta-analysis empirical prevalence posteriors compared to model output prev
## -----------------------------------------------------------------------------------

# read in the data generated in the ticks.R script
empirical<-read.csv('~/Dropbox/2014&Older/ModelFitting/FutureProof/TickPrevGLMM_MetaAnalysis_NoRacc_Results_LogitTransform.csv')
ilogit = function(x) 1/{1+exp(-x)}

# modification of post.plot for prevalence data
prev.plot<-function(param.data, param.name, y.label, color, empirical, 
                    xmax=1, alpha.lvl=0.1, add.line=NULL, add.line.color=NULL){
  hist.list<-list()
  for(chain in 1:length(param.data)){
    hist.list[[chain]]<-hist((param.data[,chain]), plot=FALSE)
  }
  hist.list[[length(param.data)+1]]<-hist(empirical, plot=FALSE)
  ymax=0
  for(hist in 1:length(hist.list)){
    if(max(hist.list[[hist]]$density)>ymax){
      ymax=max(hist.list[[hist]]$density)
    }
  }
  plot(hist.list[[1]]$mids, hist.list[[1]]$density, type='l', xlim=c(0,xmax), ylim=c(0, ymax),
       xlab=param.name, axes=FALSE, cex.lab=2, cex.axis=2, ylab='')
  for(hist in 1:length(hist.list)){
    polygon(x=c(min(hist.list[[hist]]$mids), hist.list[[hist]]$mids, max(hist.list[[hist]]$mids)),
            y=c(0, hist.list[[hist]]$density, 0), col=alpha(color, alpha.lvl))
  }
  axis(side=1, cex.axis=2, cex.lab=2)
  mtext(text=y.label, side=2, line=0.5, cex=2)
  if(!is.null(add.line)){
    abline(v=add.line, col=add.line.color, lty=2, lwd=3)
  }
  #return(hist.list)
}
 
# ticks # original color in dissertation: violetred3
ap1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/aPrev_burned_multi_chains.csv')
ap3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/aPrev_burned_multi_chains.csv')
ap5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/aPrev_burned_multi_chains.csv')
#ap7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/aPrev_burned_multi_chains.csv')
#ap9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/aPrev_burned_multi_chains.csv')
empirical.ap<-ilogit(empirical[,2])

par(mfrow=c(1,1))
prev.plot(ap1.9[,-1], 'Prevalence, Adult Ticks', 'Scenario 1', 'red', empirical.ap, 0.1)
prev.plot(ap3.7[,-1], 'Prevalence, Adult Ticks', 'Scenario 2', 'red', empirical.ap, 0.1)
prev.plot(ap5.5[,-1], 'Prevalence, Adult Ticks', 'Scenario 3', 'red', empirical.ap, 0.1)
#prev.plot(ap7.3, 'Prevalence, Adult Ticks', 'Scenario 4', 'red', empirical.ap, 0.1)
#prev.plot(ap9.1, 'Prevalence, Adult Ticks', 'Scenario 5', 'red', empirical.ap, 0.1)

ap.all<-data.frame(cbind(ap1.9[,2], ap3.7[,2], ap5.5[,2]))
prev.plot(ap.all, 'Prevalence, Adult Ticks', '', '#ef8a62', empirical.ap, 0.1, 0.5)

# deer # original color in dissertation: blue1
dp1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/dPrev_burned_multi_chains.csv')
dp3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/dPrev_burned_multi_chains.csv')
dp5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/dPrev_burned_multi_chains.csv')
#dp7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/dPrev_burned_multi_chains.csv')
#dp9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/dPrev_burned_multi_chains.csv')
empirical.dp<-ilogit(empirical[,3])

par(mfrow=c(1,1))
prev.plot(dp1.9[,-1], 'Prevalence, Deer', 'Scenario 1', 'blue', empirical.dp)
prev.plot(dp3.7[,-1], 'Prevalence, Deer', 'Scenario 2', 'blue', empirical.dp)
prev.plot(dp5.5[,-1], 'Prevalence, Deer', 'Scenario 3', 'blue', empirical.dp)
#prev.plot(dp7.3, 'Prevalence, Deer', 'Scenario 4', 'blue', empirical.dp)
#prev.plot(dp9.1, 'Prevalence, Deer', 'Scenario 5', 'blue', empirical.dp)

deer.all<-data.frame(cbind(dp1.9[,2], dp3.7[,2], dp5.5[,2]))
prev.plot(deer.all, 'Prevalence, Deer', '', '#67a9cf', empirical.dp, 1, 0.5)

# deer AB # original color in dissertation: yellow3
dabp1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/dABprev_burned_multi_chains.csv')
dabp3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/dABPrev_burned_multi_chains.csv')
dabp5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/dABPrev_burned_multi_chains.csv')
#dabp7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/dABprev_burned_multi_chains.csv')
#dabp9.1<-read.csv('~/Desktop/MV_ConstantHostPop/MV_900_100/dABprev_burned_multi_chains.csv')
empirical.dabp<-ilogit(empirical[,4])

par(mfrow=c(1,1))
prev.plot(dabp1.9, 'Antibody Prevalence, Deer', 'Scenario 1', 'green', empirical.dabp)
prev.plot(dabp3.7, 'Antibody Prevalence, Deer', 'Scenario 2', 'green', empirical.dabp)
prev.plot(dabp5.5, 'Antibody Prevalence, Deer', 'Scenario 3', 'green', empirical.dabp)
#prev.plot(dabp7.3, 'Antibody Prevalence, Deer', 'Scenario 4', 'green', empirical.dabp)
#prev.plot(dabp9.1, 'Antibody Prevalence, Deer', 'Scenario 5', 'green', empirical.dabp)

deer.ab.all<-data.frame(cbind(dabp1.9[,2], dabp3.7[,2], dabp5.5[,2]))
prev.plot(deer.ab.all, 'Antibody Prevalence, Deer', '', '#ffffbf', empirical.dabp, 1, 0.5)

#png(file='~/Desktop/Multivar_Metropolis_Prevalence_Comparison.png', height=10, width=30, units='cm', res=300)
png(file='~/Dropbox/2014&Older/ModelFitting/FutureProof/Figures for Presentations/Multivar_Metropolis_Prevalence_Comparison_RevJune2015.png', height=10, width=30, units='cm', res=300)
par(mfrow=c(1,3), mar=c(8,2,4,2))

# alt color #ef8a62
prev.plot(ap.all, '', '', 'violetred3', empirical.ap, 0.1, 0.5)
title('A', adj=0, cex.main=3)
mtext(side=1, text='Fraction Adult\nTicks Infected', line=5.5, cex=1.5)
# alt color #67a9cf
prev.plot(deer.all, '', '', 'blue1', empirical.dp, 1, 0.5)
title('B', adj=0, cex.main=3)
mtext(side=1, text='Fraction White-Tailed\nDeer Infected', line=5.5, cex=1.5)
# alt color #ffffbf
prev.plot(deer.ab.all, '', '', 'yellow3', empirical.dabp, 1, 0.5)
title('C', adj=0, cex.main=3)
mtext(side=1, text='Fraction White-Tailed\nDeer Seropositive', line=5.5, cex=1.5)

dev.off()

## -----------------------------------------------------------------------------------
## 6. Meta-analysis distributions w/original empirical data
## -----------------------------------------------------------------------------------

ticks = read.csv("~/Dropbox/2014&Older/ModelFitting/FutureProof/CompletePrevsFixed_2.csv", header=FALSE)
names(ticks) = c("species", "ncases", "ntrials")
ticks = ticks[order(ticks$species),]

tick.prev<-cbind(ticks, (ticks$ncases/ticks$ntrials))
names(tick.prev)<-c('species', 'ncases', 'ntrials', 'prev')

ilogit = function(x) 1/{1+exp(-x)}
line.colors<-c('violetred3', 'salmon4', 'blue1', 'paleturquoise3', 'yellow3')
#line.colors<-c('#ef8a62', '#67a9cf', '#ffffbf')

empirical<-read.csv('~/Dropbox/2014&Older/ModelFitting/FutureProof/TickPrevGLMM_MetaAnalysis_NoRacc_Results_LogitTransform.csv')
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

#png(file='~/Desktop/EmpiricalPrev_Posterior_wRug.png', height=10, width=30, res=300, units='cm')

## -------------------------------------------------------------------------------------
## 7. Tick Burdens
## -------------------------------------------------------------------------------------

par(mfrow=c(1,1))
deer.burden1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/deer_burden_chains.csv')
deer.burden3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/deer_burden_chains.csv')
deer.burden5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/deer_burden_chains.csv')

alt.burden1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/alt_burden_chains.csv')
alt.burden3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/alt_burden_chains.csv')
alt.burden5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/alt_burden_chains.csv')

# smoothing density plot option
smooth.hists<-function(data1, data2, pop1, pop2, colors){
  data1.hist<-hist(data1/pop1, breaks=40, plot=FALSE)
  data2.hist<-hist(data2/pop2, breaks=40, plot=FALSE)
  max.burden<-max(c(data1.hist$mids, data2.hist$mids))
  max.dens<-max(c(data1.hist$density, data2.hist$density))
  extra.points<-seq(0,max.burden,0.5)
  plot(x=NULL, y=NULL, xlim=c(0,60), ylim=c(0,max.dens), axes=FALSE, 
       xlab='Average Number of Ticks per Host', cex.lab=1.5, ylab="")
  predict1<-predict(loess(data1.hist$density~data1.hist$mids, span=0.18), 
                    newdata=extra.points)
  predict2<-predict(loess(data2.hist$density~data2.hist$mids, span=0.18), 
                    newdata=extra.points)
  d1<-cbind(extra.points, predict1)
  d2<-d1[complete.cases(d1),]
  d3<-rbind(c(min(d2[,1]),0), d2, c(max(d2[,1]), 0))
  a1<-cbind(extra.points, predict2)
  a2<-a1[complete.cases(a1),]
  a3<-rbind(c(min(a2[,1]),0), a2, c(max(a2[,1]), 0))
  polygon(d3, col=alpha(colors[1], 0.5))
  polygon(a3, col=alpha(colors[2], 0.5))
  axis(side=1, cex.axis=1.5)
}

png('~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Burdens_June2015.png', height=30, width=20, units='cm', res=300)
par(mfrow=c(3,1))
smooth.hists(deer.burden1.9[,2], alt.burden1.9[,2], 100, 900, colors=c('blue', 'red'))
title('A', adj=0, cex.main=3)
legend(x='topright', legend=c('Deer', 'Alternative Hosts'), fill=alpha(c('blue','red'),0.5),
       cex=1.5, bty='n')
smooth.hists(deer.burden3.7[,2], alt.burden3.7[,2], 300, 700, colors=c('blue', 'red'))
title('B', adj=0, cex.main=3)
smooth.hists(deer.burden5.5[,2], alt.burden5.5[,2], 500, 500, colors=c('blue', 'red'))
title('C', adj=0, cex.main=3)
dev.off()

# old figure
#png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Burdens_May2015.png', height=30, width=20, units='cm', res=300)

# reading in the 7.3 scenario, but it's not plotted above
deer.burden7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/deer_burden_chains.csv')
alt.burden7.3<-read.csv('~/Desktop/MV_ConstantHostPop/MV_700_300/alt_burden_chains.csv')

# ----------------------------------------------------------------------
# -- 8. ALL THE 95% CIs
# ----------------------------------------------------------------------

getCIs<-function(data.vector){
  m<-median(data.vector)
  q<-quantile(data.vector, probs=c(0.025, 0.975))
  out<-c(m,q)
  print(out)
  return(out)
}

pref.list<-data.frame(cbind(exp(pref1.9[,2]), exp(pref3.7[,2]), exp(pref5.5[,2]), exp(pref7.3[,2])))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(pref.list[,i])
}

rhoTD.list<-data.frame(cbind(exp(rhoTD1.9[,2]), exp(rhoTD3.7[,2]), exp(rhoTD5.5[,2]), exp(rhoTD7.3[,2])))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(rhoTD.list[,i])
}

rhoDT.list<-data.frame(cbind(exp(rhoDT1.9[,2]), exp(rhoDT3.7[,2]), exp(rhoDT5.5[,2]), exp(rhoDT7.3[,2])))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(rhoDT.list[,i])
}

deer.burden.list<-data.frame(cbind(deer.burden1.9[,2], deer.burden3.7[,2], deer.burden5.5[,2], deer.burden7.3[,2]))
deer.pop<-c(100,300,500,700)
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(deer.burden.list[,i]/deer.pop[i])
}

alt.burden.list<-data.frame(cbind(alt.burden1.9[,2], alt.burden3.7[,2], alt.burden5.5[,2], alt.burden7.3[,2]))
alt.pop<-c(900,700,500,300)
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(alt.burden.list[,i]/alt.pop[i])
}

ap.list<-data.frame(cbind(ap1.9[,2], ap3.7[,2], ap5.5[,2], ap7.3[,2]))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(ap.list[,i])
}

dp.list<-data.frame(cbind(dp1.9[,2], dp3.7[,2], dp5.5[,2], dp7.3[,2]))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(dp.list[,i])
}

dabp.list<-data.frame(cbind(dabp1.9[,2], dabp3.7[,2], dabp5.5[,2], dabp7.3[,2]))
for(i in 1:4){
  #head(pref.list[,i])
  getCIs(dabp.list[,i])
}

# ---------------------------------------------------------------------
# -- 9. Probability tick finds and feeds on deer
# ---------------------------------------------------------------------
phi<-seq(0,1,0.01)
D<-c(100,300,500)
R<-c(900,700,500)
eff.deer.1<-(phi*D[1])/(D[1]+R[1])*100
eff.deer.2<-(phi*D[2])/(D[2]+R[2])*100
eff.deer.3<-(phi*D[3])/(D[3]+R[3])*100
scenario.colors<-c('#a1dab4', '#41b6c4', '#225ea8')

png('~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/phi_behavior.png', width=20, height=20,
    units='cm', res=300)
par(mfrow=c(1,1), mar=c(5,5,3,2))
#bquote is the function for mixing text and expressions, but it's hard to use...
plot(phi, eff.deer.1, xlim=c(0,1), ylim=c(0,100), type='l', col=scenario.colors[1],
     axes=FALSE, xlab=expression(phi['D']), ylab='',
     cex.lab=1.5, lwd=3)
mtext(side=2, line=3.2, text='Pr(Tick Encounters and Feeds on Deer)', cex=1.5)
lines(phi, eff.deer.2, col=scenario.colors[2], lwd=3, lty=2)
lines(phi, eff.deer.3, col=scenario.colors[3], lwd=3, lty=3)
axis(side=1, cex.axis=1.5)
axis(side=2, cex.axis=1.5, las=1)
legend(x='topleft', col=scenario.colors[1:3], lty=c(1,2,3), lwd=3, bty='n',
       legend=c('Scenario 1 (100:900)', 'Scenario 2 (300:700)', 'Scenario 3 (500:500)'))
dev.off()

pr.encounter.feed<-function(phi, D, H){
  pr.e.f<-phi*(D/H)*100
  return(pr.e.f)
}
p1<-pr.encounter.feed(phi=1/2, D=seq(0,500,100), H=1000)
p2<-pr.encounter.feed(phi=1/3, D=seq(0,500,100), H=1000)
p3<-pr.encounter.feed(phi=1/4, D=seq(0,500,100), H=1000)
p4<-pr.encounter.feed(phi=1/5, D=seq(0,500,100), H=1000)
rel.abundance.range<-(seq(0,500,100)/1000)*100
no.pref<-(seq(0,500,100)/1000)*100
png('~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/dilution_expectations.png', width=20, height=20,
    units='cm', res=300)
par(mfrow=c(1,1), mar=c(6,6,3,2))
plot(rel.abundance.range, p1, ylim=c(0,50), las=1, ylab='', xlab='Relative Abundance, Deer',
     type='l', lty=3, lwd=2, col=scenario.colors[1], cex.lab=1.5, axes=F, 
     xlim=c(0,50))
axis(side=1, cex.axis=1.5)
axis(side=2, cex.axis=1.5)
mtext(side=2, line=3.2, text='Pr(Tick Encounters and Feeds on Deer)', cex=1.5)
lines(rel.abundance.range, p2, lty=3, lwd=2, col=scenario.colors[2])
lines(rel.abundance.range, p3, lty=3, lwd=2, col=scenario.colors[3])
lines(rel.abundance.range, p4, lty=3, lwd=2)
lines(rel.abundance.range, no.pref, lwd=2)
points(x=10, y=pr.encounter.feed(0.422, 100, 1000), pch=16, col='red')
points(x=30, y=pr.encounter.feed(0.216, 300, 1000), pch=16, col='red')
points(x=50, y=pr.encounter.feed(0.14, 500, 1000), pch=16, col='red')
legend('topleft', legend=c('No Preference Term', 'Equal Preference, 2 Hosts', 'Equal Preference, 3 Hosts', 'Equal Preference, 4 Hosts', 'Equal Preference, 5 hosts', 'Model Estimates'),
       lty=c(1, 3, 3, 3, 3, NA), col=c('black', scenario.colors, 'black', 'red'), pch=c(NA, NA, NA, NA, NA, 16), bty='n',
       lwd=c(2, 2, 2, 2, 2, NA))
dev.off()

library(asbio)
library(vegan)
deer<-c(100,300,500)
alt<-c(900,700,500)
scenarios<-cbind(deer,alt)
ds<-diversity(scenarios)
plot(deer/1000, ds, pch=16, ylim=c(0,0.7))
points(deer/1000, c(0.422, 0.216, 0.14), col='red', pch=16)
plot(c(0.422, 0.216, 0.14), ds, col='blue', pch=16, xlim=c(0,1), ylim=c(0,1))
cor.test(ds, c) # estimated pref parameter correlates with diversity
# that's not unexpected, but is it meaningful?

h3<-diversity(c(500,500,500))
h4<-diversity(c(500,500,500,500))
h5<-diversity(c(500,500,500,500,500))
h6<-diversity(c(500,500,500,500,500,500))

rel<-lm(ds~c(0.422, 0.216, 0.14))

pred.phi<-function(x){
  phi<-((coefficients(rel)[[1]])+(x*coefficients(rel)[[2]]))
  return(phi)
}

lines(c(0.422, 0.216, 0.14), pred.phi(c(0.422, 0.216, 0.14)))
pred.phi(h3)

# how does probability of tick encounter scale with different metrics of diversity?

it is a set probability of encounter + feeding that explains the disease prevalence
what would actually drive this in natural settings?
it can't be that preference really changes as a function of community composition
rather preference is likely something that is set, and different communities will have

because the model was fitted to studies drawn from a variety of different sites and with
presumably different community structures, we don't know which community assemblage makes sense



# ----------------------------------------------------------------------
# -- 10. FURTHER EXPLORATION OF PARAMETER CORRELATIONS
# ----------------------------------------------------------------------

# pairs plot for thinned chain representative of 100:900 deer:alt scenario
metrics1.9.2<-cbind(exp(pref1.9[,2]), exp(rhoDT1.9[,2]), 
                    exp(rhoTD1.9[,2]), ll1.9[,2],
                    ap1.9[,2], dp1.9[,2], dabp1.9[,2])
thin.by<-seq(1, length(metrics1.9.2[,1]), 100)
metrics1.9.2.thin<-metrics1.9.2[thin.by,]
names(metrics1.9.2.thin)<-c('pref', 'trans D->T', 'trans T->D',
                            'loglik', 'adult prev', 'deer prev', 
                            'deer AB prev')
png(file='~/Desktop/MV_ConstantHostPop/Pairs_1.9_chain1_thinned.png', height=30, width=30, units='cm', res=300)
pairs(metrics1.9.2.thin, labels=names(metrics1.9.2.thin))  
dev.off()

# pairs plot for thinned chain representative of 500:500 deer:alt scenario are consistent with the 100:900 scenario
metrics5.5.2<-cbind(exp(pref5.5[,2]), exp(rhoDT5.5[,2]), 
                    exp(rhoTD5.5[,2]), ll5.5[,2],
                    ap5.5[,2], dp5.5[,2], dabp5.5[,2])
thin.by55<-seq(1, length(metrics5.5.2[,1]), 100)
metrics5.5.2.thin<-metrics5.5.2[thin.by55,]
names(metrics5.5.2.thin)<-c('pref', 'trans D->T', 'trans T->D',
                            'loglik', 'adult prev', 'deer prev', 
                            'deer AB prev')
#png(file='~/Desktop/MV_ConstantHostPop/Pairs_1.9_chain1_thinned.png', height=30, width=30, units='cm', res=300)
pairs(metrics5.5.2.thin, labels=names(metrics5.5.2.thin))  
#dev.off()

# the deer prev/deer AB prev numbers are perfectly correlated because AB prevalence is 
# a deterministic function of infection prevalence -- the duration of infection and duration
# of immunity parameters are fixed.