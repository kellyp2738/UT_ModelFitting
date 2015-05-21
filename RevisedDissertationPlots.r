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
  for(hist in 2:5){
    polygon(x=c(min(hist.list[[hist]]$mids), hist.list[[hist]]$mids, max(hist.list[[hist]]$mids)),
            y=c(0, hist.list[[hist]]$counts, 0), col=alpha(color, 0.25))
  }
  axis(side=1, cex.axis=1.5)
  mtext(text=y.label, side=2, line=0.5, cex=1.5)
  if(!is.null(add.line)){
    abline(v=add.line, col=add.line.color, lty=2, lwd=3)
  }
  return(hist.list)
}

par(mfrow=c(3,5))
blues<-brewer.pal(5, 'Blues')
post.plot(pref1.9, expression(phi['D']), 'Scenario 1', blues[1])
post.plot(pref3.7, expression(phi['D']), 'Scenario 2', blues[2])
post.plot(pref5.5, expression(phi['D']), 'Scenario 3', blues[3])
post.plot(pref7.3, expression(phi['D']), 'Scenario 4', blues[4])
post.plot(pref9.1, expression(phi['D']), 'Scenario 5', blues[5])

#par(mfrow=c(1,5))
reds<-brewer.pal(5, 'Reds')
post.plot(rhoTD1.9, expression(rho['D'%->%'T']), 'Scenario 1', reds[1])
post.plot(rhoTD3.7, expression(rho['D'%->%'T']), 'Scenario 2', reds[2])
post.plot(rhoTD5.5, expression(rho['D'%->%'T']), 'Scenario 3', reds[3])
post.plot(rhoTD7.3, expression(rho['D'%->%'T']), 'Scenario 4', reds[4])
post.plot(rhoTD9.1, expression(rho['D'%->%'T']), 'Scenario 5', reds[5])

greens<-brewer.pal(5, 'Greens')
post.plot(rhoDT1.9, expression(rho['T'%->%'D']), 'Scenario 1', greens[1])
post.plot(rhoDT3.7, expression(rho['T'%->%'D']), 'Scenario 2', greens[2])
post.plot(rhoDT5.5, expression(rho['T'%->%'D']), 'Scenario 3', greens[3])
post.plot(rhoDT7.3, expression(rho['T'%->%'D']), 'Scenario 4', greens[4])
post.plot(rhoDT9.1, expression(rho['T'%->%'D']), 'Scenario 5', greens[5])

pref1.9.hist1<-hist(exp(pref1.9[,2]), breaks=30, plot=FALSE)
plot(pref1.9.hist1$mids, pref1.9.hist1$counts, type='l', xlim=c(0,1), ylim=c(0, 19500), xlab=expression(phi['D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(pref1.9.hist1$mids),pref1.9.hist1$mids, max(pref1.9.hist1$mids)), y=c(0,pref1.9.hist1$counts,0), col=alpha('darkolivegreen4', 0.25))
axis(side=1, cex.axis=1.5)
mtext(text='Scenario 1', side=2, line=0.5)
abline(v=0.5, col='red', lty=2, lwd=2)

pref1.9.hist2<-hist(exp(pref1.9[,3]), breaks=30, plot=FALSE)
lines(pref1.9.hist2$mids, pref1.9.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref1.9.hist2$mids),pref1.9.hist2$mids, max(pref1.9.hist2$mids)), y=c(0,pref1.9.hist2$counts,0), col=alpha('darkolivegreen4', 0.25))

pref1.9.hist3<-hist(exp(pref1.9[,4]), breaks=30, plot=FALSE)
lines(pref1.9.hist3$mids, pref1.9.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref1.9.hist3$mids),pref1.9.hist3$mids, max(pref1.9.hist3$mids)), y=c(0,pref1.9.hist3$counts,0), col=alpha('darkolivegreen4', 0.25))

pref1.9.hist4<-hist(exp(pref1.9[,5]), breaks=30, plot=FALSE)
lines(pref1.9.hist4$mids, pref1.9.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref1.9.hist4$mids),pref1.9.hist4$mids, max(pref1.9.hist4$mids)), y=c(0,pref1.9.hist4$counts,0), col=alpha('darkolivegreen4', 0.25))

pref1.9.hist5<-hist(exp(pref1.9[,6]), breaks=30, plot=FALSE)
lines(pref1.9.hist5$mids, pref1.9.hist5$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref1.9.hist5$mids),pref1.9.hist5$mids, max(pref1.9.hist5$mids)), y=c(0,pref1.9.hist5$counts,0), col=alpha('darkolivegreen4', 0.25))


#RHODT (tick to deer)
rhoDT1.9.hist1<-hist(exp(rhoDT1.9[,2]), breaks=30, plot=FALSE)
plot(rhoDT1.9.hist1$mids, rhoDT1.9.hist1$counts, type='l', xlim=c(0,0.5), xlab=expression(rho['T'%->%'D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoDT1.9.hist1$mids),rhoDT1.9.hist1$mids, max(rhoDT1.9.hist1$mids)), y=c(0,rhoDT1.9.hist1$counts,0), col=alpha('darkorange2', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.26, col='black', lty=2, lwd=2)

rhoDT1.9.hist2<-hist(exp(rhoDT1.9[,3]), breaks=30, plot=FALSE)
lines(rhoDT1.9.hist2$mids, rhoDT1.9.hist2$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT1.9.hist2$mids),rhoDT1.9.hist2$mids, max(rhoDT1.9.hist2$mids)), y=c(0,rhoDT1.9.hist2$counts,0), col=alpha('darkorange2', 0.25))

rhoDT1.9.hist3<-hist(exp(rhoDT1.9[,4]), breaks=30, plot=FALSE)
lines(rhoDT1.9.hist3$mids, rhoDT1.9.hist3$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT1.9.hist3$mids),rhoDT1.9.hist3$mids, max(rhoDT1.9.hist3$mids)), y=c(0,rhoDT1.9.hist3$counts,0), col=alpha('darkorange2', 0.25))

rhoDT1.9.hist4<-hist(exp(rhoDT1.9[,5]), breaks=30, plot=FALSE)
lines(rhoDT1.9.hist4$mids, rhoDT1.9.hist4$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT1.9.hist4$mids),rhoDT1.9.hist4$mids, max(rhoDT1.9.hist4$mids)), y=c(0,rhoDT1.9.hist4$counts,0), col=alpha('darkorange2', 0.25))

rhoDT1.9.hist5<-hist(exp(rhoDT1.9[,6]), breaks=30, plot=FALSE)
lines(rhoDT1.9.hist5$mids, rhoDT1.9.hist5$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT1.9.hist5$mids),rhoDT1.9.hist5$mids, max(rhoDT1.9.hist5$mids)), y=c(0,rhoDT1.9.hist5$counts,0), col=alpha('darkorange2', 0.25))

#RHOTD (deer to tick)
rhoTD1.9.hist1<-hist(exp(rhoTD1.9[,2]), breaks=30, plot=FALSE)
plot(rhoTD1.9.hist1$mids, rhoTD1.9.hist1$counts, type='l', xlim=c(0,1), ylim=c(0,12000), xlab=expression(rho['D'%->%'T']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoTD1.9.hist1$mids),rhoTD1.9.hist1$mids, max(rhoTD1.9.hist1$mids)), y=c(0,rhoTD1.9.hist1$counts,0), col=alpha('orchid4', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.06, col='black', lty=2, lwd=2)

rhoTD1.9.hist2<-hist(exp(rhoTD1.9[,3]), breaks=30, plot=FALSE)
lines(rhoTD1.9.hist2$mids, rhoTD1.9.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD1.9.hist2$mids),rhoTD1.9.hist2$mids, max(rhoTD1.9.hist2$mids)), y=c(0,rhoTD1.9.hist2$counts,0), col=alpha('orchid4', 0.25))

rhoTD1.9.hist3<-hist(exp(rhoTD1.9[,4]), breaks=30, plot=FALSE)
lines(rhoTD1.9.hist3$mids, rhoTD1.9.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD1.9.hist3$mids),rhoTD1.9.hist3$mids, max(rhoTD1.9.hist3$mids)), y=c(0,rhoTD1.9.hist3$counts,0), col=alpha('orchid4', 0.25))

rhoTD1.9.hist4<-hist(exp(rhoTD1.9[,5]), breaks=30, plot=FALSE)
lines(rhoTD1.9.hist4$mids, rhoTD1.9.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD1.9.hist4$mids),rhoTD1.9.hist4$mids, max(rhoTD1.9.hist4$mids)), y=c(0,rhoTD1.9.hist4$counts,0), col=alpha('orchid4', 0.25))

rhoTD1.9.hist5<-hist(exp(rhoTD1.9[,6]), breaks=30, plot=FALSE)
lines(rhoTD1.9.hist5$mids, rhoTD1.9.hist5$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD1.9.hist5$mids),rhoTD1.9.hist5$mids, max(rhoTD1.9.hist5$mids)), y=c(0,rhoTD1.9.hist5$counts,0), col=alpha('orchid4', 0.25))
#dev.off()

## --- 300/700 scenario --- 

## plots to show convergence/overlap of posteriors for different chains (for my own edification/possibly for supplement)

#PREF
pref3.7.hist1<-hist(exp(pref3.7[,2]), breaks=30, plot=FALSE)
plot(pref3.7.hist1$mids, pref3.7.hist1$counts, type='l', xlim=c(0,1), ylim=c(0,55000), xlab=expression(phi['D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(pref3.7.hist1$mids),pref3.7.hist1$mids, max(pref3.7.hist1$mids)), y=c(0,pref3.7.hist1$counts,0), col=alpha('darkolivegreen4', 0.25))
axis(side=1, cex.axis=1.5)
mtext(text='Scenario 2', side=2, line=0.5)
abline(v=0.5, col='red', lty=2, lwd=2)

pref3.7.hist2<-hist(exp(pref3.7[,3]), breaks=30, plot=FALSE)
lines(pref3.7.hist2$mids, pref3.7.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref3.7.hist2$mids),pref3.7.hist2$mids, max(pref3.7.hist2$mids)), y=c(0,pref3.7.hist2$counts,0), col=alpha('darkolivegreen4', 0.25))

pref3.7.hist3<-hist(exp(pref3.7[,4]), breaks=30, plot=FALSE)
lines(pref3.7.hist3$mids, pref3.7.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref3.7.hist3$mids),pref3.7.hist3$mids, max(pref3.7.hist3$mids)), y=c(0,pref3.7.hist3$counts,0), col=alpha('darkolivegreen4', 0.25))

pref3.7.hist4<-hist(exp(pref3.7[,5]), breaks=30, plot=FALSE)
lines(pref3.7.hist4$mids, pref3.7.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref3.7.hist4$mids),pref3.7.hist4$mids, max(pref3.7.hist4$mids)), y=c(0,pref3.7.hist4$counts,0), col=alpha('darkolivegreen4', 0.25))

pref3.7.hist5<-hist(exp(pref3.7[,6]), breaks=30, plot=FALSE)
lines(pref3.7.hist5$mids, pref3.7.hist5$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref3.7.hist5$mids),pref3.7.hist5$mids, max(pref3.7.hist5$mids)), y=c(0,pref3.7.hist5$counts,0), col=alpha('darkolivegreen4', 0.25))

#RHODT (tick to deer)
rhoDT3.7.hist1<-hist(exp(rhoDT3.7[,2]), breaks=30, plot=FALSE)
plot(rhoDT3.7.hist1$mids, rhoDT3.7.hist1$counts, type='l', xlim=c(0,0.5), ylim=c(0,40000), xlab=expression(rho['T'%->%'D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoDT3.7.hist1$mids),rhoDT3.7.hist1$mids, max(rhoDT3.7.hist1$mids)), y=c(0,rhoDT3.7.hist1$counts,0), col=alpha('darkorange2', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.26, col='black', lty=2, lwd=2)

rhoDT3.7.hist2<-hist(exp(rhoDT3.7[,3]), breaks=30, plot=FALSE)
lines(rhoDT3.7.hist2$mids, rhoDT3.7.hist2$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT3.7.hist2$mids),rhoDT3.7.hist2$mids, max(rhoDT3.7.hist2$mids)), y=c(0,rhoDT3.7.hist2$counts,0), col=alpha('darkorange2', 0.25))

rhoDT3.7.hist3<-hist(exp(rhoDT3.7[,4]), breaks=30, plot=FALSE)
lines(rhoDT3.7.hist3$mids, rhoDT3.7.hist3$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT3.7.hist3$mids),rhoDT3.7.hist3$mids, max(rhoDT3.7.hist3$mids)), y=c(0,rhoDT3.7.hist3$counts,0), col=alpha('darkorange2', 0.25))

rhoDT3.7.hist4<-hist(exp(rhoDT3.7[,5]), breaks=30, plot=FALSE)
lines(rhoDT3.7.hist4$mids, rhoDT3.7.hist4$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT3.7.hist4$mids),rhoDT3.7.hist4$mids, max(rhoDT3.7.hist4$mids)), y=c(0,rhoDT3.7.hist4$counts,0), col=alpha('darkorange2', 0.25))

rhoDT3.7.hist5<-hist(exp(rhoDT3.7[,6]), breaks=30, plot=FALSE)
lines(rhoDT3.7.hist5$mids, rhoDT3.7.hist5$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT3.7.hist5$mids),rhoDT3.7.hist5$mids, max(rhoDT3.7.hist5$mids)), y=c(0,rhoDT3.7.hist5$counts,0), col=alpha('darkorange2', 0.25))

#RHOTD (deer to tick)
rhoTD3.7.hist1<-hist(exp(rhoTD3.7[,2]), breaks=30, plot=FALSE)
plot(rhoTD3.7.hist1$mids, rhoTD3.7.hist1$counts, type='l', xlim=c(0,1), c(0,25000), xlab=expression(rho['D'%->%'T']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoTD3.7.hist1$mids),rhoTD3.7.hist1$mids, max(rhoTD3.7.hist1$mids)), y=c(0,rhoTD3.7.hist1$counts,0), col=alpha('orchid4', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.06, col='black', lty=2, lwd=2)

rhoTD3.7.hist2<-hist(exp(rhoTD3.7[,3]), breaks=30, plot=FALSE)
lines(rhoTD3.7.hist2$mids, rhoTD3.7.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD3.7.hist2$mids),rhoTD3.7.hist2$mids, max(rhoTD3.7.hist2$mids)), y=c(0,rhoTD3.7.hist2$counts,0), col=alpha('orchid4', 0.25))

rhoTD3.7.hist3<-hist(exp(rhoTD3.7[,4]), breaks=30, plot=FALSE)
lines(rhoTD3.7.hist3$mids, rhoTD3.7.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD3.7.hist3$mids),rhoTD3.7.hist3$mids, max(rhoTD3.7.hist3$mids)), y=c(0,rhoTD3.7.hist3$counts,0), col=alpha('orchid4', 0.25))

rhoTD3.7.hist4<-hist(exp(rhoTD3.7[,5]), breaks=30, plot=FALSE)
lines(rhoTD3.7.hist4$mids, rhoTD3.7.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD3.7.hist4$mids),rhoTD3.7.hist4$mids, max(rhoTD3.7.hist4$mids)), y=c(0,rhoTD3.7.hist4$counts,0), col=alpha('orchid4', 0.25))

rhoTD3.7.hist5<-hist(exp(rhoTD3.7[,6]), breaks=30, plot=FALSE)
lines(rhoTD3.7.hist5$mids, rhoTD3.7.hist5$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD3.7.hist5$mids),rhoTD3.7.hist5$mids, max(rhoTD3.7.hist5$mids)), y=c(0,rhoTD3.7.hist5$counts,0), col=alpha('orchid4', 0.25))

## --- 500/500 scenario --- 

## plots to show convergence/overlap of posteriors for different chains (for my own edification/possibly for supplement)
pref5.5.hist1<-hist(exp(pref5.5[,2]), breaks=30, plot=FALSE)
plot(pref5.5.hist1$mids, pref5.5.hist1$counts, type='l', xlim=c(0,1), xlab=expression(phi['D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(pref5.5.hist1$mids),pref5.5.hist1$mids, max(pref5.5.hist1$mids)), y=c(0,pref5.5.hist1$counts,0), col=alpha('darkolivegreen4', 0.25))
axis(side=1, cex.axis=1.5)
mtext(text='Scenario 3', side=2, line=0.5)
abline(v=0.5, col='red', lty=2, lwd=2)

pref5.5.hist2<-hist(exp(pref5.5[,3]), breaks=30, plot=FALSE)
lines(pref5.5.hist2$mids, pref5.5.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref5.5.hist2$mids),pref5.5.hist2$mids, max(pref5.5.hist2$mids)), y=c(0,pref5.5.hist2$counts,0), col=alpha('darkolivegreen4', 0.25))

pref5.5.hist3<-hist(exp(pref5.5[,4]), breaks=30, plot=FALSE)
lines(pref5.5.hist3$mids, pref5.5.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref5.5.hist3$mids),pref5.5.hist3$mids, max(pref5.5.hist3$mids)), y=c(0,pref5.5.hist3$counts,0), col=alpha('darkolivegreen4', 0.25))

pref5.5.hist4<-hist(exp(pref5.5[,5]), breaks=30, plot=FALSE)
lines(pref5.5.hist4$mids, pref5.5.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(pref5.5.hist4$mids),pref5.5.hist4$mids, max(pref5.5.hist4$mids)), y=c(0,pref5.5.hist4$counts,0), col=alpha('darkolivegreen4', 0.25))

#pref5.5.hist5<-hist(exp(pref5.5[,6]), breaks=30, plot=FALSE)
#lines(pref5.5.hist5$mids, pref5.5.hist5$counts, type='l', xlim=c(0,1))
#polygon(x=c(min(pref5.5.hist5$mids),pref5.5.hist5$mids, max(pref5.5.hist5$mids)), y=c(0,pref5.5.hist5$counts,0), col=alpha('darkolivegreen4', 0.25))

#RHODT (tick to deer)
rhoDT5.5.hist1<-hist(exp(rhoDT5.5[,2]), breaks=30, plot=FALSE)
plot(rhoDT5.5.hist1$mids, rhoDT5.5.hist1$counts, type='l', xlim=c(0,0.5), ylim=c(0,33000), xlab=expression(rho['T'%->%'D']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoDT5.5.hist1$mids),rhoDT5.5.hist1$mids, max(rhoDT5.5.hist1$mids)), y=c(0,rhoDT5.5.hist1$counts,0), col=alpha('darkorange2', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.26, col='black', lty=2, lwd=2)

rhoDT5.5.hist2<-hist(exp(rhoDT5.5[,3]), breaks=30, plot=FALSE)
lines(rhoDT5.5.hist2$mids, rhoDT5.5.hist2$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT5.5.hist2$mids),rhoDT5.5.hist2$mids, max(rhoDT5.5.hist2$mids)), y=c(0,rhoDT5.5.hist2$counts,0), col=alpha('darkorange2', 0.25))

rhoDT5.5.hist3<-hist(exp(rhoDT5.5[,4]), breaks=30, plot=FALSE)
lines(rhoDT5.5.hist3$mids, rhoDT5.5.hist3$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT5.5.hist3$mids),rhoDT5.5.hist3$mids, max(rhoDT5.5.hist3$mids)), y=c(0,rhoDT5.5.hist3$counts,0), col=alpha('darkorange2', 0.25))

rhoDT5.5.hist4<-hist(exp(rhoDT5.5[,5]), breaks=30, plot=FALSE)
lines(rhoDT5.5.hist4$mids, rhoDT5.5.hist4$counts, type='l', xlim=c(0,0.5))
polygon(x=c(min(rhoDT5.5.hist4$mids),rhoDT5.5.hist4$mids, max(rhoDT5.5.hist4$mids)), y=c(0,rhoDT5.5.hist4$counts,0), col=alpha('darkorange2', 0.25))

#rhoDT5.5.hist5<-hist(exp(rhoDT5.5[,6]), breaks=30, plot=FALSE)
#lines(rhoDT5.5.hist5$mids, rhoDT5.5.hist5$counts, type='l', xlim=c(0,0.5))
#polygon(x=c(min(rhoDT5.5.hist5$mids),rhoDT5.5.hist5$mids, max(rhoDT5.5.hist5$mids)), y=c(0,rhoDT5.5.hist5$counts,0), col=alpha('darkorange2', 0.25))

#RHOTD (deer to tick)
rhoTD5.5.hist1<-hist(exp(rhoTD5.5[,2]), breaks=30, plot=FALSE)
plot(rhoTD5.5.hist1$mids, rhoTD5.5.hist1$counts, type='l', xlim=c(0,1), c(0,20000), xlab=expression(rho['D'%->%'T']), ylab='', axes=FALSE, cex.lab=1.5)
polygon(x=c(min(rhoTD5.5.hist1$mids),rhoTD5.5.hist1$mids, max(rhoTD5.5.hist1$mids)), y=c(0,rhoTD5.5.hist1$counts,0), col=alpha('orchid4', 0.25))
axis(side=1, cex.axis=1.5)
abline(v=0.06, col='black', lty=2, lwd=2)

rhoTD5.5.hist2<-hist(exp(rhoTD5.5[,3]), breaks=30, plot=FALSE)
lines(rhoTD5.5.hist2$mids, rhoTD5.5.hist2$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD5.5.hist2$mids),rhoTD5.5.hist2$mids, max(rhoTD5.5.hist2$mids)), y=c(0,rhoTD5.5.hist2$counts,0), col=alpha('orchid4', 0.25))

rhoTD5.5.hist3<-hist(exp(rhoTD5.5[,4]), breaks=30, plot=FALSE)
lines(rhoTD5.5.hist3$mids, rhoTD5.5.hist3$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD5.5.hist3$mids),rhoTD5.5.hist3$mids, max(rhoTD5.5.hist3$mids)), y=c(0,rhoTD5.5.hist3$counts,0), col=alpha('orchid4', 0.25))

rhoTD5.5.hist4<-hist(exp(rhoTD5.5[,5]), breaks=30, plot=FALSE)
lines(rhoTD5.5.hist4$mids, rhoTD5.5.hist4$counts, type='l', xlim=c(0,1))
polygon(x=c(min(rhoTD5.5.hist4$mids),rhoTD5.5.hist4$mids, max(rhoTD5.5.hist4$mids)), y=c(0,rhoTD5.5.hist4$counts,0), col=alpha('orchid4', 0.25))

#rhoTD5.5.hist5<-hist(exp(rhoTD5.5[,6]), breaks=30, plot=FALSE)
#lines(rhoTD5.5.hist5$mids, rhoTD5.5.hist5$counts, type='l', xlim=c(0,1))
#polygon(x=c(min(rhoTD5.5.hist5$mids),rhoTD5.5.hist5$mids, max(rhoTD5.5.hist5$mids)), y=c(0,rhoTD5.5.hist5$counts,0), col=alpha('orchid4', 0.25))

dev.off()

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
post1.9colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(ll1.9[,2:length(ll1.9)])))))

ll3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/loglik_burned_multi_chains.csv')
post3.7colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(ll3.7[,2:length(ll3.7)])))))

ll5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/loglik_burned_multi_chains.csv')
post5.5colors<-colorRampPalette(c('gray20', 'yellow'))(length(1:max(abs(round(ll5.5[,2:length(ll5.5)])))))

png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Multivariate_Metropolis_ExpPairs_Colored_Legend.png', height=10, width=5, units='cm', res=300)
plot.new()
legend_image <- as.raster(matrix(post5.5colors, ncol=1))
rasterImage(legend_image, 0, 0, 1,1)
mtext(text = c('poor'), side=2, line=0, at=c(0), las=1, cex=1.5)
mtext(text = c('good'), side=2, line=0, at=c(1), las=1, cex=1.5)
mtext(text = c('Model Fit'), side=2, line=2, cex=1.5)
dev.off()

png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Multivariate_Metropolis_ExpPairs_Colored_May2015.png', height=30, width=30, units='cm', res=300)
par(mfrow=c(3,3), mar=c(5,7,2,2))

# row 1 (100:900 rel. abundance)
thin.by<-seq(1, 100001, 10)
plot(pref1.9[thin.by,2], rhoDT1.9[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post1.9colors[abs(round(ll1.9[thin.by,2]))])
points(pref1.9[thin.by,3], rhoDT1.9[thin.by,3], col=post1.9colors[abs(round(ll1.9[thin.by,3]))])
points(pref1.9[thin.by,4], rhoDT1.9[thin.by,4], col=post1.9colors[abs(round(ll1.9[thin.by,4]))])
points(pref1.9[thin.by,5], rhoDT1.9[thin.by,5], col=post1.9colors[abs(round(ll1.9[thin.by,5]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)
mtext(side=2, text='Scenario 1', cex=1.3, line=5)

plot(pref1.9[thin.by,2], rhoTD1.9[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['D'%->%'T']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post1.9colors[abs(round(ll1.9[thin.by,2]))])
points(pref1.9[thin.by,3], rhoTD1.9[thin.by,3], col=post1.9colors[abs(round(ll1.9[thin.by,3]))])
points(pref1.9[thin.by,4], rhoTD1.9[thin.by,4], col=post1.9colors[abs(round(ll1.9[thin.by,4]))])
points(pref1.9[thin.by,5], rhoTD1.9[thin.by,5], col=post1.9colors[abs(round(ll1.9[thin.by,5]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)


plot(rhoDT1.9[thin.by,2], rhoTD1.9[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(rho['D'%->%'T']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post1.9colors[abs(round(ll1.9[thin.by,2]))])
points(rhoDT1.9[thin.by,3], rhoTD1.9[thin.by,3], col=post1.9colors[abs(round(ll1.9[thin.by,3]))])
points(rhoDT1.9[thin.by,4], rhoTD1.9[thin.by,4], col=post1.9colors[abs(round(ll1.9[thin.by,4]))])
points(rhoDT1.9[thin.by,5], rhoTD1.9[thin.by,5], col=post1.9colors[abs(round(ll1.9[thin.by,5]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)


# row 2 (300:700 relative abundance)
plot(pref3.7[thin.by,2], rhoDT3.7[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post3.7colors[abs(round(ll3.7[thin.by,2]))])
points(pref3.7[thin.by,3], rhoDT3.7[thin.by,3], col=post3.7colors[abs(round(ll3.7[thin.by,3]))])
points(pref3.7[thin.by,4], rhoDT3.7[thin.by,4], col=post3.7colors[abs(round(ll3.7[thin.by,4]))])
points(pref3.7[thin.by,5], rhoDT3.7[thin.by,5], col=post3.7colors[abs(round(ll3.7[thin.by,5]))])
points(pref3.7[thin.by,6], rhoDT3.7[thin.by,6], col=post3.7colors[abs(round(ll3.7[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)
mtext(side=2, text='Scenario 2', cex=1.3, line=5)

plot(pref3.7[thin.by,2], rhoTD3.7[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['D'%->%'T']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post3.7colors[abs(round(ll3.7[thin.by,2]))])
points(pref3.7[thin.by,3], rhoTD3.7[thin.by,3], col=post3.7colors[abs(round(ll3.7[thin.by,3]))])
points(pref3.7[thin.by,4], rhoTD3.7[thin.by,4], col=post3.7colors[abs(round(ll3.7[thin.by,4]))])
points(pref3.7[thin.by,5], rhoTD3.7[thin.by,5], col=post3.7colors[abs(round(ll3.7[thin.by,5]))])
points(pref3.7[thin.by,6], rhoTD3.7[thin.by,6], col=post3.7colors[abs(round(ll3.7[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)

plot(rhoDT3.7[thin.by,2], rhoTD3.7[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(rho['D'%->%'T']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post3.7colors[abs(round(ll3.7[thin.by,2]))])
points(rhoDT3.7[thin.by,3], rhoTD3.7[thin.by,3], col=post3.7colors[abs(round(ll3.7[thin.by,3]))])
points(rhoDT3.7[thin.by,4], rhoTD3.7[thin.by,4], col=post3.7colors[abs(round(ll3.7[thin.by,4]))])
points(rhoDT3.7[thin.by,5], rhoTD3.7[thin.by,5], col=post3.7colors[abs(round(ll3.7[thin.by,5]))])
points(rhoDT3.7[thin.by,6], rhoTD3.7[thin.by,6], col=post3.7colors[abs(round(ll3.7[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)

# row 3 (500:500 relative abundance)
plot(pref5.5[thin.by,2], rhoDT5.5[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post5.5colors[abs(round(ll5.5[thin.by,2]))])
points(pref5.5[thin.by,3], rhoDT5.5[thin.by,3], col=post5.5colors[abs(round(ll5.5[thin.by,3]))])
points(pref5.5[thin.by,4], rhoDT5.5[thin.by,4], col=post5.5colors[abs(round(ll5.5[thin.by,4]))])
points(pref5.5[thin.by,5], rhoDT5.5[thin.by,5], col=post5.5colors[abs(round(ll5.5[thin.by,5]))])
points(pref5.5[thin.by,6], rhoDT5.5[thin.by,6], col=post5.5colors[abs(round(ll5.5[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)
mtext(side=2, text='Scenario 3', cex=1.3, line=5)

plot(pref5.5[thin.by,2], rhoTD5.5[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(phi['D']), ylab=expression(rho['D'%->%'T']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post5.5colors[abs(round(ll5.5[thin.by,2]))])
points(pref5.5[thin.by,3], rhoTD5.5[thin.by,3], col=post5.5colors[abs(round(ll5.5[thin.by,3]))])
points(pref5.5[thin.by,4], rhoTD5.5[thin.by,4], col=post5.5colors[abs(round(ll5.5[thin.by,4]))])
points(pref5.5[thin.by,5], rhoTD5.5[thin.by,5], col=post5.5colors[abs(round(ll5.5[thin.by,5]))])
points(pref5.5[thin.by,6], rhoTD5.5[thin.by,6], col=post5.5colors[abs(round(ll5.5[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)

plot(rhoDT5.5[thin.by,2], rhoTD5.5[thin.by,2], xlim=c(log(10^-2.5), log(1)), ylim=c(log(10^-2.5), log(1)), axes=FALSE,
     xlab=expression(rho['D'%->%'T']), ylab=expression(rho['T'%->%'D']), las=1, cex.axis=1.5, 
     cex.lab=1.5, col=post5.5colors[abs(round(ll5.5[thin.by,2]))])
points(rhoDT5.5[thin.by,3], rhoTD5.5[thin.by,3], col=post5.5colors[abs(round(ll5.5[thin.by,3]))])
points(rhoDT5.5[thin.by,4], rhoTD5.5[thin.by,4], col=post5.5colors[abs(round(ll5.5[thin.by,4]))])
points(rhoDT5.5[thin.by,5], rhoTD5.5[thin.by,5], col=post5.5colors[abs(round(ll5.5[thin.by,5]))])
points(rhoDT5.5[thin.by,6], rhoTD5.5[thin.by,6], col=post5.5colors[abs(round(ll5.5[thin.by,6]))])
magaxis(side=c(1,2), las=1, cex.axis=1.5, unlog=TRUE)

dev.off()


## -----------------------------------------------------------------------------------
## 5. Meta-analysis empirical prevalence posteriors compared to model output prev
## -----------------------------------------------------------------------------------

# read in the data generated in the ticks.R script
empirical<-read.csv('~/Dropbox/2014&Older/ModelFitting/FutureProof/TickPrevGLMM_MetaAnalysis_NoRacc_Results_LogitTransform.csv')
ilogit = function(x) 1/{1+exp(-x)}

# ticks
ap1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/aPrev_burned_multi_chains.csv')
ap1.9.hist1<-hist(ap1.9[,2], breaks=25, plot=FALSE)
ap1.9.hist2<-hist(ap1.9[,3], breaks=25, plot=FALSE)
ap1.9.hist3<-hist(ap1.9[,4], breaks=25, plot=FALSE)
ap1.9.hist4<-hist(ap1.9[,5], breaks=25, plot=FALSE)

ap3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/aPrev_burned_multi_chains.csv')
ap3.7.hist1<-hist(ap3.7[,2], breaks=25, plot=FALSE)
ap3.7.hist2<-hist(ap3.7[,3], breaks=25, plot=FALSE)
ap3.7.hist3<-hist(ap3.7[,4], breaks=25, plot=FALSE)
ap3.7.hist4<-hist(ap3.7[,5], breaks=25, plot=FALSE)
ap3.7.hist5<-hist(ap3.7[,6], breaks=25, plot=FALSE)

ap5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/aPrev_burned_multi_chains.csv')
ap5.5.hist1<-hist(ap5.5[,2], breaks=25, plot=FALSE)
ap5.5.hist2<-hist(ap5.5[,3], breaks=25, plot=FALSE)
ap5.5.hist3<-hist(ap5.5[,4], breaks=25, plot=FALSE)
ap5.5.hist4<-hist(ap5.5[,5], breaks=25, plot=FALSE)
ap5.5.hist5<-hist(ap5.5[,6], breaks=25, plot=FALSE)

empirical.ap.hist<-hist(ilogit(empirical[,2]), breaks=25, plot=FALSE)

# deer
dp1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/dPrev_burned_multi_chains.csv')
dp1.9.hist1<-hist(dp1.9[,2], breaks=25, plot=FALSE)
dp1.9.hist2<-hist(dp1.9[,3], breaks=25, plot=FALSE)
dp1.9.hist3<-hist(dp1.9[,4], breaks=25, plot=FALSE)
dp1.9.hist4<-hist(dp1.9[,5], breaks=25, plot=FALSE)
dp1.9.hist5<-hist(dp1.9[,6], breaks=25, plot=FALSE)

dp3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/dPrev_burned_multi_chains.csv')
dp3.7.hist1<-hist(dp3.7[,2], breaks=25, plot=FALSE)
dp3.7.hist2<-hist(dp3.7[,3], breaks=25, plot=FALSE)
dp3.7.hist3<-hist(dp3.7[,4], breaks=25, plot=FALSE)
dp3.7.hist4<-hist(dp3.7[,5], breaks=25, plot=FALSE)
dp3.7.hist5<-hist(dp3.7[,6], breaks=25, plot=FALSE)

dp5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/dPrev_burned_multi_chains.csv')
dp5.5.hist1<-hist(dp5.5[,2], breaks=25, plot=FALSE)
dp5.5.hist2<-hist(dp5.5[,3], breaks=25, plot=FALSE)
dp5.5.hist3<-hist(dp5.5[,4], breaks=25, plot=FALSE)
dp5.5.hist4<-hist(dp5.5[,5], breaks=25, plot=FALSE)
dp5.5.hist5<-hist(dp5.5[,6], breaks=25, plot=FALSE)

empirical.dp.hist<-hist(ilogit(empirical[,3]), breaks=25, plot=FALSE)
#ymax.d<-max(c(dp1.9.hist$density, dp3.7hist$density, dp5.5.hist$density, empirical.dp.hist$density))
ymax.d<-max(c(dp1.9.hist1$density, dp1.9.hist2$density, dp1.9.hist3$density, dp1.9.hist4$density, dp1.9.hist5$density, 
                dp3.7.hist1$density, dp3.7.hist2$density, dp3.7.hist3$density, dp3.7.hist4$density, dp3.7.hist5$density, 
                dp5.5.hist1$density, dp5.5.hist2$density, dp5.5.hist3$density, dp5.5.hist4$density, dp5.5.hist5$density,
                empirical.dp.hist$density))

# deer AB
dabp1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/dABprev_burned_multi_chains.csv')
dabp1.9.hist1<-hist(dabp1.9[,2], breaks=25, plot=FALSE)
dabp1.9.hist2<-hist(dabp1.9[,3], breaks=25, plot=FALSE)
dabp1.9.hist3<-hist(dabp1.9[,4], breaks=25, plot=FALSE)
dabp1.9.hist4<-hist(dabp1.9[,5], breaks=25, plot=FALSE)
dabp1.9.hist5<-hist(dabp1.9[,6], breaks=25, plot=FALSE)

dabp3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/dABPrev_burned_multi_chains.csv')
dabp3.7.hist1<-hist(dabp3.7[,2], breaks=25, plot=FALSE)
dabp3.7.hist2<-hist(dabp3.7[,3], breaks=25, plot=FALSE)
dabp3.7.hist3<-hist(dabp3.7[,4], breaks=25, plot=FALSE)
dabp3.7.hist4<-hist(dabp3.7[,5], breaks=25, plot=FALSE)
dabp3.7.hist5<-hist(dabp3.7[,6], breaks=25, plot=FALSE)

dabp5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/dABPrev_burned_multi_chains.csv')
dabp5.5.hist1<-hist(dabp5.5[,2], breaks=25, plot=FALSE)
dabp5.5.hist2<-hist(dabp5.5[,3], breaks=25, plot=FALSE)
dabp5.5.hist3<-hist(dabp5.5[,4], breaks=25, plot=FALSE)
dabp5.5.hist4<-hist(dabp5.5[,5], breaks=25, plot=FALSE)
dabp5.5.hist5<-hist(dabp5.5[,6], breaks=25, plot=FALSE)

empirical.dabp.hist<-hist(ilogit(empirical[,4]), breaks=25, plot=FALSE)
ymax.dab<-max(c(dabp1.9.hist1$density, dabp1.9.hist2$density, dabp1.9.hist3$density, dabp1.9.hist4$density, dabp1.9.hist5$density, 
                dabp3.7.hist1$density, dabp3.7.hist2$density, dabp3.7.hist3$density, dabp3.7.hist4$density, dabp3.7.hist5$density, 
                dabp5.5.hist1$density, dabp5.5.hist2$density, dabp5.5.hist3$density, dabp5.5.hist4$density, dabp5.5.hist5$density,
                empirical.dabp.hist$density))

png(file='~/Desktop/Multivar_Metropolis_Prevalence_Comparison.png', height=10, width=30, units='cm', res=300)
par(mfrow=c(1,3), mar=c(8,2,4,2))
plot(empirical.ap.hist$mids, empirical.ap.hist$density, type='l', ylim=c(0,250),
     axes=FALSE, xlab='', ylab='', las=1, xlim=c(0,0.04))
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nAdult Ticks', side=1, line=6, cex=1.5)
polygon(empirical.ap.hist$mids, empirical.ap.hist$density, col=alpha('violetred3',0.25))
polygon(ap1.9.hist1$mids, ap1.9.hist1$density, col=alpha('violetred3',0.1))
polygon(ap1.9.hist2$mids, ap1.9.hist2$density, col=alpha('violetred3',0.1))
polygon(ap1.9.hist3$mids, ap1.9.hist3$density, col=alpha('violetred3',0.1))
polygon(ap1.9.hist4$mids, ap1.9.hist4$density, col=alpha('violetred3',0.1))
polygon(ap3.7.hist1$mids, ap3.7.hist1$density, col=alpha('violetred3',0.1))
polygon(ap3.7.hist2$mids, ap3.7.hist2$density, col=alpha('violetred3',0.1))
polygon(ap3.7.hist3$mids, ap3.7.hist3$density, col=alpha('violetred3',0.1))
polygon(ap3.7.hist4$mids, ap3.7.hist4$density, col=alpha('violetred3',0.1))
polygon(ap3.7.hist5$mids, ap3.7.hist5$density, col=alpha('violetred3',0.1))
polygon(ap5.5.hist1$mids, ap5.5.hist1$density, col=alpha('violetred3',0.1))
polygon(ap5.5.hist2$mids, ap5.5.hist2$density, col=alpha('violetred3',0.1))
polygon(ap5.5.hist3$mids, ap5.5.hist3$density, col=alpha('violetred3',0.1))
polygon(ap5.5.hist4$mids, ap5.5.hist4$density, col=alpha('violetred3',0.1))
polygon(ap5.5.hist5$mids, ap5.5.hist5$density, col=alpha('violetred3',0.1))
title('A', adj=0, cex.main=3)

plot(c(min(empirical.dp.hist$mids), empirical.dp.hist$mids, max(empirical.dp.hist$mids)), 
     c(0,empirical.dp.hist$density,0), ylim=c(0,ymax.d), type='l', xlim=c(0,1),
     axes=FALSE, xlab='', ylab='', las=1)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Infection \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(min(empirical.dp.hist$mids), empirical.dp.hist$mids, max(empirical.dp.hist$mids)), c(0,empirical.dp.hist$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp1.9.hist1$mids),dp1.9.hist1$mids, max(dp1.9.hist1$mids)), c(0,dp1.9.hist1$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp1.9.hist2$mids),dp1.9.hist2$mids, max(dp1.9.hist2$mids)), c(0,dp1.9.hist2$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp1.9.hist3$mids),dp1.9.hist3$mids, max(dp1.9.hist3$mids)), c(0,dp1.9.hist3$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp1.9.hist4$mids),dp1.9.hist4$mids, max(dp1.9.hist4$mids)), c(0,dp1.9.hist4$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp1.9.hist5$mids),dp1.9.hist5$mids, max(dp1.9.hist5$mids)), c(0,dp1.9.hist5$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp3.7.hist1$mids),dp3.7.hist1$mids, max(dp3.7.hist1$mids)), c(0,dp3.7.hist1$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp3.7.hist2$mids),dp3.7.hist2$mids, max(dp3.7.hist2$mids)), c(0,dp3.7.hist2$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp3.7.hist3$mids),dp3.7.hist3$mids, max(dp3.7.hist3$mids)), c(0,dp3.7.hist3$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp3.7.hist4$mids),dp3.7.hist4$mids, max(dp3.7.hist4$mids)), c(0,dp3.7.hist4$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp3.7.hist5$mids),dp3.7.hist5$mids, max(dp3.7.hist5$mids)), c(0,dp3.7.hist5$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp5.5.hist1$mids),dp5.5.hist1$mids, max(dp5.5.hist1$mids)), c(0,dp5.5.hist1$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp5.5.hist2$mids),dp5.5.hist2$mids, max(dp5.5.hist2$mids)), c(0,dp5.5.hist2$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp5.5.hist3$mids),dp5.5.hist3$mids, max(dp5.5.hist3$mids)), c(0,dp5.5.hist3$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp5.5.hist4$mids),dp5.5.hist4$mids, max(dp5.5.hist4$mids)), c(0,dp5.5.hist4$density,0), col=alpha('blue1',0.25))
polygon(c(min(dp5.5.hist5$mids),dp5.5.hist5$mids, max(dp5.5.hist5$mids)), c(0,dp5.5.hist5$density,0), col=alpha('blue1',0.25))
title('B', adj=0, cex.main=3)

#par(mfrow=c(1,1))
plot(c(min(empirical.dabp.hist$mids), empirical.dabp.hist$mids, max(empirical.dabp.hist$mids)), 
     c(0,empirical.dabp.hist$density,0), ylim=c(0,ymax.dab), type='l', xlim=c(0,1),
     axes=FALSE, xlab='', ylab='', las=1)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='% Seroprevalence \nWhite-tailed Deer', side=1, line=6, cex=1.5)
polygon(c(min(empirical.dabp.hist$mids), empirical.dabp.hist$mids, max(empirical.dabp.hist$mids)), c(0,empirical.dabp.hist$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp1.9.hist1$mids),dabp1.9.hist1$mids, max(dabp1.9.hist1$mids)), c(0,dabp1.9.hist1$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp1.9.hist2$mids),dabp1.9.hist2$mids, max(dabp1.9.hist2$mids)), c(0,dabp1.9.hist2$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp1.9.hist3$mids),dabp1.9.hist3$mids, max(dabp1.9.hist3$mids)), c(0,dabp1.9.hist3$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp1.9.hist4$mids),dabp1.9.hist4$mids, max(dabp1.9.hist4$mids)), c(0,dabp1.9.hist4$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp1.9.hist5$mids),dabp1.9.hist5$mids, max(dabp1.9.hist5$mids)), c(0,dabp1.9.hist5$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp3.7.hist1$mids),dabp3.7.hist1$mids, max(dabp3.7.hist1$mids)), c(0,dabp3.7.hist1$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp3.7.hist2$mids),dabp3.7.hist2$mids, max(dabp3.7.hist2$mids)), c(0,dabp3.7.hist2$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp3.7.hist3$mids),dabp3.7.hist3$mids, max(dabp3.7.hist3$mids)), c(0,dabp3.7.hist3$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp3.7.hist4$mids),dabp3.7.hist4$mids, max(dabp3.7.hist4$mids)), c(0,dabp3.7.hist4$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp3.7.hist5$mids),dabp3.7.hist5$mids, max(dabp3.7.hist5$mids)), c(0,dabp3.7.hist5$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp5.5.hist1$mids),dabp5.5.hist1$mids, max(dabp5.5.hist1$mids)), c(0,dabp5.5.hist1$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp5.5.hist2$mids),dabp5.5.hist2$mids, max(dabp5.5.hist2$mids)), c(0,dabp5.5.hist2$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp5.5.hist3$mids),dabp5.5.hist3$mids, max(dabp5.5.hist3$mids)), c(0,dabp5.5.hist3$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp5.5.hist4$mids),dabp5.5.hist4$mids, max(dabp5.5.hist4$mids)), c(0,dabp5.5.hist4$density,0), col=alpha('yellow3',0.25))
polygon(c(min(dabp5.5.hist5$mids),dabp5.5.hist5$mids, max(dabp5.5.hist5$mids)), c(0,dabp5.5.hist5$density,0), col=alpha('yellow3',0.25))
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
#data.2.2<-read.csv('Multivar_Metropolis_200D_200R/Run1_200D_200R/1 MCMC Accepted Iteration Log 2014-10-11 17-46 839.csv')
data.2.2<-read.csv('MV_200_200_fixed/Multivar_Metropolis_200D_200R_1/1 MCMC Accepted Iteration Log 2014-10-21 22-38 486.csv')
data.2.2.2<-read.csv('MV_200_200_fixed/Multivar_Metropolis_200D_200R_3/1 MCMC Accepted Iteration Log 2014-10-21 22-38 265.csv')
data.2.2.3<-read.csv('MV_200_200_fixed/Multivar_Metropolis_200D_200R_4/1 MCMC Accepted Iteration Log 2014-10-21 22-38 809.csv')

par(mfrow=c(1,1))
deer.burden1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/deer_burden_chains.csv')
db1.9.1<-hist(deer.burden1.9[,2]/100, plot=FALSE)
db1.9.2<-hist(deer.burden1.9[,3]/100, plot=FALSE)
db1.9.3<-hist(deer.burden1.9[,4]/100, plot=FALSE)
db1.9.4<-hist(deer.burden1.9[,5]/100, plot=FALSE)
db1.9.5<-hist(deer.burden1.9[,6]/100, plot=FALSE)
maxdb1.9<-max(c(db1.9.1$density, db1.9.2$density, db1.9.3$density, db1.9.4$density, db1.9.5$density))

alt.burden1.9<-read.csv('~/Desktop/MV_ConstantHostPop/MV_100_900/alt_burden_chains.csv')
ab1.9.1<-hist(alt.burden1.9[,2]/900, plot=FALSE)
ab1.9.2<-hist(alt.burden1.9[,3]/900, plot=FALSE)
ab1.9.3<-hist(alt.burden1.9[,4]/900, plot=FALSE)
ab1.9.4<-hist(alt.burden1.9[,5]/900, plot=FALSE)
ab1.9.5<-hist(alt.burden1.9[,6]/900, plot=FALSE)
maxab1.9<-max(c(ab1.9.1$density, ab1.9.2$density, ab1.9.3$density, ab1.9.4$density, ab1.9.5$density))

png(file='~/Dropbox/DigitalLabNotebooks/TickModel/MiscFigures/Burdens_May2015.png', height=30, width=20, units='cm', res=300)
par(mfrow=c(3,2))
plot(db1.9.1$mids, db1.9.1$density, type='l', ylim=c(0,maxdb1.9), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(db1.9.1$mids), db1.9.1$mids, max(db1.9.1$mids)), 
        c(0,db1.9.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db1.9.2$mids),db1.9.2$mids,max(db1.9.2$mids)), 
        c(0,db1.9.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db1.9.3$mids),db1.9.3$mids,max(db1.9.3$mids)), 
        c(0,db1.9.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db1.9.4$mids),db1.9.4$mids,max(db1.9.4$mids)), 
        c(0,db1.9.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db1.9.5$mids),db1.9.5$mids,max(db1.9.5$mids)), 
        c(0,db1.9.5$density,0), col=alpha('blue3', 0.5))

plot(ab1.9.1$mids, ab1.9.1$density, type='l', ylim=c(0,maxab1.9), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(ab1.9.1$mids), ab1.9.1$mids, max(ab1.9.1$mids)), 
        c(0,ab1.9.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab1.9.2$mids),ab1.9.2$mids,max(ab1.9.2$mids)), 
        c(0,ab1.9.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab1.9.3$mids),ab1.9.3$mids,max(ab1.9.3$mids)), 
        c(0,ab1.9.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab1.9.4$mids),ab1.9.4$mids,max(ab1.9.4$mids)), 
        c(0,ab1.9.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab1.9.5$mids),ab1.9.5$mids,max(ab1.9.5$mids)), 
        c(0,ab1.9.5$density,0), col=alpha('blue3', 0.5))

deer.burden3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/deer_burden_chains.csv')
db3.7.1<-hist(deer.burden3.7[,2]/300, plot=FALSE)
db3.7.2<-hist(deer.burden3.7[,3]/300, plot=FALSE)
db3.7.3<-hist(deer.burden3.7[,4]/300, plot=FALSE)
db3.7.4<-hist(deer.burden3.7[,5]/300, plot=FALSE)
db3.7.5<-hist(deer.burden3.7[,6]/300, plot=FALSE)
maxdb3.7<-max(c(db3.7.1$density, db3.7.2$density, db3.7.3$density, db3.7.4$density, db3.7.5$density))

alt.burden3.7<-read.csv('~/Desktop/MV_ConstantHostPop/MV_300_700/alt_burden_chains.csv')
ab3.7.1<-hist(alt.burden3.7[,2]/700, plot=FALSE)
ab3.7.2<-hist(alt.burden3.7[,3]/700, plot=FALSE)
ab3.7.3<-hist(alt.burden3.7[,4]/700, plot=FALSE)
ab3.7.4<-hist(alt.burden3.7[,5]/700, plot=FALSE)
ab3.7.5<-hist(alt.burden3.7[,6]/700, plot=FALSE)
maxab3.7<-max(c(ab3.7.1$density, ab3.7.2$density, ab3.7.3$density, ab3.7.4$density, ab3.7.5$density))


plot(db3.7.1$mids, db3.7.1$density, type='l', ylim=c(0,maxdb3.7), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(db3.7.1$mids), db3.7.1$mids, max(db3.7.1$mids)), 
        c(0,db3.7.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db3.7.2$mids),db3.7.2$mids,max(db3.7.2$mids)), 
        c(0,db3.7.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db3.7.3$mids),db3.7.3$mids,max(db3.7.3$mids)), 
        c(0,db3.7.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db3.7.4$mids),db3.7.4$mids,max(db3.7.4$mids)), 
        c(0,db3.7.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db3.7.5$mids),db3.7.5$mids,max(db3.7.5$mids)), 
        c(0,db3.7.5$density,0), col=alpha('blue3', 0.5))

plot(ab3.7.1$mids, ab3.7.1$density, type='l', ylim=c(0,maxab3.7), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(ab3.7.1$mids), ab3.7.1$mids, max(ab3.7.1$mids)), 
        c(0,ab3.7.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab3.7.2$mids),ab3.7.2$mids,max(ab3.7.2$mids)), 
        c(0,ab3.7.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab3.7.3$mids),ab3.7.3$mids,max(ab3.7.3$mids)), 
        c(0,ab3.7.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab3.7.4$mids),ab3.7.4$mids,max(ab3.7.4$mids)), 
        c(0,ab3.7.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab3.7.5$mids),ab3.7.5$mids,max(ab3.7.5$mids)), 
        c(0,ab3.7.5$density,0), col=alpha('blue3', 0.5))

deer.burden5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/deer_burden_chains.csv')
db5.5.1<-hist(deer.burden5.5[,2]/500, plot=FALSE)
db5.5.2<-hist(deer.burden5.5[,3]/500, plot=FALSE)
db5.5.3<-hist(deer.burden5.5[,4]/500, plot=FALSE)
db5.5.4<-hist(deer.burden5.5[,5]/500, plot=FALSE)
db5.5.5<-hist(deer.burden5.5[,6]/500, plot=FALSE)
maxdb5.5<-max(c(db5.5.1$density, db5.5.2$density, db5.5.3$density, db5.5.4$density, db5.5.5$density))

alt.burden5.5<-read.csv('~/Desktop/MV_ConstantHostPop/MV_500_500/alt_burden_chains.csv')
ab5.5.1<-hist(alt.burden5.5[,2]/500, plot=FALSE)
ab5.5.2<-hist(alt.burden5.5[,3]/500, plot=FALSE)
ab5.5.3<-hist(alt.burden5.5[,4]/500, plot=FALSE)
ab5.5.4<-hist(alt.burden5.5[,5]/500, plot=FALSE)
ab5.5.5<-hist(alt.burden5.5[,6]/500, plot=FALSE)
maxab5.5<-max(c(ab5.5.1$density, ab5.5.2$density, ab5.5.3$density, ab5.5.4$density, ab5.5.5$density))

plot(db5.5.1$mids, db5.5.1$density, type='l', ylim=c(0,maxdb5.5), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(db5.5.1$mids), db5.5.1$mids, max(db5.5.1$mids)), 
        c(0,db5.5.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db5.5.2$mids),db5.5.2$mids,max(db5.5.2$mids)), 
        c(0,db5.5.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db5.5.3$mids),db5.5.3$mids,max(db5.5.3$mids)), 
        c(0,db5.5.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db5.5.4$mids),db5.5.4$mids,max(db5.5.4$mids)), 
        c(0,db5.5.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(db5.5.5$mids),db5.5.5$mids,max(db5.5.5$mids)), 
        c(0,db5.5.5$density,0), col=alpha('blue3', 0.5))

plot(ab5.5.1$mids, ab5.5.1$density, type='l', ylim=c(0,maxab5.5), axes=FALSE,
     xlab='', ylab='', las=1, cex.lab=2, cex.axis=2)
axis(side=1, cex.axis=2, cex.lab=2)
mtext(text='Average Number Ticks per Deer', side=1, line=6, cex=1.5)
polygon(c(min(ab5.5.1$mids), ab5.5.1$mids, max(ab5.5.1$mids)), 
        c(0,ab5.5.1$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab5.5.2$mids),ab5.5.2$mids,max(ab5.5.2$mids)), 
        c(0,ab5.5.2$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab5.5.3$mids),ab5.5.3$mids,max(ab5.5.3$mids)), 
        c(0,ab5.5.3$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab5.5.4$mids),ab5.5.4$mids,max(ab5.5.4$mids)), 
        c(0,ab5.5.4$density,0), col=alpha('blue3', 0.5))
polygon(c(min(ab5.5.5$mids),ab5.5.5$mids,max(ab5.5.5$mids)), 
        c(0,ab5.5.5$density,0), col=alpha('blue3', 0.5))

dev.off()




hist(deer.burden24/200, main='1:2', xlim=c(0,35))
hist(alt.burden24/400, xlim=c(0,35))
hist(deer.burden/200, main='1:1', xlim=c(0,35), breaks=25)
hist(alt.burden/200, xlim=c(0,35), breaks=25)
hist(deer.burden42/400, main='2:1', xlim=c(0,35))
hist(alt.burden42/200, xlim=c(0,35))

png(file='~/Dropbox/ModelFitting/FutureProof/Figures for Presentations/Dissertation_DeerBurden_Posteriors.png',
    height=30, width=20, units='cm', res=300)
par(mfrow=c(3,1), mar=c(5,5,2,2))

db24<-hist(deer.burden24/200, breaks=40, plot=FALSE)
db24.coords<-cbind(db24$mids, db24$density)
db24.95<-quantile(deer.burden24/200, probs=c(0.025, 0.975))
db24.95.data<-subset(db24.coords, db24.coords[,1] > db24.95[1] & db24.coords[,1] < db24.95[2])
plot(db24$mids, db24$density, xlim=c(0,30), xlab='',
     type='l', ylab='', axes=FALSE, cex.axis=1.5, cex.lab=1.5)
axis(side=1, cex.axis=1.5)
mtext(side=2, line=1.5, text='Scenario 1', cex=1.3)
polygon(db24$mids, db24$density, col=alpha('gray', 0.5))
polygon(c(min(db24.95.data[,1]),db24.95.data[,1],max(db24.95.data[,1])), 
        c(0,db24.95.data[,2],0), col=alpha('gray', 0.5))

db22<-hist(deer.burden/200, breaks=30, plot=FALSE)
db22.coords<-cbind(db22$mids, db22$density)
db22.95<-quantile(deer.burden/200, probs=c(0.025, 0.975))
db22.95.data<-subset(db22.coords, db22.coords[,1] > db22.95[1] & db22.coords[,1] < db22.95[2])
plot(db22$mids, db22$density, xlim=c(0,30), xlab='',
     type='l', ylab='', axes=FALSE, cex.axis=1.5, cex.lab=1.5)
axis(side=1, cex.axis=1.5)
mtext(side=2, line=1.5, text='Scenario 2', cex=1.3)
polygon(c(min(db22$mids),db22$mids), c(0,db22$density), col=alpha('gray', 0.5))
polygon(c(min(db22.95.data[,1]),db22.95.data[,1],max(db22.95.data[,1])), 
        c(0,db22.95.data[,2],0), col=alpha('gray', 0.5))

db42<-hist(deer.burden42/200, breaks=40, plot=FALSE)
db42.coords<-cbind(db42$mids, db42$density)
db42.95<-quantile(deer.burden42/200, probs=c(0.025, 0.975))
db42.95.data<-subset(db42.coords, db42.coords[,1] > db42.95[1] & db42.coords[,1] < db42.95[2])
plot(db42$mids, db42$density, xlim=c(0,30), xlab='Number of Ticks Per Deer',
     type='l', ylab='', axes=FALSE, cex.axis=1.5, cex.lab=1.5)
axis(side=1, cex.axis=1.5)
mtext(side=2, line=1.5, text='Scenario 3', cex=1.3)
polygon(c(min(db42$mids), db42$mids), c(0,db42$density), col=alpha('gray', 0.5))
polygon(c(min(db42.95.data[,1]),db42.95.data[,1],max(db42.95.data[,1])), 
        c(0,db42.95.data[,2],0), col=alpha('gray', 0.5))

dev.off()