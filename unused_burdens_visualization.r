
burden.plot<-function(burden.data, burden.name=NULL, y.label=NULL, color='gray', xlim=c(0,60),
                      pop.size, add.line=NULL, add.line.color=NULL, alpha.lvl=1){
  hist.list<-list()
  i=1
  for(chain in 1:length(burden.data)){
    hist.list[[chain]]<-hist((burden.data[,chain])/pop.size[i], plot=FALSE, prob=TRUE)
    if(length(pop.size)>1){
      i=i+1
    }
  }
  ymax=0
  xmax=0
  for(hist in 1:length(hist.list)){
    if(max(hist.list[[hist]]$density)>ymax){
      ymax=max(hist.list[[hist]]$density)
    }
    if(max(hist.list[[hist]]$mids)>xmax){
      xmax=max(hist.list[[hist]]$mids)
    }
  }
  #print(xmax)
  plot(hist.list[[1]]$mids, hist.list[[1]]$density, type='l', ylim=c(0, ymax),
       xlim=xlim, xlab=burden.name, axes=FALSE, cex.lab=1.5, ylab='', col='white')
  for(hist in 1:length(hist.list)){
    polygon(x=c(min(hist.list[[hist]]$mids), hist.list[[hist]]$mids, max(hist.list[[hist]]$mids)),
            y=c(0, hist.list[[hist]]$density, 0), col=alpha(color[hist], alpha.lvl))
  }
  axis(side=1, cex.axis=1.5)
  mtext(text=y.label, side=2, line=0.5, cex=1.5)
  if(!is.null(add.line)){
    abline(v=add.line, col=add.line.color, lty=2, lwd=3)
  }
  #return(hist.list)
}

burden.poly.plot<-function(burden.data, burden.name=NULL, y.label=NULL, color='gray', xlim=c(0,60),
                           pop.size, add.line=NULL, add.line.color=NULL, alpha.lvl=1){
  hist.list<-list()
  i=1
  for(chain in 1:length(data.frame(burden.data))){
    hist.list[[chain]]<-hist((burden.data[,chain])/pop.size[i], plot=FALSE)
    if(length(pop.size)>1){
      i=i+1
    }
  }
  ymax=0
  xmax=0
  for(hist in 1:length(hist.list)){
    if(max(hist.list[[hist]]$density)>ymax){
      ymax=max(hist.list[[hist]]$density)
      print(ymax)
    }
    if(max(hist.list[[hist]]$mids)>xmax){
      xmax=max(hist.list[[hist]]$mids)
    }
  }
  #print(ymax)
  plot(x=NULL, y=NULL, type='l', ylim=c(0, ymax),
       xlim=xlim, xlab=burden.name, axes=FALSE, cex.lab=1.5, ylab='', col='white')
  for(j in 1:length(burden.data)){
    polygon(density(burden.data[,j]/pop.size[j]), col=alpha(color[j], alpha.lvl))
  }
  axis(side=1, cex.axis=1.5)
  mtext(text=y.label, side=2, line=0.5, cex=1.5)
  if(!is.null(add.line)){
    abline(v=add.line, col=add.line.color, lty=2, lwd=3)
  }
  #return(hist.list)
}


burden.plot(deer.burden1.9[-1], 'Number of Ticks per Deer', '', 'slategray3', pop.size=100, alpha.lvl=0.1)
burden.plot(deer.burden3.7[-1], 'Number of Ticks per Deer', '', 'slategray3', pop.size=300, alpha.lvl=0.1)
burden.plot(deer.burden5.5[-1], 'Number of Ticks per Deer', '', 'slategray3', pop.size=500, alpha.lvl=0.1)
burden.plot(alt.burden1.9[-1], 'Number of Ticks per Alt Host', '', 'slategray3', pop.size=900, alpha.lvl=0.1)
burden.plot(alt.burden3.7[-1], 'Number of Ticks per Alt Host', '', 'slategray3', pop.size=700, alpha.lvl=0.1)
burden.plot(alt.burden5.5[-1], 'Number of Ticks per Alt Host', '', 'slategray3', pop.size=500, alpha.lvl=0.1)

## barplot option
aggregate.burdens<-cbind(deer.burden1.9[,2]/100, alt.burden1.9[,2]/900,
                         deer.burden3.7[,2]/300, alt.burden3.7[,2]/700,
                         deer.burden5.5[,2]/700, alt.burden5.5[,2]/300)
boxplot(aggregate.burdens, col=c(rep(scenario.colors[1], 2), rep(scenario.colors[2], 2), rep(scenario.colors[3], 2)),
        las=1, names=rep(c('Deer', 'Alternate'),3), ylab='Average Number Ticks per Animal')

## violin plot option
library(vioplot)
source('~/Desktop/UT_ModelFitting/multicolor_vioplot.r')
multicolor.vioplot(deer.burden1.9[,2]/100, alt.burden1.9[,2]/900, 
                   deer.burden3.7[,2]/300, alt.burden3.7[,2]/700,
                   deer.burden5.5[,2]/500, alt.burden5.5[,2]/500,
                   col=c(rep(scenario.colors[1], 2), rep(scenario.colors[2], 2), rep(scenario.colors[3], 2)))

## histogram/density plot option
burden1.9<-data.frame(cbind(deer.burden1.9[,2], alt.burden1.9[,2]))
burden3.7<-data.frame(cbind(deer.burden3.7[,2], alt.burden3.7[,2]))
burden5.5<-data.frame(cbind(deer.burden5.5[,2], alt.burden5.5[,2]))
burden.name=c('Average Number of Ticks per Host')
par(mfrow=c(4,1), mar=c(3,2,3,2))
burden.poly.plot(burden1.9, '', color=c('blue', 'red'), pop.size=c(100,900), alpha.lvl=0.25)
burden.poly.plot(burden3.7, '', color=c('blue', 'red'), pop.size=c(300,700), alpha.lvl=0.25)
burden.poly.plot(burden5.5, burden.name, color=c('blue', 'red'), pop.size=c(500,500), alpha.lvl=0.25)
