## Visualize results from the toy model for tick borne disease transmission

setwd('~/Desktop/UT_ModelFitting/')

ticks<-read.csv('toy_model_scenarios_tick_out.csv')
deer<-read.csv('toy_model_scenarios_deer_out.csv')
deerAB<-read.csv('toy_model_scenarios_deerAB_out.csv')

#par(mfrow=c(3,1), mar=c(5,12,2,2))
pdf(file='~/Dropbox/EEID2015/Toy_Model_TickPrev_Motivation.pdf', height=15/2.54, width=15/2.54)
par(mfrow=c(1,1), mar=c(5,6,2,2))
plot(ticks$prefs, (ticks$Scenario_1)*100, type='l', xlab='Fraction of Ticks Feeding on Reservoir', ylab='',
     #lty=3, 
     col='#54278f', lwd=4,
     cex.lab=1.5, cex.axis=1.5, las=1, ylim=c(0, 3), axes=FALSE)
axis(side=1, cex.axis=1.5)
axis(side=2, las=1, cex.axis=1.5, at=c(0, 1, 2, 3))
#mtext(text='Infection \nPrevalence \nin Ticks', side=2, line=5, cex=1.5)
mtext(text='Infection Prevalence in Ticks (%)', side=2, line=2, cex=1.5)
lines(ticks$prefs, (ticks$Scenario_2)*100, 
      #lty=4)
      col='#807dba', lwd=4)
lines(ticks$prefs, (ticks$Scenario_3)*100, 
      #lty=1)
      col='#bcbddc', lwd=4)
legend(x='topleft', legend=c('1:9', '3:7', '5:5'), 
       cex=1.2, bty='n', title='Reservoir:Non-Reservoir',
       #lty=c(1,4,3))
       fill=c('#bcbddc', '#807dba', '#54278f'))
dev.off()

plot(deer$prefs, deer$Scenario_1, type='l', xlab='', ylab='', axes=FALSE,
     lty=3, cex.lab=1.5, cex.axis=1.5, las=1, ylim=c(0, 0.2))
axis(side=1, cex.axis=1.5)
axis(side=2, las=1, cex.axis=1.5, at=c(0, 0.05, 0.10, 0.15, 0.2))
mtext(text='Infection \nPrevalence \nin Deer', side=2, line=5, cex=1.5)
lines(deer$prefs, deer$Scenario_2, lty=4)
lines(deer$prefs, deer$Scenario_3, lty=1)

plot(deerAB$prefs, deerAB$Scenario_1, type='l', xlab='Tick Preference for Deer', ylab='',
     lty=3, cex.lab=1.5, cex.axis=1.5, las=1, ylim=c(0, 1), axes=FALSE)
axis(side=1, cex.axis=1.5)
axis(side=2, las=1, cex.axis=1.5)
mtext(text='Antibody \nPrevalence \nin Deer', side=2, line=5, cex=1.5)
lines(deerAB$prefs, deerAB$Scenario_2, lty=4)
lines(deerAB$prefs, deerAB$Scenario_3, lty=1)
legend(x='bottomright', legend=c('Scenario 1', 'Scenario 2', 'Scenario 3'), 
       cex=1.5, lty=c(3,4,1), bty='n')


