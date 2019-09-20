#Rethink cat predation model (based on Loss et al 2013)
rm(list=ls())
#owned cat predation rates
pred=read.csv("~/kim/reThink/catpredrecent/ownedpred.csv")
xpred=pred[pred$exclude!="X",]

library("fitdistrplus")
abirds=xpred$annbirds
fit=fitdist(abirds, "lnorm")
fit2=fitdist(abirds, "norm")
fit3=fitdist(abirds, "unif")
par(mfrow=c(1,2))
denscomp(list(fit,fit2,fit3), addlegend = FALSE, main = "", xlab = "birds per owned cat per year", fitcol = c("blue", "red", "green"), 
         xlim = c(0,12), 
         ylim=c(0, .3),lwd=2, datacol="gray", cex.lab=1.5)

mams=xpred$annmammals[!is.na(xpred$annmammals)]
fitm=fitdist(mams, "lnorm")
fitm2=fitdist(mams, "norm")
fitm3=fitdist(mams, "unif")
denscomp(list(fitm,fitm2, fitm3), addlegend = FALSE, main = "", xlab = "mammals per owned cat per year", fitcol = c("blue", "red", "green"), 
         xlim = c(0,35), ylim=c(0, .08),
        lwd=2, datacol="gray", cex.lab=1.5)
legend(15,0.08, c("log norm", "norm", "unif"),lwd=2,cex=1.5,bty="n",lty=c(1,2), col=c("blue", "red", "green"))

#cat population
n=10000 # number of simulations
UScatsurvey=c(54*10^6, 54*10^6, 94*10^6)
cmean=mean(UScatsurvey)
csd=sd(UScatsurvey)
Ccat=rnorm(n, mean=8.3*10^6, sd=.25*10^6) #Owned cats
UScat=rnorm(n, mean=cmean, sd=csd)
NAcat=Ccat+UScat
Nmean=mean(NAcat)
Nsd=sd(NAcat)

#our model
cat=NAcat#Owned cats
out=runif(n, min=0.25, max=0.5) #proporition of owned cats that go outdoors
cathunt=runif(n,min=0.5, max=0.8) #proportion of outdoor owned cats that hunt
crt=runif(n,min=1.2, max=3.3) # Correction for owned cats not returning prey
catus=(cat*out)*cathunt #number of owned outdoor hunting US cats

#bird predation

predb=rlnorm(n, meanlog = log(exp(fit$estimate[1])*1), sdlog = fit$estimate[2])
birds=catus*(predb*crt) #annual birds per owned outdoor hunting US cats


#mammal predation
predm=rlnorm(n, meanlog = log(exp(fitm$estimate[1])*1), sdlog = fitm$estimate[2])
mamms=catus*(predm*crt) #annual mammals per owned outdoor hunting US cats

#median bird kills
median(birds)/10^6
lossdiff=684*10^6-median(birds)#Difference between my median and Loss reported
lossdiff/10^6
ql=quantile(birds, c(0.05, 0.95))/10^6
ql

#median mammal kills
median(mamms)/10^6
lossdiff=1249*10^6-median(mamms)#Difference between my median and Loss reported
lossdiff/10^6
ql=quantile(mamms, c(0.05, 0.95))/10^6
ql

#cost effect of confine
percat=(birds+mamms+2.4+6.3)/(cat*out)
ancon=percat*1000/30
median(ancon)
qlc=quantile(ancon, c(0.05, 0.95))
qlc



#sensitivity analysis for birds
dmb=data.frame(cat, out, cathunt, crt, predb)
dmb2=cbind(birds, dmb)
lrel=lm(birds~., data=dmb2)
af <- anova(lrel)
afss <- af$"Sum Sq"
afr2=(cbind(af,PctExp=afss/sum(afss)))
afr2$adj=1-((1-afr2$PctExp)*(n-1)/(n-ncol(dmb)-1))


dmm=data.frame(cat, out, cathunt, crt, predm)
dmm2=cbind(mamms, dmm)
mrel=lm(mamms~., data=dmm2)
mf <- anova(mrel)
mfss <- mf$"Sum Sq"
mfr2=(cbind(af,PctExp=mfss/sum(mfss)))
mfr2$adj=1-((1-mfr2$PctExp)*(n-1)/(n-ncol(dmb)-1))
par(mfrow=c(1,1))
factors<-c("cat", "out", "cathunt", "return", "pred")
barplot(rbind(mfr2$adj[1:5],afr2$adj[1:5]), names=factors, cex.lab=1.5,
        beside=TRUE,ylab="% variance explained", ylim=c(0, 0.65))
legend("topleft", c("birds", "mammals"), fill=c("black", "gray"),
       cex=1.5, bty = "n")

library(sensitivity)
sob=sensitivity::src(dmb, birds, nboot=50)
sob2<-sob$SRC
pob<-sensitivity::pcc(birds, dmb, rank=TRUE, nboot=50)
pob2<-pob$PRCC


plot(pob2$original, ylab="regression and correlation coefficients", bty="l", 
     col="red", xlab="", xaxt = "n",cex.lab=1.5, cex=2, pch=15, 
     #ylim=c(-0.1, 1), 
     xlim=c(1,6))
bfacts=c("pop", "out", "hunters", "prey return", "pred")
axis(1, at=1:5, labels=bfacts)
arrows( 1:5,pob$`min. c.i.`,  1:5, pob$`max. c.i.`, 
        length=0.05, angle=90, code=3)
xseq=seq(from=1.25, to=5.25, by=1)
points(xseq,sob2$original, 
     col="blue", cex=2, pch=19)
arrows( xseq,sob2$`min. c.i.`,  xseq, sob2$`max. c.i.`, 
        length=0.05, angle=90, code=3)
abline(h=0, lty=3)
legend("bottomright", c("prcc", "src"), pch=c(15,19), col=c("red", "blue"),
       bty="n", cex=1.5)
