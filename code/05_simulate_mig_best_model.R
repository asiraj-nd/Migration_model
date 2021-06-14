##############################################################################################
#### ##### ### # ##### ###          ### ## ######     ######## ######## ### ######### # # # #
### ## ########      ######### ######## ### ######### # # # ##### ##### ### # ##### ###    #### ###  
######                          simulate                  ############
######                        best model                ############ 
##### # # # ##### ##### ### # ##### ###    #### ###  ######### ######## ### ########
#############################################################################################

#rm(list=ls())
#########################################################################################
setwd("~/Dropbox/tempx/COL_mig_model/code")
#setwd("C:/Users/asiraj/Dropbox/tempx/COL_mig_model/code")

if(!require(inline)){
  install.packages('inline'); library(inline)
}
if(!require(BayesianTools)){
  install.packages('BayesianTools'); library(BayesianTools)
}


citation("BayesianTools")
#set.seed(123)
sessionInfo()


ncols = 18

rm(list=ls())
#########################################################################################
setwd("~/Dropbox/tempx/COL_mig_model/code")
#setwd("C:/Users/asiraj/Dropbox/tempx/COL_mig_model/code")

if(!require(inline)){
  install.packages('inline'); library(inline)
}
if(!require(BayesianTools)){
  install.packages('BayesianTools'); library(BayesianTools)
}




citation("BayesianTools")

#set.seed(123)
sessionInfo()

ncols = 18


#########################################################################################
###### SET UP THE DATA MATRICES 
#########################################################################################

deptused  = read.table("./COL_dept_mig.dat", sep=" " ,skip=0, header=FALSE)
munidata  = read.table(paste("./COL_dist_feed.dat",sep=""), sep=" " ,skip=0, header=FALSE)
head(munidata)
thisexcl<- c(4:14)
bauexcl= thisexcl

muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))

row.num=1
deptused<- matrix(as.vector(unlist(deptused)),ncol=2)

source("./00_numfunctions.R")
# -----------------------------------------
given.params = 
  c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, 0.375, -9.7)[-bauexcl]

pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = c(rep(-10, (length(given.params)-1)),-15), upper = rep(10, length(given.params)))

mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)

deptused<- matrix(as.vector(unlist(deptused)),ncol=2)

row.num=1
deptused<- matrix(as.vector(unlist(deptused)),ncol=2)


models_1 = c(3,16,37,41,57,67,75,84,94,101)
summaryexcls<- matrix(as.vector(unlist(read.csv("../model_inputs//summary_excl.csv", stringsAsFactors = F))),nrow=12)

### multiple runs
# -----------------------------------------
iter = 300000
#head(dataused)

### adaptive mcmc prior optimization
settings <- list(iterations = iter, adapt = T, nrChains = 1, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = TRUE)

#bestmodel = model 84, c(3,16,37,41,57,67,75,84,94,101)[8]

  kk=8
  model =   models_1[kk]
  run= 1
  
  bauexcl= summaryexcls[kk,] [which(!is.na(summaryexcls[kk,]))]
  
  dim(deptused)
  #bauexcl= which(bauexcl==99)
  
  muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))
  head(muniused)
  source("./00_numfunctions.R")
  
  #########################################################################################
  ###### Multi-variate MCMC 
  #########################################################################################
  
  #given.params = c(-3.25,0.02,0.11,-7.6)
  given.params = 
    c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, -0.600,0.375, -9.7)[-bauexcl]
  
  # -----------------------------------------
  pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = c(rep(-10, (length(given.params)-1)),-15), upper = rep(10, length(given.params)))
  
  mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
  bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)
  
  #set.seed(123)
  
  ### multiple runs
  # -----------------------------------------
  iter = 300000
  #head(dataused)
  models_1
  ### adaptive mcmc prior optimization
  settings <- list(iterations = iter, adapt = T, nrChains = 1, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = TRUE)
  load(paste('../output/mcmc_v',model,'_run',run,'.RData',sep=''))
  summary(out.1)
  
  postSamp<- getSample(out.1, start = 16650, end = NULL, thin = 10, whichParameters = 1:(16-length(bauexcl)))
  head(postSamp) 
  
  ########## simulate based on sample districts (n=1122)
  
  simulateSample<- NULL
  
  for (jj in 1: nrow(postSamp)){
    if (jj%%10==0) {
      simulateSample<- cbind(simulateSample, simulagg.func(postSamp[jj,]))
      print(jj)
    }
  }
  
dim(simulateSample)

write.csv(simulateSample, "../processed_data/simul_2_5k_agg_dept.csv", row.names=F, quote=F)



########################### plot scatter
deptdata = read.csv("../model_inputs/COL_dept_datamaster.csv", stringsAsFactors = F)
head(deptdata)
dim(deptdata)

simulateSample<- read.csv("../processed_data/simul_2_5k_agg_dept.csv", stringsAsFactors = F)
medPostSamp<- apply(postSamp,2,median)  
write.csv(medPostSamp, "../processed_data/model_84_median_parames.csv", row.names=F, quote=F)
head(simulateSample)

######################## simulate at district level
#munidata<- read.csv("../model_inputs/COL_dist_1122_datamaster.csv", stringsAsFactors = F)
simulateMuniSampleMed<-simul.func(medPostSamp)
#modelmig_muni<- simulateMuniSampleMed * munidata$POPi
#write.csv(simulateMuniSampleMed, "../processed_data/simul_median_muni_mig_population.csv", row.names=F, quote=F)
sum(modelmig_muni)
sum(modelmig)
sum(simulateMuniSampleMed)

length(simulateMuniSampleMed)

head(muniused)

###### simulate median migration probability
medPostSamp<- as.vector(unlist(medPostSamp))
simulateSampleMed<- simulagg.func(medPostSamp)
write.csv(simulateSampleMed, "../processed_data/simul_median_agg_dept.csv", row.names=F, quote=F)

simulateSampleMed<- read.csv("../processed_data/simul_median_agg_dept.csv", stringsAsFactors = F)

####### convert to mig population
modelmig<- (simulateSampleMed * deptdata$POPI)[,1]
head(modelmig)
#write.csv(modelmig, "../processed_data/simul_median_agg_dept_mig_population.csv", row.names=F, quote=F)
#modelmig<- read.csv("../processed_data/simul_median_agg_dept_mig_population.csv", stringsAsFactors = F)[,1]
sum(modelmig)
datamig<- deptdata$MIGIJ/ deptdata$IPUMPSPOP * deptdata$POPI
which.max(modelmig)
mxmig<- round(vec_to_mx_byrow(modelmig,78),0)
head(mxmig[,1:10])

plot(datamig, modelmig)  
cor(datamig, modelmig)^2

dd<- datamig
mm<- modelmig
ylm=c(0,12)

png(filename = paste("../plot/fig3_a_v3.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
  layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()

sqrt(.56)
0.744546^2
########## posterior samples

vecid<- 0
layout(matrix(1:10,nrow=5, byrow=TRUE))
par(mar=c(2.5,.5,1,.5))
allqts<- NULL
prob0<- NULL
#layout(1)
#j=2
dim(postSamp)
postSamp1<- postSamp

ttl0<- c('DIST','POP','POP','CONT','UrbanProp','UrbanProp','Perc','Perc','Tiny','Tiny',
         'MajCen','MajCen','GECON','GECON','InREG', 'Intercept')

vsuball<- c('ij','i','j','ij','i','j','i','j','i','j',
            'i','j','i','j','ij',"")
vnames<- ttl0[-bauexcl]
vsub<- vsuball[-bauexcl]

pdf(file = paste0('../plot/mrg_like_v',model,'_run',run,'.pdf'), height=11, width=6.5)
  marginalPlot(out.1)
dev.off()
  
summary(out.1)

summaryexcls

xlms<- c(-2.5,1)
xlms<- rbind(xlms,c(-0.68,1))
xlms<- rbind(xlms,c(-0.136,0.6))
xlms<- rbind(xlms,c(-2,3.5))
xlms<- rbind(xlms,c(-0.6,0.6))
xlms<- rbind(xlms,c(-1,1.5))
xlms<- rbind(xlms,c(-0.8,.5))
xlms<- rbind(xlms,c(-0.8,.5))
xlms<- rbind(xlms,c(-1,1.75))
xlms<- rbind(xlms,c(-1.25,.5))
xlms<- rbind(xlms,c(-0.5,0.6))
xlms<- rbind(xlms,c(-11,-10))

j=1
allqts<- NULL
prob0<- NULL

png(filename = paste('../plot/Figure2_mx_new.png',sep=""), width = 7.25, height = 9.25, units = 'in', res = 400)
  layout(matrix(1:12,nrow=6, byrow=TRUE))
  par(mar=c(2.75,.5,1,.5))
  for(j in 1:(16-length(bauexcl))){
    ylm<- c(0,max(density(postSamp1[,j])$y))
    plot(density(postSamp1[,j]), main=bquote(paste(.(vnames[j])[.(vsub[j])])), xlim = xlms[j,],xlab="", ylim= ylm, yaxt='n')
    pr<-mle.norm(cbind(density(postSamp1[,j])$y,density(postSamp1[,j])$x))
    allqts<- rbind(allqts, as.vector(unlist(quantile(postSamp1[,j], probs=c(.025, 0.5, 0.975)))))
    prob0<- c(prob0,1-pnorm(0, mean=pr[1], sd=pr[2]))
    abline(v=0, lty=2, lwd=.5)
    polygon(
      c(density(postSamp1[,j])$x,0,-100),
      c(density(postSamp1[,j])$y,0,0),
      col=rgb(1,0,0,.5),border=NA)
    axis(1,at=seq(0,20,5))
  }

dev.off()







##############################################################
### characterizing routes

zoneattrib<-read.csv("../data/COL_dept_collector.csv", stringsAsFactors = F)
head(zoneattrib)

##### routes to urban in red
urbprop<- zoneattrib$urb_Pop07/ (zoneattrib$urb_Pop07+zoneattrib$rur_Pop07)
tourb<- which(urbprop[deptdata$CODEJ]>0.5)


ddurb<- datamig[tourb]
mmurb<- modelmig[tourb]

png(filename = paste("../plot/fig3_a_urb.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
  layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()


##### routes to urban in red
urbprop<- zoneattrib$urb_Pop07/ (zoneattrib$urb_Pop07+zoneattrib$rur_Pop07)
tourb<- which(urbprop[deptdata$CODEJ]>0.5)
tourb<- tourb[-which(tourb %in% 2182)]

ddurb<- datamig[tourb]
mmurb<- modelmig[tourb]

png(filename = paste("../plot/fig3_a_urb.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
  layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()


##### routes to rich zones
highinc<- which(zoneattrib$gecon[deptdata$CODEJ]> quantile(zoneattrib$gecon[deptdata$CODEJ], p=0.9))
highinc<- highinc[-which(highinc %in% 2182)]

ddurb<- datamig[highinc]
mmurb<- modelmig[highinc]

png(filename = paste("../plot/fig3_a_highinc.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()


##### routes to poor zones
lowinc<- which(zoneattrib$gecon[deptdata$CODEJ]> quantile(zoneattrib$gecon[deptdata$CODEJ], p=0.1))
lowinc<- lowinc[-which(lowinc %in% 2182)]

ddurb<- datamig[lowinc]
mmurb<- modelmig[lowinc]

png(filename = paste("../plot/fig3_a_lowinc.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
  layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()

head(deptdata)

##### routes to within region
inreg<- which(deptdata$InREGIJ ==1)
#inreg<- inreg[-which(inreg %in% 2182)]

ddurb<- datamig[inreg]
mmurb<- modelmig[inreg]

png(filename = paste("../plot/fig3_a_inreg.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()

##### routes to highpop
highpop<- which(zoneattrib$pop_2013[deptdata$CODEJ]> quantile(zoneattrib$gecon[deptdata$CODEJ], p=0.9))
highpop<- highpop

ddurb<- datamig[highpop]
mmurb<- modelmig[highpop]

png(filename = paste("../plot/fig3_a_highpop.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()

deptdata[deptdata$CODEJ,][2182,]

##### routes to contiguous
contreg<- which(deptdata$CONTIJ ==1)
contreg<- contreg[-which(contreg %in% 2182)]

ddurb<- datamig[contreg]
mmurb<- modelmig[contreg]

png(filename = paste("../plot/fig3_a_contig.png",sep=""), width = 5, height = 5, units = 'in', res = 400)
layout(1)
  plot(log1p(dd), log1p(mm), col=rgb(0,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  abline(a=0,b=1, col='red', ylim= ylm, xlim=ylm)
  points(log1p(ddurb), log1p(mmurb), col=rgb(1,0,0,.3), ylim= ylm, xlim=ylm, pch =20, ylab="predicted (log)", xlab="observed (log)")
  text(2,12,paste("Cor = ", round(cor(mm,dd),2),sep=""))
  mtext("2009-2013 Predicted vs Observed, Pop-weighted Migration\n aggregated at 78 zones x 78 zones", side=3, line=2)
dev.off()

##########
cor(datamig, modelmig)







###################################################################
############# Misc


pdf(file = paste0('../plot/mrg_like_v',model,'_run',run,'.pdf'), height=11, width=6.5)
  marginalPlot(out.1)
dev.off()

########################### verify data is correct

####### continuous variables should have 0 mean and 0.5 SD
mean(muniused[,1])
sd(muniused[,1])

###### dummy variables should have 0 mean and 1 range
mean(muniused[,4])
diff(range(muniused[,4]))
table(muniused[,4])

mean(muniused[,10])
diff(range(muniused[,10]))
table(muniused[,10])


pdf(file = paste0('../plot/chain__v',model,'_run',run,'.pdf'), height=11, width=6.5)
  plot(out.1)
dev.off()

