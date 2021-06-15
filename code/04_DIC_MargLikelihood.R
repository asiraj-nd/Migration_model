##############################################################################################
#### ##### ### # ##### ###          ### ## ######     ######## ######## ### ######### # # # #
### ## ########      ######### ######## ### ######### # # # ##### ##### ### # ##### ###    #### ###  

######                          Analysis                  ############
######                        Generate DIC                ############ 

######  ##  ###            # Amir Siraj, Alex Perkins
######        ##        # asiraj@nd.edu, taperkins@nd.edu
######                       #   May 2020
##### # # # ##### ##### ### # ##### ###    #### ###  ######### ######## ### ########
#############################################################################################

setwd("~/Dropbox/tempx/COL_mig_model/code")
#setwd("C:/Users/asiraj/Dropbox/tempx/COL_mig_model/code")

if(!require(inline)){
  install.packages('inline'); library(inline)
}
if(!require(BayesianTools)){
  install.packages('BayesianTools'); library(BayesianTools)
}

if(!require(mvtnorm)){
  install.packages('mvtnorm'); library(mvtnorm)
}



citation("BayesianTools")

#set.seed(123)
sessionInfo()


ncols = 18


#########################################################################################
###### SET UP the input feeder for a particular model
#########################################################################################

deptused  = read.table("./COL_dept_mig.dat", sep=" " ,skip=0, header=FALSE)
munidata  = read.table(paste("./COL_dist_feed.dat",sep=""), sep=" " ,skip=0, header=FALSE)
thisexcl<- c(4:14)
bauexcl= thisexcl

head(deptused)

muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))
head(muniused)
dim(deptused)

row.num=1
deptused<- matrix(as.vector(unlist(deptused)),ncol=2)

source("./00_numfunctions.R")

######### for each round find the model with the 
######## lowest marginal_likelihood star
#######  good for same number of independent variables 

######### generate sumamry indices
summaryexcls<- matrix(NA,nrow=12,ncol=12)
initrnd<- c(20,seq(40,100,10))
thisbest<- c(10)
models_1<- NULL
rnd=1
allbest<- 1:3

for (rnd in 1:10){
  vardomain_init<- (1:15)[-(allbest)]
  allmlx<- read.csv(paste0("../output/round",rnd,"_mll.csv"),  stringsAsFactors = F)[,1]
  thisbest<- vardomain_init[which.max(allmlx)]
  allbest<- c(allbest,thisbest)
  summaryexcls[rnd, 1:(15-rnd-3)]<- (1:15)[which(!(1:15)%in%allbest)]
  models_1<- c(models_1,as.vector(verinits[rnd]+which.max(allmlx)))
}

######################## Generate DIC
alldic<- list()
for (kk in 1:10){
  
  model = models_1[kk]
  run= 1
  
  bauexcl= summaryexcls[kk,] [which(!is.na(summaryexcls[kk,]))]
  
  head(deptused)
  #bauexcl= which(bauexcl==99)
  
  muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))
  source("./00_numfunctions.R")
  
  #########################################################################################
  ###### Multi-variate MCMC 
  #########################################################################################
  
  #given.params = c(-3.25,0.02,0.11,-7.6)
  given.params = 
    c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, 0.375,-9.7)[-bauexcl]
  
  # -----------------------------------------
  pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = c(rep(-10, (length(given.params)-1)),-15), upper = rep(10, length(given.params)))
  
  mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
  bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)
  
  #set.seed(123)
  
  ### multiple runs
  # -----------------------------------------
  iter = 300000
  #head(dataused)
  
  ### adaptive mcmc prior optimization
  settings <- list(iterations = iter, adapt = T, nrChains = 1, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = TRUE)
  
  load(paste0('../output/mcmc_v',model,'_run',run,'.RData'))
  #summary(out.1)  
  #dim(out.1)
  #postSamp<- getSample(out.1, start = 1, end = NULL, thin = 1, whichParameters = 1:(16-length(bauexcl)))
  postSamp<- getSample(out.1, start = 16650, end = NULL, thin = 10, whichParameters = 1:(16-length(bauexcl)))
  #if (kk==1) postSamp<- getSample(out.1, start = 150000, end = NULL, thin = 10, whichParameters = 1:(16-length(bauexcl)))
  tail(postSamp)  
  dim(postSamp)  
  thisll<- NULL
  for (jj in 1: nrow(postSamp)){
    thisll<- c(thisll, log.like.func(postSamp[jj,]))
    print(c(kk,jj))
  }
  x<- postSamp
  lik<- thisll
  dicer<- calc.dic(x,lik,log.like.func)
  alldic [[kk]]<- dicer
}

allrw<- NULL
for (kk in 1:10){
  thisrw<- NULL
  for (jj in 1:6) {
    thisrw<- c(thisrw, alldic[[kk]][jj])
  }
  allrw<- rbind(allrw, thisrw )
}

which(diff(unlist(allrw[,1]))>0)[1]

#write.csv(allrw, "../output/dic_summary.csv", quote=F, row.names=F)


######################## Generate Marginal LL

allmargll<- list()

for (kk in 1:10){
  
  model = models_1[kk]
  run= 1
  
  bauexcl= summaryexcls[kk,] [which(!is.na(summaryexcls[kk,]))]
  
  head(deptused)
  #bauexcl= which(bauexcl==99)
  
  muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))
  source("./00_numfunctions.R")
  
  #########################################################################################
  ###### Multi-variate MCMC 
  #########################################################################################
  
  #given.params = c(-3.25,0.02,0.11,-7.6)
  given.params = 
    c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, 0.375, -9.7)[-bauexcl]
  
  # -----------------------------------------
  pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = c(rep(-10, (length(given.params)-1)),-15), upper = rep(10, length(given.params)))
  
  mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
  bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)
  
  #set.seed(123)
  
  ### multiple runs
  # -----------------------------------------
  iter = 300000
  #head(dataused)
  
  ### adaptive mcmc prior optimization
  settings <- list(iterations = iter, adapt = T, nrChains = 1, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = TRUE)
  
  load(paste0('../output/mcmc_v',model,'_run',run,'.RData'))
  #summary(out.1)  
  #dim(out.1)
  #postSamp<- getSample(out.1, start = 1, end = NULL, thin = 1, whichParameters = 1:(16-length(bauexcl)))
  postSamp<- getSample(out.1, start = 16650, end = NULL, thin = 10, whichParameters = 1:(16-length(bauexcl)))
  #if (kk==1) postSamp<- getSample(out.1, start = 150000, end = NULL, thin = 10, whichParameters = 1:(16-length(bauexcl)))
  tail(postSamp)  
  dim(postSamp)  
  thisll<- NULL
  for (jj in 1: nrow(postSamp)){
    thisll<- c(thisll, log.like.func(postSamp[jj,]))
    print(c(kk,jj))
  }
  x<- postSamp
  lik<- thisll
  marger<- calc.marginal.ll(x,lik,log.like.func, prior.like.chy, ,thisll,)
  allmargll [[kk]]<- marger
}

allmrw<- NULL
for (kk in 1:10){
  thisrw<- NULL
  for (jj in 1:4) {
    thisrw<- c(thisrw, allmargll[[kk]][jj])
  }
  allmrw<- rbind(allmrw, thisrw )
}

#write.csv(allmrw, "../output/margll_summary.csv", quote=F, row.names=F)

