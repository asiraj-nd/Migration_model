##############################################################################################
#### ##### ### # ##### ###          ### ## ######     ######## ######## ### ######### # # # #
### ## ########      ######### ######## ### ######### # # # ##### ##### ### # ##### ###    #### ###  

######                          fitting models                  ############
######                        identify the best            ############ 

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

  deptused  = read.table("./COL_zone_mig.dat", sep=" " ,skip=0, header=FALSE)
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
  
  #########################################################################################
  ###### Multi-variate MCMC 
  #########################################################################################
  
  #given.params = c(-3.25,0.02,0.11,-7.6)
  given.params = 
    c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, 0.375,-9.7)[-bauexcl]
  
  prior.like.chy(given.params)
  print(sampler.chy())
  
  thist<- Sys.time()
  log.like.func(given.params)
  print(Sys.time() - thist)
  
  # -----------------------------------------
  pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = rep(-10, length(given.params)), upper = rep(10, length(given.params)))
  
  mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
  bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)
  
  
########################################################################################################################  

  ######### iterative run for 10 rounds each with kk models

########################################################################################################################  
  
  # -----------------------------------------
  iter = 300000
  #head(dataused)
  
  #####iteration_version model name initialer
  verinits<- c(0,15,seq(30,100,10)) 

#################################################################################
allbest<- 1:3   # default variable list (dist, popi, popj)
  
for (rnd in 1:10) {  

    thisroundexcls<-   matrix(NA,nrow=12,ncol=12)
    initexcl<- NULL
    vardomain_init<- (1:15)[-(allbest)]
    
    ##### variable exclusion rows
    for (jj in 1:length(vardomain_init)){
        thisexcl<- vardomain_init[!(vardomain_init %in% vardomain_init[jj])]
        thisroundexcls[jj,1:length(thisexcl)]<- thisexcl
      }
  
    ##### remove blank rows
    if (length(which(apply(thisroundexcls,1,function(tr){all(is.na(tr))})))>0) {
      thisroundexcls<- thisroundexcls[-which(apply(thisroundexcls,1,function(tr){all(is.na(tr))})),]
    }
  
    for (kk in 1:nrow(thisroundexcls)) {
      model = verinits[rnd] + kk
      run= row.num
      bauexcl= thisroundexcls[kk,] [which(!is.na(thisroundexcls[kk,]))]
    
      muniused = matrix( as.vector(unlist(munidata[,-bauexcl])),ncol=(ncols-length(bauexcl)))
      source("./00_numfunctions.R")
    
      #########################################################################################
      ###### Multi-variate MCMC 
      #########################################################################################
    
      given.params = 
        c(-1.504,0.032,0.271,3.122,-0.374,0.986,0.221,0.0, -0.423, 0.0, -0.640, 1.035, 0.080, 0.375, -9.7)[-bauexcl]
    
      # -----------------------------------------
      pll<- createPrior(density = prior.like.chy, sampler = sampler.chy, lower = c(rep(-10, (length(given.params)-1)),-15), upper = rep(10, length(given.params)))
    
      mylikelihood <- createLikelihood(log.like.func, names=NULL, parallel=F, catchDuplicates = T, sampler=NULL, parallelOptions = NULL)
      bayesianSetup <- createBayesianSetup (likelihood=mylikelihood,  prior=pll)
    
      ### multiple runs
      # -----------------------------------------
      iter = 300000

      ### adaptive mcmc prior optimization
      settings <- list(iterations = iter, adapt = T, nrChains = 1, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = T,  message = TRUE)
      out.1 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
      save(out.1,file=paste('../output/mcmc_v',model,'_run',run,'.RData',sep=''))
      print(kk)
    }
  
    ######## Compare marginal likelihood in this round
    allmlx<- rep(NA, nrow(thisroundexcls))
    for (i in 1:nrow(thisroundexcls)) {
      model<- verinits[rnd] + i
      load(paste('../output/mcmc_v',model,'_run',run,'.RData',sep=''))
      allmlx[i]<- marginalLikelihood(out.1)$ln.lik.star
      print(i)
    }

    thisbest<- vardomain_init[which.max(allmlx)]
    allbest<- c(allbest,thisbest)
    write.csv(allmlx, paste0("../output/round",rnd,"_mll.csv"),  row.names=F, quote=F)
  #allmlx<- read.csv(paste0("../output/round",rnd,"_mll.csv"),  stringsAsFactors = F)[,1]
  
} #### loop through rounds 


