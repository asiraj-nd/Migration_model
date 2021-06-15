##############################################################################################
#### ##### ### # ##### ###          ### ## ######     ######## ######## ### ######### # # # #
### ## ########      ######### ######## ### ######### # # # ##### ##### ### # ##### ###    #### ###  

######     ######        generate centered covariates      ############
######      ######          fed into C++ snippets            ############ 

######  ##  ###            # Amir Siraj, Alex Perkins
######        ##        # asiraj@nd.edu, taperkins@nd.edu
######                       #   May 2020
##### # # # ##### ##### ### # ##### ###    #### ###  ######### ######## ### ########
#############################################################################################

setwd("~/Dropbox/tempx/COL_mig_model/code")
#setwd("C:/Users/asiraj/Dropbox/tempx/COL_mig_model/code")

library(inline)
library(BayesianTools)


source("./00_numfunctions.R")

deptdata = read.csv("../model_inputs/COL_dept_datamaster.csv", stringsAsFactors = F)
head(deptdata)
dim(deptdata)/33

if(!require(mgcv)){install.packages('mgcv');library(mgcv)}
if(!require(mvtnorm)){install.packages('mvtnorm');library(mvtnorm)}
if(!require(sm)){install.packages('sm');library(sm)}


##### variable name labels
ttl<- c('DISTIJ','POPI','POPJ','CONTIJ','UrbanPropI','UrbanPropJ','PercI','PercJ','TinyI','TinyJ',
        'MajCenI','MajCenJ','GECONI','GECONJ', 'Intercept')

migdata<- read.csv("../model_inputs/dept_data_long.csv", stringsAsFactors = F)


##########################################################
####### DEVELOP departent level centered covariate

  ddx<- matricx.cust(migdata[,2],33)
  ddx[which(is.na(ddx))]<- rowSums(ddx, na.rm=T)/32
  head(ddx[,1:10])
  ddy<- matricx.cust(migdata[,1],33)
  dim(migdata)
  head(ddy[,1:10])
  ddy[which(is.na(ddy))]<- ddx[which(is.na(ddy))] - rowSums(ddy, na.rm=T)

  head(cbind(as.vector(unlist(ddy)), as.vector(unlist(ddx))))

  #write.table(round(migdata[,c(15,5)],0),"COL_zone_mig.dat",sep= " ", row.names=F, quote=F, col.names=F)
  #write.table(round(migdata,0),"COL_zone_mig.dat",sep= " ", row.names=F, quote=F, col.names=F)

  ### save non-centered  
  write.table(round(cbind(as.vector(unlist(t(ddy))), as.vector(unlist(t(ddx)))),0),"COL_zone_mig.dat",sep= " ", row.names=F, quote=F, col.names=F)

  ###### center continuous and dummy variables separately 

  varsel = c('DEPTSNI','DEPTSNJ',ttl)[-(2+length(ttl))]
  deptused= data.frame(deptdata[,varsel],Intercept=1)
  head(deptused)

  deptused[,-(1:5)] = round(deptused[,-(1:5)],5)  # round to the nearest e-05
  deptused[,15:16] = round(deptused[,15:16],2)  # GECon
  head(deptused)
  
  bbn<- c(6,11:14)  # dummy variables

  ## center covariates
  deptused[,-c(1,2,bbn,ncol(deptused))] = apply(deptused[,-c(1,2,bbn,ncol(deptused))],2,scale)/2
  deptused[,bbn] = apply(deptused[,bbn],2, function(tr) {tr - optim(0,fn= function(kink){abs(mean(tr-kink))}, method="Brent", lower=0, upper=1)$par})
  deptusedx = deptused[,c(3:(dim(deptused)[2]),1,2)]
  deptusedx[,(dim(deptused)[2]-1):(dim(deptused)[2])] = deptusedx[,(dim(deptused)[2]-1):(dim(deptused)[2])]-1  #c++ index
  deptusedx = deptusedx[,-(dim(deptused)[2]-2)] # drop intercept as data
  names(deptusedx)
  head(deptusedx)
  head(deptused)
  write.table(deptusedx,paste("./COL_zone_feed.dat",sep=""), sep= " ", row.names=F, quote=F, col.names=F) 

  ### save parameter for de-centering variable values
  bbn2<- bbn-2
  allpars<-NULL
  for (i in 1:(ncol(deptusedx)-2)){
    clx<-which(names(deptdata)==names(deptusedx)[i] ) 
    bfr_mean<- mean(as.vector(unlist(deptdata[,clx])))
    bfr_sd<- sd(as.vector(unlist(deptdata[,clx])))
    aft_cnt<- deptusedx[1,i]
    bfr_cnt<- deptdata[1,clx]
    kink=0
    if (i %in% bbn2) {
      kink<-optim(0,fn= function(kink){abs(mean(deptdata[,clx]-kink))}, method="Brent", lower=0, upper=1)$par
      print(round(bfr_cnt,1) == round(aft_cnt+kink,1))
      allpars<- rbind(allpars,c(0,kink))
    } else{  
      print(round(bfr_cnt,1) == round((aft_cnt*2*bfr_sd+bfr_mean)[1],1))
      allpars<- rbind(allpars,c(bfr_mean,bfr_sd))
    }
  ((200-bfr_mean)/bfr_sd)/2
  }

write.table(allpars, paste("./centdept.dat",sep=""), sep= " ", row.names=F, quote=F, col.names=F) 


##########################################################
####### DEVELOP municipal level centered covariate

munidata<- read.csv("../model_inputs/COL_dist_1122_datamaster.csv", stringsAsFactors = F)

length(munidata[,8])

ttl<- c('DISTij','POPi','POPj','CONTij','UrbanPropi','UrbanPropj','Perci','Percj','Tinyi','Tinyj',
        'MajCeni','MajCenj','GECONi','GECONj','Intercept')

head(munidata)

#### rearange variable order
  varsel = c('DEPTSNi','DEPTSNj',ttl)[-(2+length(ttl))]
  varsel<- c(varsel, "CountDisti","CountDistj")
  muniused= data.frame(munidata[,varsel],Intercept=1)

  head(muniused)
  muniused[,-(1:5)] = round(muniused[,-(1:5)],5)  # round to the nearest e-05
  muniused[,15:16] = round(muniused[,15:16],2)  # GECon
  head(muniused)
  dim(muniused)

  bbn<- c(6,11:14) # dummy variables

  ## center covariates
  muniused[,-c(1,2,bbn,(ncol(muniused)-2): ncol(muniused))] = apply(muniused[,-c(1,2,bbn,(ncol(muniused)-2): ncol(muniused))],2,scale)/2
  muniused[,bbn] = apply(muniused[,bbn],2, function(tr) {tr - optim(0,fn= function(kink){abs(mean(tr-kink))}, method="Brent", lower=0, upper=1)$par})
  muniusedx = muniused[,c(3:(dim(muniused)[2]),1,2)]
  muniusedx[,(dim(muniused)[2]-1):(dim(muniused)[2])] = muniusedx[,(dim(muniused)[2]-1):(dim(muniused)[2])]-1  #c++ index
  muniusedx = muniusedx[,-(dim(muniused)[2]-2)] # drop intercept as data
  names(muniusedx)
  head(muniusedx)
  max(muniusedx$DEPTSNi)
  write.table(muniusedx,paste("../code/COL_dist_feed.dat",sep=""), sep= " ", row.names=F, quote=F, col.names=F) 
  dim(muniusedx)/752

  ### save de-centering parameters for each variable
  allpars<-NULL
  bbn2<- bbn+2
  for (i in 1:(ncol(muniusedx)-4)){
    clx<-which( names(munidata)== names(muniusedx)[i] ) 
    bfr_mean<- mean(as.vector(unlist(munidata[,clx])))
    bfr_sd<- sd(as.vector(unlist(munidata[,clx])))
    aft_cnt<- muniusedx[1,i]
    bfr_cnt<- munidata[1,clx]
    kink=0
    if (i %in% bbn2) {
      kink<-optim(0,fn= function(kink){abs(mean(munidata[,clx]-kink))}, method="Brent", lower=0, upper=1)$par
      print(round(bfr_cnt,1) == round(aft_cnt+kink,1))
      allpars<- rbind(allpars,c(0,kink))
    } else{  
      print(round(bfr_cnt,1) == round((aft_cnt*2*bfr_sd+bfr_mean)[1],1))
      allpars<- rbind(allpars,c(bfr_mean,bfr_sd))
    }
    ((200-bfr_mean)/bfr_sd)/2
  }

  write.table(allpars, paste("./centmuni.dat",sep=""), sep= " ", row.names=F, quote=F, col.names=F) 

######### Generate Log-likelihood K values 

  deptdata  = read.table("./COL_zone_mig.dat", sep=" " ,skip=0, header=FALSE)
  n=deptdata[,2]
  x=deptdata[,1]

  findMK(n,x,33)

