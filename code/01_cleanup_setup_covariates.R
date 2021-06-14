##############################################################################################
#### ##### ### # ##### ###          ### ## ######     ######## ######## ### ######### # # # #
### ## ########      ######### ######## ### ######### # # # ##### ##### ### # ##### ###    #### ###  

######     ######   setup  large (pair-wise) covarate matrix     ############

######  ##  ###            # Amir Siraj, Alex Perkins
######        ##        # asiraj@nd.edu, taperkins@nd.edu
######                       #   May 2020
##### # # # ##### ##### ### # ##### ###    #### ###  ######### ######## ### ########
#############################################################################################

#setwd("~/Dropbox/tempx/COL_mig_model/code")
setwd("C:/Users/asiraj/Dropbox/tempx/COL_mig_model/code")

source("./00_numfunctions.R")

deptcollect<- read.csv("../data/COL_dept_collector.csv", stringsAsFactors = F)    
dist1122collect<- read.csv("../data/COL_dist_1112_collector.csv", stringsAsFactors = F)    
lcx<- read.csv("../data/COL_lookup_all.csv", stringsAsFactors = F) 

deptmigmx<- read.csv("../data/COL_dept_migmatrix.csv", stringsAsFactors = F)
deptdistmx<- read.csv("../data/COL_dept_distmatrix.csv", stringsAsFactors = F)
deptcontmx<- read.csv("../data/dept_contig_matrix.csv", stringsAsFactors = F)

df = deptmigmx
allposts=NULL
IJserial = matrix(1:(nrow(df)*ncol(df)), nrow=nrow(df)[1], byrow=T)
zpop<- deptcollect$pop_2013
additional_miss
for (thispair in 1:(nrow(df) * ncol(df))) {
  thisrw = which (thispair==IJserial, arr.ind=T)[1]
  thiscol = which (thispair==IJserial, arr.ind=T)[2]
  allposts=rbind(allposts,c(deptcollect$id[thisrw],deptcollect$id[thiscol],
                            which (thispair==IJserial, arr.ind=T), zpop[thisrw],zpop[thiscol],
                            as.vector(apply(df,1,sum))[thisrw], urbprop[thisrw],
                            urbprop[thiscol], deptcollect$gecon[thisrw],
                            deptcollect$gecon[thiscol], deptcollect$coordx[thisrw], 
                            deptcollect$coordy[thisrw],deptcollect$coordx[thiscol], 
                            deptcollect$coordy[thiscol]))
}

head(allposts)
dim(allposts)

deptdata = data.frame(allposts,as.vector(unlist(deptdistmx)),as.vector(unlist(deptcontmx)))
names(deptdata) = c("CODEI","CODEJ","DEPTSNI","DEPTSNJ","POPI","POPJ","IPUMPSPOP","UrbanPropI","UrbanPropJ", "GECONI","GECONJ", "LONI","LATI","LONJ","LATJ", "DISTIJ", "CONTIJ")

#### add pop-based additional covariates
  deptdata$MIGIJ = as.vector(unlist(t(df)))
  deptdata= deptdata[which(deptdata[,1]!=deptdata[,2]),] ## remove diagonal
  deptdata$MajCenI=with(deptdata, ifelse(POPI>quantile(POPI,.9), 1,0))
  deptdata$MajCenJ=with(deptdata, ifelse(POPJ>quantile(POPJ,.9), 1,0))
  deptdata$TinyI=with(deptdata, ifelse(POPI>quantile(POPI,.10), 1,0))
  deptdata$TinyJ=with(deptdata, ifelse(POPJ>quantile(POPJ,.10), 1,0))
  deptdata$PercI=with(deptdata, rank(deptdata$POPI)/length(POPI)) ## ask Ale why they used order instead of rank
  deptdata$PercJ=with(deptdata, rank(deptdata$POPJ)/length(POPJ))
  deptdata$TotPop=with(deptdata, sum(POPI))
  head(deptdata)

  length(rowSums(df))
  n=matrix(rep(rowSums(df),78),nrow=78)
  y=df
  #write.csv(cbind(mx_to_vec(y), mx_to_vec(n)), "../model_inputs/dept_data_long.csv", row.names=F, quote=F)

##############################################################  
##### 1122 districts level #creat dist-dist list

dist1122collect<- read.csv("../data/COL_dist_1122_collector.csv", stringsAsFactors = F)
distmx<- read.csv("../data/COL_dist_1122_distmatrix.csv", stringsAsFactors = F)
contgmx<- read.csv("../data/dist_1122_contig_matrix.csv", stringsAsFactors = F)



########## setup 1122 x 1122 (excluding diagonal) (long dataset)

nsel = 1122

cl_1=NULL  # DCODE_i  (repeating)
cl_2 = rep(dist1122collect$deptID,nsel) # DCODE_j (alternating) 

cl_3=NULL # MCODE_i (repeating)
cl_4 =  rep(dist1122collect[,1],nsel) # MCODE_j (alternating)

cl_5=NULL   # DeptSN_j  (repeating)
cl_6 = rep( sapply(dist1122collect$deptID,function(tr) {which(deptcollect$id==tr)}),nsel) # DeptSN_j  (alternating)

cl_7=NULL # MunSN_j (repeating)
cl_8= rep( sapply(dist1122collect[,1],function(tr) {which(dist1122collect[,1]==tr)}),nsel) # MunSN_j (alternating)

cl_9=NULL # POPi (repeating)
cl_10 = rep(dist1122collect$pop_2010,nsel) # POPj (alternating)

cl_11=NULL # AREAi (repeating)  # not used
cl_12=NULL  # AREAj (alternating)  # not used

cl_13=NULL # URBANPROPi (repeating)
cl_14 = rep(dist1122collect$urbprop,nsel) # URBANPROPj (alternating)

cl_15=NULL # GECONi (repeating)
cl_16 = rep(dist1122collect$gecon,nsel) # GECONj (alternating)

cl_17=NULL # LAT_i (repeating)
cl_18 = rep(dist1122collect$lat,nsel) #LAT_j (alternating)

cl_19=NULL #LON_i (repeating)
cl_20 = rep(dist1122collect$lon,nsel) #LON_j (alternating)

cl_25 = NULL # COUNTDist_i
cl_26 = rep(dist1122collect$CountDistinDept,nsel) #COUNTDIST_j (alternating)

clsel = (0:(nsel-1))*nsel

for (kk in 1:nsel) {
  cl_1 = c(cl_1,rep(dist1122collect$deptID,nsel) [clsel+kk])
  cl_3 = c(cl_3,rep(dist1122collect[,1],nsel) [clsel+kk])
  cl_5=  c(cl_5,rep(sapply(dist1122collect$deptID,function(tr) {which(deptcollect$id==tr)}),nsel) [clsel+kk])
  cl_7=  c(cl_7,rep (sapply(dist1122collect[,1],function(tr) {which(dist1122collect[,1]==tr)}),nsel) [clsel+kk])
  cl_9 =  c(cl_9,rep(dist1122collect$pop_2010,nsel) [clsel+kk])
  #cl_11 =  c(cl_11,rep(munilist[,5],nsel) [clsel+kk])
  cl_13 =  c(cl_13,rep(dist1122collect$urbprop,nsel) [clsel+kk])
  cl_15 =  c(cl_15,rep(dist1122collect$gecon,nsel) [clsel+kk])
  cl_17 =  c(cl_17,rep(dist1122collect$lat,nsel) [clsel+kk])
  cl_19 =  c(cl_19,rep(dist1122collect$lon,nsel) [clsel+kk])
  cl_25 =  c(cl_25,rep(dist1122collect$CountDistinDept,nsel) [clsel+kk])
  print(kk)
}

munidata = data.frame(cbind(cl_1,cl_2,cl_3,cl_4,cl_5,cl_6,cl_7,cl_8,as.vector(unlist(distmx[,-msx])),as.vector(unlist(contgmx[,-msx])),cl_9,cl_10,NA,NA,cl_13,cl_14,cl_15,cl_16,cl_17,cl_18,cl_19,cl_20,cl_25,cl_26))

names(munidata) = c("DCODEi","DCODEj","MCODEi","MCODEj","DEPTSNi","DEPTSNj","MUNSNi","MUNSNj","DISTij","CONTij","POPi","POPj","AREAi","AREAj","UrbanPropi","UrbanPropj","GECONi","GECONj","LATFR","LATTO","LONFR","LONTO", "CountDisti","CountDistj")

##### delete diagonal
excl = which(unlist(munidata$MUNSNi)==unlist(munidata$MUNSNj))

munidata= munidata[-excl,] ## remove diagonal
munidata$MajCeni= ifelse(unlist(munidata$POPi)>quantile(unlist(munidata$POPi),.9), 1,0)
munidata$MajCenj= ifelse( unlist(munidata$POPj)>quantile(unlist(munidata$POPj),.9), 1,0)
munidata$Tinyi=   ifelse(unlist(munidata$POPi)>quantile(unlist(munidata$POPi),.10), 1,0)
munidata$Tinyj=   ifelse(unlist(munidata$POPj)>quantile(unlist(munidata$POPj),.10), 1,0)
munidata$Perci=  rank(unlist(munidata$POPi))/ length(munidata$POPi)
munidata$Percj=   rank(unlist(munidata$POPj))/ length(munidata$POPj)

names(munidata) = c("DCODEi","DCODEj","MCODEi","MCODEj","DEPTSNi","DEPTSNj","MUNSNi","MUNSNj","DISTij","CONTij","POPi","POPj","AREAi","AREAj","UrbanPropi","UrbanPropj","GECONi","GECONj","LATFR","LATTO","LONFR","LONTO","CountDisti", "CountDistj",
                    "MajCeni","MajCenj","Tinyi","Tinyj","Perci","Percj")

#fwrite(munidata, file="../model_inputs/COL_dist_752_datamaster.csv")


