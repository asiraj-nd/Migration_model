

################## 

####### numeric functions


loggamma <- cxxfunction(signature(xin = "numeric"),
                        plugin="Rcpp",
                        body=
                          '
                        Rcpp::NumericVector tx(xin);
                        Rcpp::NumericVector xout(1);
                        double x = sum(tx);
                        double result, y, xnum, xden;
                        int i;
                        const double d1 = -5.772156649015328605195174e-1;
                        const double p1[] = { 
                        4.945235359296727046734888e0, 2.018112620856775083915565e2, 
                        2.290838373831346393026739e3, 1.131967205903380828685045e4, 
                        2.855724635671635335736389e4, 3.848496228443793359990269e4, 
                        2.637748787624195437963534e4, 7.225813979700288197698961e3 
                        };
                        const double q1[] = {
                        6.748212550303777196073036e1, 1.113332393857199323513008e3, 
                        7.738757056935398733233834e3, 2.763987074403340708898585e4, 
                        5.499310206226157329794414e4, 6.161122180066002127833352e4, 
                        3.635127591501940507276287e4, 8.785536302431013170870835e3
                        };
                        const double d2 = 4.227843350984671393993777e-1;
                        const double p2[] = {
                        4.974607845568932035012064e0, 5.424138599891070494101986e2, 
                        1.550693864978364947665077e4, 1.847932904445632425417223e5, 
                        1.088204769468828767498470e6, 3.338152967987029735917223e6, 
                        5.106661678927352456275255e6, 3.074109054850539556250927e6
                        };
                        const double q2[] = {
                        1.830328399370592604055942e2, 7.765049321445005871323047e3, 
                        1.331903827966074194402448e5, 1.136705821321969608938755e6, 
                        5.267964117437946917577538e6, 1.346701454311101692290052e7, 
                        1.782736530353274213975932e7, 9.533095591844353613395747e6
                        };
                        const double d4 = 1.791759469228055000094023e0;
                        const double p4[] = {
                        1.474502166059939948905062e4, 2.426813369486704502836312e6, 
                        1.214755574045093227939592e8, 2.663432449630976949898078e9, 
                        2.940378956634553899906876e10, 1.702665737765398868392998e11, 
                        4.926125793377430887588120e11, 5.606251856223951465078242e11
                        };
                        const double q4[] = {
                        2.690530175870899333379843e3, 6.393885654300092398984238e5, 
                        4.135599930241388052042842e7, 1.120872109616147941376570e9, 
                        1.488613728678813811542398e10, 1.016803586272438228077304e11, 
                        3.417476345507377132798597e11, 4.463158187419713286462081e11
                        };
                        const double c[] = {
                        -1.910444077728e-03, 8.4171387781295e-04, 
                        -5.952379913043012e-04, 7.93650793500350248e-04, 
                        -2.777777777777681622553e-03, 8.333333333333333331554247e-02, 
                        5.7083835261e-03
                        };
                        const double a = 0.6796875;
                        
                        if((x <= 0.5) || ((x > a) && (x <= 1.5))) {
                        if(x <= 0.5) {
                        result = -log(x);
                        /*  Test whether X < machine epsilon. */
                        if(x+1 == 1) {
                        xout[0] = result;
                        return (xout);
                        }
                        }
                        else {
                        result = 0;
                        x = (x - 0.5) - 0.5;
                        }
                        xnum = 0;
                        xden = 1;
                        for(i=0;i<8;i++) {
                        xnum = xnum * x + p1[i];
                        xden = xden * x + q1[i];
                        }
                        result += x*(d1 + x*(xnum/xden));
                        }
                        else if((x <= a) || ((x > 1.5) && (x <= 4))) {
                        if(x <= a) {
                        result = -log(x);
                        x = (x - 0.5) - 0.5;
                        }
                        else {
                        result = 0;
                        x -= 2;
                        }
                        xnum = 0;
                        xden = 1;
                        for(i=0;i<8;i++) {
                        xnum = xnum * x + p2[i];
                        xden = xden * x + q2[i];
                        }
                        result += x*(d2 + x*(xnum/xden));
                        }
                        else if(x <= 12) {
                        x -= 4;
                        xnum = 0;
                        xden = -1;
                        for(i=0;i<8;i++) {
                        xnum = xnum * x + p4[i];
                        xden = xden * x + q4[i];
                        }
                        result = d4 + x*(xnum/xden);
                        }
                        /*  X > 12  */
                        else {
                        y = log(x);
                        result = x*(y - 1) - y*0.5 + .9189385332046727417803297;
                        x = 1/x;
                        y = x*x;
                        xnum = c[6];
                        for(i=0;i<6;i++) {
                        xnum = xnum * y + c[i];
                        }
                        xnum *= x;
                        result += xnum;
                        }
                        //System.out.println(x);
                        xout[0] = result;
                        return (xout);
                        '
)


####################################################
generate.Header<-function(nx,ny,minlon,maxlat,txt,cellWidth) {
  A1=c("ENVI description", "=", "{File Resize Result, x resize factor: 1.000000, y resize factor: 1.000000. [Sun Nov 05 11:08:48 2012]}")
  A2=c("samples", "=",nx)
  A3=c("lines", "=",ny)
  A4=c("bands",  "=",1)
  A5=c("header offset",  "=",0)
  A6=c("file type", "=","ENVI Standard")
  A7=c("data type", "=",4)
  A8=c("interleave", "=","bsq")
  A9=c("sensor type", "=","Unknown")
  A10=c("byte order ", "=",0)
  A11=c("x start", "=",1)
  A12=c("y start", "=",1)
  A13=c("map info", "=",paste("{Geographic Lat/Lon, 1.0000, 1.0000,", minlon, ",",maxlat,",", cellWidth,",", cellWidth,", WGS-84, units=Meters}",sep=""))
  A14=c("wavelength units", "=","Unknown")
  A15=c("data ignore value", "=",-9.99900000e+003)
  A16=c("band names", "=","{Resize (Band 1:tmin_1.bil)}")
  write.table(rbind(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16), txt, sep=" ",col.names=F,row.names=F, quote=F)
  return(0)
}

## raster definition function
world_2_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-8640
  ny<-4320
  xllCorner<- -180
  maxlat<-90
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.04166667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

colo2_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-432
  ny<-480
  xllCorner<- -83
  maxlat<-15.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.04166667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

#### raster definition fuctions
eth2_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-365
  ny<-279
  xllCorner<- 32.875
  maxlat<- 14.95833333
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.04166667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

raster_collapse<- function(gmx,rt,aggrFun) {
  gx<-dim(gmx)/rt
  shrinkmx<-matrix(NA,gx[1],gx[2])
  for (i in 0:(gx[1]-1))
    for (j in 0:(gx[2]-1)) { 
      if (aggrFun==1 & any(!is.na(gmx[(1:rt)+i*rt,(1:rt)+j*rt]))) shrinkmx[i+1, j+1]<- mean(gmx[(1:rt)+i*rt,(1:rt)+j*rt],na.rm=T)
      if (aggrFun==2 & any(!is.na(gmx[(1:rt)+i*rt,(1:rt)+j*rt]))) shrinkmx[i+1, j+1]<- sum (gmx[(1:rt)+i*rt,(1:rt)+j*rt],na.rm=T)
    }
  return(shrinkmx)
}


## linear interpolate
interpolateCols <- function(vec, nSubsteps)
{
  matInterp <- rep(0,length(vec)+(length(vec)-1)*nSubsteps)
  
  for(i in 1:(length(vec)-1)) {
    interpVals <- approx(vec[i:(i+1)], n=(nSubsteps+2))
    lilloc<-i+nSubsteps*(i-1)
    matInterp[lilloc:(lilloc+nSubsteps+1)] <- interpVals$y
  }
  return(matInterp)
}

mx_to_vec<- function(mx,rw=1) {
  if (rw==1) {
    mmx<- as.vector(unlist(t(mx)))
  }
  
  if (rw==2) {
    mmx<- as.vector(unlist(mx))
  }

  dx<- diag( matrix(1: length(mmx), nrow= nrow(mx), byrow=TRUE)) 
#  print(  dx)
  #  print(  seq(1,length(mmx),nrow(mx)))
  
  return (mmx[-(dx)])
}
    

vec_to_mx_byrow<- function(sim,dx){
  
  allmx<- matrix(1,dx,dx)
  for (kk in 1:dx) {
    allmx[kk,kk]<- NA    
  }
  
  allmx[-which(is.na(allmx))]<- sim
  
  return(t(allmx))
}



eth_trim<- function(worldmx){
  
  return(worldmx[-(1:1801),-(1:5109)][1:279,1:365])
}

pivoter<- function (x,y, vals) {
  mx<- matrix(NA, nrow= length(unique(y)), ncol= length(unique(x)))
  for (i in 1: length(unique(y))) {
    for (j in 1: length(unique(x))){
      vi<- unique(y)[i]
      vj<- unique(x)[j]
      aggij<- which(x==vj & y==vi)
      mx[i,j] <- mean(vals[aggij])
    }
  }
  
  mx<- data.frame(mx)
  names(mx)<- unique(x)
  row.names(mx)<- unique(y)
  return(mx)  
}


prior.like.chy<- function (params) {
  aprior = dcauchy(params[-length(params)], scale=2.5, location=0, log=T)
  bprior = dcauchy(params[length(params)], scale=10, location=0, log=T)
  return(sum(c(aprior,bprior)))
}

sampler.chy<- function (){
  s = rcauchy(length(given.params)-1, scale=2.5, location=0)
  s= c(s, rcauchy(1, scale=10, location=0))
}


##  C++ version of multinomial likelihood function
logllMx <- cxxfunction(signature(xin = "numeric", mxin = "numeric", ddin="numeric"),
                      plugin="Rcpp",
                      body=
                        '
                      int dptsrc;
                      int dptdest;
                      double logll = 0;
                      Rcpp::NumericVector params(xin);
                      Rcpp::NumericMatrix mrec(mxin);
                      Rcpp::NumericMatrix dept(ddin);
                      Rcpp::NumericVector migij(1257762);                                              
                      Rcpp::NumericMatrix aggMig(33,33);
                      Rcpp::NumericMatrix srcPop(33,33);
                      Rcpp::NumericVector probMig(1056);
                      Rcpp::NumericVector lglv(1089);
                      Rcpp::NumericVector sump(33);
                      int nrows = mrec.nrow();
                      int ncols = mrec.ncol();
                      
                      Rcpp::NumericVector xout(5);
                      for (int j= 0; j<nrows; j++) {
                      
                        dptsrc = mrec[nrows*(ncols-2)+j];
                        dptdest = mrec[nrows*(ncols-1)+j];
                        migij[j] = params[(ncols-4)];
                        for (int k= 0; k<(ncols-4); k++) {
                            migij[j] += params[k] * mrec[nrows*k + j];
                        }
                            migij[j] = exp(migij[j]) / (1+ exp(migij[j])) * (2 * 2.535538e5 * mrec[nrows*1+j] + 4.144954e4);
                           aggMig[33 * dptsrc + dptdest] += (double) migij[j];
                           srcPop[33 * dptsrc + dptdest ] += (double)(2 * 2.535538e5 * mrec[nrows*1+j] + 4.144954e4)/ mrec[nrows*(ncols-3)+j];
                      }  // j loop

                      int idr = 0;
                      for (int ll=0; ll<33; ll++){
                        sump[ll] = 0;
                        for (int kk=0; kk<33; kk++){
                          if (ll!=kk) {   // l then kk (because orig matrix and this are similar = by row)
                             probMig[idr] = (double) (aggMig[33 * ll + kk]/ srcPop[33 * ll + kk]);  //comment - DEFUNCT weighted by number of municipalities in dest dept
                             sump[ll] +=  probMig[idr];
                            idr++;
                          }	
                        }	
                      }
                      
                      // cap at cummulative prob of 0.51, at or above that- deflate LL so likely gets rejected
                      for (int i=0; i<33; i++){
                        if (sump[i]>0.5) sump[i]=0.51;
                      }  
                      xout[1] = 0;  
                      int odr = 0;
                      idr = 0;
                      logll = 0;
                     for (int ss=0; ss<33; ss++){
                          for (int dd=0; dd<33; dd++) {
                            if (ss==dd) {
                                 //             y   *  log (p)  + K
                                  lglv[odr] = dept[1089* 0 + odr] * log(1- sump[ss]);  // missing local mobility on purpose
                                  logll  += lglv[odr];
                                  if (sump[ss]==0.51) logll += -3.0e5;   //deflate LL
                            }
                            if (ss!=dd) {
                                lglv[odr] = dept[1089* 0 + odr] * log(probMig[idr]);
                                logll  += lglv[odr];
                                idr++;
                              } 
                            odr++;
                          }
                        }
                        logll += 790786.3; // 790k = sum(logFactorial(n[i]) -  logFactorial(y[i]) - logFactorial(n[i]-y[i]))
                        
                        xout[0]= logll;
                        return(xout);
                      '
)

####### Likelihood function caller
log.like.func<- function(given.params) {
  logllMx(given.params, muniused,deptused)[1]
}


############# simulate aggregate migration estimates at department level
simulxagg <- cxxfunction(signature(xin = "numeric", mxin = "numeric"),
                       plugin="Rcpp",
                       body=
                         '
                      int dptsrc;
                      int dptdest;
                      Rcpp::NumericVector params(xin);
                      Rcpp::NumericMatrix mrec(mxin);
                      Rcpp::NumericVector migij(1257762);                                              
                      Rcpp::NumericMatrix aggMig(33,33);
                      Rcpp::NumericMatrix srcPop(33,33);
                      Rcpp::NumericVector probMig(1056);
                      Rcpp::NumericVector lglv(1089);
                      Rcpp::NumericVector sump(33);
                      int nrows = mrec.nrow();
                      int ncols = mrec.ncol();
                      
                      Rcpp::NumericVector xout(5);
                      for (int j= 0; j<nrows; j++) {
                      
                        dptsrc = mrec[nrows*(ncols-2)+j];
                        dptdest = mrec[nrows*(ncols-1)+j];
                        migij[j] = params[(ncols-4)];
                        for (int k= 0; k<(ncols-4); k++) {
                            migij[j] += params[k] * mrec[nrows*k + j];
                        }
                            migij[j] = exp(migij[j]) / (1+ exp(migij[j])) * (2 * 2.535538e5 * mrec[nrows*1+j] + 4.144954e4);
                           aggMig[33 * dptsrc + dptdest] += (double) migij[j];
                           srcPop[33 * dptsrc + dptdest ] += (double)(2 * 2.535538e5 * mrec[nrows*1+j] + 4.144954e4)/ mrec[nrows*(ncols-3)+j];
                      }  // j loop

                      int idr = 0;
                      for (int ll=0; ll<33; ll++){
                        sump[ll] = 0;
                        for (int kk=0; kk<33; kk++){
                          if (ll!=kk) {   // l then kk (because orig matrix and this are similar = by row)
                             probMig[idr] = (double) (aggMig[33 * ll + kk]/ srcPop[33 * ll + kk]);  //comment - DEFUNCT weighted by number of municipalities in dest dept
                            idr++;
                          }	
                        }	
                      }
                      return(probMig);
                      '
)

######## aggregate simulatation function caller
simulagg.func<- function(given.params) {
  simulxagg(given.params, muniused)
}

############# simulate migration estimates at district level

simulx <- cxxfunction(signature(xin = "numeric", mxin = "numeric"),
                           plugin="Rcpp",
                           body=
                             '
                      Rcpp::NumericVector params(xin);
                      Rcpp::NumericMatrix mrec(mxin);
                      Rcpp::NumericVector migij(1257762);                                              
                      int nrows = mrec.nrow();
                      int ncols = mrec.ncol();
                      for (int j= 0; j<nrows; j++) {
                        migij[j] = params[(ncols-4)];
                        for (int k= 0; k<(ncols-4); k++) {
                            migij[j] += params[k] * mrec[nrows*k + j];
                        }
                        migij[j] = exp(migij[j]) / (1+ exp(migij[j])) * (2 * 2.535538e5 * mrec[nrows*1+j] + 112243.0);
                      }
                      return(migij);
                    '
)


############# simulate migration function caller
simul.func<- function(given.params) {
  simulx(given.params, muniused)
}


####### factorial function caller
logFactorial<- function(x){
  return(loggamma(x+1))
}

###### generate multinomial function likelihood contant term
findMK<- function (n,x,t=33){
  allk<- 0
  t2<- length(n)/t
  for (i in 1:t) {
    dpti = (i-1)*t2
    moved = 0
    thisk <- logFactorial(n[dpti+1])  
    for (j in 1: t2) {
      thisk<- thisk - logFactorial(x[dpti+j])
      #   print(c(j,n[dpti+1], dpti+j))
      moved <- moved + x[dpti+j]
    }    
    thisk <- thisk - logFactorial (n[dpti+1] - moved)
    # print( c(i,n[dpti+1] - moved))
    allk<- allk + thisk
  }
  return(allk)  
}

calc.dic <- function(x,lik,lik.fun,...) {
  D.bar <- -2*mean(lik)
  theta.bar <- as.vector(apply(x,2,mean))
  D.hat <- -2*lik.fun(theta.bar)
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}

calc.marginal.ll <- function(x,lik,lik.fun,prior.fun,num.samples=1000,lik.j,log=TRUE) {
  
  y <- summary(x)
  theta.star <- as.vector(apply(x,2,mean))
  
  lik.star <- lik.fun(theta.star)
  V = cov(x,x)
  g <- sample(1:nrow(x),num.samples,replace=TRUE)
  theta.g <- x[g,]
  q.g <- dmvnorm(theta.g,mean=theta.star,sigma=V,log=FALSE)
  
  lik.g <- lik[g]
  alpha.g <- sapply(lik.g,function(tr) min(1,exp(lik.star-tr)))
  
  alpha.j <- sapply(lik.j,function(tr) min(1,exp(tr-lik.star)))
  
  pi.hat <- mean(alpha.g*q.g)/mean(alpha.j)
  pi.star <- 0
  if (!is.null(prior.fun)) pi.star <- prior.fun(theta.star)
  ln.m <- lik.star + pi.star - log(pi.hat)
  if (! log) ln.m <- exp(ln.m)
  list(ln.m=ln.m,ln.lik.star=lik.star,ln.pi.star=pi.star,ln.pi.hat=log(pi.hat))
}

mle.norm<- function(barer) {
  return (optim(c(0,1),function(par){
    -sum(barer[,1] / sum(barer[,1]) * log(dnorm(barer[,2], par[1], par[2])))
  })$par)
}



  
  
