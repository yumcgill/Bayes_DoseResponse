library(MASS)
library(nlme)
library(lme4)
library(mnormt)
library(mgcv)
library(gamm4)
library(boot)
library(splines)
library(parallel)
library(np)
library(MCMCpack)
library(geepack)
options(mc.cores = 23)

S<-20
N<-1000
mux<-0.2
muu<-1.0
muz<-0.2
varx<-0.1
varu<-0.6
varz<-0.1
dose1<-3.00
dose2<-4.00
dose3<-5.00

# function to lag longitudinal data
longlag <- function(x) {
  y <- ts(x$y, start = x$time[1]) #Time series objects, y;
  d <- ts(x$d, start = x$time[1]) #Time series objects, d;
  idx <- seq(length = length(y))
  yl <- cbind(y, stats::lag(y, -1))[idx,]
  dl <- cbind(d, stats::lag(d, -1))[idx,]
  cbind(x[,1:2], yl, dl)
}

# function to difference longitudinal data
longdiff <- function(x) {
  y <- ts(x$y, start = x$time[1])
  logy <- ts(log(x$y), start = x$time[1])
  d <- ts(x$d, start = x$time[1])
  idx <- seq(length = length(y))
  yd <- cbind(y, y-stats::lag(y, -1))[idx,] #absolute difference from previous y;
  dlogy <- cbind(log(y), log(y)-stats::lag(log(y), -1))[idx,] #absolute difference on log-scale from previous y;
  dd <- cbind(d, d-stats::lag(d, -1))[idx,] #absolute difference from previous d;
  cbind(x[,1:2], yd, dlogy, dd)
}


  seed <- 123
  set.seed(seed)
  x<-rnorm(1000,mux,sqrt(varx)) 
  v1<-rnorm(1000,muu,sqrt(varu)) 
  v2<-rnorm(100,muz,sqrt(varz)) 
  
  #u<-rep(v1,10)
  u <- v1
  z<-rep(v2,10)
  count<-c(1:100) #subjects;
  obs<-c(1:1000) #observations;
  id<-rep(count,10) #100 subjects each measured 10 times = 1000 observations;
  time<-c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),rep(7,100),rep(8,100),rep(9,100),rep(10,100))
  
  # treatment assignment mechanism;
  bx<-1.0
  bu<-2.0
  bz<-4.0
  vared<-1
  ed<-rnorm(1000,0,sqrt(vared))
  d<-(1.0+bx*x+bu*u+bz*z+ed) 
  #normal dose; place holder for non-symmetric dose dist;
  # x and u are time-varying covariates, z are fixed baseline variable;
  # dose is time-varying not single point treatment assignment!;
  
  ey<-rnorm(1000,0,1)
  y<-(20*exp((d+x-0.25*u + 4*z)))*exp(ey) 
  #continuous y; fairly large Y with exponentiated variance;

  
  yi0<-rep(y[1:100],10) #we are not using the rest 900 simulated Y;
  di0<-rep(d[1:100],10)
  xi0<-rep(x[1:100],10)
  lyi0<-log(yi0)
  eyi0<-mean(y[1:100])
  elyi0<-mean(log(y[1:100]))
  dy<-y-yi0
  dd<-d-di0
  dx<-x-xi0
  dlogy<-log(y)-log(yi0) #logging this because count data?;
  
  A<-data.frame(obs,id,time,y,log(y),d,dy,dlogy,dd,dx,di0)
  sort.A <-A[order(id) , ]
  llA<-do.call("rbind", by(sort.A, sort.A$id, longlag))
  ldA<- do.call("rbind", by(sort.A, sort.A$id, longdiff))
  lid<-as.numeric(llA$"id")
  ldAy<-as.numeric(llA$"stats::lag(y, -1)")
  ldAd<-as.numeric(llA$"stats::lag(d, -1)")
  lobs<-as.numeric(llA$"obs")
  procA<-data.frame(ly=ldAy,ld=ldAd,lobs=lobs,lid=lid,time=sort.A$time)
  lagdat<-procA[order(lobs) , ]
  lagdat$id<-id
  lagdat$y<-y
  lagdat$logy<-log(y)
  lagdat$x<-x
  lagdat$u<-u
  lagdat$z<-z
  lagdat$d<-d
  lagdat$dy<-dy
  lagdat$dd<-dd
  lagdat$dx<-dx
  lagdat$di0<-di0
  lagdat$dlogy<-dlogy #numerically assign dy dx dd with 0 for first visit;
  
  mud1<-log(20)+mean(dose1+x-0.25*u+ 4*z) #setting up dose change on mean of logy;
  mud2<-log(20)+mean(dose2+x-0.25*u+ 4*z) 
  mud3<-log(20)+mean(dose3+x-0.25*u+ 4*z)
  
  # y<-(20*exp((d+x-0.25*u + 0.5*z)))*exp(ey) #original dose-response curve;
  
  BB<-1000 #for faster run;
  Boot_resp_or<-array(0,c(BB,100))
  Boot_resp_weigthed<-array(0,c(BB,100))
  BB_or_dose1 <- rep(NA,BB)
  BB_or_dose2 <- rep(NA,BB)
  BB_or_dose3 <- rep(NA,BB)
  BB_orw_dose1 <- rep(NA,BB)
  BB_orw_dose2 <- rep(NA,BB)
  BB_orw_dose3 <- rep(NA,BB)
  
  for(i in 1:BB){
    ww<-rdirichlet(1,rep(1,100)) #weights on 100 patients;
    w<-rep(ww,each=10) #10 time points; weights on sampling patients (and their entire data history);
    
    # time-varying dose, predicted by time-varying and time-fixed covariates;
    ps.bb<-geeglm(d ~ x+u,id=id, data = lagdat,weights=w, family = 'gaussian',
                  scale.fix=FALSE, corstr="independence")
    # checked estimates, pretty accurate to simulated truth;
    
    est.sig_bb<-as.numeric(sqrt(summary(ps.bb)$geese$scale[1]))
    ps_bb<-dnorm(lagdat$d,mean=ps.bb$fitted.values,est.sig_bb)
    
    ### GPS as a covariate
    or_bb<-geeglm(logy ~ ns(ps_bb,3) + d,id=id, 
                  data = lagdat,weights=w,
                  family = 'gaussian', scale.fix=FALSE, corstr="independence")
    
    BB_or_dose1[i]<- mean(predict(or_bb, newdata = data.frame(d = rep(dose1,1000)), type = "response"))
    BB_or_dose2[i]<-mean(predict(or_bb, newdata = data.frame(d = rep(dose2,1000)), type = "response"))
    BB_or_dose3[i]<-mean(predict(or_bb, newdata = data.frame(d =rep(dose3,1000)), type = "response"))
    
    
    ### GPS in weighted GEE
    d_density <- npudens(~ d, data = data.frame(d = lagdat$d)) #density of d;
    
    # Evaluate the marginal density at each observed T
    f_d <- predict(d_density, newdata = data.frame(d = lagdat$d))
    
    ## SWGPS
    lagdat$sw <- f_d / ps_bb  
    #density weights; smart and simple!;
    #weights nolonger subject-specific, now subject-visit-specific;
    
    or_weighted_bb<-geeglm(logy ~ ns(d,3),id=id, data = lagdat,weights=w * sw,
                           family = 'gaussian', scale.fix=FALSE, corstr="independence")
    
    # summary(or_weighted_bb) #weighted estimates are quite different from regression ones, and lower-rank;
  
    # adding checks for three fixed dose levelsl
    BB_orw_dose1[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose1,1000)), type = "response"))
    BB_orw_dose2[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose2,1000)), type = "response"))
    BB_orw_dose3[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose3,1000)), type = "response"))
    
    ### Dose response curve
    d_vals <- seq(min(lagdat$d), max(lagdat$d), length.out = 100) #not affected by BB but do change over data iterations;
    
    #at each dose level we use the observed GPS to calculate the predicted outcome;
    predicted_d_or <- numeric(length(d_vals))
    predicted_d_weighted <- numeric(length(d_vals))
    
    # Define margin for d;
    dose_margin <- 1 #not sure about this, a minimum threshold?;
    for (j in 1:length(d_vals)) {
      # meet the margin
      idx <- which(abs(lagdat$d - d_vals[j]) <= dose_margin)
      if (length(idx) > 0) {
        
        # Extract observed GPS values for these observations
        observed_gps <- ps_bb[idx]
        predicted_or_Y <- predict(or_bb, newdata = data.frame(d = d_vals[j], ps_bb = observed_gps))
        predicted_weighted_Y <- predict(or_weighted_bb, newdata = data.frame(d = d_vals[j], ps_bb = observed_gps))
        
        # average over the same does level;
        predicted_d_or[j] <- mean(predicted_or_Y)
        predicted_d_weighted[j] <- mean(predicted_weighted_Y)
      } else {
        predicted_d_or[j] <- NA
        predicted_d_weighted[j] <- NA
      }
    }
    
    Boot_resp_or[i,] <- predicted_d_or
    Boot_resp_weigthed[i,] <- predicted_d_weighted


  }
  

# looking at density plots

# Create the plot
par(mfrow=c(2,2))
plot(NULL, ylab = 'Response', xlab = 'Dose', ylim=c(-5,20),main="COV BB",
     xlim=c(-3,12))
for (i in 1:dim(Boot_resp_or)[2]){
  lines(d_vals, Boot_resp_or[i,],type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
}
points(lagdat$d,lagdat$logy,pch=16,cex=0.8,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))

plot(NULL, ylab = 'Response', xlab = 'Dose', ylim=c(-5,20),main="WOR BB",
     xlim=c(-3,12))
for (i in 1:dim(Boot_resp_weigthed)[2]){
  lines(d_vals, Boot_resp_weigthed[i,],type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
}
points(lagdat$d,lagdat$logy,pch=16,cex=0.8,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))


