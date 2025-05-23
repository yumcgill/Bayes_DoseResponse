library(splines)
library(np)
library(MCMCpack)
library(geepack)

###simulate the data
S<-20
N<-1000
mux<-0.2
muu<-1.0
muz<-0.2
varx<-0.1
varu<-0.6
varz<-0.1
dose1<-1.00
dose2<-4.00
dose3<-7.00


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


# function to lag longitudinal data
x<-rnorm(1000,mux,sqrt(varx)) 
v1<-rnorm(1000,muu,sqrt(varu)) 
v2<-rnorm(100,muz,sqrt(varz)) 

#u<-rep(v1,10)
u <- v1
z<-rep(v2,10)
count<-c(1:100) #subjects;
obs<-c(1:1000) #observations;
id<-rep(count,10) #100 subjects each measured 10 times = 1000 observations;
time<-c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100),rep(6,100),
        rep(7,100),rep(8,100),rep(9,100),rep(10,100))

# treatment assignment mechanism;
bx<-4.0
bu<-2.0
bz<-1.0
vared<-1
ed<-rnorm(1000,0,sqrt(vared))
d<-(1.0+bx*x+bu*u+bz*z+ed) 
#normal dose; place holder for non-symmetric dose dist;
# x and u are time-varying covariates, z are fixed baseline variable;
# dose is time-varying not single point treatment assignment!;
expit <- function(x) {1/(1+exp(-x))}
y <- rbinom(1000, 1, p=expit(-2 + 0.5*d - 0.4*x + 0.2*u - 0.4*z))
table(y)


A<-data.frame(obs,id,time,y,x,u,z,d)
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
lagdat$x<-x
lagdat$u<-u
lagdat$z<-z
lagdat$d<-d


mud1<-mean(expit(-2 +0.5*dose1-0.4*x + 0.2*u - 0.4*z)) 
mud2<-mean(expit(-2 +0.5*dose2-0.4*x + 0.2*u - 0.4*z)) 
mud3<-mean(expit(-2 +0.5*dose3-0.4*x + 0.2*u - 0.4*z))

##sampling the stick breaking weights 
stick.breaking<-function(av,Nv){
  u<-rbeta(Nv,1,av)	
  v<-u
  w<-c(1,cumprod(1-u[-Nv]))
  return(v*w)
}

##alpha
al = 5


### split the data as a list
split.data<-function(data){
  split.data <- split(data, data$id)
  split.data <- 
    lapply(split.data, function(x) data.frame(x, obs.num = 1:nrow(x)))
  return(split.data)
}

sim.data<-lagdat
ps_gee <- geeglm(d ~ x+u+z,id=id, data = lagdat, family = 'gaussian', scale.fix=FALSE, corstr="independence")
est.sig<-as.numeric(sqrt(summary(ps_gee)$geese$scale[1]))
ps_tran<-dnorm(sim.data$d,ps_gee$fitted.values,est.sig)
pred<-predict(geeglm(logy ~ ns(ps_tran,3) + d,id=id, data = lagdat,family = 'gaussian', scale.fix=FALSE, corstr="independence"))
sim.data$pred<-pred
split.data_dp <- split.data(sim.data)

BB<-1000
## sample size for the new dataset
Nv<-500



Boot_resp_or_dp<-array(0,c(BB,100))
Boot_resp_weigthed_dp<-array(0,c(BB,100))
DP_or_dose1 <- rep(NA,BB)
DP_or_dose2 <- rep(NA,BB)
DP_or_dose3 <- rep(NA,BB)
DP_orw_dose1 <- rep(NA,BB)
DP_orw_dose2 <- rep(NA,BB)
DP_orw_dose3 <- rep(NA,BB)
for (i in 1:BB){ 
  ind<-as.numeric(sample(unique(sim.data$id),size=Nv,replace = TRUE,prob = table(sim.data$id)/length(sim.data$id)))
  u<-runif(Nv)
  res_ind<-which(u< al/(al+length(unique(sim.data$id))))
  w_dp1<-stick.breaking(al+length(unique(sim.data$id)),Nv)
  newdata<-NULL
  ## resample the new data
  for (j in 1:Nv){
    newdata1<-split.data_dp[[ind[j]]]
    newdata1$w_dp<-w_dp1[j]
    if (j %in% res_ind){
      newdata1$logy <- (rnorm(dim(newdata1)[1],newdata1$pred,0.5))}
    newdata <- rbind(newdata,newdata1) 
  }
  
  resamp.data<-newdata
  
  ps.bb<-geeglm(d ~ x+u,id=id, data = resamp.data,weights=w_dp, family = 'gaussian', scale.fix=FALSE, corstr="independence")
  est.sig_bb<-as.numeric(sqrt(summary(ps.bb)$geese$scale[1]))
  ps_bb<-dnorm(resamp.data$d,ps.bb$fitted.values,est.sig_bb)
  
  or_bb<-geeglm(logy ~ ns(ps_bb,3) + d,id=id, data = resamp.data,weights=w_dp,family = 'gaussian', scale.fix=FALSE, corstr="independence")
  
  DP_or_dose1[i]<- mean(predict(or_bb, newdata = data.frame(d = rep(dose1,1000)), type = "response"))
  DP_or_dose2[i]<-mean(predict(or_bb, newdata = data.frame(d = rep(dose2,1000)), type = "response"))
  DP_or_dose3[i]<-mean(predict(or_bb, newdata = data.frame(d =rep(dose3,1000)), type = "response"))
  
  ## SWGPS
  ### GPS in weighted GEE
  d_density <- npudens(~ resamp.data$d) #density of d;
  
  # Evaluate the marginal density at each observed T
  f_d <- predict(d_density, newdata = data.frame(d = resamp.data$d))
  resamp.data$sw <- f_d / ps_bb  
  #density weights; smart and simple!;
 
  or_weighted_bb<-geeglm(logy ~ ns(d,3),id=id, data = resamp.data, weights=w_dp * sw,
                         family = 'gaussian', scale.fix=FALSE, corstr="independence")
  DP_orw_dose1[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose1,1000)), type = "response"))
  DP_orw_dose2[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose2,1000)), type = "response"))
  DP_orw_dose3[i]<-mean(predict(or_weighted_bb, newdata = data.frame(d = rep(dose3,1000)), type = "response"))
  
  #at each dose level we use the observed GPS to calculate the predicted outcome;
  predicted_d_or_dp <- numeric(length(d_vals))
  predicted_d_weighted_dp <- numeric(length(d_vals))
  
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
      predicted_d_or_dp[j] <- mean(predicted_or_Y)
      predicted_d_weighted_dp[j] <- mean(predicted_weighted_Y,na.rm = TRUE)
    } else {
      predicted_d_or_dp[j] <- NA
      predicted_d_weighted_dp[j] <- NA
    }
  }
  
  Boot_resp_or_dp[i,] <- predicted_d_or_dp
  Boot_resp_weigthed_dp[i,] <- predicted_d_weighted_dp
  
  
}


par(mfrow=c(2,2))
plot(NULL, ylab = 'Response', xlab = 'Dose', ylim=c(-5,20),main="COV DP",
     xlim=c(-3,12))
for (i in 1:dim(Boot_resp_or_dp)[2]){
  lines(d_vals, Boot_resp_or_dp[i,],type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
}
points(lagdat$d,lagdat$logy,pch=16,cex=0.8,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))

plot(NULL, ylab = 'Response', xlab = 'Dose', ylim=c(-5,20),main="WOR DP",
     xlim=c(-3,12))
for (i in 1:dim(Boot_resp_weigthed_dp)[2]){
  lines(d_vals, Boot_resp_weigthed_dp[i,],type="l", col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
}
points(lagdat$d,lagdat$logy,pch=16,cex=0.8,col=rgb(red = 1, green = 0, blue = 0, alpha = 0.2))
