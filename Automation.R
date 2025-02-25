rm(list=ls())
##RUnning the experiment using the above tau leap function
##---------------------------------------------------------

library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=10, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 50 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)

## epidemic parameters

betavalues <- c(5,10,20)##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta <- 40
b <- 1
area.host<-1
infbegin<-1
iter<-5

##################################add a timer##############################################################

ts<-proc.time()

###########################################################################################################
##Concatenating a list of metric values
##-----------------------------------------

dim<-1000
hosts<-900
radiusCluster<-50
lambdaParent<-.05
lambdaDaughter<-25
randmod<-0

tempbind<-c()

sim_par <- function(i=NULL){
  for (j in betavalues){
    #for (l in 1:(length(betavalues)*iter)){
    # print(l)
    #l<-1  
    
    rExt=radiusCluster; #extension parameter -- use cluster radius
    xDeltaExt=dim+rExt;
    yDeltaExt=dim+rExt;
    numbparents<-rpois(1,xDeltaExt*lambdaParent)
    
    xxParent<-runif(numbparents,0-rExt,xDeltaExt)
    yyParent<-runif(numbparents,0-rExt,yDeltaExt)
    
    
    numbdaughter<-rpois(numbparents,(lambdaDaughter))
    sumdaughter<-sum(numbdaughter)
    
    
    
    thetaLandscape<-2*pi*runif(sumdaughter)
    
    rho<-radiusCluster*sqrt(runif(sumdaughter))
    
    
    
    xx0=rho*cos(thetaLandscape)
    yy0=rho*sin(thetaLandscape)
    
    
    xx<-rep(xxParent,numbdaughter)
    yy<-rep(yyParent,numbdaughter)
    
    
    
    xx<-xx+xx0
    
    yy<-yy+yy0
    
    booleInside=((xx>=0)&(xx<=dim)&(yy>=0)&(yy<=dim));
    #retain points inside simulation window
    xx=xx[booleInside]; 
    yy=yy[booleInside]; 
    
    landscape3<-ppp(x=xx,y=yy,window=owin(xrange=c(0,1000),yrange=c(0,1000)))
    
    dfland<-data.frame(landscape3)
    while (nrow(dfland)<hosts){
      dif<-hosts-nrow(dfland)
      extraparentxx<-sample(xxParent,dif,replace = TRUE)
      extraparentyy<-sample(yyParent,dif,replace = TRUE)
      extrathetaLandscape<-2*pi*runif(dif)
      extrarho<-radiusCluster*sqrt(runif(dif))
      newextracoodsxx<-extrarho*cos(extrathetaLandscape)
      newextracoodsyy<-extrarho*sin(extrathetaLandscape)
      extraxx<-extraparentxx+newextracoodsxx
      extrayy<-extraparentyy+newextracoodsyy
      dflandextra<-data.frame(x=extraxx,y=extrayy)
      dfland<-rbind(dfland,dflandextra)
    }
    
    equalhosts<-sample(1:nrow(dfland),hosts,replace=F)
    dfequal<-dfland[equalhosts,]
    
    
    landscape1<-(gridcenters(window=owin(xrange=c(0,dim),yrange=c(0,dim)),30,30))
    landscape1<-data.frame(landscape1)
    randfunction<-function(x){
      x<-replace()
    }
    randselect<-sample(1:nrow(dfland),floor(hosts*randmod),replace=F)
    landscape1[randselect,]<-dfequal[randselect,]
    landscape2<-ppp(x=landscape1$x,y=landscape1$y,owin(xrange=c(0,dim),yrange=c(0,dim)))
    
    data <- data.frame(x=landscape2$x, y=landscape2$y, id=1:hosts)
    
    ## design a function that will be called
    
    
    set.seed(seed=NULL)
    marks(landscape2)<- sample(c(rep(TRUE,infbegin), rep(FALSE, hosts-infbegin)))
    
    
    output <- tauLeapG(beta = j, theta = theta, b = b,
                       sigma = 0, delta.t = delta.t,
                       ppp = landscape2)
    
    temp <- output[[2]][,1:2][order(output[[2]][,2]),]
    temp<-cbind(temp,beta=j)#,sim=l)
    tempbind<-rbind(tempbind,temp)
    #l<-l+1
  }
  
  datatest1<-data.frame(time=tempbind$time, who=tempbind$who, x=landscape2$x[tempbind$who], y=landscape2$y[tempbind$who],beta=tempbind$beta)#,sim=tempbind$sim)
  
}


#}



## create a cluster with the set number of cores, say nmax-1
cl <- makeCluster(mc <- getOption("cl.cores", 3))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("tidyverse"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1:iter, fun=sim_par, cl=cl) ## test with 10 first, but then replace 10 by 1000
clusterEvalQ(cl,sim_par)
#simtest<-clusterSplit(cl,seq=1:(iter*length(betavalues)))
#par_results<-cbind(par_results,simtest)
## stop the cluster
stopCluster(cl)
## call cbind on your list of lines to find the matrix you expect
data <- do.call("rbind", par_results)

simtest2<-rep((1:(iter*length(betavalues))),hosts)
simtest2<-simtest2[order(simtest2)]

data<-cbind(data,sim=simtest2)



##################################add a timer############################################################
proc.end<-proc.time()-ts
proc.end
beep()

t2<- proc.time()

head(data)
data<-data.frame(data)
times <- sort(unique(data$time))



data_logistic <- function(i=NULL){
  data  %>% group_by(sim) %>%
    do(data.frame(beta=sample(.$beta,size = length(times)), time=times, infected=sapply(times, function(x) sum(.$time <= x))))
}
## make a logistic df from this data
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("data","times"),envir = environment())
par_data_logistic<-parLapply(1,fun=data_logistic,cl=cl)
stopCluster(cl)
data_log<-data.frame(par_data_logistic)



## prepare a logistic function of r to fit
temp <- filter(data_log, infected <= (hosts/4))
#temp$simdigit<-as.numeric(temp$sim)

###############################linear regression#####################################################################


intloop<-c()
rloop<-c()
betaloop<-c()
simloop<-c()

r_lnreg<-function(i=NULL){
  for (g in unique(temp$sim)){
    betaVal <- unique(temp$beta[temp$sim==g])
    lmoutput<-lm(formula=log(infected)~time,data=filter(temp,sim==g))
    lmoutput_int<-as.numeric(exp(lmoutput$coefficients[1]))
    lmoutput_r<-as.numeric(lmoutput$coefficients[2])
    #lmoutputtransform<-exp(lmoutput1)
    #rloop<-c(rloop,lmoutputtransform)
    intloop<-c(intloop,lmoutput_int)
    rloop<-c(rloop,lmoutput_r)
    betaloop<-c(betaloop,betaVal)
    simloop<-c(simloop,g)
  }
  rdata<-data.frame(int=intloop,r=rloop,beta=betaloop,sim=simloop)
  
  temp$predExp <- NA
  temp$predLog <- NA
  for(g in unique(temp$sim)){
    temp$predExp[temp$sim==g] <- rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g])
    temp$predLog[temp$sim==g] <- (hosts*rdata$int[rdata$sim==g]*exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))/(hosts + rdata$int[rdata$sim==g]*((exp(rdata$r[rdata$sim==g]*temp$time[temp$sim==g]))-1))
  }
 rdata 
}
#another cluster
cl <- makeCluster(mc <- getOption("cl.cores", 3))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("temp","hosts","rloop","betaloop","simloop","intloop"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg,cl=cl)
stopCluster(cl)
rdataset <- do.call("rbind", par_r1)


rdatamean<-rdataset%>%group_by(beta)%>%summarise_at(vars(r),list(r_mean = mean))
#library(tidyverse)
#tempLongPreds <- pivot_longer(temp, cols=c(infected,predExp,predLog), names_to="estimate", values_to="inf")

#ggplot(tempLongPreds, aes(x=time)) + 
# geom_line(aes(y=inf,group=interaction(beta,sim,estimate),colour=as.factor(estimate))) + 
#facet_wrap(facets="beta") +
# ylim(c(0,max(temp$infected)))

logis <- function(t, r, K=1, s=0, q0){
  pmin(
    K*q0*exp(r*(t+s)) / (K + q0*(exp(r*(t+s)) - 1)),
    K) # numerical errors can happen for high r and sigma
}

predstore<-c()
for (s in unique(rdatamean$r_mean)){
  calstore<-logis(r=s, t=times, K=hosts, q0=1)
  predstore<-c(predstore,calstore)
}  
pred_data <- data.frame(time=rep(times,length(rdatamean$beta)), infected=predstore,beta=rep(rdatamean$beta,each=length(times)))

ggplot(data_log) + geom_line(aes(x=time, y=(infected/hosts), group=interaction(beta,sim),colour=as.factor(beta)), size=.2) +
  geom_line(data=filter(pred_data, infected<(hosts)), aes(x=time, y=(infected/hosts), group=beta, colour=as.factor(beta)), size=2)+
  #ggtitle(paste0("Figure check"))+
  #geom_line(data=filter(pred_data1, infected<hosts), aes(x=time, y=(infected/hosts), group=pred_data1$beta, colour=as.factor(pred_data$beta)), size=5)+
  theme_tufte()+
  labs(x="Time",
       y="Prevalence",
       colour="Beta")
##################################################################################################################################

####################################Surveillance structure##################################################


detectionfunction<-function(data,sim,n,s){
res<-NULL
for (sim in unique(data$sim)){ ## looping over simulations
  n <- 10 ## we sample 20 hosts
  stti<-sample(1:10,1)
  infmax<-max(data$time)
  samp.time <- seq(from = stti, to = infmax, by = 10) ## at sampling time 1 and 4
  ## for ease we make a temporary dataframe of the current simulation
  temp <- data[data$sim==sim,]
  ## we sample n=20 hosts
  ## loop over the sampling times
  for (t in samp.time){##here you see that the hosts sampled is after the loop detecting
    ##time that hosts are detected as infected, thus changing every time you update the sampling
    test <- sample(temp$who, n,replace=FALSE)
    
    ## get those hosts infection time
    inf.time <- temp[temp$who %in% test, "time"]
    #print(paste("at time",inf.time))
    ## if an infection time is anterior to the sampling time, we have a detection event
    m <- sum(inf.time <= t) ## we sum to know how many host are seen as infected
    ## we can also measure the true incidence at sampling time
    q <- mean(temp$time<=t)
    ## and increment the result table in the loop
    res <- rbind(res, data.frame(sim, m, q, t=t, n))
    if (m>=1){
      break
    }
  }
}
  return(res)
}

res<-detectionfunction(data=data,sim=sim,n=10,s=10)

avg.sur<-res%>% distinct(res$m,res$sim,.keep_all=TRUE)
avg.sur2<-subset(avg.sur, !m==0)
sum.q<-mean(avg.sur2$q)
anq<-(mean_r*10/20)

ggplot(avg.sur2,aes(x=q))+geom_histogram(color="black",fill="lightblue",aes(y=..count../sum(..count..)),binwidth=0.002)+
  geom_vline(aes(xintercept=mean(avg.sur2$q)),color="blue",linetype="dashed")+
  geom_vline(aes(xintercept=(anq)),color="blue")+
  ylab("Probability")+
  xlim(NA,1)+
  ylim(NA,0.2)

absdif<-abs(anq-sum.q)
reldif<-absdif/sum.q 

#################################Calculating 95% confidence limits#########################################
sum.q #mean discovery prevalence for simulations 
sd95<-sd(avg.sur2$q) #standard deviation of simulations
sample95<-sqrt(1000) #sample size for simulations
zvalue95<-1.96 #selected confidence interval

perct2.5<-zvalue95*sd95/sample95
Upper<-sum.q+perct2.5
Lower<-sum.q-perct2.5

##########################Adding the cluster metric#######################################################

kk<-Kest(landscape)
plot(kk)
kk_iso<-kk$iso
kk_pois<-kk$theo

kk_div_na<-kk_iso/kk_pois
kk_div_0<-replace_na(kk_div_na,0)
kk_mean<-round(mean(kk_div_0),3)

######################################creating a data table from all the simulations#######################

Figure2<-data.frame(sum.q,Upper,Lower,anq,mean_r)
Figure2namefile