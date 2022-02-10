  res<-NULL
detectionfunction<-function(data,sim,n,s){
  for (sim in unique(data$sim)){ ## looping over simulations
    ## we sample 20 hosts
    stti<-sample(1:s,1)
    infmax<-max(data$time)
    samp.time <- seq(from = stti, to = infmax, by = s) ## at sampling time 1 and 4
    ## for ease we make a temporary dataframe of the current simulation
    #temp <- data[data$sim==sim,]
    ## we sample n=20 hosts
    ## loop over the sampling times
    for (t in samp.time){##here you see that the hosts sampled is after the loop detecting
      ##time that hosts are detected as infected, thus changing every time you update the sampling
      test <- sample(data$who, n,replace=FALSE)
      
      ## get those hosts infection time
      inf.time <- data[data$who %in% test, "time"]
      #print(paste("at time",inf.time))
      ## if an infection time is anterior to the sampling time, we have a detection event
      m <- sum(inf.time <= t) ## we sum to know how many host are seen as infected
      ## we can also measure the true incidence at sampling time
      q <- mean(data$time<=t)
      ## and increment the result table in the loop
      res <- rbind(res, data.frame(sim, m, q, t=t, n,beta=head((data$beta[data$sim==sim]),n=1)))
      if (m>=1){
        break
      }
      print(data$sim)
    }
  }
  return(res)
}
n=10
s=10

dataoutput<-function(i=NULL){
  res<-detectionfunction(data=data,sim=sim,n=n,s=s)
  
  
  #predictionstore<-c()
  #for(b in unique(res$beta)){
}


cl <- makeCluster(mc <- getOption("cl.cores", 3))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("tidyverse"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1, fun=dataoutput, cl=cl) ## test with 10 first, but then replace 10 by 1000
#simtest<-clusterSplit(cl,seq=1:(iter*length(betavalues)))
#par_results<-cbind(par_results,simtest)
## stop the cluster
stopCluster(cl)
comparison_values <- do.call("rbind", par_results)

  avg.sur<-res%>% distinct(res$m,res$sim,.keep_all=TRUE)
  avg.sur2<-subset(avg.sur, !m==0)
  sum.q<-avg.sur2%>%group_by(beta)%>%summarise_at(vars(q),list(r_mean = mean))
  anq<-((rdatamean$r_mean)*n/s)
  
  absdif<-abs(anq-sum.q$r_mean)
  reldif<-absdif/sum.q$r_mean 
  
  datatemp<-data.frame(predicted=anq,mean=sum.q$r_mean,absdif=absdif,reldif=reldif,beta=rdatamean$beta)
