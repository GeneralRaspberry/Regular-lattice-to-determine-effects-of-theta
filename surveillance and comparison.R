####################################Surveillance structure##################################################

df<-data
res<-NULL
for (sim in unique(df$sim)){ ## looping over simulations
  n <- 10 ## we sample 20 hosts
  stti<-sample(1:20,1)
  infmax<-max(df$time)
  samp.time <- seq(from = stti, to = infmax, by = 10) ## at sampling time 1 and 4
  ## for ease we make a temporary dataframe of the current simulation
  temp <- df[df$sim==sim,]
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