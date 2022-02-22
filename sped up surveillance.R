
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
library("purrr")

n<-30
fre<-30

s<-split(data,data$sim)


test1<-function(x){
  f<-max(pluck(x,"time"))
  #stti<-sample(1:g)
  #de<-c(f,stti)
}

maxtimes<-map(s,test1)
stti<-sample(1:fre,length(s),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = max(unlist(maxtimes)), by = fre))


test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=30,replace=FALSE),]
    print(d)
    print(match(g,y))
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      mylist<-list(q,m,theta=head(r5$theta,1),beta=head(r5$beta,1),sim=head(r5$sim,1))
      return(mylist)
    }
  }
}

dftest4<-map2(s,samp.time,test4)

#dftest4unlist<-data.frame(unlist(dftest4))

dftestlistdocall<-data.frame(do.call(rbind,dftest4))
dftestlistdocall$q<-as.numeric(dftestlistdocall$V1)
sum.q<-dftestlistdocall%>%group_by(beta,theta)%>%summarise_at(vars(q),list(q_mean = mean))

for(i in dftestlistdocall$theta){
  for (j in dftestlistdocall$beta){
ggplot(data=filter(dftestlistdocall,theta==40 ,beta==100))+geom_histogram(aes(x=q))
  }
}