
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
library("purrr")


detectionfunction<-function(y){
  colmax(y)
  return(y[,1])
}

colmax<-function(data)lapply(data,max)

s<-split(data,data$sim)




  testfuncfunc<-function(f,time){
  stti<-sample(1:f,1)
  infmax<-max(time)
  samp.time <- seq(from = stti, to = infmax, by = s)
  return(samp.time)
  
}

testfunc<-map(s,testfuncfunc(f=30,time="time"))

nested_lapply<-function(data,fun){
  lapply(data,function(sublist){lapply(sublist,fun)})
}

nested_lapply(s,testfuncfunc(f=30,time="time"))


anotherfuckingfunction<-function(tm,n,f){
  stti<-sample(1:f,1)
  infmax<-max(tm)
  samp.time <- seq(from = stti, to = infmax, by = f)
}

tm<-"time"
n<-30
f<-30
map(s,anotherfuckingfunction())

test1<-function(x,g){
  f<-max(pluck(x,"time"))
  stti<-sample(1:g)
  de<-c(f,stti)
}

map(s,test1(s,1000))
