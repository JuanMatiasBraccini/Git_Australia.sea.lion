# EXPLORING MORTALITY RATE AND TOTAL MORTALITY ESTIMATIONS (based on Goldworthy et al 2010)

n.mortalities=12
n.shots=234
Obs.km.hrs=5794 #observed km hours
foraging.effort=sample(1:43,n.shots,replace=T) #dummy foraging effort (seal days/yr) for each fishing shot

# 1. Observed mortalities
#Trial of different ways of calculating mortality rate

num.seals=sample(c(1,0),n.shots, replace=TRUE,prob=c(0.1,0.9))
num.seals=c(rep(1,n.mortalities),rep(0,(n.shots-n.mortalities)))
effort=c(rep(5794/n.shots,n.shots)) #Km net hours, no specific info
effort=effort*exp(runif(n.shots,-0.1,0.1))#add variability

data=data.frame(num.seals,effort)

weight=effort/sum(effort) #if using this weight, then mean.weighted=total
#weight=c(.2,.3.,4.,) #assuming relative area (smaller area is more representative)
  
total=sum(data$num.seals)/sum(data$effort)
mean=mean(data$num.seals/data$effort)

mean.weighted=weighted.mean(data$num.seals/data$effort,weight)


#2. Bycatch rate estimation

#add foraging effort
data$forg.eff=foraging.effort

data=data[order(data$forg.eff),]#sort by foraging effort

#create bins of foraging effort
bins=5  #number of bins for sorting data
Q=quantile(data$forg.eff,probs=seq(0,1,1/bins))

data$Bins=factor(ifelse(data$forg.eff<Q[2],1,ifelse(data$forg.eff>=Q[2] & data$forg.eff<Q[3],2,ifelse(data$forg.eff>=Q[3] 
          & data$forg.eff<Q[4],3,ifelse(data$forg.eff>=Q[4] & data$forg.eff<Q[5],4,5)))))


#calculate observed total mortality and fishing effort per bin
Mortality.by.bin=by(data[, 1], data[,"Bins"], sum)
Mortality.by.bin=do.call(rbind,as.list(Mortality.by.bin))
Effort.by.bin=by(data[, 2], data[,"Bins"], sum)
Effort.by.bin=do.call(rbind,as.list(Effort.by.bin))

mean.bin.forg.eff=by(data[, 2], data[,"Bins"], median)
mean.bin.forg.eff=do.call(rbind,as.list(mean.bin.forg.eff))

Observed.mort.rate.per.bin=Mortality.by.bin/Effort.by.bin


# compare with means for each bin approach
data$ratio=data$num.seals/data$effort
#mean.bin.mort=by(data[, 5], data[,"Bins"], median)
#mean.bin.mort=do.call(rbind,as.list(mean.bin.mort))
mean.bin.mort=NULL
for(i in 1:bins)
  {
  a=subset(data,Bins==i)
  b=mean(a$ratio)
  mean.bin.mort=rbind(mean.bin.mort,b)
  }


# compare with simple means approach
Observed.mort.rate.means=mean(data$num.seals/data$effort)


#regressions
regression=lm(Observed.mort.rate.per.bin~mean.bin.forg.eff)
ANOVA=anova(regression)
plot(mean.bin.forg.eff,Observed.mort.rate.per.bin)
lines(mean.bin.forg.eff,predict(regression),col=2)

plot(data$forg.eff,data$num.seals)
