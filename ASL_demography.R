setwd("H:/Matias WA Fisheries/Analyses/Sea lions")
rm(list=ls(all=TRUE))

library(popbio)   #matrix calculations

#---DATA SECTION---



#---PARAMETERS SECTION---

#source: Goldsworthy et al 2010
A=27 #maximum age
step=1.5 #1.5 year cycle

#ASL
#age specific survival (as proportion surviving from previous age class)

#survivorship  (lx)
surv.Simon=c(1,.355,.298,.282,.267,.253,.239,.226,.214,.203,.189,.174,.159,.143,.115,.081,.049,.024,0)/1.0026

#probability of surviving (px)
px.Simon=c(1,.354,.837,rep(.944,7),.93,.92,.91,.9,.8,.7,.6,.5,0)

#how calculate probability of surviving from lx (Caswell 2001)
#1.for birth-flow reproduction
#a=NULL;for(i in 2:(length(surv.Simon)-1)) a=c(a,(surv.Simon[i]+surv.Simon[i+1])/(surv.Simon[i-1]+surv.Simon[i])) 
#a=c(a,0)

#2.for birth-pulse, post-breeding census reproduction
#b=NULL;for(i in 2:(length(surv.Simon)-1)) b=c(b,(surv.Simon[i])/(surv.Simon[i-1])) 
#b=c(b,0)

#3.for birth-pulse, pre-breeding census reproduction
#p=NULL;for(i in 2:(length(surv.Simon)-1)) p=c(p,(surv.Simon[i+1])/(surv.Simon[i])) 
#p=c(p,0)

# plot(a,pch=19,xlab="age class",ylab="prob of survival")
# points(b,col=2,pch=19)
# points(p,col=3,pch=19)
# legend("center",c("birth-flow", "birth-pulse post bred", "birth-pulse pre bree"),bty='n',col=1:3,pch=19)

#age specific fecundity
Rep.rate.Simon=c(0,0,0,.2,.27,.34,.38,.41,.42,.42,.42,.41,.4,.375,.34,.29,.2,.1,0) #reproductive rate
fec.Simon=Rep.rate.Simon/px.Simon
fec.Simon[length(fec.Simon)]=0

fec.alt=c(rep(0,3),rep(.5,15),0)  #set fecundity to .5 (i.e. 1 pup per cycle, 0.5 sex ratio, cycle of 1.5 year already)

#how to calculate fertilities from fecundity
#mx=fec.alt  #average number of female offspring per female

#1.for birth-flow reproduction
#l_0.5=surv.Simon[1]*(surv.Simon[2])^0.5
#a=NULL;for(i in 1:(length(surv.Simon)-1)) a=c(a,l_0.5*((mx[i]+px.Simon[i]*mx[i+1])/2)) 
#a=c(0,a)

#2.for birth-pulse, post-breeding census reproduction
#b=mx*px.Simon

#3.for birth-pulse, pre-breeding census reproduction
#f=mx*surv.Simon[2]

#  plot(a,pch=19,ylim=c(0,.5),xlab="age class",ylab="fertility")
#  points(b,col=2,pch=19)
#  points(f,col=3,pch=19)
#  legend("topright",c("birth-flow", "birth-pulse post bred", "birth-pulse pre bree"),bty='n',col=1:3,pch=19)



#Northern fur seal (Barlow and Boveng 1991)
north.furseal.px=c(1,0.539,0.702,0.814,0.885,0.928,0.951,0.964,0.969,0.969,0.964,0.956,0.943,0.923,0.894,
                   0.854,0.799,0.726,0.631,0.515,0.384)
north.furseal.lx.cum=cumprod(north.furseal.px)

north.furseal.birth.rate=c(0,0.000,0.000,0.015,0.015,0.220,0.395,0.395,0.425,0.460,0.445,0.450,0.440,0.430,0.420,
                   0.410,0.390,0.330,0.355,0.265,0.255)
north.furseal.fec=north.furseal.birth.rate/north.furseal.px
north.furseal.fec[length(north.furseal.fec)]=0
Amax.NFS=20
step.NFS=1


#Stellar sea lion (birth pulse, post breeding census, Holmes and York 2003)
stellar.surv=c(.782,.782,.93,.909,.895,.884,.875,.867,.859,.853,.847,.841,.841,.836,.831,.827,.822,.818,.814,.81,.807,
               .803,.80,.797,.794,.791,.788,.785,.782,.78,.77)
stellar.fec=c(0,0,0,.1008,.1795,.2614,rep(.315,25))
stellar.fertility= stellar.surv*stellar.fec 
  
#Estimable pars
log.r=log(0.01)



#---PROCEDURE SECTION---

#------------LIFE TABLE ----------

  #-Estimation of population parameters                  
demographic=function(log.r,surv,fec.at.age)
{
  r=exp(log.r)
 # f=log.f   
  
  #survivorship
 # S=exp(-m)    
#  lx=S^(age-1)			
  lx=surv
  
  #reproductive schedules
  mx=fec.at.age		                            #proportion of female pups fer female
  lx_mx=lx*mx					                                            #reproductive rate
  lx_mx_X=lx_mx*age			                                        	#reproductive rate times age
  
  e_rx=exp(-r*age)
  lx_mx_erx=lx_mx*e_rx               
  EulerLotka=sum(lx_mx_erx)         #Euler-Lotka equation
  Ro=sum(lx_mx)                     #net reproductive rate
  G=sum(lx_mx_X)/Ro                 #mean generation length
  t2=log(2)/r                       #population doubling time

  epsilon=(1-EulerLotka)^2	 				#objective function
		
  return(list(EulerLotka=EulerLotka,epsilon=epsilon,r=r,Ro=Ro,G=G,t2=t2))
}



#---MAIN SECTION---
#Estimate population parameters   
fn=function(log.r,surv,fec.at.age)demographic(log.r,surv,fec.at.age)$epsilon

EDAD=LISTA=r.estim=list()
LISTA[[1]]=rbind(surv.Simon,fec.Simon)
LISTA[[2]]=rbind(surv.Simon,fec.alt)
LISTA[[3]]=rbind(north.furseal.lx.cum,north.furseal.fec)
EDAD[[1]]=c(0,A,step)
EDAD[[2]]=c(0,A,step)
EDAD[[3]]=c(0,Amax.NFS,step.NFS)  
  
for (i in 1: length(LISTA))
{
  start.age=EDAD[[i]][1]
  max.age=EDAD[[i]][2]
  STEP=EDAD[[i]][3]
  age=seq(start.age,max.age,STEP)  
  surv=LISTA[[i]][1,]
  fec=LISTA[[i]][2,]
  fit=optimize(fn,lower=log(0.0001),upper=log(0.99999),surv=surv,fec.at.age=fec)   
  r.estim[[i]]=exp(fit$minimum)
  
}

       

#--REPORT SECTION---

#fecundity rates and r
par(mfcol=c(2,1),mai=c(1.05,1.05,.2,.2))
plot(seq(0,A,step),fec.Simon,col=3,type='l',lwd=2,ylab="Fecundity",ylim=c(0,1))
lines(seq(0,A,step),fec.alt,col=2,lwd=2)
lines(seq(0,Amax.NFS,step.NFS),north.furseal.fec,col=1,lwd=2)
legend(0,1,paste("simon fec, r=",round(r.estim[[1]],5), sep=""),bty='n',lty=1,col=3,lwd=2)
legend(0,.925,paste("0.5 fec, r=",round(r.estim[[2]],5),sep=""),bty='n',lty=1,col=2,lwd=2)
legend(0,.85,paste("NFS, r=",round(r.estim[[3]],5),sep=""),bty='n',lty=1,col=1,lwd=2)

plot(seq(0,A,step),surv.Simon,col=3,type='l',lwd=2,ylab="Fecundity",ylim=c(0,1))
lines(seq(0,Amax.NFS,step.NFS),north.furseal.lx.cum,col=1,lwd=2)




#-----------LESLIE MATRIX-------------------
#put data on different species in a list
px=list();px[[1]]=px.Simon[-1];px[[2]]=north.furseal.px[-1];px[[3]]= stellar.surv
bx=list();bx[[1]]=Rep.rate.Simon[-1]; bx[[2]]=north.furseal.birth.rate[-1];bx[[3]]=stellar.fertility
type=c("ASL","NFS","SSL") 
demographic.data=list()
for (i in 1:length(px))
{
  PX=px[[i]]
  PX=PX[-length(PX)]
  BX=bx[[i]]
  n=length(BX)
  Data=matrix(0,nrow=n,ncol=n)
  diag(Data)[-nrow(Data)]=PX
  Data=rbind(matrix(0,nrow=1,ncol=n),Data)
  Data=Data[-(n+1),]
  Data[1,]=BX
  rownames(Data)=colnames(Data)=0:length(PX)
  demographic.data[[type[i]]]=Data
}

#Apply leslie matrix analyses to a list of different species
Net.rep.value=sapply(demographic.data, net.reproductive.rate)
LAMBDA=sapply(demographic.data, lambda)
r=log(LAMBDA)
Rep.value=sapply(demographic.data,reproductive.value)



#Reproductive value
v=reproductive.value(ASL.data)
dotchart(v, pch=16, xlab="Reproductive value")


#Stable age distribution
fullon<-eigen.analysis(ASL.data)

barplot(fullon$stable.stage, col="green", ylim=c(0,1), 
       ylab="Stable stage proportion", xlab="Stage class", main="Teasel")
box()

#op<-par(mfrow=c(2,2))
#projection matrix
image2(ASL.data, cex=.8, mar=c(0.5,3,4,1) )
title("ASL projection matrix", line=3)

#elasticities
image2(fullon$elasticities, cex=.8, mar=c(0.5,3,4,1) )
title("Elasticity matrix", line=3)

#sensititivies
## default is sensitivity for non-zero elements in matrix
image2(fullon$sensitivities, cex=.8, mar=c(0.5,3,4,1) )
title("Sensitivity matrix 1", line=3)

## use zero=FALSE to get sensitivities of all elements
image2(eigen.analysis(ASL.data, zero=FALSE)$sensitivities, cex=.8, mar=c(0.5,3,4,1) )
title("Sensitivity matrix 2", line=3)
#par(op)

#Net reproductive rate
Rep.rate=net.reproductive.rate(ASL.data)

#Generation time
Gen.time=generation.time(ASL.data)


#ACA: apply to a range of seals!!!!!
#Apply demographic analysis to a list of species
data(calathea)
sapply(calathea[9:12], net.reproductive.rate)

