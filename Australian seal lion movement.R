 # MOVEMENT MODELLING OF SEA LIONS BASED ON GOLDWORTHY ET AL 2009

 library(MASS) #for fitting distributions
 library(trip) #for tracking data analysis
 library(geosphere)      #for spatial statistics and trigonometry in a sphere
 library(sp)
 
#-- PARAMETER SECTION
n.seals=10 
n=1000  #number of hits per seal
colony.lat=-35.789
colony.long=137.297
hotspot=c(141.282,-33.95372)   #foraging hotspot
 dispersion=20
 
#-- PROCEDURE SECTION
 
#1. Simulate tracking data
 
 walk.2d<-function(x)
{
  rw <- matrix(0, ncol = 2, nrow = n)

  # generate the indices to set the deltas
  indx <- cbind(seq(n), sample(c(1, 2), n, TRUE))

  # now set the values
  rw[indx] <- sample(c(-.1, .1), n, TRUE)


  # start walk from colony
  rw[1,]=cbind(colony.long,colony.lat)

  # cumsum the columns
  rw[,1] <- cumsum(rw[, 1])
  rw[,2] <- cumsum(rw[, 2])

  #focus walk to foraging hot spot
  steps=(nrow(rw)*.6):(nrow(rw)*.8)
  rw[steps,]=gcIntermediate(rw[length(steps)-2,], hotspot, n=length(steps), addStartEnd=F)
  rw[steps,]=jitter(rw[steps,],dispersion)
  
  steps=(nrow(rw)*.8):nrow(rw)
  rw[steps,]= matrix(rep(hotspot,length(steps)),length(steps),2,byrow=T)
  

  
  rw  # return value
}


rw <- lapply(1:n.seals,walk.2d)

colores=1:n.seals
plot(0, type="n",xlab="Long",ylab="Lat",main="Seal Tracking Simulation",xlim=c(132,142),ylim=c(-38,-30))
for(i in 1:length(rw))
{
  lines(rw[[i]],col=colores[i])
}
points(rw[[1]][1,1],rw[[1]][1,2],cex=3,col=1,pch=19)
 
 
rw=do.call(rbind,rw) #put as dataframe
 

 #convert to trip object and plot trajectories
rw=data.frame(rw,tms=Sys.time()+1:n,id=gl(n.seals,n))
 colnames(rw)[1:2]=c("Long","Lat")
coordinates(rw)=~Long+Lat 
rw=trip(rw,c("tms","id")) 
plot(rw,axes=T)
lines(rw)


#2. Construct map of time (in seconds) spent in cell
gt <- makeGridTopology(rw, cells.dim = c(10, 10)) #grid of 10 by 10 cells
 trg=tripGrid.interp(rw,method="count",grid=gt,dur=600) #10 mins. Straight counts, use "kde" for kernel density
 boundaries=bbox(trg)  #polygon boundaries
 count=trg[[1]]         #counts per cell
 
#plot
# image(trg, col = bpy.colors(200))
 image(trg, col = heat.colors(200))
 contour(trg,add=T)
summary(trg)  ## it's a SpatialGridDataFrame

 
#2. Construct probability distributions of distance from colony and depth


#3. Estimate probability of being in a given cell