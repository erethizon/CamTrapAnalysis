### Chapter 6

source("TEAM library 1.7.R")
library(chron) 
library(reshape)
library(ggplot2)
library(vegan)
library(unmarked)
library(AICcmodavg)
library(MuMIn)
library(plyr)
library(R2jags)

### loading data
team_data<-read.csv(file="teamexample.csv", sep=",",h=T,stringsAsFactors=F)
iucn.full<-read.csv("IUCN.csv", sep=",",h=T)
iucn<-iucn.full[,c("Class","Order","Family","Genus","Species")]
team<-merge(iucn, team_data, all.y=T)
fd<-fix.dta(team)
yr2009<-fd[fd$Sampling.Event =="2009.01" & fd$Class=="MAMMALIA",]
 

### load covariate data
cov<-read.table("covariates.txt", header=TRUE)
workingcam<-which(cov$Sampling.Unit.Name %in% unique(yr2009$Sampling.Unit.Name)) # removing cameras that did not work
cov.or<-cov[workingcam, ] # retain only working cameras in 2009
cov.num<-cov.or[,sapply(cov.or,is.numeric)]
cov.std<-decostand(cov.num,method="standardize")
cov.fac<-cov.or[,sapply(cov.or,is.factor)]  # extract factor variables
covs<-data.frame(cov.fac, cov.std)
covs


## create matrices for each species
mat.udz.09<-f.matrix.creator(yr2009)
names(mat.udz.09) 
naivetable<-naive(mat.udz.09) 
naivetable


#======================================#
#  Cercocebus sanjei; Sanje mangabey   #
#======================================#
Cs<-shrink(mat.udz.09[["Cercocebus sanjei"]],5)
umCs<-unmarkedFrameOccu(y=Cs,siteCovs= covs)

m0<- occu(~1~1,umCs)
d1<- occu(~edge~1,umCs)
d2<- occu(~border~1,umCs)
d3<- occu(~edge+border~1,umCs)
o1<- occu(~1~border,umCs)
o2<- occu(~1~habitat,umCs)
o3<- occu(~1~habitat+border,umCs)
m1<- occu(~edge~border,umCs)         
m2<- occu(~border~border,umCs)
m3<- occu(~edge+border~border,umCs)  
m4<- occu(~edge~habitat,umCs)
m5<- occu(~border~habitat,umCs)
m6<- occu(~edge+border~habitat,umCs)
m7<- occu(~edge+border~habitat+border,umCs)

## examine model m1 as an example:
m1
backTransform(linearComb(m1, coefficients = c(1, 0), type = "det")) 
backTransform(linearComb(m1, coefficients = c(1, 0), type = "state"))


## model selection
dlist<-fitList(Nullo = m0,d1=d1,d2=d2,d3=d3,o1=o1,o2=o2,o3=o3,m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7)
selmod<-modSel(dlist,nullmod="Nullo")
selmod   

newhab<-data.frame(habitat=c("Deciduous","Montane"))
pred<-predict(o2,type="state",newdata=newhab,appendData=T)

ggplot(pred,aes(x=habitat,y=Predicted))+
  geom_point(size=4) +
  ylab("Predicted Psi Cercocebus sanjei") +
  theme_bw()+
  geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE), width=.2)

# prepare the raster matrix with standardized covariates
map<-read.table("covs100x100.txt",h=T)  # covs100x100.txt is a matrix with the covariate values on a grid of points
mapst<-data.frame(x=map$x,y=map$y, habitat=map$habitat,   # standardize the matrix using mean and sd of covs measured at camera points
                  edge=(map$edge-mean(cov.or$edge))/sd(cov.or$edge),
                  border=(map$border-mean(cov.or$border))/sd(cov.or$border),
                  river=(map$river-mean(cov.or$river))/sd(cov.or$river))
# map
predmap<-predict(o2,type="state",newdata=mapst,appendData=T) # it takes ~ 30 sec
levelplot(predmap$Predicted ~ x + y, map, aspect="iso", xlab="Easting (m)", ylab="Northing (m)", col.regions=terrain.colors(100))



#===============================================#
# Rhynchocyon udzungwensis; Grey-faced sengi    #
#===============================================#
Ru<-shrink(mat.udz.09[["Rhynchocyon udzungwensis"]],5)
umRu<-unmarkedFrameOccu(y=Ru,siteCovs=covs)

m0<- occu(~1~1,umRu)
d1<- occu(~edge~1,umRu)
d2<- occu(~border~1,umRu)
o1<- occu(~1~edge,umRu)
o2<- occu(~1~border,umRu)
o3<- occu(~1~habitat,umRu)
o4<- occu(~1~river,umRu)
o5<- occu(~1~edge+habitat+border,umRu)
o6<- occu(~1~habitat+border,umRu)
o7<- occu(~1~edge+habitat,umRu)

dlist<-fitList(Nullo=m0,d1=d1,d2=d2,o1=o1,o2=o2,o3=o3,o4=o4,o5=o5,o6=o6,o7=o7)
sel<-modSel(dlist,nullmod="Nullo")
sel
best<-list(o7,o3,o5)
avgmod <- model.avg(best, fit=T)
summary(avgmod)
modnames<-as.character(c(o5,o3,o6))
modavgpred(best,modnames=modnames,newdata=site.cov,parm.type="psi")

# map

predmap<-predict(avgmod,type="state",newdata=mapst,appendData=T) # it takes about 253 sec
levelplot(predmap$fit ~ x + y, map, aspect="iso", xlab="Easting (m)", ylab="Northing (m)",col.regions=terrain.colors(100))


#########################
# Multiseason analyses  #
#########################

## reload data
teamc<-read.csv(file="team.yr2009_2013.csv", sep=",",h=T,stringsAsFactors=F)
fd<-fix.dta(teamc) 

#Separate data by year
samp<-unique(fd$Sampling.Event)
res<-numeric()
for(i in 1:length(samp)){
  temp<-which(fd$Sampling.Event==samp[i])
  fd2<-fd[temp,]
  res<-c(res,list(fd2))
}
data.byYear<-res
names(data.byYear)<-samp

data.byYear<-lapply(X=data.byYear,f.minusBirds)

mats<-sapply(data.byYear,f.matrix.creator)

#======================================#
#  Cercocebus sanjei; Sanje mangabey   #
#======================================#
spl<-sapply(mats, "[[", "Cercocebus sanjei")  # this is the list with all the years for this species
sl<-llply(spl,.fun=function (x) shrink(x,5))  # apply shrink to the list
mat.yrs<-adjMult(sl)                          # transform the list in a dataframe adjusting the number of working cameras and sampling events

colnames(mat.yrs) # check periods after shrinking

#JAGS
J<-array(as.matrix(mat.yrs),dim=c(60,25,5)) #60 sites, 25 events, 5 years

# the model
modJ.bug<- function () {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   p[k] ~ dunif(0, 1)
   }
p[nyear] ~ dunif(0, 1)

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1)
   for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i

# Observation model
for (i in 1:nsite){
   for (j in 1:nrep){
      for (k in 1:nyear){
         muy[i,j,k] <- z[i,k]*p[k]
         y[i,j,k] ~ dbern(muy[i,j,k])
         } #k
      } #j
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1]<-sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k] <- psi[k]/psi[k-1]
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
}

# Initial values
zst <- apply(J, c(1, 3), sum, na.rm=T)	# Observed occurrence as inits for z
zsti<-ifelse(zst>0,1,0)
jags.inits <- function(){ list(z = zsti)}
jags.params <- c("psi", "phi", "gamma", "p", "n.occ", "turnover")

# call jags
jagsfit<-jags(data=list(y = J, nsite = dim(J)[1], nrep = dim(J)[2], nyear = dim(J)[3]), inits=jags.inits, jags.params,
n.iter=2000, n.chains=3, model.file=modJ.bug)

print(jagsfit, dig=2)

# checking convergence and diagnostic
library(coda)
summary(jagsfit)
codaout.jags <- as.mcmc(jagsfit)
plot(codaout.jags , ask=TRUE)
densityplot(codaout.jags)

# naive for Cercocebus sanjei
nv<-function (x) sum(ifelse(rowSums(x, na.rm=T)>=1,1,0))/nrow(x)
naive<-aaply(J,3,.fun=nv)  # naive occupancy

# Summarize posteriors and plot
meanB<-jagsfit$BUGSoutput$mean$psi
CRI2.5<-jagsfit$BUGSoutput$summary[20:24,3] # credible interval
CRI97.5<-jagsfit$BUGSoutput$summary[20:24,7]

plot(2009:2013,meanB, type="b",xlab = "Year", ylab = "Occupancy Sanje mangabey", ylim = c(0,1))
segments(2009:2013, CRI2.5, 2009:2013,CRI97.5, lwd = 1)
points(2009:2013, naive, type = "b", col = "blue", pch=16) # add this line to compare real occupancy with naive (blue line)
legend(x=2011, y=0.2, legend=c("estimated occupancy","naive occupancy"), col=c("black","blue"), pch=c(1,16), lty=c(1,1), bty="n") 



