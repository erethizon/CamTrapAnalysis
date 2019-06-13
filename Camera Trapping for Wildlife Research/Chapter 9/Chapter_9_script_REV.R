######################################
# Chapter 9 - R and JAGS model codes
# Simone Tenan
#
# last modified Wednesday 17 June 2015 
######################################


###############################
# Extract the data 
###############################

# set your working directory
setwd("/path/to/your/folder/")

### extract 
source("TEAM library 1.7.R")
library(chron)
library(reshape)

#load the data and fix them (see Chapter 5)
team_data<-read.csv(file="teamexample.csv", sep=",",h=T,stringsAsFactors=F)
iucn.full<-read.csv("IUCN.csv", sep=",",h=T)
iucn<-iucn.full[,c("Class","Order","Family","Genus","Species")]
team<-merge(iucn, team_data, all.y=T)
data<-fix.dta(team)
data<- droplevels(data[data$bin!="Homo sapiens", ]) # remove Homo sapiens from data set

# extract the mammals   
mam<-data[data$Class=="MAMMALIA",]

# select year by year, extract events by species and camera trap (example for 2009)
ev2009<-event.sp(dtaframe=mam, year=2009.01, thresh=1440)  
rownames(ev2009)<-ev2009[["Sampling.Unit.Name"]]

# transpose matrix and save it
sp<-t(ev2009)                              
write.table(sp,"spudz_2009.txt")
#write.csv(sp, "spudz_2009.csv")

# calculate camera days (as done in Chapter 5)
camera_days<-cam.days(data,2009.01)
summary(camera_days[,2:4]) # check the maximum number of days
write.table(camera_days, file="camera_days_2009.txt",quote=F, sep="\t",row.names = F)



###############################
# Static multi-species occupancy model (Section 9.3.1) 
###############################
# set your working directory
setwd("/path/to/your/folder/")

# load detection data
Y <- as.matrix(read.table(file="spudz_2009.txt",header=T,sep=" "))
# load effort
effort <- read.table(file="camera_days_2009.txt",header=T,sep="\t")

# load libraries
library(snow)
library(rjags)
library(dclone)

# set seed
set.seed(1980)

# BUGS model
modelFilename = "smsom.txt"

cat("
model {

# Priors for community-level parameters
omega ~ dunif(0,1)
psi.mean ~ dunif(0,1)
beta <- log(psi.mean) - log(1-psi.mean)	# logit(psi.mean)
p.mean ~ dunif(0,1)
alpha <- log(p.mean) - log(1-p.mean)	# logit(p.mean)
sigma.psi ~ dunif(0,10)
sigma.p ~ dunif(0,10)
tau.psi <- pow(sigma.psi,-2)
tau.p <- pow(sigma.p,-2)

# Likelihood
for (i in 1:M) {
	w[i] ~ dbern(omega)

	# occupancy
	phi[i] ~ dnorm(beta, tau.psi)
	# detectability
	eta[i] ~ dnorm(alpha, tau.p)
	
	logit(psi[i]) <- phi[i]
	mu.psi[i] <- psi[i]*w[i]
	logit(p[i]) <- eta[i]
	
	for (j in 1:n.site) {
		Z[i,j] ~ dbern(mu.psi[i])

		mu.p[i,j] <- p[i]*Z[i,j]
		Y[i,j] ~ dbin(mu.p[i,j], K[j])
	}
}

# compute species richness
	N <- sum(w[])

}
", fill=TRUE, file=modelFilename)



# number of sampling occasions for each trap
K <- effort$ndays

# number of traps
nsites <- dim(Y)[2]

# number of observed species
nspecies <- dim(Y)[1]

# Augment data set
nzeros <- 100
Y_aug <- rbind(Y,matrix(0,nrow=nzeros,ncol=nsites))

# Latent states
w <- c(rep(1, nspecies), rep(NA, nzeros))

# Parameters monitored
parameters <- c("omega","psi.mean","sigma.psi","p.mean","sigma.p","N")

# data
bugs.data <- list(M=(nspecies+nzeros),n.site=nsites,K=K,Y=Y_aug, Z=(Y_aug>0)*1,w=w)

# initial values
inits <- function() { list(omega=runif(1),
                           psi.mean=runif(1), p.mean=runif(1),
                           sigma.psi=runif(1,0,4), sigma.p=runif(1,0,4))}

#mcmc settings
n.adapt <- 5000		#pre-burnin
n.update <- 10000	#burnin
n.iter <- 30000		#iterations post-burnin
thin <- 10
chains<-3

# run the model and record the run time
cl <- makeCluster(chains, type = "SOCK")
start.time = Sys.time()
out <- jags.parfit(cl, data = bugs.data, 
                  params = parameters, 
                  model = "smsom.txt",
                  inits = inits,
                  n.adapt = n.adapt,	 
                  n.update = n.update,
                  n.iter = n.iter,
                  thin = thin, n.chains = chains)

end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes\n', sep=''))
stopCluster(cl)

#Posterior computed in 3.37720005909602 minutes

# Summarize posteriors
summary(out)

# check output
xyplot(out[,c("omega","psi.mean","sigma.psi","p.mean","sigma.p"), drop=F])
densityplot(out[,c("omega","psi.mean","sigma.psi","p.mean","sigma.p"), drop=F])
acfplot(out[,c("omega","psi.mean","sigma.psi","p.mean","sigma.p"), drop=F])

### plot species richness (Fig. 9.1)

# bind chains for the plot
out2 <- mcmc(do.call(rbind, out))

# plot
par(mar = c(5,4.5,4,1.5)+.1, cex.axis=1.5, cex.lab=2, tcl = 0.25)
hist(out2[,"N"],breaks=seq(20,100,by=1),main="",xlab="Species richness",ylab="Posterior probability",freq=F)
abline(v=nspecies,col="red",lwd=3)


###############################
# Dynamic multi-species occupancy model (Section 9.4.1) 
###############################
# set your working directory
setwd("/path/to/your/folder/")

# load detection data and effort for the period 2009-2013
# get files' names
filesY <- list.files(pattern="spudz_")
names_objY <- paste("Y_",substring(filesY, first=9, last=10),sep="")
fileseff <- list.files(pattern="camera_days_")
names_objeff <- paste("effort_",substring(fileseff, first=15, last=16),sep="")

# open files and assign them names
for(i in seq(along=filesY)) {
	assign(names_objY[i], as.matrix(read.table(filesY[i],header=T,sep=" ")))
	assign(names_objeff[i], read.table(fileseff[i],header=T,sep="\t"))
}

# check dimensions of each data set
dimcheck <- sapply(mget(names_objY),FUN=dim)
dimcheck
# the number of detected species ranges from 24 to 28
# the number of active cameras ranges from 58 to 60

### let's keep only sites active every year for at least two sampling occasions
# get trap names for each year
Yl <- list()
for (i in 1:length(names_objeff)){
	Yl[[i]] <- get(names_objeff[i])[,1]
}

# homogenize trap names between detection matricies and effort data sets
Yl2 <- lapply(Yl,function(x) gsub("-",".",x))
for (i in 1:length(names_objeff)){
	tmp <- get(names_objeff[i])
	tmp[,1] <- Yl2[[i]]
	assign(paste(names_objeff[i],"b",sep=""),tmp)
}

# check which traps in a specific year (Yl2[[n]]) were not active in other years (unlist(Yl2[-n]))
traptodel <- lapply(1:length(Yl2), function(n) setdiff(unlist(Yl2[-n]),Yl2[[n]]))
traptodel2 <- unlist(traptodel)
traptodel3 <- traptodel2[-duplicated(unlist(traptodel2))]
traptodel3


# delete traps not active every year
names_objeff_b <- objects(pattern="b$")

for(i in 1:length(names_objeff_b)){
	col_todel_Y <- which(colnames(get(names_objY[i])) %in% traptodel3)
	row_todel_eff <- which(Yl2[[i]] %in% traptodel3)
	assign(paste("Y",i,sep=""),get(names_objY[i])[,-col_todel_Y])
	assign(paste("effort",i,sep=""),get(names_objeff_b[i])[-row_todel_eff,])
}

### insert all-zeros rows for species not detected in a certain year
names_objY_b <- objects(pattern="^Y[1-5]")

# get species names for each year-specific dataset and a list of species names detected at least once
colnamesYs <- data.frame(matrix(nrow=28, ncol=5))
allnames <- NULL
for(i in 1:(length(names_objY_b))){
	tmp <- rownames(get(names_objY_b[i]))
	colnamesYs[(1:length(rownames(get(names_objY_b[i])))),i] <- tmp
	allnames <- c(allnames,tmp)
}
allnames2 <- allnames[-(which(duplicated(allnames)==T))]
allnames3 <- allnames2[order(allnames2)]

# total number of observed species during the period 2009-2013
allnames3


# see which species are not present in the different years
for(i in 1:(length(names_objY_b))){
	print(which(sapply(allnames3,"%in%",colnamesYs[,i])==FALSE))
	print("---------------------------")
}



# edit the detection matrices in order to incorporate undetected species
not_detected <- rep(0,dim(Y1)[2])
Y1b <- rbind(Y1[1:9,],not_detected,Y1[10:12,],not_detected,Y1[13:26,],not_detected,not_detected)
Y2b <- rbind(Y2[1:6,],not_detected,Y2[7:12,],not_detected,Y2[13:18,],not_detected,Y2[19:25,],not_detected,Y2[26,])
Y3b <- rbind(Y3[1:6,],not_detected,Y3[7:11,],not_detected,Y3[12:28,])
Y4b <- rbind(Y4[1:6,],not_detected,Y4[7:8,],not_detected,Y4[9:11,],not_detected,Y4[12:13,],not_detected,Y4[14:16,],not_detected,Y4[17:20,],not_detected,Y4[21:24,])
Y5b <- rbind(Y5[1:6,],not_detected,Y5[7:8,],not_detected,Y5[9:10,],not_detected,not_detected,Y5[11:26,])


#### augment detection data
names_objY_c <- objects(pattern="^Y[1-5]b")
nzeros <- 100
for(i in 1:(length(names_objY_c))){
	nc <- dim(get(names_objY_c[i]))[2]
	assign(paste(names_objY_c[i],"_aug",sep=""),
	       rbind(get(names_objY_c[i]), matrix(0, nrow=nzeros, ncol=nc))
	      )
}

# bind augmeted matricies
library(abind)
Yaug_tot <- abind(Y1b_aug,Y2b_aug,Y3b_aug,Y4b_aug,Y5b_aug, along=3)
attr(Yaug_tot, "dimnames") <- NULL

# number of sampling occasions for each camera in each year
names_objeff_c <- objects(pattern="^effort[1-5]")
K_tot <- matrix(NA, nrow=dim(Yaug_tot)[2], ncol=length(names_objeff_c))
for(i in 1:length(names_objeff_c)){
	K_tot[,i] <- get(names_objeff_c[i])[,"ndays"]
}


# Target species for the WPI
speciesWPI <- c("Atilax paludinosus",
"Bdeogale crassicauda",
"Cephalophus harveyi",
"Cephalophus spadix",
"Genetta servalina",
"Loxodonta africana",
"Mellivora capensis",
"Nandinia binotata",
"Nesotragus moschatus",
"Panthera pardus")

nspeciesWPI <- length(speciesWPI)


# find out to which rows of the detection matrix these target species correspond
id_speciesWPI <- which(allnames3 %in% speciesWPI)


# load libraries
library(snow)
library(rjags)
library(dclone)

# set seed
set.seed(1980)

# BUGS model
modelFilename = "dmsom.txt"

cat("
model {
# priors
omega ~ dunif(0,1)
psiMean ~ dunif(0,1)
for (t in 1:T) {
	pMean[t] ~ dunif(0,1)
}
for (t in 1:(T-1)) {
	phiMean[t] ~ dunif(0,1)
	gamMean[t] ~ dunif(0,1)
}
lpsiMean <- log(psiMean) - log(1-psiMean)
for (t in 1:T) {
	lpMean[t] <- log(pMean[t]) - log(1-pMean[t])
}
for (t in 1:(T-1)) {
	lphiMean[t] <- log(phiMean[t]) - log(1-phiMean[t])
	lgamMean[t] <- log(gamMean[t]) - log(1-gamMean[t])
}
lpsiSD ~ dunif(0,10)
lpsiPrec <- pow(lpsiSD,-2)
lpSD ~ dunif(0,10)
lphiSD ~ dunif(0,10)
lgamSD ~ dunif(0,10)
for (t in 1:T) {
	lpPrec[t] <- pow(lpSD,-2)
}
for (t in 1:(T-1)) {
	lphiPrec[t] <- pow(lphiSD,-2)
	lgamPrec[t] <- pow(lgamSD,-2)
}

# likelihood
for (i in 1:M) {
	# initial occupancy state at t=1
	w[i] ~ dbern(omega)
	b0[i] ~ dnorm(lpsiMean, lpsiPrec)T(-12,12)
	lp[i,1] ~ dnorm(lpMean[1], lpPrec[1])T(-12,12)
	p[i,1] <- 1/(1+exp(-lp[i,1]))
	for (j in 1:n.site) {
		lpsi[i,j,1] <- b0[i]
		psi[i,j,1] <- 1/(1 + exp(-lpsi[i,j,1]))
		mu.z[i,j,1] <- w[i] * psi[i,j,1]
		Z[i,j,1] ~ dbern(mu.z[i,j,1])
		mu.y[i,j,1] <- p[i,1]*Z[i,j,1]
		Y[i,j,1] ~ dbin(mu.y[i,j,1], K_tot[j,1])
	}
	# model of changes in occupancy state for t=2, ..., T
	for (t in 1:(T-1)) {
		lp[i,t+1] ~ dnorm(lpMean[t+1], lpPrec[t+1])T(-12,12)
		p[i,t+1] <- 1/(1+exp(-lp[i,t+1]))
		c0[i,t] ~ dnorm(lgamMean[t], lgamPrec[t])T(-12,12)
		d0[i,t] ~ dnorm(lphiMean[t], lphiPrec[t])T(-12,12)
		for (j in 1:n.site) {
			lgam[i,j,t] <- c0[i,t]
			gam[i,j,t] <- 1/(1+exp(-lgam[i,j,t]))
			lphi[i,j,t] <- d0[i,t]
			phi[i,j,t] <- 1/(1+exp(-lphi[i,j,t]))
			psi[i,j,t+1] <- phi[i,j,t]*psi[i,j,t] + gam[i,j,t]*(1-psi[i,j,t])
			mu.z[i,j,t+1] <- w[i] * (phi[i,j,t]*Z[i,j,t] + gam[i,j,t]*(1-Z[i,j,t]))
			Z[i,j,t+1] ~ dbern(mu.z[i,j,t+1])
			mu.y[i,j,t+1] <- p[i,t+1]*Z[i,j,t+1]
			Y[i,j,t+1] ~ dbin(mu.y[i,j,t+1], K_tot[j,t+1])
		}
	}
}

# Derive total species richness for the metacommunity
N_tot <- sum(w[])

# Derive yearly species richness
for (i in 1:M) {
	for (t in 1:T) {
		tmp[i,t] <- sum(Z[i,,t])
		tmp2[i,t] <- ifelse(tmp[i,t]==0,0,1)
	}
}
for (t in 1:T) {
	N[t] <- sum(tmp2[,t])
}

# Derive WPI
for (i in 1:nspeciesWPI) {
	for (t in 1:T) {
		psi_ratio[i,t] <- psi[id_speciesWPI[i],1,t]/psi[id_speciesWPI[i],1,1] 
		log_psi_ratio[i,t] <- log(psi_ratio[i,t])
	}
}
for (t in 1:T) {
	WPI[t] <-  exp((1/(nspeciesWPI))*sum(log_psi_ratio[1:nspeciesWPI,t]))
}

# rate of change in WPI
for (t in 1:(T-1)) {
	lambda_WPI[t] <- WPI[t+1]/WPI[t]	
}

}
", fill=TRUE, file=modelFilename)


# number of traps
nsites <- dim(Yaug_tot)[2]

# number of primary periods (years)
T <- dim(Yaug_tot)[3]

# latent states
w <- c(rep(1,length(allnames3)),rep(NA,nzeros))

# Parameters monitored
parameters <- c("omega","psiMean","pMean","phiMean","gamMean",
                "lpsiSD","lpSD","lphiSD","lgamSD",
                "N","N_tot","WPI","lambda_WPI")

# data
bugs.data <- list(M=dim(Yaug_tot)[1],n.site=nsites,K_tot=K_tot,Y=Yaug_tot,T=T,
                  nspeciesWPI=nspeciesWPI,id_speciesWPI=id_speciesWPI)

# Initial values
inits <- function() { list(omega=runif(1),Z=(Yaug_tot>0)*1,w=w,
                           psiMean=runif(1),pMean=runif(T),phiMean=runif(T-1),gamMean=runif(T-1),
                           lpsiSD=runif(1,0,4),lpSD=runif(1,0,4),lphiSD=runif(1,0,4),lgamSD=runif(1,0,4))}

#mcmc settings
n.adapt <- 5000		#pre-burnin
n.update <- 10000	#burnin
n.iter <- 50000		#iterations post-burnin
thin <- 10
chains<-3

# run the model and record the run time
cl <- makeCluster(chains, type = "SOCK")
start.time = Sys.time()
out <- jags.parfit(cl, data = bugs.data, 
                  params = parameters, 
                  model = "dmsom.txt",
                  inits = inits,
                  n.adapt = n.adapt,	 
                  n.update = n.update,
                  n.iter = n.iter,
                  thin = thin, n.chains = chains)

end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes\n', sep=''))
stopCluster(cl)

#Posterior computed in 137.004419155916 minutes

# Summarize posteriors
summary(out)


# check output
xyplot(out[,c("omega","psiMean","lpsiSD","lpSD","lphiSD","lgamSD"), drop=F])
xyplot(out[,c("pMean[1]","pMean[2]","pMean[3]","pMean[4]","pMean[5]"), drop=F])
xyplot(out[,c("gamMean[1]","gamMean[2]","gamMean[3]","gamMean[4]"), drop=F])
xyplot(out[,c("phiMean[1]","phiMean[2]","phiMean[3]","phiMean[4]"), drop=F])


### plot year-specific species richness
# bind chains for the plot
out2 <- mcmc(do.call(rbind, out))

median_N <- lower_N <- upper_N <- NULL
for (i in 1:T){
	lab <- paste("N[",i,"]", sep="")
	median_N[i] <- quantile(out2[,lab], 0.500)
	lower_N[i] <- quantile(out2[,lab], 0.025)
	upper_N[i] <- quantile(out2[,lab], 0.975)
}

par(mar = c(5,4.5,4,1.5)+.1, cex.axis=1.5, cex.lab=2, tcl = 0.25)	# Note: many setting here... (tcl is for thick marks inside)
plot(x=1:T, y=median_N, type = "b", pch = 16, ylab = "Number of species", 
     xlab = "", bty = "n", ylim=c(20,35), 
     cex = 1.5, axes = FALSE, main="", lwd=2, cex.axis=1.5, cex.lab=1.5)
axis(1, las = 1, at = 1:5, labels = seq(2009,2013,1), cex = 1.5, lwd=2, cex.axis=1.5, cex.lab=1.5)
axis(2, cex = 1.5, lwd=2, cex.axis=1.5, cex.lab=1.5)
segments((1:T), lower_N, (1:T), upper_N, lwd=2)
points((1:T),dimcheck[1,])

### plot WPI
median_WPI <- lower_WPI <- upper_WPI <- NULL
for (i in 1:T){
	lab <- paste("WPI[",i,"]", sep="")
	median_WPI[i] <- quantile(out2[,lab], 0.500)
	lower_WPI[i] <- quantile(out2[,lab], 0.025)
	upper_WPI[i] <- quantile(out2[,lab], 0.975)
}

par(mar = c(5,4.5,4,1.5)+.1, cex.axis=1.5, cex.lab=2, tcl = 0.25)	# Note: many setting here... (tcl is for thick marks inside)
plot(0, 0, ylim=c(0,2), xlim = c(1,T), ylab = "WPI", xlab = " ",
col = "red", type = "l", lty = 5, lwd = 2, axes = F)
axis(2, las = 1)
axis(1, at = 1:T, labels = c("2009","2010","2011","2012","2013"))
polygon(x = c(1:T, T:1), y = c(lower_WPI, upper_WPI[T:1]), col = "grey90",
border = "grey90")
points(median_WPI,type = "l",lwd = 2.5)
abline(h=1,lty=2,col="red")
