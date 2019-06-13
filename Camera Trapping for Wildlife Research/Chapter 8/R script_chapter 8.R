####################################################################
# Chapter 8 - Comparison of activity patterns
# 
# Adapted from Meredith and Ridout (2014)
#
# Scripting Fridolin Zimmermann and Danilo Foresti
#
# August 2015
####################################################################
# To set the working directory type the following command line
# This will open an explorer window for easy browsing
setwd(choose.dir())
# To get the current working directory and verify that it is correct
getwd()
# Import the .csv table into R when the data are separated by ","
events <- read.delim(file.choose(),header=TRUE,sep=',')
# Choose the file containing the data (true_events_NWA2013_14.csv)
# Display the header of the table
# to verify that the file has been correctly imported
# This will recall you the exact name of the columns in your dataset
# This is crucial to successfully run the functions presented hereafter
head(events)
# Before the analyses, it is good to explore the dataset
# To see the different study areas the data come from type
table(events$area)
# To see how many species and events per species are available type
summary(events$species)
# To check if the time format is correct (0-1 range) type
range(events$time)
# Package overlap works entirely in radian units
# The conversion is straightforward
timeRad<-events$time*2*pi

#################################################################### 
# Fitting the kernel density
# Load the package overlap (Meredith and Ridout 2014)
# Extract the data for the Eurasian lynx, plot a kernel density curve
library(overlap)
lynx<-timeRad[events$area==1 & events$species== 'Lynx lynx']
densityPlot(lynx, rug=T)
densityPlot(lynx, rug=T, adjust=2)
densityPlot(lynx, rug=T, adjust=0.5)

####################################################################
# Coefficient of overlap
# Practical example with the north-western Swiss Alps dataset
# Extract the data for the Eurasian lynx and its prey the roe deer
lynx<-timeRad[events$area==1 & events$species=='Lynx lynx']
roe<-timeRad[events$area==1 & events$species=='Capreolus capreolus']
# Get the size of the smaller of the two samples type
min(length(lynx), length(roe))
# Calculate the overlap with the three estimators
lynxroeest<-overlapEst(lynx,roe)
lynxroeest
# Plot the curves
overlapPlot(lynx,roe, rug=T)
legend('topright', c("Eurasian lynx", "Roe deer"), lty=c(1,2), col=c(1,4), bty='n')

####################################################################
# Bootstrap analysis to estimate the CI of the coefficient of overlapping
# Generate 10000 smoothed bootstrap for Eurasian lynx and roe deer
lynxboot<-resample(lynx,10000)
roeboot<-resample(roe, 10000)
dim(lynxboot)
dim(roeboot)
# To generate estimates of overlap from each pair of samples
# these two matrices are passed to the function bootEst()
lynxroe<-bootEst(lynxboot, roeboot, adjust=c(NA,1,NA))
dim(lynxroe)
BSmean<-colMeans(lynxroe)
BSmean
tmp<-lynxroe[,2]
bootCI(lynxroeest[2],tmp)
tmp<-lynxroe[,2]
bootCIlogit(lynxroeest[2],tmp)

####################################################################
# Overlap for lynx and chamois including estimation of the CI
# estimated following the same procedure
chamois<-timeRad[events$area==1 & events$species=='Rupicapra rupicapra']
# Get the size of the smaller of the two samples type
min(length(lynx), length(chamois))
lynxchamoisest<-overlapEst(lynx,chamois)
lynxchamoisest
overlapPlot(lynx,chamois, rug=T)
legend('topright', c("Eurasian lynx", "Chamois"), lty=c(1,2), col=c(1,4), bty='n')
# Bootstrap analysis to estimate the CI of the coefficient of overlapping
chamoisboot<-resample(chamois, 10000)
lynxchamois<-bootEst(lynxboot, chamoisboot, adjust=c(NA,1,NA))
BSMean<-colMeans(lynxchamois)
BSMean
tmp<-lynxchamois[,2]
bootCI(lynxchamoisest[2],tmp)
bootCIlogit(lynxchamoisest[2],tmp)

####################################################################
# Function compareCkern () of package activity (Rowcliffe 2015) to test
# that two sets of circular observations come from the same distribution

library(activity)
compareCkern(lynx, roe, reps = 10000)
compareCkern(lynx, chamois, reps = 10000)
compareCkern(roe, chamois, reps = 10000)

####################################################################
# Wald test using the function compareAct(fits) to estimate
# the significance of pairwise comparisons between overall activity levels
library(activity)
# A circular kernel density on the original dataset is fitted 
# The coverage of the confidence interval seems to be better estimated
# with sample = "model" when the sample size is greater than 100-200, 
# whereas smaller sample sizes should be investigated using sample = "data"
# The bootstrapping on a model basis for lynx and roe deer was chosen 
# given their sample size (> 100)
lynxactmod<-fitact(lynx,adj=1, sample="model", reps=10000)
lynxactmod
roeactmod<-fitact(roe,adj=1, sample="model", reps=10000)
roeactmod
compareAct(list(lynxactmod,roeactmod))

# Same analysis for the chamois
# As the sample size for chamois is less than 100, 
# the argument sample of the function was set to data

chamoisactmod<-fitact(chamois,adj=1, sample="data", reps=10000)
chamoisactmod
compareAct(list(lynxactmod,chamoisactmod))
compareAct(list(roeactmod,chamoisactmod))
