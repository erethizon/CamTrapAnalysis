############################################################################
######                            Chapter 7                           ######
######        Spatially explicit capture-recapture with secr          ######
######                               ***                              ######
######     Scripting by Danilo Foresti and Fridolin Zimmermann        ######
######                          (November 2015)                       ######
############################################################################

############################################################################
## Preparation
############################################################################

## Define working directory
setwd(choose.dir())

## Load secr package
library(secr)

############################################################################
## Read dataset & prepare for analysis
############################################################################

## Read Lynx_data.txt: All lynx photos of the session of interest:
lynx_data <- read.delim("Lynx_data.txt", header=T, stringsAsFactors=F)
lynx_data[,"Time"]<-as.POSIXct(lynx_data$Time) # Adopt a workable format for Time column 
head(lynx_data)
dim(lynx_data)

## Import Sites.txt table:
sites<-read.delim("Sites.txt",stringsAsFactors=F);head(sites);dim(sites)
sites<-data.frame(LOC_ID=1:dim(sites)[1],sites) # Add a column with a single numeric identifier for each site
head(sites)

## Import Trapnights.txt calendar:
trap_nights<-read.delim("Trapnights.txt",header=T, stringsAsFactors=F);head(trap_nights)

############################################################################
## Convert the juveniles' detections in mother's detections
## If they have to be considered, otherwise delete them before

for (i in 1:dim(lynx_data)[1]){
	if(lynx_data$Mother[i]!="" & lynx_data$Year_Born[i]==2013){	# If known mother and lynx is defined as juvenile for the specified year
		 lynx_data[i,"Lynx_Name"]<-lynx_data$Mother[i]		# Replace juvenile's name with the mother's name
		 lynx_data[i,"Sex"]<-"f"					# and insert "f" as the corresponding sex 
	}};rm(i)

## Delete column with descendance information, since no more consistent
lynx_data<-lynx_data[,-which(names(lynx_data)=="Mother")]
lynx_data<-lynx_data[,-which(names(lynx_data)=="Year_Born")]
head(lynx_data)

############################################################################
## Define period:
 
## Start and stop of the sampling period:
start_date <- as.POSIXlt("2013-11-29 12:00:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Berlin")
end_date <- as.POSIXlt("2014-01-28 12:00:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Berlin")
## Length of sampling occasion IN DAYS:
length_occasion<-5
## Calculate number of sampling occasions in the defined sampling period 
max_occasions<-as.numeric(difftime(end_date,start_date,unit='days')/length_occasion)
max_occasions

############################################################################
## Compute the trap nights (method is fraction):

traps_table<-matrix(nrow=dim(sites)[1], ncol=max_occasions)
for (i in 1:dim(traps_table)[1]){
	for (j in 1:max_occasions){
		traps_table[i,j]<-sum(trap_nights[i,(2+length_occasion*(j-1)):(1+length_occasion*j)])/length_occasion
		}
	};rm(i,j)
head(traps_table)

############################################################################
## Calculate capture events for all lynx detections

captures<-data.frame(SESSION=integer(), ANIMAL_ID=integer(),SO=integer(), LOC_ID=integer(), SEX=integer())
for (i in 1:dim(lynx_data)[1]){                                               # For every detection...
	captures[i,"SESSION"]<-1								# ...assign a constant value for session name, then...
	captures[i,"LOC_ID"]<-sites[sites[,"Site"]==lynx_data$Site[i],"LOC_ID"] # ...find the LOC_ID...
	captures[i,"ANIMAL_ID"]<-lynx_data$Lynx_Name[i]                         # ...the lynx ID...
	captures[i,"SO"]<-floor(difftime(lynx_data$Time[i],start_date,units='days')/length_occasion)+1	# ...and assign detections to a specific sampling occasion
	captures[i,"SEX"]<-lynx_data$Sex[i]							# Assign the sex of the animal
	if (captures[i,"SEX"]=="") captures[i,"SEX"]<-NA				# If sex not known specify NA
	};rm(i)
head(captures)
dim(captures)[1]

## Remove repeated observations from captures (from "count" to "proximity" dataset)
captures_prox<-unique(captures)
rownames(captures_prox)<-1:length(rownames(captures_prox))

head(captures_prox)
dim(captures)[1] 		# Amount of detections before computation
dim(captures_prox)[1] 	# Amount of detections after compression
head(captures_prox)
dim(captures_prox)[1]


############################################################################
## Create SECR 'proximity' objects:

## Create SECR traps "proximity" object
traps_analysis_nousage<-read.traps(data=sites[,c(1,3:4)], detector="proximity", binary.usage=FALSE)
summary(traps_analysis_nousage)
## Add usage info
traps_analysis<-traps_analysis_nousage
usage(traps_analysis)<-traps_table
summary(traps_analysis)

## Create SECR capthist "proximity" object
lynx_capthist<-make.capthist(captures_prox, traps_analysis, fmt="trapID", noccasions=max_occasions, covnames="SEX")
summary(lynx_capthist)

############################################################################
## Analysis
############################################################################

## Create masks for mask inspection
mask_1 <-make.mask(traps_analysis, buffer=1000, type="traprec", spacing=1000)
mask_3 <-make.mask(traps_analysis, buffer=3000, type="traprec", spacing=1000)
mask_5 <-make.mask(traps_analysis, buffer=5000, type="traprec", spacing=1000)
mask_7 <-make.mask(traps_analysis, buffer=7000, type="traprec", spacing=1000)
mask_9 <-make.mask(traps_analysis, buffer=9000, type="traprec", spacing=1000)
mask_11 <-make.mask(traps_analysis, buffer=11000, type="traprec", spacing=1000)
mask_13 <-make.mask(traps_analysis, buffer=13000, type="traprec", spacing=1000)
mask_15 <-make.mask(traps_analysis, buffer=15000, type="traprec", spacing=1000)
mask_17 <-make.mask(traps_analysis, buffer=17000, type="traprec", spacing=1000)
mask_20 <-make.mask(traps_analysis, buffer=20000, type="traprec", spacing=1000)
mask_25 <-make.mask(traps_analysis, buffer=25000, type="traprec", spacing=1000)
mask_30 <-make.mask(traps_analysis, buffer=30000, type="traprec", spacing=1000)

## Run SECR analysis for different masks - find correct buffer
M0_1 <-secr.fit(lynx_capthist, mask=mask_1, model = list(D~1, g0~1, sigma~1))
M0_3 <-secr.fit(lynx_capthist, mask=mask_3, model = list(D~1, g0~1, sigma~1))
M0_5 <-secr.fit(lynx_capthist, mask=mask_5, model = list(D~1, g0~1, sigma~1))
M0_7 <-secr.fit(lynx_capthist, mask=mask_7, model = list(D~1, g0~1, sigma~1))
M0_9 <-secr.fit(lynx_capthist, mask=mask_9, model = list(D~1, g0~1, sigma~1))
M0_11 <-secr.fit(lynx_capthist, mask=mask_11, model = list(D~1, g0~1, sigma~1))
M0_13 <-secr.fit(lynx_capthist, mask=mask_13, model = list(D~1, g0~1, sigma~1))
M0_15 <-secr.fit(lynx_capthist, mask=mask_15, model = list(D~1, g0~1, sigma~1))
M0_17 <-secr.fit(lynx_capthist, mask=mask_17, model = list(D~1, g0~1, sigma~1))
M0_20 <-secr.fit(lynx_capthist, mask=mask_20, model = list(D~1, g0~1, sigma~1))
M0_25 <-secr.fit(lynx_capthist, mask=mask_25, model = list(D~1, g0~1, sigma~1))
M0_30 <-secr.fit(lynx_capthist, mask=mask_30, model = list(D~1, g0~1, sigma~1))

densities<-c(1.347605e-04,1.203294e-04,1.111139e-04,1.061034e-04,1.039526e-04,1.032479e-04,1.030708e-04,1.030367e-04,1.030318e-04,1.030314e-04,1.030314e-04,1.030314e-04)
SE<-c(2.850448e-05,2.550388e-05,2.363025e-05,2.266154e-05,2.228379e-05,2.217866e-05,2.215876e-05,2.215678e-05,2.215694e-05,2.215712e-05,2.215713e-05,2.215714e-05)

library (Hmisc)
errbar(c(1,3,5,7,9,11,13,15,17,20,25,30),densities*10000,(densities+SE)*10000,(densities-SE)*10000,type='b', 
	col='black',pch=19,xlab='Buffer size (km)',errbar.col='black',ylim=c(0,max(densities+SE)*10000),
	ylab='Fitted (real) parameter for D (individuals/100 km2)')

## Choose good buffer: 13km
mask_prox<-mask_13
summary(mask_prox)

########################################################################################

## Run several predefined models (without mixture on sex)

M0 <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~1, sigma~1))
Mb <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~b, sigma~1))
Mt <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~t, sigma~1))
Mbt <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~b+t, sigma~1))
MT <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~T, sigma~1))
MB <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~B, sigma~1))
Mbk <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~bk, sigma~1))
MBk <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~Bk, sigma~1))
MK <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~K, sigma~1))
MK_NM <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~K, sigma~1),method="Nelder-Mead")
Mk <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~k, sigma~1))
Mh2 <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~h2, sigma~1))
Mh2_NM <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~h2, sigma~1),method="Nelder-Mead")
Mh2_NM_SL <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~h2, sigma~1),method="Nelder-Mead",start=list(D=0.00012,g0=0.24,sigma=3200))
Mh2_NM_SM0 <-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~h2, sigma~1),method="Nelder-Mead",start=M0)


# AIC analysis
AIC(M0,Mb,Mt,Mbt,MT,MB,Mbk,MBk,MK_NM,Mk,Mh2_NM_SM0, criterion="AIC")
# Best model is Mbk

########################################################################################

# Compare spatial versus non-spatial models with SECR
# region.N to transform SECR model into Npop estimations
region.N(Mbk, region=mask_prox) 

########################################################################################

# Analysis of sex differences 
M_g0bk_sigma<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~bk, sigma~1), hcov="SEX")
M_g0_sigmasex<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~1, sigma~h2), hcov="SEX")
M_g0sex_sigmasex<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~h2, sigma~h2), hcov="SEX")
M_g0bk_sigmasex<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~bk, sigma~h2), hcov="SEX")
M_g0bksex_sigma<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~(bk+h2), sigma~1), hcov="SEX")
M_g0bksex_sigmasex<-secr.fit(lynx_capthist, mask=mask_prox, model = list(D~1, g0~(bk+h2), sigma~h2), hcov="SEX")

# AIC analysis
AIC(M_g0_sigmasex,M_g0sex_sigmasex,M_g0bk_sigma,M_g0bk_sigmasex,M_g0bksex_sigma,M_g0bksex_sigmasex, criterion="AIC")
# There is no outstanding model
# We thus apply model averaging to the two competing models (delta AIC < 2; Burnham and Anderson 2002) to get unbiased parameter estimates
model.average(M_g0bk_sigmasex,M_g0bksex_sigmasex,criterion="AIC")
########################################################################################

## Including the habitat mask

## Write a text file with the mask_prox for further analysis in GIS
write.table(mask_prox,'mask_prox_GIS.txt',sep='\t',row.names=F)

## Load the library maptools
library (maptools)

## Create the polygon object in R
HabitatPol <- readShapePoly("suitable_habitat_secr")

## Create the habitat mask
habitat_mask <- make.mask(traps_analysis, spacing=1000, type="polygon", poly = HabitatPol)
summary(habitat_mask)

## Draw a map to see if the sites and the shapefile are aligned to the same area
plot(HabitatPol, axes=F)
plot(habitat_mask, pch=20, add=T)
points(traps_analysis[,"x"], traps_analysis[,"y"], pch=20, bg=grey(1), col='black')

## Run an analysis with a habitat mask on the previous models (with and without sex)
M_g0bk_sigma_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~bk, sigma~1), hcov="SEX")
M_g0sex_sigmasex_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~h2, sigma~h2), hcov="SEX")
M_g0_sigmasex_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~1, sigma~h2), hcov="SEX")
M_g0bk_sigmasex_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~bk, sigma~h2), hcov="SEX")
M_g0bksex_sigma_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~(bk+h2), sigma~1), hcov="SEX")
M_g0bksex_sigmasex_habitat<-secr.fit(lynx_capthist, mask=habitat_mask, model = list(D~1, g0~(bk+h2), sigma~h2), hcov="SEX")
# AIC analysis
AIC(M_g0bk_sigma_habitat, M_g0sex_sigmasex_habitat, M_g0_sigmasex_habitat, M_g0bk_sigmasex_habitat, M_g0bksex_sigma_habitat, M_g0bksex_sigmasex_habitat, criterion="AIC")
# There is no outstanding model as in the previous case without the habitat mask
# We thus apply model averaging to the two competing models (delta AIC < 2; Burnham and Anderson 2002) to get unbiased parameter estimates
model.average(M_g0bk_sigmasex_habitat,M_g0bksex_sigmasex_habitat,criterion="AIC")
save.image("Chapter_7.RData")
