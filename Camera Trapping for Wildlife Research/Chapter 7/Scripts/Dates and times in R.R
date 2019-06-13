############################################################################
######                           Section 7.4.2                        ######
######                            for Chapter 7                       ######
######                                ***                             ######
######             Scripting by Danilo Foresti  (August 2015)         ######
############################################################################

as.POSIXlt("2015-03-29 01:30:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Berlin")
as.POSIXct("2015-03-29 01:30:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Berlin")
day<-"29.3.2015 00:00:00"
clocktime<-"30.12.1899 01:30:00"
strsplit(day, split=" ")
day<-strsplit(day, split=" ")[[1]][1]
clocktime<-strsplit(clocktime, split=" ")[[1]][2]
day; clocktime
x<-paste(day,clocktime)
x
xPOSIXlt<-as.POSIXlt(x, format="%d.%m.%Y %H:%M:%S", tz="Europe/Berlin")
xPOSIXlt
xPOSIXlt2<-as.POSIXlt("2015-03-29 03:30:00", tz="Europe/Berlin")
xPOSIXlt3<-as.POSIXlt("2015-03-29 03:30:00", tz="Australia/Melbourne")
xPOSIXlt; xPOSIXlt2; xPOSIXlt3
difftime(xPOSIXlt2,xPOSIXlt)
difftime(xPOSIXlt2,xPOSIXlt3)
xUTC<-format(as.POSIXct(xPOSIXlt),tz="UTC")
xUTC2<-format(as.POSIXct(xPOSIXlt2),tz="UTC")
xUTC;xUTC2
xUTClocal<-as.POSIXlt(xUTC,tz="UTC") + (60*60)
xUTClocal2<-as.POSIXlt(xUTC2,tz="UTC") + (60*60)
xUTClocal;xUTClocal2
