# Misc exploratory scripts - Mark needs to go through and clean up.

# Prelims
setwd("/Users/Mark/Dropbox (BOSTON UNIVERSITY)/Main/Rwork/2Stream_Canopy")
source('Rscripts/functionDefs.R')
source("Rscripts/LeafDat.R")

# set overpass time
ovp.time <- '15:15:00'

# First, read data
hfdat <- read.csv('Data/HF-LAI-fPAR-Data/hf_fpar_lai_11_1.csv',header=T)[,c(-1)]

# Convert date to DOY, add to dateframe, rename cals
doy <- as.integer(strftime(hfdat[,'Date'], format = "%j"))
hfdat <- data.frame(doy,hfdat)
colnames(hfdat) <- c('doy','date','glai','pai','fpar','red','nir')

# compute VIs, add to dataframe
ndvi <- (hfdat[,'nir']-hfdat[,'red'])/(hfdat[,'nir']+hfdat[,'red'])
evi2 <- 2.5*(hfdat[,'nir']-hfdat[,'red'])/(hfdat[,'nir']+2.4*hfdat[,'red']+1)
nirv <- hfdat[,'nir']*ndvi

# compute SZAs at 10:05 local time (compromise L8overpass~9:50, S2~10:20 )
solar.info <- solar(as.POSIXct(paste(hfdat[,'date'],ovp.time),tz="GMT"))
hfmu <- cos(zenith(solar.info,-71.1,42.3)*pi/180)

# make data frame
hfdat <- data.frame(hfdat,ndvi,evi2,nirv,hfmu)

# remove outlier HLS reflectances
indx <- !(hfdat[,'red']>0.08 | hfdat[,'nir']>0.43)
indx[is.na(indx)] <- TRUE
hfdat <- hfdat[indx,]

# get range of dates with valid green LAI values
valid.lai <- !is.na(hfdat[,'glai'])
doy.min <- min(hfdat[valid.lai,'doy'])-1
doy.max <- max(hfdat[valid.lai,'doy'])+1

# estimate leaf-off fpar (i.e., for LAI < 0.1), and replace all cases with fpar < base w/base value
glai.indx=hfdat[,'glai'] < 0.1
base.fpar <- median(hfdat[glai.indx,'fpar'],na.rm=T)
base.indx <- hfdat[,'fpar']<base.fpar & !is.na(hfdat[,'fpar'])
hfdat[base.indx,'fpar'] <- base.fpar

# get red, nir pre-leaf emergence
red.min <- median(subset(hfdat[,'red'],hfdat[,'doy']>100 & hfdat[,'doy']<120),na.rm=T)
nir.min <- median(subset(hfdat[,'nir'],hfdat[,'doy']>100 & hfdat[,'doy']<120),na.rm=T)

# now subset growing season data
hfdat.gs <- subset(hfdat,hfdat[,'doy'] > doy.min & hfdat[,'doy'] < doy.max)
apar.gs <- 1-exp(-0.5*hfdat.gs[,'glai']/hfdat.gs[,'hfmu'])

# rescale LAI and fPAR to 0-1 for plotting to see better comparison of timing
min.lai <- min(hfdat.gs[,'glai'],na.rm=T)
max.lai <- max(hfdat.gs[,'glai'],na.rm=T)
lai.rs <- (hfdat.gs[,'glai']-min.lai)/(max.lai-min.lai)
min.fpar <- min(hfdat.gs[,'fpar'],na.rm=T)
max.fpar <- max(hfdat.gs[,'fpar'],na.rm=T)
fpar.rs <- (hfdat.gs[,'fpar']-min.fpar)/(max.fpar-min.fpar)
min.apar <- min(apar.gs,na.rm=T)
max.apar <- max(apar.gs,na.rm=T)
apar.rs <- (apar.gs-min.apar)/(max.apar-min.apar)

# filter outlier fPAR values
fpar.rs[fpar.rs > 0.1 & hfdat.gs[,'doy']<120] = 0

# now plot
plot(hfdat.gs[,'doy'],
     fpar.rs,
     pch=16,
     col='black',
     xlab='DOY',
     ylab='Normalized LAI/fpar')
     
points(hfdat.gs[,'doy'],
     lai.rs,
     pch=16, col='red')

points(hfdat.gs[,'doy'],apar.rs,pch=16,col='green')

#######################################
# ok - let's compare against 2-Stream
#######################################

# First set up optics - run 'LeafDat.R
setCanopyParms(case='VIs')

# now get dates for which LAI != NA
good.dat <- subset(hfdat.gs,!is.na(hfdat.gs[,'glai']))
lai.range <- good.dat[,'glai']
hfmu.range <- good.dat[,'hfmu']

# run 2-stream 
ts.vis <- twostr_canopy(L=lai.range, 
                        chiL=chi.canopy, 
                        mu=hfmu.range,
                        alpha=alpha.vis,
                        tau=tau.vis+ 0.02, #        +0.02 to tune to HLS
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=lai.range, 
                        chiL=chi.canopy, 
                        mu=hfmu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)   # -0.1 to Tune to HLS

# compute VIs
VIs <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=lai.range,varname='LAI')

# now plot modeled vs observed reflectances and EVI2
plot(hfdat[,'doy'],hfdat[,'red'],xlab='doy',ylab='Red', pch=16,col='red',ylim=c(0.01,0.08))
points(good.dat[,'doy'],VIs$vis.rho,pch=16,col='black')
plot(hfdat[,'doy'],hfdat[,'nir'],xlab='DOY',ylab='NIR',pch=16,col='blue',ylim=c(0.12,0.45))
points(good.dat[,'doy'],VIs$nir.rho,pch=16,col='black')
plot(hfdat[,'doy'],hfdat[,'evi2'],xlab='DOY',ylab='EVI2',pch=16,col='green',ylim=c(0.15,0.7))
points(good.dat[,'doy'],VIs$evi2,pch=16,col='black')

######################################
# linear mixture of shaded vs sunlit #
######################################

# # compute sunlit vs shaded for entire range
# # first get SZA at Landsat overpass, Jun 21-Oct 1 (172-274)
# jun21 <- solar(as.POSIXct("2017-06-21 15:00:00",tz="GMT"))
# mu.jun21 <- cos(zenith(jun21,-71.1,42.3)*pi/180)
# oct1 <- solar(as.POSIXct("2017-10-01 15:00:00","GMT"))
# mu.oct1 <- cos(zenith(oct1,-71.1,42.3)*pi/180)
# 
# #set max LAI
# lai.max <- max(hfdat.gs[,'glai'],na.rm=T)
# 
# # get K
# K.jun21 <- twostr_canopy(L=lai.max, 
#                    chiL=chi.canopy, 
#                    mu=mu.jun21,
#                    alpha=alpha.vis,
#                    tau=tau.vis,
#                    rhos=soil.rho.vis)$K
# 
# K.oct1 <- twostr_canopy(L=lai.max, 
#                          chiL=chi.canopy, 
#                          mu=mu.oct1,
#                          alpha=alpha.vis,
#                          tau=tau.vis,
#                          rhos=soil.rho.vis)$K
# 
# K.sza0 <- twostr_canopy(L=lai.max, 
#                        chiL=chi.canopy, 
#                        mu=1,
#                        alpha=alpha.vis,
#                        tau=tau.vis,
#                        rhos=soil.rho.vis)$K
# 
# 
# L.sun.sza0 <- (1-exp(-K.sza0*lai.max))/K.sza0
# L.sun.jun21 <- (1-exp(-K.jun21*lai.max))/K.jun21
# L.sun.oct1 <- (1-exp(-K.oct1*lai.max))/K.oct1
# 
# # now assume changes in L.sun only approximate changes in sunlit at top of canopy in sensor FOV
# # and assume L.sun.0 = fully illuminated (no shade)
# p.shaded.jun21 <- (L.sun.sza0-L.sun.jun21)/L.sun.sza0
# p.shaded.oct1 <- (L.sun.sza0-L.sun.oct1)/L.sun.sza0
# 
# # set up range of variation from solstice to equinox
# p.shade.range <- seq(p.shaded.jun21,p.shaded.oct1,length.out=2)
# 
# # set up SZA range
# mu.range <- c(mu.jun21,mu.oct1)
# 
# # get canopy refleflances: NOTE USING FIXED SZA FOR NOW!!!!
# ts.vis <- twostr_canopy(L=lai.max, 
#                         chiL=chi.canopy, 
#                         mu=mu.range,
#                         alpha=alpha.vis,
#                         tau=tau.vis,
#                         rhos=soil.rho.vis)
# 
# ts.nir <- twostr_canopy(L=lai.max, 
#                         chiL=chi.canopy, 
#                         mu=mu.range,
#                         alpha=alpha.nir,
#                         tau=tau.nir,
#                         rhos=soil.rho.nir)
# 
# # now compute adjusted reflectances - note need to adjustment to proportion diffuse in red vs NIR
# ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=0.02,p.shaded=p.shade.range,mu=mu.range)
# ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=0.15,p.shaded=p.shade.range,mu=mu.range)
# 
# ts.noshadow <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=mu.range,varname='Mu')
# ts.noshadow

# peek at range 
# range(2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1))

#############################################################################################
# Simulate growing season, incoporating effects of (1) shadows and (2) changing leaf optics #
#############################################################################################

# lai and cos(sza) for all days with valid glai values
lai.range <- good.dat[,'glai']
mu.range <- good.dat[,'hfmu']

# include variation in vis rho over season; assume nir constant
alpha.leaf.vis=matrix(0,365)    # rho red for every DOY
alpha.leaf.vis[210] <- 0.0399   # values from Dillen measurements
alpha.leaf.vis[244] <- 0.0469
alpha.leaf.vis[272] <- 0.0510
alpha.leaf.vis[282] <- 0.0575
alpha.leaf.vis[293] <- 0.0925

# set up vector of DOYs and red reflectance
dillen.days=c(1:209,210,244,272,282,293,303:365)
dillen.rhos=c(rep(0.0399,209),vis.med,rep(0.1,length(303:365)))

# interpolate to daily values
dillen.rhos=spline(x=dillen.days,y=dillen.rhos,method='nat',n=365)$y
plot(dillen.rhos,ylim=c(0,0.13))
abline(v=c(210,244,272,282,293))
abline(h=c(0.04,0.047,0.051,0.0575,0.0925))

# subset vis alpha and tau from daily interpolated values.
alpha.good.dat <- dillen.rhos[good.dat[,'doy']]
tau.good.dat <- alpha.good.dat

# get canopy refleflances using time varying LAI, SZA, and alpha.vis interpolated from Dillen daat 
ts.vis <- twostr_canopy(L=lai.range, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.good.dat,
                        tau=alpha.good.dat,      #tau=alpha.good.dat+0.02,
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=lai.range, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)

# compute K SZA = 0 and K for SZA at overpass times
K.sza0 <- twostr_canopy(L=lai.range, 
                        chiL=chi.canopy,
                        mu=1,
                        alpha=alpha.vis,
                        tau=tau.vis,
                        rhos=soil.rho.vis)$K

K.mu <- twostr_canopy(L=lai.range,
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)$K

# Now compute proportion sunlit LAI for SZA=0 and for SZA at overpass time
L.sun <- (1-exp(-K.mu*lai.range))/K.mu
L.sun.sza0 <- (1-exp(-K.sza0*lai.range))/K.sza0

# compute fraction shaded based on ratio of sunlit LAI for SZA at overpass vs SZA = 0
p.shaded.gs <- (L.sun.sza0-L.sun)/L.sun.sza0
p.shaded.gs[is.na(p.shaded.gs)] <- 0      # replace NAs for LAI = 0 with p.shaded = 0
p.shaded.gs <- p.shaded.gs#*1.25/mu.range

# check on this
plot(good.dat[,'doy'],p.shaded.gs,xlab='DOY',ylab="Proportion Shade",pch=16)

# now do linear mixing for shadow area
ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=0.075,p.shaded=p.shaded.gs,mu=mu.range)
ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=0.1,p.shaded=p.shaded.gs,mu=mu.range)

# and compute EVI w/shadow and no shadow
shadow.evi2 <- 2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1)
shadow.ndvi<- (ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj)

ts.noshadow <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=mu.range,varname='Mu')

# plot results for EVI
plot(good.dat[,'doy'],ts.noshadow$evi2,pch=16,col='gray',xlab='DOY',ylab='EVI2',ylim=c(0.1,0.8))
points(good.dat[,'doy'],shadow.evi2,pch=16,col='red')
points(hfdat[,'doy'],hfdat[,'evi2'],pch=16,col='green')

# plot results for individual bands
plot(hfdat[,'doy'],hfdat[,'red'],xlab='doy',ylab='Red', pch=16,col='red')
points(good.dat[,'doy'],ts.vis.adj,pch=16,col='black')
points(good.dat[,'doy'],ts.vis$Iupdir,col='green')
points(good.dat[,'doy'],ts.vis$Iupdif,col='blue')

plot(hfdat[,'doy'],hfdat[,'nir'],xlab='DOY',ylab='NIR',pch=16,col='blue')
points(good.dat[,'doy'],ts.nir.adj,pch=16,col='black')
points(good.dat[,'doy'],ts.nir$Iupdir,col='green')
points(good.dat[,'doy'],ts.nir$Iupdif,col='blue')

# splined results
doy.order <- order(good.dat[,'doy'])
good.doys <- good.dat[doy.order,'doy']
start.doy <- min(good.dat[,'doy'],na.rm=T)
end.doy <- max(good.dat[,'doy'],na.rm=T)

evi2.ns.sp <- loess.smooth(x=good.doys,y=ts.noshadow$evi2[doy.order],
                           span=0.15,
                           evaluation=length(start.doy:end.doy))$y
evi2.shd.sp <- loess.smooth(x=good.doys,y=shadow.evi2[doy.order],
                            span=0.15,
                            evaluation=length(start.doy:end.doy))$y

plot(start.doy:end.doy,evi2.ns.sp,type='l',ylim=c(0.1,0.8),xlab='DOY',ylab='EVI2')
lines(start.doy:end.doy,evi2.shd.sp,col='blue',lwd=3)
points(good.doys,shadow.evi2[doy.order],col='blue',pch=16,cex=0.5)

doy.order <- order(hfdat[,'doy'])
good.doys <- hfdat[doy.order,'doy']
start.doy <- min(hfdat[,'doy'],na.rm=T)
end.doy <- max(hfdat[,'doy'],na.rm=T)
evi2.hf <- loess.smooth(x=good.doys,y=hfdat[doy.order,'evi2'],
                            span=0.075,
                            evaluation=length(start.doy:end.doy))$y
lines(start.doy:end.doy,evi2.hf,col='green',lwd=3)
points(good.doys,hfdat[doy.order,'evi2'],pch=16,col='green',cex=0.5)



###############################################################################
# Do deep canopy for VI.inf
year.o.dates <- seq(as.Date("2019/1/1"), as.Date("2019/12/31"), "days")
doy <- as.integer(strftime(year.o.dates, format = "%j"))
solar.info <- solar(as.POSIXct(paste(year.o.dates,'15:30:00'),tz="GMT"))
mu.range<- cos(zenith(solar.info,-71.1,42.3)*pi/180)
alpha.vis.fullyear <- dillen.rhos
tau.vis.fullyear <- dillen.rhos

ts.vis <- twostr_canopy(L=25, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.vis.fullyear,
                        tau=tau.vis.fullyear,
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=25, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)

# now compute adjusted reflectances - note need to adjustment to proportion diffuse in red vs NIR
K.sza0 <- twostr_canopy(L=25, 
                        chiL=chi.canopy,
                        mu=1,
                        alpha=alpha.vis,
                        tau=tau.vis,
                        rhos=soil.rho.vis)$K
K <- ts.nir$K
L.sun <- (1-exp(-K*25))/ts.nir$K
L.sun.sza0 <- (1-exp(-K.sza0*25))/K.sza0

# now assume changes in L.sun only approximate changes in sunlit at top of canopy in sensor FOV
# and assume L.sun.0 = fully illuminated (no shade)
p.shaded.fullyear <- (L.sun.sza0-L.sun)/L.sun.sza0

# plot p.shaded
plot(p.shaded.fullyear,xlab='DOY',ylab="Proportion Shade",pch=16)


# now do linear mixing for shadow area
ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=0.02,p.shaded=p.shaded.fullyear,mu=mu.range)
ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=0.15,p.shaded=p.shaded.fullyear,mu=mu.range)

ts.noshadow <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=mu.range,varname='Mu')

shadow.evi2 <- 2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1)
shadow.ndvi <- (ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj)
shadow.nirv <- ts.nir.adj*shadow.ndvi

# 


dailyVIinf <- data.frame(doy,
                         mu.range,
                         p.shaded.fullyear,
                         ts.noshadow$ndvi,
                         shadow.ndvi,
                         ts.noshadow$evi2,
                         shadow.evi2,
                         ts.noshadow$nirv,
                         shadow.nirv)

colnames(dailyVIinf) <- c('DOY','cos(SZA)','Prop.Shaded','NDVI.NS','NDVI.SHDW','EVI2.NS','EVI2.SHDW','NIRv.NS','NIRv.SHDW')

plot(1:365,dailyVIinf[,2],type='l',pch=16,col='black',lwd=2,xlab='DOY',ylab='VI',ylim=c(0,1))#xlim=c(100,330),ylim=c(0,1))
points(1:365,dailyVIinf[,3],pch=16,cex=0.7)
lines(1:365,dailyVIinf[,4],col='green',lwd=2)
points(1:365,dailyVIinf[,5],pch=16,cex=0.7,col='green')
lines(1:365,dailyVIinf[,6],col='red',lwd=2)
points(1:365,dailyVIinf[,7],pch=16,cex=0.7,col='red')
points(1:365, p.shaded.fullyear)


###### quick check on letiticia results

lai.doys <- good.dat[,'doy']
evi.lai.doys <- evi2.hf[lai.doys]
VIinf.lai.days <-dailyVIinf[lai.doys,'EVI2.SHDW']
# make data frame
all.dat <- data.frame(lai.doys,good.dat[,'glai'],evi.lai.doys,VIinf.lai.days,evi2.hf[100])
colnames(all.dat) <- c('DOY','LAI','EVI2','VI.inf','VI.g')

#all.dat <- subset(all.dat,all.dat[,'LAI']>1.0)
all.dat[,'VI.inf'] <- all.dat[,'VI.inf']

# K.bg <- getK(LAI=all.dat[,'LAI'],VI=all.dat[,'EVI2'],VI.inf=all.dat[,'VI.inf'],VI.g = all.dat[,'VI.g'])
# K.mn <- median(K.bg)
# K.mn
# K.ts <- ts.nir$K[lai.doys]

# compute K.vi from 2-stream - use LAI = 3
ts.vis <- twostr_canopy(L=3, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.vis.fullyear,
                        tau=tau.vis.fullyear,
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=3, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)

ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=0.02,p.shaded=p.shaded.fullyear,mu=mu.range)
ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=0.15,p.shaded=p.shaded.fullyear,mu=mu.range)
shadow.evi2 <- 2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1)
#K.ts <- getK(LAI=3,VI=shadow.evi2,VI.inf=dailyVIinf[,'EVI2.SHDW'],VI.g = evi2.hf[100])
K.ts <- K.vi[lai.doys]  # from B+GModel.R

lai.bg <- -log((all.dat[,'EVI2']-all.dat[,'VI.inf'])/(all.dat[,'VI.g']-all.dat[,'VI.inf']))/K.ts
plot(all.dat[,'LAI'],lai.bg)
abline(0,1)


################################33



p <- ggplot(data=hfdat,
            mapping=aes(x=doy,y=glai),
            color=year,fill=year)
p+geom_point(mapping=aes(color=as.factor(year))) +
        geom_smooth(method='gam')

p <- ggplot(data=hfdat,
            mapping=aes(x=doy,y=glai))
p+geom_point(mapping=aes(color=as.factor(year))) +
        geom_smooth(method='gam') + facet_wrap(~year)

p <- ggplot(data=hfdat,
            mapping=aes(x=doy,y=evi2))
p+geom_point(mapping=aes(color=as.factor(year))) +
        geom_smooth(method='gam') + facet_wrap(~year)

lai.inv <- rep(NA,dim(hfdat)[1])
test.df <- data.frame(hfdat,lai.inv)

head(test.df)
indx=!(is.na(hfdat[,'evi2']))
test.df[indx,'lai.inv'] <- predict(evi.gam)

p <- ggplot(data=test.df,
            mapping=aes(x=doy,y=lai.inv),
            color=year,fill=year)
p+geom_point(mapping=aes(color=as.factor(year))) +
        geom_smooth(method='gam')

p <- ggplot(data=test.df,
            mapping=aes(x=doy,y=lai.inv))
p+geom_point(mapping=aes(color=as.factor(year))) +
        geom_smooth(method='gam') + facet_wrap(~year)


rug.fac=1.25
p.shade.self <- 1-(tmp$L.v/tmp$L.v.sza0)
p.shade.rug <- (1-mu.val^2)^rug.fac
p.shade.tot <- p.shade.rug+(1-p.shade.rug)*p.shade.self


#####
# create data frame for plotting w/ggplot
# hf.df <- data.frame(hfdat[,c('doy','lai')],rep('Measured',dim(hfdat)[1]))
# inverte.df
# colnames(final.df) <- c('doy','LAI.m','LAI.inv','LAI.smo')
# final.df[evi.doy.indx,'LAI.inv'] <- lai.invert.hf
# final.df[evi.doy.indx,'LAI.smo'] <- lai.invert.gam
# 
# p <- ggplot(data=final.df,
#             mapping=aes(x=doy,y=lai.inv),
#             color=year,fill=year)
# p+geom_point(mapping=aes(color=as.factor(year))) +
#   geom_smooth(method='gam')
# 

# Framework for computing interception by stems vs leaves in spring vs fall.
# needs to be cleaned up!

X=matrix(rnorm(600),nrow=6)
y=matrix(rnorm(100),ncol=1)
X%*%y
Z=matrix(rnorm(36),ncol=6)

k=0.5
S=1.5
L=seq(0,5,0.1)

# during leaf emergence, leaves progressively intercept/absorb 
# more radiation, which decreases relative role of stems
fpar.l=1-exp(-k*L)
fpar.s=(1-exp(-k*S))*(1-fpar.l)
fpar.t=fpar.l+fpar.s

plot(L,fpar.t,ylim=c(0,1))
lines(L,fpar.l)
lines(L,fpar.s)

# to separate role of canopy remove effect of canopy subtract fpar.s from measured fpar




# in fall, leaves color, and decrease absorption, 
# fraction intercepted PAR does not change, but fpar does
L=c(rep(5,30),seq(4.7,0,-0.3)) # 30 days of leaf coloration; 16 days of leaf drop
# larger k = stronger extinction - set k to decrease linearly w/time during senescence
k.s=seq(0.5,0.3,length.out=30)    # k for leaves during senescence
k.ld=rep(0.3,16)
k.fall=c(k.s,k.ld)
fpar.l=1-exp(-k.fall*L)
plot(1:46,fpar.l)

# total fpar is then fpar.l + fpar.s
fpar.s=(1-exp(-0.3*S))*(1-fpar.l)
fpar.t=fpar.l+fpar.s










