########################################################
### Script to Generate Deep Canopy VIs for B+G Model ###
########################################################

#################
### 0. Setup  ###
#################

# Scripts to compute deep canopy VI (VI.inf) in Baret & Guyot using 2-Stream
setwd("~/Dropbox (BOSTON UNIVERSITY)/Main/Rwork/2Stream_Canopy")

# define constants used in script
ovp.time <- '15:10:00' # overpass time: midpoint btwn L8 and S2 (LT = GMT - 5)
rug.coef <- 1.25       # for deep canopy, assume uniform (baseline = 1.25)
L.deep <- 10           # deep canopy LAI
p.diff.red <- 0.15     # original = 0.15
p.diff.nir <- 0.10     # original = 0.10
tau.adjust.vis <- 0.02   # adjust tau for red based on fits with data
tau.adjust.nir <- -0.03  # adjust tau for vis based on fits with data

# note: for p.diff in red, nir I used 6S to estimate p.diff in each band 
# http://www-loa.univ-lille1.fr/documents/LOA/informatique/logiciels/Wsixs/
# Landsat8 bands, mid-latitude summer, continental atmosphere; tau=0.2 at 0.5 micrometers

# define functions
source('Rscripts/functionDefs.R')

# Ingest Sophie Dillen's leaf optics data
source("Rscripts/LeafDat.R")

#####################################################
### 1. Read in field data, and set up data frames ###
#####################################################

# Read data from Harvard Forest
hfdat <- read.csv('Data/HF-LAI-fPAR-Data/hf_fpar_lai_11_1.csv',header=T)[,c(-1)]

# Convert date to DOY, add DOY, year to dateframe, rename cals
doy <- as.integer(strftime(hfdat[,'Date'], format = "%j"))
year <- as.integer(strftime(hfdat[,'Date'], format = "%y"))+2000
hfdat <- data.frame(doy,year,hfdat)
colnames(hfdat) <- c('doy','year','date','glai','pai','fpar','red','nir')

# compute VIs
ndvi <- (hfdat[,'nir']-hfdat[,'red'])/(hfdat[,'nir']+hfdat[,'red'])
evi2 <- 2.5*(hfdat[,'nir']-hfdat[,'red'])/(hfdat[,'nir']+2.4*hfdat[,'red']+1)
nirv <- hfdat[,'nir']*ndvi

# compute SZAs at overpass time
solar.info <- solar(as.POSIXct(paste(hfdat[,'date'],ovp.time),tz="GMT"))
hf.mu <- cos(zenith(solar.info,-71.1,42.3)*pi/180)

# make data frame with Harvard Forest Data
hfdat <- data.frame(hfdat,ndvi,evi2,nirv,hf.mu)

#########################################################
### 2. Set up data to model deep canopy for DOY 1:365 ###
#########################################################

# first get leaf reflectances from Sophie Dillen data, interpolate to DOYs 1:365
dillen.days=c(1,210,244,272,282,294,365)
dillen.rhos=c(vis.med[1],vis.med,vis.med[5])
alpha.vis.fullyear <- approx(dillen.days,dillen.rhos,n=365)$y

# assume alpha and tau are equal 
tau.vis.fullyear <- alpha.vis.fullyear

# create vector of DOYs and cos(SZAs) at overpass time for each DOY
year.o.dates <- seq(as.Date("2019/1/1"), as.Date("2019/12/31"), "days")
doy <- as.integer(strftime(year.o.dates, format = "%j"))
solar.info <- solar(as.POSIXct(paste(year.o.dates,ovp.time),tz="GMT"))
mu.range<- cos(zenith(solar.info,-72.171478,42.537755)*pi/180)

#######################################################
### 3. Run two-stream for each DOY for LAI = L.deep ###
#######################################################

# initialize generic canopy parameters for Temperate Decid Forest
setCanopyParms(case='VIs')

# 2-stream for deep canopy for each day of year in VIS and NIR
ts.vis <- twostr_canopy(L=L.deep, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.vis.fullyear,
                        tau=tau.vis.fullyear + tau.adjust.vis,
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=L.deep, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir + tau.adjust.nir,
                        rhos=soil.rho.nir)

# compute K for SZA = 0 
K.sza0 <- twostr_canopy(L=L.deep, 
                        chiL=chi.canopy,
                        mu=1,
                        alpha=alpha.vis,
                        tau=tau.vis,
                        rhos=soil.rho.vis)$K

# compute K for SZA at overpass time for each DOY
K <- ts.nir$K

# estimate shade fraction; uniform canopy changes in illuminated LAI approximate shadow fraction
p.shaded.fullyear <- getPsunlit(lai.val=L.deep,
                               chiL=chi.canopy,
                               mu.val=mu.range,
                               rug.fac=0)$p.shade.tot

# compute VIs ignoring impact of shadows (don't pllot them)
VIs.noshadow.deep <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=mu.range,varname='Mu',plotVIs=FALSE)

# now do linear mixing for to adjust vis, nir reflectances for seasonally varying shadows
ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,
                             p.diff=p.diff.nir,
                             p.shaded=p.shaded.fullyear,
                             mu=mu.range)  
ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,
                             p.diff=p.diff.red,
                             p.shaded=p.shaded.fullyear,
                             mu=mu.range)

# create dataframe with modeled deep canopy EVI2, NDVI, & NIRv, including effects of shading
VIs.shadow.deep <- data.frame(2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1),
                         (ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj),
                         ts.nir.adj*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj))
colnames(VIs.shadow.deep)=c('evi2','ndvi','nirv')

# put everything in a single data frame
dailyVIinf <- data.frame(doy,
                         mu.range,
                         p.shaded.fullyear,
                         VIs.noshadow.deep$ndvi,
                         VIs.shadow.deep$ndvi,
                         VIs.noshadow.deep$evi2,
                         VIs.shadow.deep$evi2,
                         VIs.noshadow.deep$nirv,
                         VIs.shadow.deep$nirv)
colnames(dailyVIinf) <- c('DOY','cos(SZA)','Prop.Shaded','NDVI.NS','NDVI.SHDW','EVI2.NS','EVI2.SHDW','NIRv.NS','NIRv.SHDW')

# write results to CSV
write.csv(dailyVIinf,file='DailyDeepCanopyVIs.csv',row.names=FALSE)


