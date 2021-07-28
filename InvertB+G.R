#############################################
#       Script to invert B+G Model          #
#############################################

# set up data sets and deep canopy VI's
source('Rscripts/DeepCanopy.R')
library(mgcv)

# define growing season
grow.season <- 121:304     # define growing season as May 1 - Nov 1

# fudge factor!
dailyVIinf$EVI2.SHDW <- dailyVIinf$EVI2.SHDW+0.2

# linear adjustment of glai for woody fraction
W.frac <- 1-(hfdat[,'pai']-min(hfdat[,'pai'],na.rm=T))/max(hfdat[,'pai'],na.rm=T)
hfdat[,'glai'] <- hfdat[,'pai']*(1-W.frac)

# create daily time series of LAI using LAI data from all years vs DOY
lai.loess <- loess(hfdat[,'glai']~hfdat[,'doy'],span=0.2)           # smooth time series across years
lai.days <- c(1,min(lai.loess$x)-1,lai.loess$x,max(lai.loess$x)+1,365)  # LAI = 0 outside of growing season
lai.smo <- c(0,0,lai.loess$fitted,0,0)
lai.ave <- approx(lai.days,lai.smo,n=365)$y                         # interpolate to daily; ignore warning

# HF data to exclude data outside of GS
hfdat <- subset(hfdat,hfdat[,'doy'] %in% grow.season)

# no do 2-stream for canopy using 'lai.ave'
ts.vis <- twostr_canopy(L=lai.ave, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.vis.fullyear,
                        tau=tau.vis.fullyear,
                        rhos=soil.rho.vis)

ts.nir <- twostr_canopy(L=lai.ave, 
                        chiL=chi.canopy, 
                        mu=mu.range,
                        alpha=alpha.nir,
                        tau=tau.nir,
                        rhos=soil.rho.nir)

# Get shade fraction: first, compute K for SZA = 0 
K.sza0 <- twostr_canopy(L=lai.ave, 
                        chiL=chi.canopy,
                        mu=1,
                        alpha=alpha.vis,
                        tau=tau.vis,
                        rhos=soil.rho.vis)$K

# estimate shade fraction; changes in illuminated LAI approximate shadow fraction in sensor FOV
p.shaded.fullyear <- getPsunlit(lai.val=lai.ave,
                                chiL=chi.canopy, 
                                mu.val=mu.range,
                                rug.fac=rug.coef)$p.shade.tot

# compute VIs ignoring impact of shadows
VIs.noshadow <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=mu.range,varname='Mu')

# now do linear mixing for to account for shadow area
ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=p.diff.nir,p.shaded=p.shaded.fullyear,mu=mu.range)
ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=p.diff.red,p.shaded=p.shaded.fullyear,mu=mu.range)

# and put in a single data frame
VIs.shadow.ave <- data.frame(2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1),
                             (ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj),
                             ts.nir.adj*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+ts.vis.adj))
colnames(VIs.shadow.ave)=c('evi2','ndvi','nirv')

# plot modeled and observered to compare
plot(hfdat[,'doy'],hfdat[,'evi2'],
     ylab='DOY',xlab='EVI2',
     pch=16,ylim=c(0.2,0.8),col='green')
lines(grow.season,VIs.shadow.ave$evi2[grow.season],
      col='lightgreen',lwd=3)
lines(grow.season,dailyVIinf$EVI2.SHDW[grow.season],
      col='darkgreen',lwd=3)

# clean up EVI data
# evi.doy.rws <- rownames(hfdat[!is.na(hfdat[,'evi2']),])   # rows where HLS data are available
# evi.doys <- hfdat[evi.doy.rws,'doy']             # doys when HLS are available
# high.evi.indx <- hfdat[evi.doy.rws,'evi2']>dailyVIinf$EVI2.SHDW[evi.doys] # indx: T when EVI2 > deep canopy for that DOY
# 
# # now loop through and replace/adjust extreme values
# for (i in 1:length(evi.doys)) {
#   # first test for anomalously high (> deep canopy!)
#   if (hfdat[evi.doy.rws[i],'evi2']>(dailyVIinf$EVI2.SHDW[evi.doys[i]]-0.02)) { 
#     # replace with value just greater than than Deep Canopy
#     # hfdat[evi.doy.rws[i],'evi2'] <- dailyVIinf$EVI2.SHDW[evi.doys[i]]-0.005
#     hfdat[evi.doy.rws[i],'evi2'] <- NA
#   # now screen for anomaloustly low (cloud)
#   } else if (VIs.shadow.ave$evi2[evi.doys[i]] - hfdat[evi.doy.rws[i],'evi2']>0.075) {
#     # cloudy - get rid of!
#     hfdat[evi.doy.rws[i],'evi2'] <- NA
#   }
# }


# now estimate K from B+G  - start by setting up basic inputs for inversion: 
VIg <- VIs.shadow.ave$evi2[110]                   # Set VI based on background May 1         
VIinf <- dailyVIinf$EVI2.SHDW[grow.season]        # VI.in varies with SZA and leaf optics May 1 - Nov 1
lai.gs <- lai.ave[grow.season]                    # representative lai for growing season
VI.doy=VIs.shadow.ave$evi2[grow.season]

# estimate k for each day in growing season by inverting B+G
k.doy <- (-log((VI.doy-VIinf)/(VIg-VIinf))/lai.gs)

# add DOY as names to k.doy and VI.inf 
names(k.doy) <- grow.season
names(VIinf) <- grow.season

# Set up index/list of days with HLS data at HF
evi.doy.indx <- !is.na(hfdat[,'evi2'])
hf.evi.doys <- as.character(hfdat[evi.doy.indx,'doy'])  # list of days, with HLS data

# and invert lai from HLS
lai.invert.hf <- log((hfdat[evi.doy.indx,'evi2']-VIinf[hf.evi.doys])/(VIg-VIinf[hf.evi.doys]))/-(k.doy[hf.evi.doys])

# plot results
plot(hfdat[,'doy'],hfdat[,'glai'],pch=16,ylim=c(0,6),ylab='LAI',xlab='DOY')
points(hfdat[evi.doy.indx,'doy'],lai.invert.hf,pch=16,col='red')

# smooth evi, re-invert lai, and plot
evidoy.df <- na.omit(data.frame(hfdat[evi.doy.indx,'evi2'],hfdat[evi.doy.indx,'doy']))
colnames(evidoy.df)=c('evi','doy')
evi.gam <- gam(evi~s(doy),data=evidoy.df)
lai.invert.gam <- log((predict(evi.gam)-VIinf[hf.evi.doys])/(VIg-VIinf[hf.evi.doys]))/-(k.doy[hf.evi.doys])
points(hfdat[evi.doy.indx,'doy'],lai.invert.gam,pch=16,col='blue')


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



