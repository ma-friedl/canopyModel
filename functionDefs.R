# Definitions for various functions used for canopy reflectance modeling based on 
# combination of two-stream and Baret and Guyot Models.

#source('~/Dropbox/Main/Courses/GE503/RTools/EE503-ToolBox/Two_Stream_Canopy.R')
library(GeoLight)

# function to read sophie dillen leaf data
getLeafDat <- function(path=NULL, bandmin=400, bandmax=700){
  setwd(path)
  fils <- list.files()
  nfils <- length(fils)
  leafdat <- NULL
  bands <- read.csv(fils[1])[,1]
  rng <- bands>bandmin & bands<bandmax
  for (fil in fils[(-nfils)]) {
    leafdat <- c(leafdat,read.csv(fil)[rng,2])
  }
  return(as.numeric(leafdat))
}

doVIplots <- function(vis,nir,dep.var,varname, p.diff.vis=0.15, p.diff.nir=0.02,plotVIs=FALSE) {
  vis.rho <- (1-p.diff.vis)*vis$Iupdir + p.diff.vis*vis$Iupdif
  nir.rho <- (1-p.diff.nir)*nir$Iupdir + p.diff.nir*nir$Iupdif
  evi2 <- 2.5*(nir.rho-vis.rho)/(nir.rho+2.4*vis.rho+1)
  ndvi <- (nir.rho-vis.rho)/(nir.rho+vis.rho)
  nirv <- nir.rho*ndvi
  if (plotVIs) {
    plot(dep.var,vis.rho,xlab=varname,ylab='Rho.VIS',type='l',col='blue')
    plot(dep.var,nir.rho,xlab=varname,ylab='Rho.NIR',type='l',col='blue')
    plot(dep.var,evi2,xlab=varname,ylab='EVI2',type='l',col='blue')
    plot(dep.var,ndvi,xlab=varname,ylab='NDVI',type='l',col='blue')
    plot(dep.var,nirv,xlab=varname,ylab='NIRv',type='l',col='blue')
  }
  return(data.frame(dep.var,vis.rho,nir.rho,ndvi,evi2,nirv))
}

fARPlot <- function(ts.mod,dep.var,varname) {
  plot(dep.var,ts.mod$fAR.dir,
       type='l',
       col='black',
       xlab=varname,
       ylab='fAR',
       ylim=c(0,1))
  lines(dep.var,ts.mod$fAR.dif,col='blue')
  lines(dep.var,ts.mod$fAR.net,col='green')
}

getK <- function(LAI,VI,VI.inf,VI.g) {
  K.vi <- -log((VI-VI.inf)/(VI.g-VI.inf))/LAI
  return(K.vi)
}

shadowMixModel <- function(twost.flx,p.diff=0.1,p.shaded=0.1,mu=0.1) {
  
  # CRUDE adjustment of p.diff for SZA
  p.diff <- p.diff/mu
  
  # partition downwelling fluxes into beam, diffuse
  Idwn.dif  <-  p.diff
  Idwn.dir <- (1-p.diff)
  
  # diffuse, beam albedoes
  rho.dif <- twost.flx$Iupdif
  rho.dir <- twost.flx$Iupdir
  
  # upwelling exitance in shaded area
  Iup.shaded <- rho.dif*p.diff 
  
  # upwelling exitance in sunlit area
  Iup.sunlit <- rho.dir*Idwn.dir + rho.dif*Idwn.dif
  
  # total upwelling
  Iup <- p.shaded*Iup.shaded + (1-p.shaded)*Iup.sunlit

  return(Iup)
}

fourScale <- function(fs.sl=0.25,fs.sh=0.25,fv.sl=0.25,fv.sh=0.25,
                      rhos.sl=0.25,rhos.sh=0.25,rhoc.sl=0.25,rhoc.sh=0.25) {
  rho.surf <- fs.sl*rhos.sl + 
    fs.sh*rhos.sh +
    fv.sl*rhoc.sl +
    fv.sh*rhoc.sh
  return(rho.surf)
}

setCanopyParms <- function(case='PAR') {
  # default SZA, LAI and CHiL
  mu.canopy <<- 0.853          # SZA = 30.5; ~mid-way btwn solstice and equinox
  lai.canopy <<- 3             # leaf area index
  
  # now set up optical properties
  if (case=='CLM') {
    
    chi.canopy <<- 0.25           # plagiophile for harvard forest
    
    # use default values for broadlead decid forest from CLM
    alpha.vis <<- 0.10           # leaf reflectance in vis from CLM5.0
    tau.vis <<- 0.05             # leaf transmittance in vis from CLM5.0
    alpha.nir <<- 0.45           # leaf reflectance in nir from CLM5.0
    tau.nir <<- 0.25             # leaf transmittance in nir from CLM
    
    soil.rho.vis <<- 0.20        # background reflectance, in vis from ~CLM5
    soil.rho.nir <<- 0.40        # background reflectance, in nir from ~CLM5
    
  } else if (case=='VIs') {
    
    chi.canopy <<- 0.5           # planophile for Harvard Forest
    
    # Set up narrow-band leaf alpha & tau based on Dillen measurements
    alpha.vis <<- 0.04            # median vis reflectance, July 28-30, 2010 from Dillen data 
    alpha.nir <<- 0.51            # median nir reflectance, July 28-30, 2010 from Dillen data
    
    tau.vis <<- alpha.vis         # assume leaf transmittance in vis same as reflectance
    tau.nir <<- 1-alpha.nir-0.1   # leaf transmittance in nir, assume absorption 0.1
    #soil.rho.vis <<- 0.063        # background reflectance in vis from HLS
    soil.rho.vis <<- 0.080
    soil.rho.nir <<- 0.210        # background reflectance in nir from HLS        
  } else {
    print('No Action: Prescribe Case: PAR or VIs')
  }
  
}


##########################################################################
#                                                                        #
#         Two Stream Model of Canopy Radiative Transfer                  #
#                                                                        #
#  Based on basic approach described by Dickinson 1983 and Sellers 1985, #
#     and as implemented by Bonan 1996 in the CLM Land Surface model     #
#                                                                        #
##########################################################################

twostr_canopy <- function(L=3, chiL=0.001, mu=0.5,alpha=0.1,tau=0.1,rhos=0.1,p.diff=0.1) {
  # L = LAI
  # chiL:  Ross (1975) parameterization for LAD
  # mu = cosine solar zenith angle
  # alpha = leaf reflectance coef
  # tau = leaf transmission coef
  # rhos = soil reflection coef
  # p.diff = proportion diffuse to total irradiance
  
  # Trap for chiL = 0
  chiL[chiL==0]=0.00001  
  
  # trap for chiL out of range
  if (sum(chiL > 0.59) + sum(chiL < -0.4) > 0) {
    print('chiL out of range')
    return(NULL)
  }
  
  # leaf scattering coefficient 
  omega <- alpha+tau
  
  # parameters for Goudriaan LAD function; from Sellers 1985; Bonan 1996 
  phi1 <- 0.5 - 0.633*chiL - 0.33*chiL^2 
  phi2 <- 0.877*(1-2*phi1) 
  
  # Ross-Goudriaan projected area of leaves in direction acos(mu)
  Gmu <- phi1+phi2*mu 
  
  # Optical depth per unit LAI
  K <- Gmu/mu 
  
  # single scattering albedo; from Bonan 1996, page 20
  asmu <- (omega/2)*(Gmu/(mu*phi2+Gmu))*(1-((mu*phi1)/(mu*phi2+Gmu))*log((mu*phi1+mu*phi2+Gmu)/(mu*phi1))) 
  
  # ave inverse diffuse optical depth/leaf area; from Bonan, 1996, page 19
  mubar <- (1/phi2)*(1-(phi1/phi2)*log((phi1+phi2)/phi1))
  
  # p-scatter coef for diffuse  radiation and beam radiation; from Bonan ,1996; page 20
  beta <- 0.5*(alpha+tau+(alpha-tau)*((1+chiL)/2)^2)/(alpha+tau)  
  beta0 <- (1+mubar*K)*asmu/(omega*mubar*K) 
  
  # and the rest of the terms required for solution
  b <- 1-omega+omega*beta
  c <- omega*beta 
  d <- omega*mubar*K*beta0 
  f <- omega*mubar*K*(1-beta0) 
  h <- (b^2-c^2)^0.5/mubar
  
  sigma <- (mubar*K)^2 + c^2 - b^2 
  
  u1 <- b-c/rhos 
  u2 <- b-c*rhos 
  u3 <- f+c*rhos 
  
  S1 <- exp(-h*L) 
  S2 <- exp(-K*L) 
  
  p1 <- b+mubar*h 
  p2 <- b-mubar*h 
  p3 <- b+mubar*K 
  p4 <- b-mubar*K 
  
  D1 <- p1*(u1-mubar*h)/S1 - p2*(u1+mubar*h)*S1 
  D2 <- (u2+mubar*h)/S1 - (u2-mubar*h)*S1 
  
  h1 <- -d*p4-c*f 
  h2 <- (1/D1)*((d-(h1/sigma)*p3)*(u1-mubar*h)/S1 - p2*(d-c-(h1/sigma)*(u1+mubar*K))*S2) 
  h3 <- -(1/D1)*((d-(h1/sigma)*p3)*(u1+mubar*h)*S1 - p1*(d-c-(h1/sigma)*(u1+mubar*K))*S2) 
  h4 <- -f*p3-c*d 
  h5 <- -(1/D2) * (h4*(u2+mubar*h)/(S1*sigma) + (u3-h4*(u2-mubar*K)/sigma)*S2) 
  h6 <- (1/D2) * (h4*(u2-mubar*h)*S1/sigma + (u3-h4*(u2-mubar*K)/sigma)*S2) 
  h7 <- (c/D1)*(u1-mubar*h)/S1 
  h8 <- (-c/D1)*(u1+mubar*h)*S1 
  h9 <- (1/D2)*(u2+mubar*h)/S1 
  h10 <- (-1/D2)*(u2-mubar*h)*S1 
  
  # compute fluxes
  Iupdir <- h1/sigma + h2+ h3 
  Iupdif <- h7 + h8 
  Idwndir <- (h4*exp(-K*L))/sigma + h5*S1 + h6/S1
  Idwndif <- h9*S1 + h10/S1
  
  # compute fraction absorbed radiation
  fAR.dir <- 1-Iupdir-(1-rhos)*(exp(-K*L) + Idwndir)
  fAR.dif <- 1-Iupdif-(1-rhos)*Idwndif
  fAR.net <- (1-p.diff)*fAR.dir + p.diff*fAR.dif
  
  return(data.frame(Iupdir,Iupdif,Idwndir,Idwndif,fAR.dir,fAR.dif,fAR.net,K))
}

getKval <- function(laival=3.0,chival=0.0001,muval=0.5,
                    alphav=0.1,alphan=0.4,tauv=0.1,taun=0.1,
                    rhosv=0.1,rhosn=0.2,
                    VIinfval=0.7,VIgval=0.1,
                    p.shade.range=-99) {

  #3## NOTE CURRENTLY SET UP TO DO ONLY EVI2!!!!!   #####
  ts.vis <- twostr_canopy(L=laival, 
                          chiL=chival, 
                          mu=muval,
                          alpha=alphav,
                          tau=tauv,
                          rhos=rhosv)
  
  ts.nir <- twostr_canopy(L=laival, 
                          chiL=chival, 
                          mu=muval,
                          alpha=alphan,
                          tau=taun,
                          rhos=rhosn)
  
  if (p.shade.range[1]==-99) {
    VIs <- doVIplots(vis=ts.vis,nir=ts.nir,dep.var=muval,varname='cos(sza)')[,'evi2']
  } else {
    ts.nir.adj <- shadowMixModel(twost.flx=ts.nir,p.diff=0.02,p.shaded=p.shade.range,mu=muval)
    ts.vis.adj <- shadowMixModel(twost.flx=ts.vis,p.diff=0.15,p.shaded=p.shade.range,mu=muval)
    VIs <- 2.5*(ts.nir.adj-ts.vis.adj)/(ts.nir.adj+2.4*ts.vis.adj+1)
  }
  
  # Extract K
  
  K.vi <- getK(LAI=laival,VI=VIs,VI.inf=VIinfval,VI.g=VIgval)
  return(K.vi)
}

# Function to get p.shaded, for given LAI, and SZA
getPshade <- function(laival=3.0,chival=0.0001,muval=0.5) {
  # designed for fixed LAI and range of cos(SZA)
  # can prescribe optics, no effect on result
  alphav <- 0.1
  tauv <- 0.1
  rhosv <- 0.1
  
  # get K for range of SZA
  K.ts <- twostr_canopy(L=laival,
                      chiL=chival, 
                      mu=muval,
                      alpha=alphav,
                      tau=tauv,
                      rhos=rhosv)$K
  
  # Get K for SZA = 0
  K.ts.mu0 <- twostr_canopy(L=laival,
                          chiL=chival, 
                          mu=1,
                          alpha=alphav,
                          tau=tauv,
                          rhos=rhosv)$K
  
  # compute sunlit LAI when SZA is 0
  L.snlt.sza0 <- (1-exp(-K.ts.mu0*laival))/K.ts.mu0
  
  # compute sunlit LAI for range of SZA
  L.mu <- (1-exp(-K.ts*laival))/K.ts
  
  # compute shade fraction
  p.shaded <- (L.snlt.sza0-L.mu)/L.snlt.sza0
  #p.shaded <- exp(-K.ts.mu0*laival) + (1-exp(-K.ts.mu0*laival))*p.shaded
  return(data.frame(p.shaded,L.mu,L.snlt.sza0,K.ts,K.ts.mu0))
  
}


 getPsunlit <- function(lai.val=1.0,chiL=0.000001, mu.val=1,rug.fac=1.0){
   # compute sunlight & shaded canopy fractions and background in FOV; For spherical LAD
   # See pg 222 in Jones and Vaughn and page 266 in Norman + Campbell
   # This version from Jones and Vaughn - seems to be wrong!

   # Trap for chiL = 0
   chiL[chiL==0]=0.000001

   # get K
   phi1 <- 0.5 - 0.633*chiL - 0.33*chiL^2
   phi2 <- 0.877*(1-2*phi1)

   # Ross-Goudriaan projected area of leaves in direction acos(mu)
   Gmu <- phi1+phi2*mu.val

   # Optical depth per unit LAI
   K <- Gmu/mu.val

   # Do K for SZA = 0 (nadir view); cos(0)=1, so:
   K.sza0 <- phi1+phi2

   # sunlit leaf area with sun directly overhead
   L.v.sza0 <- (1-exp(-K.sza0*lai.val))/K.sza0

   # sunlit leaf area for sun and sza given by mu
   L.v <- (1-exp(-K*lai.val))/K

   # fraction soil exposed
   f.soil <- exp(-K.sza0*lai.val)

   # percent canopy self-shadowed
   p.shade.self <- 1-(L.v/L.v.sza0)
   
   # p.shade estimates prop self-shading in illuminated canopy
   # add shading due to canopy roughness/rugosity - varies non-linearly w/SZA
   if(rug.fac==0) {
     p.shade.tot <- p.shade.self
   } else {
     p.shade.rug <- (1-mu.val)^rug.fac
     p.shade.tot <- p.shade.rug+(1-p.shade.rug)*p.shade.self
   }

   return(data.frame(mu.val,Gmu,f.soil,K,K.sza0,L.v.sza0,L.v,p.shade.tot))
 }


BG.lai <- function(vi=0.5,vi.g=0.2,vi.inf=0.8,k.vi=0.5) {
  # invert VI using Baret and Guyot to get LAI
  lai <- log((vi-vi.inf)/(vi.g-vi.inf))/-k.vi
  return(lai)
}

BG.fpar <- function(vi=0.5,vi.g=0.2,vi.inf=0.8,k.p=1,k.vi=0.5,fpar.inf=0.95) {
  # invert VI using Baret and Guyot to get fPAR
  fpar <- fpar.inf*(1-((vi.inf-vi)/(vi.inf-vi.g))^(k.p/k.vi))
  return(fpar)
}

G.fun <- function(x=1,sza=pi/3) {
  # compute G
  G <- sqrt((x^2)*(cos(sza)^2) + sin(sza)^2)/(x+1.774*(x+1.182)^(-0.733))
  return(G)
}

G.trans <- function(x=1,sza=pi/6,lai=3) {
  # compute transmission through canopy
  beta <- pi/2 - sza
  G <- G.fun(x,sza)
  par.tr <- exp(-G*lai/sin(beta))
  return(par.tr)
}

####################################################################