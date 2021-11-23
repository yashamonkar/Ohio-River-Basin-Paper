#________________________________________________________________________________#
###Code for PC-Wavelet Plots
#To get the PC-Wavelet Plots for the annual maximum data





#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")

#Load Dependencies
library(dplyr)
library(NbClust)
library(biwavelet)
library(plotrix)
library(maps)
library(cowplot)



#________________________________________________________________________________#
####Reading the Data####

#Read the Annual Maximum Data
input_data <- read.table("data/Annual_Maximum_Streamflow.txt", 
                         sep="", header = TRUE)

#Read the Site Information
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)


#Data Cleaning
Years <- input_data$Year
input_data$Year <- NULL


#________________________________________________________________________________#
####Functions to compute the PC-Wavelets

###Input
#1. Data Frame of the Annual Maximum values(Dat)
#2. Number of PCs to be analyzed (npcs)
#3. Lat-Lon coordinates of the Sites (coord)
#4. Years (yrs)

Dat <- input_data
npcs <- 3
site_loc <- data.frame(Lat = site_info$dec_lat_va, 
                       Lon = site_info$dec_long_va)
yrs <- Years

#________________________________________________________________________________#

#####PCwav Function#####
PCwav=function(x,npcs,coords,nam,yr1){
  
  #Load Dependencies
  require("biwavelet")
  require("plotrix")
  require("maps")
  
  #Load the Map Parameters
  world <- map_data("world")
  us <- map_data("state")

  #Set up data parameters
  Dat=as.matrix(Dat)
  nr=nrow(Dat)
  nc=ncol(Dat)
  cx=cor(Dat)

  #Principal Component Analysis(PCA)
  #supplying correlation matrix to the princomp so that nr<nc case does not lead to a problem
  pcs=prcomp(Dat, scale = TRUE)
  var <- cumsum(pcs$sdev^2)
  
  pc_var <- data.frame(NP = 1:length(var),
                            Var = var/max(var))
  
  pvar <- ggplot(pc_var) +
    geom_point(aes(x = NP, y = Var)) +
    geom_line(aes(x = NP, y = Var)) +
    geom_vline(xintercept = 3, linetype = 'dashed') +
    scale_x_continuous(name = "Number of PCs") +
    annotate("text", x = 12, y = var[3]/max(var), 
             label = paste0("Variance explained by 3 PCs - ",
                            100*round(var[3]/max(var),2),"%")) +
    scale_y_continuous(name = "Variance Explained") +
    labs(title = "Variance explained by the PCs") + 
    theme_bw() +
    theme(plot.title = element_text(size=15),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10)) 
  print(pvar)
  
  pdf("Trial.pdf")
  #--------------------------------------------------------------------------------#
  #PC-Wavelet Analysis
  for(i in 1:npcs){
    
    #Get the PC
    pc=pcs$x[,i]
    
    #--------------------------------------------------------------------------------#
    #Plot the PC-Score
    plt_dataset <- data.frame(PC = pc,
                              Year = yrs,
                              Loess = lowess(yrs,pc,f=1/9)$y)
    p1 <- ggplot(plt_dataset) +
      geom_line(aes(x = Year, y = Loess), size = 1.2, color ='red') +
      geom_point(aes(x = Year, y = PC), size = 0.1) +
      geom_line(aes(x = Year, y = PC), size = 0.1) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "PC-Score") +
      labs(title = paste0("PC-",i," Score")) + 
      theme_bw() +
      theme(plot.title = element_text(size=12),
            axis.text=element_text(size=5),
            axis.title=element_text(size=10)) 
    
    #--------------------------------------------------------------------------------#
    #Plot the PC Eigenvectors
    eigen_mean <- mean(pcs$rotation[,i])
    eigen_min <- min(pcs$rotation[,i])
    eigen_max <- max(pcs$rotation[,i])
    
    
    p2 <- ggplot() +
      geom_map(data=us, map=us,
               aes(x=long, y=lat, map_id=region),
               fill="#D3D3D3", color="#000000", size=0.15) +
      scale_x_continuous(name = "lon", limits = c(-91, -78)) +
      scale_y_continuous(name = "lat", limits = c(36.5, 43)) +
      geom_point(data = site_loc, aes(x= Lon, y = Lat, 
                                     color = pcs$rotation[,i])) +
      scale_color_gradient2(low="blue", high="red") +
      labs(title = paste0("PC-",i," Eigenvectors")) + 
      labs(color="Eigenvectors")  +
      theme_bw() +
      theme(legend.text=element_text(size=7),
            legend.title=element_text(size=5),
            axis.text=element_text(size=0),
            axis.title=element_text(size=0),
            axis.ticks = element_blank(),
            plot.title = element_text(size=12),
            legend.key.height  = unit(0.75, "cm"))
    
    
    #Wavelet Analysis 
    wlt=wavelet(pc)
    Cw=CI(0.9,pc,"w")
    C=CI(0.9,pc,"r")
    
    #--------------------------------------------------------------------------------#
    #Global Wavelet Spectrum
    plt_dataset <- data.frame(Period = wlt$period,
                              Power = wlt$p.avg,
                              W_noise = Cw$sig,
                              R_noise = C$sig)
    
    p3 <- ggplot(plt_dataset) +
      geom_point(aes(x = Power, y = Period)) +
      geom_path(aes(x = Power, y = Period)) +
      geom_line(aes(x = W_noise, y = Period)) +
      geom_line(aes(x = R_noise, y = Period), color ='red') +
      ggtitle("Global Wavelet Spectrum") + 
      scale_x_continuous(name = "Variance") +
      scale_y_continuous(name = "Period (Years)",
                         limits = c(0, 30)) +
      theme_bw() +
      theme(legend.text=element_text(size=15),
            legend.title=element_text(size=0),
            axis.text=element_text(size=10),
            axis.title=element_text(size=10),
            plot.title = element_text(size=12))
    
    #Wavelet Power Spectrum
    wt1=wt(cbind(Years,pc))
    x <- wt1
    
    
    print(plot_grid(p1,p2,p3,
              nrow = 2,
              labels = c("A", "B", "C"),
              label_size = 12))
    
  }
  dev.off()  
    
    #--------------------------------------------------------------------------------#
    #Plotting the Wavelet Power Spectrum - Continuous Wavelet Transform
    wt1=wt(cbind(Years,pc))
    
    
    
    
    nam1=paste("PC",as.character(i),nam)
    nam2=paste("Wavelet -PC",as.character(i))
    par(mfrow=c(2,2));
    par(mar=c(4,2,3,0.5))
    plot(yr1:yr2,p,xlab="Year",ylab=nam,main=nam1)
    lines(lowess(yr1:yr2,p,f=1/9),lwd=2,col="red")
    #above plots PC time series and lowess through it
    plot(wlt$period,wlt$p.avg,xlim=c(0,32),main="Global Wavelet Spectrum",xlab="Period",ylab="Variance"); lines(wlt$period,wlt$p.avg);
    lines(wlt$period,Cw$sig)
    lines(wlt$period,C$sig,col="red")
    #above plots PC loadings
    wt1=wt(cbind(yr1:yr2,p));plot(wt1, type="power.corr.norm", xlab="Year",main=nam2)
    #above plots Global wavelet and conf limits
    #par(mfrow=c(1,1))
    ju=abs(pcs$rotation[,i])
    par(mar=c(0,0,0,0))
    map('state', region = c("Ohio","Indiana", "Illinois","West Virginia","Kentucky","Pennsylvania","Virginia"))
    points(coords$lon,coords$lat,pch=19,cex=1,col=color.scale(ju,c(1,0.5,0),c(0,0.5,0),c(0,0,1),color.spec="rgb"))
    #points(dams_dataset$Longitude,dams_dataset$Latitude,pch=9,cex=1)
    
    
  }
  #above requires wavelets package from R and plots regular wavelets plot 
  par(mfrow=c(1,1))
}

### Wavelet Function

#WAVELET  1D Wavelet transform with optional singificance testing
#
#   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
#
# INPUTS:
#
#    Y = the time series of length N.
#    DT = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)^2.
#    Its units are sigma^2 (the time series variance).
#
#
# OPTIONAL INPUTS:
# 
# *** Note *** setting any of the following to -1 will cause the default
#               value to be used.
#
#    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    DJ = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
#        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
#
#    MOTHER = the mother wavelet function.
#             The choices are 'MORLET', 'PAUL', or 'DOG'
#
#    PARAM = the mother wavelet parameter.
#            For 'MORLET' this is k0 (wavenumber), default is 6.
#            For 'PAUL' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#
# OPTIONAL OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
#            where J1+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot COI lines on a contour plot by doing:
#
#              contour(time,log(period),log(power))
#              plot(time,log(coi),'k')
#
#----------------------------------------------------------------------------
#   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made. This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------
wavelet=function(Y,dj=0.025){
  
  #Y is time series to be analyzed
  DT=1# is timestep for annual data, 1
  pad=1
  #dj=0.025
  param=6
  #pad data ? 0=F, 1=T
  #dj= spacing between discrete scales (.025)
  #param = wavenumber (6)
  
  s0=2*DT
  
  n1 = length(Y)
  J1=floor((log2(n1*DT/s0))/dj)
  
  
  #....construct time series to analyze, pad if necessary
  x = Y - mean(Y)
  
  
  if (pad == 1){
    base2 = trunc(log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
    x = c(x, rep(0, 2^(base2 + 1) - n1))
  }
  n = length(x)
  
  #....construct wavenumber array used in transform [Eqn(5)]
  k = (1:trunc(n/2))
  k = k*((2*pi)/(n*DT))
  k = c(0, k, -rev(k[1:floor((n-1)/2)]))
  
  #....compute FFT of the (padded) time series
  f = fft(x)    # [Eqn(3)]
  
  #....construct SCALE array & empty PERIOD & WAVE arrays
  scale = s0*2^((0:J1)*dj)
  period = scale;
  wave = matrix(data=0, ncol=n, nrow=J1+1)  # define the wavelet array
  wave = as.complex(wave)  # make it complex
  wave=matrix(data=wave, ncol=n, nrow=J1+1)
  
  # loop through all scales and compute transform
  for(a1 in 1:(J1+1)){
    scl=scale[a1]		
    
    nn = length(k);
    k0 = param
    expnt = -(scl*k - k0)^2/(2*(k > 0))
    norm = sqrt(scl*k[2])*(pi^(-0.25))*sqrt(nn)    # total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt)
    daughter = daughter*(k > 0)    # Heaviside step function
    fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)) # Scale-->Fourier [Sec.3h]
    coi = fourier_factor/sqrt(2)                  # Cone-of-influence [Sec.3g]
    dofmin = 2                                   # Degrees of freedom
    
    out <- list(daughter=daughter, fourier_factor=fourier_factor,coi=coi,dofmin=dofmin)
    
    daughter=out$daughter
    fourier_factor=out$fourier_factor
    coi=out$coi
    dofmin=out$dofmin	
    wave[a1,] = fft((f*daughter), inverse = TRUE)/(length(f*daughter))  # wavelet transform[Eqn(4)]
  }
  
  period = fourier_factor*scale
  
  coi = coi*c(1:(floor(n1 + 1)/2), rev(1:floor(n1/2))) * DT
  
  wave = wave[,1:n1]  # get rid of padding before returning
  power=abs(wave)^2
  ncol=length(power[1,])
  nrow=length(scale)
  avg.power=apply(power,1,mean)
  result=list(wave=wave, period=period, scale=scale, power=power, coi=coi,nc=ncol,nr=nrow,p.avg=avg.power)
  return(result)
}



### Confidence level function

CI=function(conf, dat,type){
  
  #### enter confidence as decimal 0-1
  #### two types of tests available 1) red noise enter: "r" , white noise enter: "w"
  # requires the wavelet function
  
  na=length(dat)
  wlt=wavelet(dat)
  
  if(type=="r"){
    
    zz=arima(dat/10^10, order = c(1, 0, 0))
    alpha=zz$coef[1]
    print(alpha)
  } else{
    alpha=0
  }
  
  ps=wlt$period
  LP= length(ps)
  
  freq = 1/ps
  
  CI=1:LP    ## confidence interval
  
  for(i in 1:LP){
    
    P=(1-(alpha^2))/(1+(alpha^2)-(2*alpha*cos(2*pi*freq[i])))    # descrete fourier power spectrum page 69 [qn 16] ( torrence and compo)... null hypothesis test
    df=2*sqrt(1+((na/(2.32*ps[i]))^2))
    CI[i] =P*(qchisq(conf, df)/df)*var(dat)          #divide by 2 removes power of 2.....for mean no chi dist.[ eqn 17]
  }
  
  
  list(sig=CI)
  
}


####Calling PC-Wav####
x <- input_data
npcs <- 3
coords <- as.data.frame(cbind(site_info$dec_lat_va,site_info$dec_long_va))
colnames(coords) <- c("lat","lon")
nam <- c("Ann_Max")
yr1 <- 1934

#pdf(file = 'plots/PcWav Plots.pdf')
PCwav(x=x, npcs = npcs, coords = coords, nam = nam, yr1=yr1 )
#dev.off()