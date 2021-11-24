#________________________________________________________________________________#
###Code for Clust-PC-Wavelets
#To get the Clust-PC-Wavelets Plots for the annual maximum data





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
library(ggplot2)

#Load Functions
source('functions/Get_Power_Spectrum.R')
source('functions/Get_Wavelets.R')


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
####Functions_Needed###########

#Function 1- Clust-PC-Wavelet Analysis
###Input
#1. Data Frame of the Annual Maximum values(Dat)
#2. Lat-Lon coordinates of the Sites (coord)
#3. Years (yrs)

###Output (TBD)
#1. For each selected Cluster
#1(a) - PC-1 Score
#1(b) - PC-1 Eigenvectors
#1(c) - PC-1 Global Wavelet Spectrum
#1(d) - PC-1 Wavelet Power Spectrum


wavclust=function(Dat,coord,yrs){
  
  #Load Dependencies
  require("biwavelet")
  require("plotrix")
  require("maps")
  require("stats")
  require('ggplot2')
  
  #Load the Map Parameters
  world <- map_data("world")
  us <- map_data("state")
  
  #Set up data parameters
  nr=nrow(Dat)
  nc=ncol(Dat)
  
  #Setting up storage for Wavelet Frequencies
  np=length(wavelet(Dat[,1])$p.avg) #Number of wavelet frequencies
  w.arr=array(NA,dim=c(nc,np)) #Place holder global wavelet for each series
  x=Dat #Standardize
  for(i in 1:nc) {
    x[,i]=(Dat[,i]-mean(Dat[,i]))/sd(Dat[,i])
    } 
  for(i in 1:nc) {
    w.arr[i,]=wavelet(x[,i])$p.avg
  }
  
  #Plotting Hierarchical Clustering
  par(mfrow=c(1,1))
  clus=hclust(dist(w.arr),method="ward.D2")
  plot(clus)
  (cls <- identify(clus))
  nclus=length(cls)
  dev.off()
  
  #Setting up the PDF
  pdf("figures/Waveclust_Analysis.pdf")
  
  #Plotting CLuster Membership
  ju=rep(NA,nc)
  for (i in 1:nclus)ju[cls[[i]]]=rep(i,length(cls[[i]]))
  coord$Cluster <- as.factor(ju)
  
  pclust <- ggplot() +
    geom_map(data=us, map=us,
             aes(x=long, y=lat, map_id=region),
             fill="#D3D3D3", color="#000000", size=0.15) +
    scale_x_continuous(name = "lon", limits = c(-91, -78)) +
    scale_y_continuous(name = "lat", limits = c(36.5, 43)) +
    geom_point(data = coord, aes(x= Lon, y = Lat, color = Cluster),
               size = 1.5) +
    labs(title = paste0("Cluster Membership")) + 
    labs(color="Clusters")  +
    theme_bw() +
    theme(legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          axis.text=element_text(size=0),
          axis.title=element_text(size=0),
          axis.ticks = element_blank(),
          plot.title = element_text(size=15),
          legend.key.height  = unit(1.75, "cm"))
  
  pclust
  
  #PC-Wavelet Analysis on each Cluster
  for (i in 1:nclus){
    
    #Select the current cluster members
    tpd=Dat[,cls[[i]]]
    
    #PCA to get PC1
    pcw=prcomp(tpd, scale = TRUE)
    pc1=pcw$x[,1]
    
    #--------------------------------------------------------------------------------#
    #Plot the PC-Score
    plt_dataset <- data.frame(PC = pc1,
                              Year = yrs,
                              Loess = lowess(yrs,pc1,f=1/9)$y)
    p1 <- ggplot(plt_dataset) +
      geom_line(aes(x = Year, y = Loess), size = 1.2, color ='red') +
      geom_point(aes(x = Year, y = PC), size = 0.1) +
      geom_line(aes(x = Year, y = PC), size = 0.1) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "PC-Score") +
      labs(title = paste0("Cluster-",i," PC-1 Score")) + 
      theme_bw() +
      theme(plot.title = element_text(size=12),
            axis.text=element_text(size=5),
            axis.title=element_text(size=10)) 
    
    #--------------------------------------------------------------------------------#
    #Plot the Cluster Memberships and Eigenvectors
    
    current_members <- coord %>% filter(Cluster == i)
    eigen_mean <- mean(pcw$rotation[,i])
    eigen_min <- min(pcw$rotation[,i])
    eigen_max <- max(pcw$rotation[,i])
    
    
    p2 <- ggplot() +
      geom_map(data=us, map=us,
               aes(x=long, y=lat, map_id=region),
               fill="#D3D3D3", color="#000000", size=0.15) +
      scale_x_continuous(name = "lon", limits = c(-91, -78)) +
      scale_y_continuous(name = "lat", limits = c(36.5, 43)) +
      geom_point(data = coord, aes(x= Lon, y = Lat), shape = 1) +
      geom_point(data = current_members, aes(x= Lon, y = Lat, 
                                             color = pcw$rotation[,i])) +
      scale_color_gradient2(low="blue", high="red") +
      labs(title = paste0("Cluster-",i," PC-1 Eigenvectors")) + 
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
    wlt=wavelet(pc1)
    Cw=CI(0.9,pc1,"w")
    C=CI(0.9,pc1,"r")
    
    #--------------------------------------------------------------------------------#
    #Global Power Spectrum
    wt1=wt(cbind(Years,pc1))
    p4 <- get_power_spectrum(wt1)
    
    #--------------------------------------------------------------------------------#
    #Global Wavelet Spectrum
    plt_dataset <- data.frame(Period = wlt$period,
                              Power = wlt$p.avg,
                              W_noise = Cw$sig,
                              R_noise = C$sig)
    
    plt_dataset <- plt_dataset %>% filter(Period > 2^p4$ylim[1] &
                                            Period < 2^p4$ylim[2])
    
    #Transformation Function
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv, 
                log_breaks(base = base), 
                domain = c(1e-100, Inf))
    }
    
    p3 <- ggplot(plt_dataset) +
      geom_point(aes(x = Power, y = Period)) +
      geom_path(aes(x = Power, y = Period)) +
      geom_line(aes(x = W_noise, y = Period)) +
      geom_line(aes(x = R_noise, y = Period), color ='red') +
      ggtitle("Global Wavelet Spectrum") + 
      scale_x_continuous(name = "Variance") +
      scale_y_continuous(trans = reverselog_trans(2),
                         name = "Period (Years)") +
      theme_bw() +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10),
            plot.title = element_text(size=12))
    
    
    
    
    print(plot_grid(p1,p2,p3,p4$ps,
                    nrow = 2,
                    labels = c("A", "B", "C", "D"),
                    label_size = 12))
    
    plot(wt1)
    
  }
  
  dev.off()
}


#________________________________________________________________________________#
###Running the functions
site_info <- data.frame(Lat = site_info$dec_lat_va, 
                    Lon = site_info$dec_long_va)

wavclust(Dat = input_data,
         coord = site_info,
         yrs = Years)
