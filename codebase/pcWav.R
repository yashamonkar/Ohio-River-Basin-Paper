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
####Functions to compute the PC-Wavelets

###Input
#1. Data Frame of the Annual Maximum values(Dat)
#2. Number of PCs to be analyzed (npcs)
#3. Lat-Lon coordinates of the Sites (coords)
#4. Years (yrs)



#________________________________________________________________________________#

#####PCwav Function#####
PCwav=function(Dat,npcs,coords,yrs){
  
  #Load Dependencies
  require("biwavelet")
  require("plotrix")
  require("maps")
  require('ggplot2')
  
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
    #Global Power Spectrum
    wt1=wt(cbind(Years,pc))
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
 
}



#________________________________________________________________________________#
#Running the PC-Wavelet Analysis
site_loc <- data.frame(Lat = site_info$dec_lat_va, 
                       Lon = site_info$dec_long_va)

pdf("figures/PC_Wavelet_Analysis.pdf")
PCwav(Dat = input_data, 
      npcs = 3,
      coords = site_loc,
      yrs = Years)
dev.off()
