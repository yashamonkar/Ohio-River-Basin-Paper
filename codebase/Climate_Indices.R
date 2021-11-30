#________________________________________________________________________________#
###Code to analyze relationship between climate indices and AM Rainfall.

###Output Needed - For each Climate Index
#1. Global Wavelet Spectrum (Climate Index)
#2. Wavelet Power Spectrum (Climate Index)
#3. Wavelet Coherence Spectrum (Climate Index and PC-1 AM Rainfall)


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
source('functions/Get_Wavelet_Coherence.R')

#________________________________________________________________________________#
####Reading the Data####

#Read the Annual Maximum Data
input_data <- read.table("data/Annual_Maximum_Streamflow.txt", 
                         sep="", header = TRUE)

#Read the Site Information
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)

#Read Climate Indices
load("data/Climate_Indices.RData")

#Data Cleaning
Years <- input_data$Year
input_data$Year <- NULL


#________________________________________________________________________________#
####Data-Wrangling Climate Indices

#Select the Months
Season <- c(1,2,3)

enso <- climate_indices$ENSO %>% 
  group_by(Year) %>%
  filter(Month == Season & Year > (Years[1]-1)) %>%
  summarise(ENSO = mean(ENSO))

nao <- climate_indices$NAO %>% 
  group_by(Year) %>%
  filter(Month == Season & Year > (Years[1]-1)) %>%
  summarise(NAO = mean(NAO))

pdo <- climate_indices$PDO %>% 
  group_by(Year) %>%
  filter(Month == Season & Year > (Years[1]-1)) %>%
  summarise(PDO = mean(PDO))

amo <- climate_indices$AMO %>% 
  group_by(Year) %>%
  filter(Month == Season & Year > (Years[1]-1)) %>%
  summarise(AMO = mean(AMO))


#________________________________________________________________________________#
###Function to analyze the Climate Connections
###Input
#1. Annual Maximum Data(AM_Data)
#2. Climate Index(Climate_Index)
#3. Climate Index Name(Field)
#4. Years(Yrs)

get_clim_connections <- function(AM_Data, Climate_Index,
                                 Yrs, Field){
  
  #Load Dependencies
  library(ggplot2)
  library(biwavelet)
  
  #PCA on AM Rainfall
  pcs <- prcomp(AM_Data, scale = TRUE)
  pc1 <- pcs$x[,1]
  
  #--------------------------------------------------------------------------------#
  ###Wavelet Analysis on Climate Index
  wlt=wavelet(Climate_Index)
  Cw=CI(0.9,Climate_Index,"w")
  C=CI(0.9,Climate_Index,"r")
  
  ###Global Power Spectrum
  wt1=wt(cbind(Yrs,Climate_Index))
  p2 <- get_power_spectrum(wt1)
  
  ###Global Wavelet Spectrum
  plt_dataset <- data.frame(Period = wlt$period,
                            Power = wlt$p.avg,
                            W_noise = Cw$sig,
                            R_noise = C$sig)
  
  plt_dataset <- plt_dataset %>% filter(Period > 2^p2$ylim[1] &
                                          Period < 2^p2$ylim[2])
  
  #Transformation Function
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  
  p1 <- ggplot(plt_dataset) +
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
  
  #--------------------------------------------------------------------------------#
  #Cross Wavelet Spectrum
  t1 = cbind(Yrs, Climate_Index)
  t2 = cbind(Yrs, pc1)
  wtc.AB = wtc(t1, t2)
  
  p3 <- get_wavelet_coherence(x=wtc.AB,
                              Fields = paste0("PC-1 & ", Field))
  
  
  ###Plotting the results
  print(plot_grid(p1,p2$ps,p3,
                  nrow = 1,
                  labels = c("A", "B", "C"),
                  label_size = 12))
  
}


#________________________________________________________________________________#
###Function to analyze the results

pdf("figures/Climate_Connections.pdf")
get_clim_connections(AM_Data = input_data,
                     Climate_Index = enso$ENSO,
                     Yrs = enso$Year,
                     Field = "ENSO")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = nao$NAO,
                     Yrs = nao$Year,
                     Field = "NAO")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = amo$AMO,
                     Yrs = amo$Year,
                     Field = "AMO")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = pdo$PDO,
                     Yrs = pdo$Year,
                     Field = "PDO")

dev.off()

