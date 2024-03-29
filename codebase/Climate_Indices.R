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
library(scales)
library(pracma)


#Load Functions
source('functions/Get_Power_Spectrum.R')
source('functions/Get_Wavelets.R')
source('functions/Get_Wavelet_Coherence.R')
source('functions/Get_Water_Year.R')

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

#Remove the sites with visible flow regulation. 
site_info <- site_info[-c(11,12,15),]
input_data <- input_data[,-c(11,12,15)]



#________________________________________________________________________________#
####Data-Wrangling Climate Indices


###ENSO
Season <- c(1,2,12) #Select the Months

climate_indices$ENSO$Water_Year <- get_water_year(Yrs = climate_indices$ENSO$Year,
                                                  Mns = climate_indices$ENSO$Month)
enso <- climate_indices$ENSO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(ENSO = mean(ENSO))
enso <- enso[complete.cases(enso), ]
#enso$ENSO <- detrend(enso$ENSO, 'linear')

###NAO
Season <- c(1,2,12) #Select the Months

climate_indices$NAO$Water_Year <- get_water_year(Yrs = climate_indices$NAO$Year,
                                                 Mns = climate_indices$NAO$Month)
nao <- climate_indices$NAO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(NAO = mean(NAO))
nao <- nao[complete.cases(nao),]
#nao$NAO <- detrend(nao$NAO, 'linear')

###PDO
Season <- c(1,2,12) #Select the Months

climate_indices$PDO$Water_Year <- get_water_year(Yrs = climate_indices$PDO$Year,
                                                 Mns = climate_indices$PDO$Month)
pdo <- climate_indices$PDO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(PDO = mean(PDO))
pdo <- pdo[complete.cases(pdo),]
#pdo$PDO <- detrend(pdo$PDO, 'linear')

###AMO
Season <- c(1,2,12) #Select the Months

climate_indices$AMO$Water_Year <- get_water_year(Yrs = climate_indices$AMO$Year,
                                                 Mns = climate_indices$AMO$Month)
amo <- climate_indices$AMO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(AMO = mean(AMO))
amo <- amo[complete.cases(amo),]
#amo$AMO <- detrend(amo$AMO, 'linear')


###Interactions

ENSO_PDO = (enso$ENSO - mean(enso$ENSO))*(pdo$PDO - mean(pdo$PDO))
#ENSO_PDO = detrend(ENSO_PDO, 'linear')

ENSO_AMO = (enso$ENSO - mean(enso$ENSO))*(amo$AMO - mean(amo$AMO))
#ENSO_AMO = detrend(ENSO_AMO, 'linear')

ENSO_NAO = (enso$ENSO - mean(enso$ENSO))*(nao$NAO - mean(nao$NAO))
#ENSO_NAO = detrend(ENSO_NAO, 'linear')

PDO_AMO = (pdo$PDO - mean(pdo$PDO))*(amo$AMO - mean(amo$AMO))
#PDO_AMO = detrend(PDO_AMO, 'linear')

PDO_NAO = (pdo$PDO - mean(pdo$PDO))*(nao$NAO - mean(nao$NAO))
#PDO_NAO = detrend(PDO_NAO, 'linear')

AMO_NAO = (amo$AMO - mean(amo$AMO))*(nao$NAO - mean(nao$NAO))
#AMO_NAO = detrend(AMO_NAO, 'linear')


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
  
  #Plotting the Climate Index
  plt_dataset <- data.frame(Clim = Climate_Index,
                            Year = Yrs,
                            Loess = lowess(Yrs,Climate_Index,f=1/9)$y)
  p1 <- ggplot(plt_dataset) +
    geom_line(aes(x = Year, y = Loess), size = 1.2, color ='red') +
    geom_point(aes(x = Year, y = Clim), size = 0.1) +
    geom_line(aes(x = Year, y = Clim), size = 0.1) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Index Value") +
    labs(title = paste0(Field)) + 
    theme_bw() +
    theme(plot.title = element_text(size=12),
          axis.text=element_text(size=5),
          axis.title=element_text(size=10)) 
  
  
  #PCA on AM Rainfall
  pcs <- prcomp(AM_Data, scale = TRUE)
  pc1 <- pcs$x[,1]
  pc1 = detrend(pc1, 'linear')
  
  #--------------------------------------------------------------------------------#
  ###Wavelet Analysis on Climate Index
  wlt=wavelet(Climate_Index)
  Cw=CI(0.9,Climate_Index,"w")
  C=CI(0.9,Climate_Index,"r")
  
  ###Global Power Spectrum
  wt1=wt(cbind(Yrs,Climate_Index))
  p4 <- get_power_spectrum(wt1)
  
  ###Global Wavelet Spectrum
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
    ggtitle(paste0(Field)) + 
    scale_x_continuous(name = "Variance") +
    scale_y_continuous(trans = reverselog_trans(2),
                       name = "Period (Years)") +
    theme_bw() +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=15),
          plot.title = element_text(size=15))
  
  #--------------------------------------------------------------------------------#
  #Wavelet Coherence
  t1 = cbind(Yrs, Climate_Index)
  t2 = cbind(Yrs, pc1)
  wtc.AB = wtc(t1, t2)
  
  p2 <- get_wavelet_coherence(x=wtc.AB,
                              Fields = paste0("with PC-1"))
  
  
  ###Plotting the results
  #print(plot_grid(p1,p2,p3,p4$ps,
  #                nrow = 2,
  #                labels = c("A", "B", "C", "D"),
  #                label_size = 12))
  
  print(plot_grid(p3,p4$ps,p2,
                  nrow = 1,
                  labels = c("A", "B", "C"),
                  label_size = 12))
  
  plot(wtc.AB, plot.phase = TRUE)
  
}


#________________________________________________________________________________#
###Function to analyze the results

pdf("figures/Climate_Connections.pdf", height=1850/300, width=5000/300)
#pdf("Trial.pdf", height=1850/300, width=5000/300)
get_clim_connections(AM_Data = input_data,
                     Climate_Index = enso$ENSO,
                     Yrs = enso$Water_Year,
                     Field = "ENSO - DJF")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = nao$NAO,
                     Yrs = nao$Water_Year,
                     Field = "NAO - DJF")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = amo$AMO,
                     Yrs = amo$Water_Year,
                     Field = "AMO - DJF")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = pdo$PDO,
                     Yrs = pdo$Water_Year,
                     Field = "PDO - DJF")


get_clim_connections(AM_Data = input_data,
                     Climate_Index = ENSO_PDO,
                     Yrs = pdo$Water_Year,
                     Field = "ENSO-PDO")

get_clim_connections(AM_Data = input_data,
                     Climate_Index = ENSO_AMO,
                     Yrs = pdo$Water_Year,
                     Field = "ENSO-AMO")

get_clim_connections(AM_Data = input_data,
                     Climate_Index = ENSO_NAO,
                     Yrs = pdo$Water_Year,
                     Field = "ENSO-NAO")

get_clim_connections(AM_Data = input_data,
                     Climate_Index = PDO_AMO,
                     Yrs = pdo$Water_Year,
                     Field = "PDO-AMO")

get_clim_connections(AM_Data = input_data,
                     Climate_Index = PDO_NAO,
                     Yrs = pdo$Water_Year,
                     Field = "PDO-NAO")

get_clim_connections(AM_Data = input_data,
                     Climate_Index = AMO_NAO,
                     Yrs = pdo$Water_Year,
                     Field = "AMO_NAO")

dev.off()

