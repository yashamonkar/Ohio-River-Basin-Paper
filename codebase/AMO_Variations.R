#________________________________________________________________________________#
###Code to analyze relationship between AMO and its multidecadal variations.

###Output Needed - 
#1. Plot of AMO and Count Aggregated Rainfall



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

pdf("figures/AMO_Exceedances.pdf")

#________________________________________________________________________________#
####Data-Wrangling Climate Indices



#Count Agggregate Exceedances
thresh <- apply(input_data, 2, quantile , probs = 0.9 , na.rm = TRUE )
exceedances <- input_data

for(i in 1:ncol(exceedances)){
  exceedances[,i] <- ifelse(input_data[,i] < thresh[i], 0, 1)
}

frac_exceedance <- rowSums(exceedances)/ncol(exceedances)
count_exceedances <- rowSums(exceedances)



###AMO
Season <- c(1,2,12) #Select the Months

climate_indices$AMO$Water_Year <- get_water_year(Yrs = climate_indices$AMO$Year,
                                                 Mns = climate_indices$AMO$Month)
amo <- climate_indices$AMO %>% 
  group_by(Water_Year) %>%
  summarise(AMO = mean(AMO))
amo <- amo[complete.cases(amo),]
amo$AMO <- detrend(amo$AMO, 'linear')
amo$AMO <- scale(amo$AMO)



#Plotting the Climate Index
plt_dataset <- data.frame(Clim = amo$AMO,
                          Year = amo$Water_Year,
                          Loess = lowess(amo$Water_Year,amo$AMO,f=1/9)$y)
p1 <- ggplot(plt_dataset) +
  geom_line(aes(x = Year, y = Loess), size = 1.2, color ='red') +
  geom_point(aes(x = Year, y = Clim), size = 0.1) +
  geom_line(aes(x = Year, y = Clim), size = 0.1) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Index Value (Scaled)") +
  labs(title = paste0("Atlantic Multidecadal Oscillation")) + 
  theme_bw() +
  theme(plot.title = element_text(size=18),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15)) 
p1


###Wavelet Analysis on Climate Index
wlt=wavelet(amo$AMO)
Cw=CI(0.9,amo$AMO,"w")
C=CI(0.9,amo$AMO,"r")

###Global Power Spectrum
wt1=wt(cbind(amo$Water_Year,amo$AMO))
p2 <- get_power_spectrum(wt1)

p2$ps



#--------------------------------------------------------------------------------#
#Scatter Plot

amo <- climate_indices$AMO %>% 
  group_by(Water_Year) %>%
  filter(Water_Year > (Years[1]-1)) %>%
  summarise(AMO = mean(AMO, na.rm=TRUE))
amo <- amo[complete.cases(amo),]
amo$AMO <- detrend(amo$AMO, 'linear')



plt_dataset <- data.frame(AMO = amo$AMO,
                          Year = amo$Water_Year,
                          Exd = count_exceedances,
                          Smooth_AMO = lowess(amo$Water_Year,amo$AMO,f=1/9)$y)





#--------------------------------------------------------------------------------#
#Joint Plot

d1 = data.frame(Year = amo$Water_Year, 
                y=amo$AMO, 
                Loess = lowess(amo$Water_Year,amo$AMO,f=1/6)$y) 
d2 = data.frame(Year = amo$Water_Year, 
                y=frac_exceedance, 
                Loess = lowess(amo$Water_Year,frac_exceedance,f=1/6)$y)
d1$z <- "AMO"
d2$z <- "Exceedances"
d3 <- within(d2, { y = y})
d3 <- within(d3, { Loess = Loess})
d4 <- rbind(d1, d3)

mycolors <- c("AMO"="blue", "Exceedances"="red")

p2 <- ggplot(d4, aes(x=Year, y=y, group=z, color=z)) +
  geom_path() +
  geom_point() +
  geom_line(d4, mapping = aes(x=Year, y=Loess , color = z, group = z), size = 1.25) +
  scale_y_continuous(name="Atlantic Multidecadal Oscillation",
                     sec.axis = sec_axis(~ ., name="Fractional Exceedances")) +
  scale_color_manual(name=" ", values = mycolors) +
  labs(title = paste0("AMO and Exceedances across the Ohio River Basin")) + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(color = mycolors["AMO"]),
    axis.text.y = element_text(color = mycolors["AMO"]),
    axis.title.y.right = element_text(color = mycolors["Exceedances"]),
    axis.text.y.right = element_text(color = mycolors["Exceedances"]),
    plot.title = element_text(size=18),
    axis.text=element_text(size=12),
    axis.title=element_text(size=15),
    legend.text=element_text(size=15)
  )

p2

dev.off()