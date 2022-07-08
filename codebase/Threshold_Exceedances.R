#________________________________________________________________________________#
###Code to plot the threshold exceedances###

###Output Needed - 
#1. Plot of Threhold Exceedances



#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")

#Load Dependencies
library(dplyr)
library(plotrix)
library(maps)
library(cowplot)
library(ggplot2)
library(scales)

#Load Functions
source('functions/Get_Water_Year.R')

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

pdf("figures/Exceedances.pdf")

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

plt_dataset <- data.frame(Exceedances = count_exceedances,
                          Year = Years)

ggplot(plt_dataset, aes(x=Years, y = count_exceedances)) +
  geom_point() +
  geom_line() +
  ylab("Count Exceedances") +
  labs(title = "Exceedances across the Ohio River Basin") + 
  theme_bw() +
  theme(plot.title = element_text(size=18),
    axis.text=element_text(size=12),
    axis.title=element_text(size=15),
    legend.text=element_text(size=15))


dev.off()