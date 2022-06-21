#________________________________________________________________________________#
###Code to analyze seasonality in the Annual Maximum Streamflow occurrent.

###Output Needed - 
#1. Rose plot


#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")

#Load Dependencies
library(dplyr)
library(corrplot)


#Load Functions


#________________________________________________________________________________#
####Reading the Data####

#Read the Annual Maximum Data
input_data <- read.table("data/Annual_Maximum_Dates.txt", 
                         sep="", header = TRUE)

#Read the Site Information
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)


#Data Cleaning
Years <- input_data$Year
input_data$Year <- NULL

#PDF to save plots
pdf("figures/Occurrence_Seasonality.pdf")


#________________________________________________________________________________#
####Rose plot of the entire data

#Compute the Day of Year
input_doy <- list()
for(i in 1:ncol(input_data)){
  input_doy[[i]] <- as.numeric(strftime(as.Date(input_data[,i], "%Y-%m-%d"), format = "%j"))
}
input_doy <- unlist(input_doy)

#Plot the histogram
hist(input_doy, xlim = c(0,366), breaks = 366,
     main = "Annual Maximum Streamflow Occurrence", xaxt='n',
     xlab = "Day of Year")
axis(side=1, at=c(16, 47, 76, 107, 137, 168, 198, 229, 260, 290, 321, 350),
     labels=c("Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug", "Sep", "Oct",
              "Nov", "Dec"))




#________________________________________________________________________________#
####Rose plots for the site level data
par(mfrow = c(3,3))

for(s in 1:ncol(input_data)){

  #Compute the Day of Year
  input_doy <- as.numeric(strftime(as.Date(input_data[,s], "%Y-%m-%d"), format = "%j"))

  
  #Plot the histogram
  hist(input_doy, xlim = c(0,366), breaks = 30,
       main = paste0("Site - ", s), xaxt='n', yaxt='n',
       xlab = "Day of Year", ylab = NULL)
  axis(side=1, at=c(16, 47, 76, 107, 137, 168, 198, 229, 260, 290, 321, 350),
       labels=c("Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug", "Sep", "Oct",
                "Nov", "Dec"))
  
  

}

dev.off()
