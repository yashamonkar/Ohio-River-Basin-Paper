#________________________________________________________________________________#
###Code to analyze relationships between annual maximum rainfall and climate indices.
###RANDOM FORESTS

###Output Needed - For each Climate Index
#1. Just correlation plots


#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")

#Load Dependencies
library(dplyr)
library(corrplot)
library(pracma)
library(randomForest)


#Load Functions
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

###Monthly Temperatures
monthly_temp <- read.table("data/Monthly_Temperatures.txt", sep="")

#Data Cleaning
Years <- input_data$Year
input_data$Year <- NULL

#Remove the sites with visible flow regulation. 
site_info <- site_info[-c(11,12),]
input_data <- input_data[,-c(11,12)]




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


####PCA on Rainfall data
pcs <- prcomp(input_data, scale = TRUE)
pc1 <- pcs$x[,1]
pc2 <- pcs$x[,2]



###Annual Temperatures
annual_temp <- head(monthly_temp,-1)
annual_temp <- data.frame(Year=annual_temp[,1], avg_temp=rowMeans(annual_temp[,-1]))
annual_temp <- annual_temp %>% filter(Year > 1934)


#_______________________________________________________________________________#
###----------------------------------Random Forests-------------------------###

#Create the PDF
pdf("figures/Random_Forest.pdf")

#Set the grid design
par(mfrow = c(1,1))


#_______________________________________________________________________________#
#PC-1
#Create the dataset
pc_type <- 1
reg_data <- data.frame(PC = scale(pc1),
                       ENSO = enso$ENSO,
                       PDO = pdo$PDO,
                       AMO = amo$AMO,
                       NAO = nao$NAO,
                       ENSO_PDO = enso$ENSO*pdo$PDO,
                       ENSO_AMO = enso$ENSO*amo$AMO,
                       ENSO_NAO = enso$ENSO*nao$NAO,
                       PDO_AMO = pdo$PDO*amo$AMO,
                       PDO_NAO = pdo$PDO*nao$NAO,
                       AMO_NAO = amo$AMO*nao$NAO,
                       Temp = annual_temp$avg_temp)



###Fitting the RANDOM FOREST
rf.fit <- randomForest(PC ~ ., data=reg_data, ntree=10000,
                       keep.forest=FALSE, importance=TRUE, mtry = 4)
#print(rf.fit)
#plot(rf.fit)

#Plot the Variable Importance
varImpPlot(rf.fit, sort = TRUE,  type='2',
           main = paste0("Variable Importance \n PC-",pc_type))




#_______________________________________________________________________________#
#PC-2
#Create the dataset
pc_type <- 2
reg_data <- data.frame(PC = scale(pc2),
                       ENSO = enso$ENSO,
                       PDO = pdo$PDO,
                       AMO = amo$AMO,
                       NAO = nao$NAO,
                       ENSO_PDO = enso$ENSO*pdo$PDO,
                       ENSO_AMO = enso$ENSO*amo$AMO,
                       ENSO_NAO = enso$ENSO*nao$NAO,
                       PDO_AMO = pdo$PDO*amo$AMO,
                       PDO_NAO = pdo$PDO*nao$NAO,
                       AMO_NAO = amo$AMO*nao$NAO,
                       Temp = annual_temp$avg_temp)



###Fitting the RANDOM FOREST
rf.fit <- randomForest(PC ~ ., data=reg_data, ntree=10000,
                       keep.forest=FALSE, importance=TRUE, mtry = 4)
#print(rf.fit)
#plot(rf.fit)

#Plot the Variable Importance
varImpPlot(rf.fit, sort = TRUE,  type='2',
           main = paste0("Variable Importance \n PC-",pc_type))


#Set the grid design
par(mfrow = c(1,1))

dev.off()