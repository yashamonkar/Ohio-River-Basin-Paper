#________________________________________________________________________________#
#Code to clean climate indices from KNMI and save clean-standard version (Rdata)
#Source:- https://climexp.knmi.nl/selectindex.cgi?id=someone@somewhere
###Input
#1. iersst_nino3.4a_rel.dat.txt --> ENSO
#2. inao.dat.txt --> NAO
#3. iamo_ersst.dat --> AMO
#4. ipdo_ersst.dat --> PDO

###Output
#List of lists (saved as Rdata)
#Format data - Time Series


#________________________________________________________________________________#
#Set up the Path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Libraries
library(dplyr)

#Reading in the climate indices. 
enso <- read.table("iersst_nino3.4a_rel.dat.txt",  sep ="")
pdo <- read.table("ipdo_ersst.dat.txt", header = TRUE, sep ="", dec=".")
nao <- read.table("inao.dat.txt", header = TRUE, sep ="", dec=".")
amo <- read.table("iamo_ersst.dat.txt", header = TRUE, sep ="", dec=".")


#________________________________________________________________________________#
#--------------------------------------------------------------------------------#
###ENSO
nino_ts <- ts(enso$V2, start = c(1854,1), end = c(2021,9), frequency = 12)
plot(nino_ts)

#--------------------------------------------------------------------------------#
###NAO
st_yr <- nao$X1821[1]
nao$X1821 <- NULL #Remove the Year
nao[,13] <- NULL #Remove the average
nao[nao == -999.900] <- NA
nao_ts <- c(t(nao))
nao_ts <- ts(nao_ts, start = c(st_yr,1), frequency = 12)
plot(nao_ts)

#--------------------------------------------------------------------------------#
##AMO
st_yr <- amo$X1880[1]
amo$X1880 <- NULL
amo[amo == -999.900] <- NA
amo_ts <- c(t(amo))
amo_ts <- ts(amo_ts, start = c(st_yr,1), frequency = 12)
plot(amo_ts)

#--------------------------------------------------------------------------------#
##PDO
st_yr <- pdo$X1880[1]
pdo$X1880 <- NULL
pdo[pdo == -999.900] <- NA
pdo_ts <- c(t(pdo))
pdo_ts <- ts(pdo_ts, start = c(st_yr,1), frequency = 12)
plot(pdo_ts)

#________________________________________________________________________________#
#Saving the data
climate_indices <- list(ENSO = nino_ts,
                        NAO = nao_ts,
                        PDO = pdo_ts,
                        AMO = amo_ts)

save(climate_indices, file="Climate_Indices.RData")
