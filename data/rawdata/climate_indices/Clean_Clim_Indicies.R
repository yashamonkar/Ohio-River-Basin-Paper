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
#Format data - Data Frame - With Month and Year


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
enso[2014:2016,] <- NA
enso$Year <- rep(1854:2021, each = 12)
enso$Month <- rep(1:12, 168)
enso$V1 <- NULL
colnames(enso)[1] <- "ENSO"
plot(enso$ENSO, type='l')

#--------------------------------------------------------------------------------#
###NAO

st_yr <- nao$X1821[1]
nao$X1821 <- NULL #Remove the Year
nao[,13] <- NULL #Remove the average
nao[nao == -999.900] <- NA
nao_ts <- c(t(nao))
nao <- data.frame(NAO = nao_ts,
                  Year = rep(st_yr:2021, each = 12),
                  Month = rep(1:12, length(st_yr:2021)))
plot(nao$NAO, type='l')

#--------------------------------------------------------------------------------#
##AMO
st_yr <- amo$X1880[1]
amo$X1880 <- NULL
amo[amo == -999.900] <- NA
amo_ts <- c(t(amo))
amo <- data.frame(AMO = amo_ts,
                  Year = rep(st_yr:2021, each = 12),
                  Month = rep(1:12, length(st_yr:2021)))
plot(amo$AMO, type='l')

#--------------------------------------------------------------------------------#
##PDO
st_yr <- pdo$X1880[1]
pdo$X1880 <- NULL
pdo[pdo == -999.900] <- NA
pdo_ts <- c(t(pdo))
pdo <- data.frame(PDO = pdo_ts,
                  Year = rep(st_yr:2021, each = 12),
                  Month = rep(1:12, length(st_yr:2021)))
plot(pdo$PDO, type='l')


#________________________________________________________________________________#
#Saving the data
climate_indices <- list(ENSO = enso,
                        NAO = nao,
                        PDO = pdo,
                        AMO = amo)

save(climate_indices, file="Climate_Indices.RData")
