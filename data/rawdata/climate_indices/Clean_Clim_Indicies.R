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
Months <- c("Jan","Feb","Mar","April","May","June", "Jul","Aug","Sept","Oct","Nov","Dec")


#________________________________________________________________________________#
###ENSO

nino <- ts(enso$V2, start = c(1854,1), frequency = 12)


temp <- do.call(rbind, strsplit(as.character(nino_34[,1]),"\\."))
nino_34$V1 <- NULL
colnames(nino_34) <- c("ENSO")
nino_34$Year <- nino_34$Month <- temp[,1] 
nino_34$Month[1:1788] <- rep(Months,149)
nino_34$Month[1789:1798] <- Months[1:10]
nino_34$Water_Year <- NA
for(i in 1:nrow(nino_34)){
  if(nino_34$Month[i] %in% c("Oct","Nov","Dec")){
    nino_34$Water_Year[i] <- as.numeric(nino_34$Year[i])+1
  }
  else{
    nino_34$Water_Year[i] <- nino_34$Year[i]
  }
}
write.table(nino_34, "Nino_34.txt", sep = " ")

#NAO
nao_temp <- list()
for(i in 1:nrow(nao)){
  nao_temp[[i]] <- nao[i,2:13]}
nao <- data.frame(Year = rep(1822:2019,each = 12), Month = rep(Months, 198), NAO=unlist(nao_temp))
nao <- head(nao,-3) #Removing the months where we do not have data. 
nao$Water_Year <- NA
for(i in 1:nrow(nao)){
  if(nao$Month[i] %in% c("Oct","Nov","Dec")){
    nao$Water_Year[i] <- as.numeric(nao$Year[i])+1
  }
  else{
    nao$Water_Year[i] <- nao$Year[i]
  }
}
nao$NAO[nao$NAO == -999.9] <- NA
nao <- nao[complete.cases(nao), ]
write.table(nao, "NAO.txt", sep = " ")

#PDO
pdo_temp <- list()
for(i in 1:nrow(pdo)){
  pdo_temp[[i]] <- pdo[i,2:13]}
pdo <- data.frame(Year = rep(1881:2020,each=12), Month = rep(Months, 140), PDO = unlist(pdo_temp))
pdo$Water_Year <- NA
for(i in 1:nrow(pdo)){
  if(pdo$Month[i] %in% c("Oct","Nov","Dec")){
    pdo$Water_Year[i] <- as.numeric(pdo$Year[i])+1
  }
  else{
    pdo$Water_Year[i] <- pdo$Year[i]
  }
}
pdo$PDO[pdo$PDO == -999.9] <- NA
pdo <- pdo[complete.cases(pdo), ]
write.table(pdo, "PDO.txt", sep = " ")


#AMO
amo_temp <- list()
for(i in 1:nrow(amo)){
  amo_temp[[i]] <- amo[i,2:13]}
amo <- data.frame(Year = rep(1855:2020,each=12), Month = rep(Months, 166), AMO = unlist(amo_temp))
amo$Water_Year <- NA
for(i in 1:nrow(amo)){
  if(amo$Month[i] %in% c("Oct","Nov","Dec")){
    amo$Water_Year[i] <- as.numeric(amo$Year[i])+1
  }
  else{
    amo$Water_Year[i] <- amo$Year[i]
  }
}
write.table(amo, "AMO.txt", sep = " ")