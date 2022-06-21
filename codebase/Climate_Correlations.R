#________________________________________________________________________________#
###Code to analyze correlations between annual maximum rainfall and climate indices.

###Output Needed - For each Climate Index
#1. Just correlation plots


#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")

#Load Dependencies
library(dplyr)
library(corrplot)


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

#Data Cleaning
Years <- input_data$Year
input_data$Year <- NULL


#________________________________________________________________________________#
####Data-Wrangling Climate Indices

#Select the Months
Season <- c(2,3,4)

#ENSO
climate_indices$ENSO$Water_Year <- get_water_year(Yrs = climate_indices$ENSO$Year,
                                                  Mns = climate_indices$ENSO$Month)
enso <- climate_indices$ENSO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(ENSO = mean(ENSO))
enso <- enso[complete.cases(enso), ]

#NAO
climate_indices$NAO$Water_Year <- get_water_year(Yrs = climate_indices$NAO$Year,
                                            Mns = climate_indices$NAO$Month)
nao <- climate_indices$NAO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(NAO = mean(NAO))
nao <- nao[complete.cases(nao),]

#PDO
climate_indices$PDO$Water_Year <- get_water_year(Yrs = climate_indices$PDO$Year,
                                           Mns = climate_indices$PDO$Month)
pdo <- climate_indices$PDO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(PDO = mean(PDO))
pdo <- pdo[complete.cases(pdo),]

#AMO
climate_indices$AMO$Water_Year <- get_water_year(Yrs = climate_indices$AMO$Year,
                                           Mns = climate_indices$AMO$Month)
amo <- climate_indices$AMO %>% 
  group_by(Water_Year) %>%
  filter(Month %in% Season & Water_Year > (Years[1]-1)) %>%
  summarise(AMO = mean(AMO))
amo <- amo[complete.cases(amo),]

#________________________________________________________________________________#
###Analyze the correlations###

#PCA on Rainfall data
pcs <- prcomp(input_data, scale = TRUE)
pc1 <- pcs$x[,1]



#----------------Correlation on Indices-----------------------------------------#
#Code for significance
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}








#----------------Correlation on Interactions------------------------------------#
cor_data <- data.frame(PC1 = pc1,
                       ENSO = enso$ENSO,
                       PDO = pdo$PDO,
                       AMO = amo$AMO,
                       NAO = nao$NAO,
                       ENSO_PDO = enso$ENSO*pdo$PDO,
                       ENSO_AMO = enso$ENSO*amo$AMO,
                       ENSO_NAO = enso$ENSO*nao$NAO,
                       PDO_AMO = pdo$PDO*amo$AMO,
                       PDO_NAO = pdo$PDO*nao$NAO,
                       AMO_NAO = amo$AMO*nao$NAO)


# matrix of the p-value of the correlation
cor.mat <- cor(cor_data)
p.mat <- cor.mtest(cor_data)


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

#pdf("figures/Climate_Correlations.pdf")
corrplot(cor.mat, 
         method="color", 
         col=col(200),  
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.8, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.1, 
         insig = "blank")
#dev.off()