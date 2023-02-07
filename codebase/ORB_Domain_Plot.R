


#________________________________________________________________________________#
#Set-up Directory
setwd("~/GitHub/Ohio-River-Basin-Paper")


# Load library
require(ggplot2)
require(sf)
require(rgeos)
require(sp)
require(rgdal)



#________________________________________________________________________________#
####Reading the Data####

#Read the Site Information
site_info <- read.table("data/site_information.txt", sep="", header = TRUE)


# Get HydroRivers Dataset
hydrorivers <- sf::st_read(paste0("data/rawdata/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp"))


#Get the ORB Shape File
orb=st_read("data/rawdata/ORBShape/ORBShape/WBDHU6.shp")%>%st_transform(crs="+proj=longlat +datum=NAD83")


#________________________________________________________________________________#
####Initial Data Cleaning

#Remove the sites with visible flow regulation. 
site_info <- site_info[-c(11,12,15),]


####Subset Rivers based on drainage area
threshold = 250 # Set threshold
shape_river_small <- hydrorivers[as.numeric(hydrorivers$UPLAND_SKM) > threshold,]


###Plotting the Data PC-1
world <- map_data("world")
us <- map_data("state")


#Getting the Intersection -- TO BE DONE IN SP 
ORB.map <- readOGR("data/rawdata/ORBShape/ORBShape/WBDHU6.shp") 
Rivers.map  <- as_Spatial(shape_river_small)
proj4string(Rivers.map) <- proj4string(ORB.map) #Getting them to the same projections

ORB_Rivers <- over(Rivers.map, ORB.map)
shape_river_small$Presence <- ORB_Rivers$TNMID
shape_river_small <- shape_river_small[!is.na(shape_river_small$Presence), ]

#________________________________________________________________________________#

#Add the Rivers to the Ohio Basin Data

pdf("figures/ORB_Domain.pdf")

ORB_gauges <- ggplot() + 
  geom_map(data=us, map=us,
           aes(x=long, y=lat, map_id=region),
           fill="#D3D3D3", color="#000000", size=0.15) +
  geom_sf(data = orb, size = 0.1, color = "cyan1", fill = "cyan1") + 
  geom_sf(data =  shape_river_small, color = 'black', fill = 'black') +
  geom_map(data=us, map=us,
           aes(x=long, y=lat, map_id=region),
           fill=NA, color="#000000", size=0.15) +
  ggtitle("Ohio River Basin") + 
  geom_point(data=site_info, aes(x=dec_long_va,y=dec_lat_va), color = 'red', 
             size = 3) +
  coord_sf() +
  xlab("Longitude")+
  ylab("Latitude") +
  scale_x_continuous(limits = c(-90, -78)) +
  scale_y_continuous(limits = c(35, 43)) +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank())

plot(ORB_gauges)


dev.off()