#______________________________________________________________________________#
###Script to get Streamflow Data and Annual Maximum Streamflow###


###Set up the Path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Get Water Year Annual Mximum Streamflow Data
source("functions/Get_WaterYear_Streamflow.R")
pdf("Streamflow_API.pdf")
get_max_wateryear_streamflow(lat_min=35, 
                   lat_max=42, 
                   lon_min=-90, 
                   lon_max=-78, 
                   drain_area=5791*0.25,
                   start_date= "1934-10-01",
                   end_date = "2021-09-30",
                   missing_data = 0.01 #Percent
)
dev.off()
