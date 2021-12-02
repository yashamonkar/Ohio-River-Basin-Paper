#________________________________________________________________________________#
###Code to compute the Water Year
#A water year is defined from Oct-Sept instead of Jan-Dec.

###Input
#1. Actual Year
#2. Current Month

###Output
#1. Water Year

#________________________________________________________________________________#
get_water_year <- function(Yrs,Mns){
  wt_yrs <- rep(NA,length(Yrs))
  for(i in 1:length(Mns)){
    if(Mns[i] > 9){
      wt_yrs[i] = Yrs[i]+1
    } else{ 
        wt_yrs[i] = Yrs[i]
    }
  }
  
  return(wt_yrs)
}
