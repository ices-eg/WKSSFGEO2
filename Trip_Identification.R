# File: Trip_Identification.R
# Author: Guillermo Martin 
# Template Created: Tue Mar 14 12:12:14 2023
# ---------------------------------------------------------------------------
# Description:
# Take geoposition data, identify points in harbour, and define individual trips by vessel
# when interval between two consecutive pings > 2Hr (7200sec).

# Most of the code extracted from the HMM classification function developed for the Marine Institute by:

# Programmer: Sophie de Grissac - DioMed&a Science
# COMPANY N° (SIRET): 882 532 245 00017
# last update: 29/03/2020
# email: s.degrissac@gmail.com
# Contact me for any update, question, problem related to this function
# ---------------------------------------------------------------------------

library(geosphere)


####    SOURCES

# burst function

## dtime is the time between one location and the NEXT
## this function return "burtsid", a numeric vector corresponding to burst numbers, i.e. blocks of successive locations forming a unique trip.
## threshold is the minimum time period between 2 locs after which you consider that these 2 locs are 2 different trips.
## In the case of fishing vessel, fishing trips are tipically separated by at least 12h. The treshold is arbitrary set to 2 hours based on the assumption that vms do not stop transmitting while at sea.

create_burst <- function(dtime, treshold)
{
  burstid=0
  j=0
  for (i in 1:(length(dtime)-1))
  {
    if (dtime[i] < treshold)
    {
      burstid <- c(burstid,j)
    }
    else
    {
      burstid <- c(burstid,j+1)
      j <- j+1
    }
  }
  return(burstid)
}


#__________________________________________________________________________________________________________________________


####    DATA INPUT

## open vessel data
# GPS location of .csv format
vessel <- read.csv(file=file.path(path_vessel)) # open GPS locations file already filter with harbour locations and with ping rate information

# harbour location in .csv format (only center lat/longs required)
harbour <- read.csv(path_harbour) 

#_________________________________________________________________________________________________________________________________________


#####     FORMATING

# names and date formating
# Modify accordingly based on your data to select VESSEL ID, Date_TIME, lon and lat 

vessel <- vessel[,1:4]
names(vessel) <- c('ID', 'Date', 'longitude', 'latitude')
vessel$Date <- as.POSIXct(vessel[,2], format = "%d/%m/%Y %H:%M:%S", tz='UTC') # format the date
#reorder by ID and Date
vessel <- vessel[order(vessel$ID, vessel$Date),]
vessel$ID <- as.factor(vessel$ID)


## remove harbour locations
v <- vessel
pb <- txtProgressBar(1, nrow(harbour), style = 3) #GM_edit: Progress bar 

message("Defining Harbour locations...") #GM_edit
for (i in 1:nrow(harbour)) {
  setTxtProgressBar(pb, i) #GM_edit:
  dist_harbour <- distHaversine(c(harbour[i,"Longitude_Center"],harbour[i,"Latitude_Center"]),v[,c("longitude","latitude")],r=6378.160) # distance to harbour in km
  v <- v[which(dist_harbour > 1),]
}
vessel <- v
rm(v, dist_harbour)


#### calculate dtime (time between successive locations) and split the different daily trips of each vessel
#### & create burst, i.e. separate trips when time between pings is 2h or more (arbitrary)
#### This uses the function "create burst" sourced in the first part of this script

pb <- txtProgressBar(0, max = length(unique(vessel$ID)), style = 3) #GM_edit: Progress bar 
p=0

message("Defining individual vessel trips...") #GM_edit
v <- NULL
#vessel$ID_trip <- NA
for (j in 1:length(unique(vessel$ID))) {
  p<-p+1 #GM_edit
  setTxtProgressBar(pb, p) #GM_edit
  
  d <- vessel[vessel$ID == unique(vessel$ID)[j],]
  d$dtime <- NA
  for (i in 1:nrow(d)-1){
    d$dtime[i] <- difftime(d$Date[i+1],d$Date[i],units="sec")
    # add a column with trip ID ("vesselID_TripNumber")
  }
  b <- create_burst(d$dtime, 7200)
  d$ID_trip <- paste(d$ID,b,sep="_")  
  v <- rbind(v, d)
}

vessel <- v
rm(d,b, v)

vessel$ID_trip <- as.factor(vessel$ID_trip)

#set NA dtime for the last point of each trip
pb <- txtProgressBar(0, max = length(unique(vessel$ID_trip)), style = 3) #GM_edit: Progress bar

message("set NA dtime for the last point of each trip") #GM_edit
for (i in 1:length(unique(vessel$ID_trip))) {
  setTxtProgressBar(pb, i)
  vessel$dtime[vessel$ID_trip == unique(vessel$ID_trip)[i]][length(vessel$dtime[vessel$ID_trip == unique(vessel$ID_trip)[i]])] <- NA
}

#reorder lines by ID and Date
vessel <- vessel[order(vessel$ID, vessel$Date),]

