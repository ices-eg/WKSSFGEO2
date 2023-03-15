## Function originally developed by: 


# Programmer: Sophie de Grissac - DioMed&a Science
# COMPANY N° (SIRET): 882 532 245 00017
# last update: 29/03/2020
# email: s.degrissac@gmail.com
# Contact me for any update, question, problem related to this function.

# Modified slightly for the purpose of WKSSFGEO2 (column names, allowance of 10min ping rates,
# and definition of initial parameters of HMM models using the EM algorithm)

# COLUMN NAMES AND DATE_TIME column needs to be in the format defined during WKSSFGEO2!!! 

# c("Source","BoatID","TripID","DATE_TIME","LATITUDE","LONGITUDE","Gear","FsihingOperation",
#"speed","bearing","Quality_index") 

#Date_Time format: 
# "%Y-%m-%d %H:%M:%S"



###### CLASSIFICATION OF FISHING VESSEL ACTIVITY #######
# ______________ START OF PROGRAM _______________________________________________________________________________________
#__________________________________________________________________________________________________________________________

rm(list = ls())

#Path folder for results
path_folder_result<-""

#Data on the sharepoint
dat<-read.csv(file=file.path(""))

#Renaming columns if necessary
if ( ! any(grepl("BoatID", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                           c("ID_VESSEL","VesselID"))] <- "BoatID"
}

if ( ! any(grepl("TripID", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                           c("ID_TRIP","IDTRIP"))] <- "TripID"
}

if ( ! any(grepl("Longitude", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                                     c("LONGITUDE"))] <- "Longitude"
}

if ( ! any(grepl("Latitude", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                           c("LATITUDE"))] <- "Latitude"
}

if ( ! any(grepl("Latitude", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                           c("LATITUDE"))] <- "Latitude"
}

if ( ! any(grepl("Time", colnames(dat)))) {
  colnames(dat)[which(toupper(colnames(dat)) %in%
                           c("DATE_TIME","DATE","TIME","TIMESTAMP"))] <- "Time"
}

#Formatting columns
dat$TripID <- as.factor(dat$TripID)
dat$BoatID <- as.factor(dat$BoatID)
dat$Time<-as.POSIXct(dat$Time,format="%Y-%m-%d %H:%M:%S",tz="UTC")


HMM_classification_function <- function (vessel,path_folder_result) {
  

t1 <- Sys.time()

requiredPackages = c('moveHMM','geosphere','adehabitatHR','sp','CircStats',"mixtools")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, dependencies = TRUE)
  library(p,character.only = TRUE)
}


#__________________________________________________________________________________________________________________________


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



#### calculate dtime (time between successive locations) 

pb <- txtProgressBar(0, max = length(unique(vessel$BoatID)), style = 3) #GM_edit: Progress bar 
p=0

message("Defining time difference") #GM_edit
v <- NULL

for (j in 1:length(unique(vessel$BoatID))) { 
  p<-p+1 #GM_edit
  setTxtProgressBar(pb, p) #GM_edit
  
  d <- vessel[vessel$ BoatID == unique(vessel$ BoatID)[j],]
  d$dtime <- NA
  for (i in 1:nrow(d)-1){
    d$dtime[i] <- difftime(d$Time[i+1],d$Time[i],units="sec")
    # add a column with trip ID ("vesselTripIDNumber")
  }
  v <- rbind(v, d)
}
vessel <- v
rm(d, v)

#set NA dtime for the last point of each trip
pb <- txtProgressBar(0, max = length(unique(vessel$TripID)), style = 3) #GM_edit: Progress bar

message("set NA dtime for the last point of each trip") #GM_edit
for (i in 1:length(unique(vessel$TripID))) {
  setTxtProgressBar(pb, i)
  vessel$dtime[vessel$TripID == unique(vessel$TripID)[i]][length(vessel$dtime[vessel$TripID == unique(vessel$TripID)[i]])] <- NA
}

#reorder lines by ID and Date
vessel <- vessel[order(vessel$BoatID, vessel$Time),]

ntot_trip <- length(unique(vessel$TripID))

#______________________________________________________________________________________________________________________________________

#### DIFFERENTIATE TRIPS WITH DIFFERENT PING RATES (1, 3, 5 or 10 min)

## vessel data = separate tracks by ping rates
pb <- txtProgressBar(0, max = length(unique(vessel$TripID)), style = 3) #GM_edit: Progress bar
message("Separate tracks by ping rates") #GM_edit

vessel_ping <- NULL
for (i in 1:length(unique(vessel$TripID))) {
  setTxtProgressBar(pb, i) #GM_edit
  d <- vessel[vessel$TripID == unique(vessel$TripID)[i],]
  d$ping_rate_min <- ifelse(median(d$dtime, na.rm=TRUE) > 500, 10,
                            ifelse(median(d$dtime, na.rm=TRUE) > 200, 5,
                                   ifelse(median(d$dtime, na.rm=TRUE) > 90, 3, 1)))
  vessel_ping <- rbind(vessel_ping, d)
}

#reorder
vessel <- vessel_ping[order(vessel_ping$BoatID, vessel_ping$Time),]


#### Exclude short trips that have to few pings for proper classification and may lead to errors in the model
del_trips_short <- droplevels(unique(vessel$TripID[vessel$TripID %in% names(which(tapply(vessel$BoatID, vessel$TripID, length) <= 50))]))
vessel <- droplevels(vessel[vessel$TripID %in% names(which(tapply(vessel$BoatID, vessel$TripID, length) > 50)),]) #keeps only trips with > 10 pings


### Exclude trips with too irregular ping rate
vessel_filter <- NULL
del_trips_irreg <- NULL

pb <- txtProgressBar(0, max = length(unique(vessel$ping_rate_min)), style = 3) #GM_edit: Progress bar
message("Exclude trips with too irregular ping rate") #GM_edit

for (i in 1:length(unique(vessel$ping_rate_min))) {
  setTxtProgressBar(pb, i) #
  #print(i)
  n <- unique(vessel$ping_rate_min)[i]*60
  v <- droplevels(vessel[vessel$ping_rate_min == unique(vessel$ping_rate_min)[i],])
  keep <- row.names(as.data.frame(tapply(v$dtime, v$TripID, sd, na.rm=T)/n))[which(tapply(v$dtime, v$TripID, sd, na.rm=T)/n <= 1)]
  v2 <- droplevels(v[v$TripID %in% keep,])
  del <- row.names(as.data.frame(tapply(v$dtime, v$TripID, sd, na.rm=T)/n))[which(tapply(v$dtime, v$TripID, sd, na.rm=T)/n > 1)]
  vessel_filter <- rbind(vessel_filter, v2)
  del_trips_irreg <- c(del_trips_irreg, del)
}
#_______________________________________________________________________________________________________________________________________________


### Tell user how many trips have been removed from analyse

print(paste0('In addition, there was ', length(del_trips_irreg),'/', ntot_trip,' trips with ping rates too irregular to be analysed'))
print(paste0('In addition, there was ', length(del_trips_short),'/', ntot_trip,' trips with too few locations to be analysed'))
print(paste0('A total of ', length(unique(vessel$BoatID)),' vessels (', length(unique(vessel$TripID)),' daily trips) will be analysed'))

#reorder
vessel <- vessel_filter[order(vessel_filter$BoatID, vessel_filter$Time),]


#_____________________________________________________________________________________________________________________________________________

### RUN HMM MODELS

vessel_result_all <- NULL

for (j in 1:length(unique(vessel$ping_rate_min))) {
  list_mod <- list()
  vessel_result_ping <- NULL
  vessel_pp <- droplevels(vessel[vessel$ping_rate_min == unique(vessel$ping_rate_min)[j],])
  
  n <- unique(vessel_pp$ping_rate_min)
  
  #delete some duplicated data
  if (length(which(duplicated(vessel_pp[,c("Time","TripID")])==TRUE)) > 0) {
    vessel_pp <- droplevels(vessel_pp[-which(duplicated(vessel_pp[,c("Time","TripID")])==TRUE),])
  }
  
  #create a folder to store and check the results
  
  if (dir.exists(file.path(path_folder_result,paste0("Ping_",n,"_min_results")))==FALSE) {
    dir.create(file.path(path_folder_result,paste0("Ping_",n,"_min_results")))
  }
  
  setwd(file.path(path_folder_result,paste0("Ping_",n,"_min_results")))
  cwd <- getwd()
  
  print(paste0('Processing vessels with ping rate ', j, ' min.'))
  
  for (i in 1:length(unique(vessel_pp$BoatID))) {
    print(paste0('Processing vessel number ', unique(vessel_pp$BoatID)[i]))
    
    # create a specific folder to store results for each vessel
    if (dir.exists(paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i]))==FALSE) {
      dir.create(paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i]))
    }
    path_res <- paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i])
    
    if (dir.exists(paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i],'/plot_tracks'))==FALSE) {
      dir.create(paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i],'/plot_tracks'))
    }
    path_plot <- paste0(cwd,"/Vessel_",unique(vessel_pp$BoatID)[i],'/plot_tracks')
    
    
    ## select vessel i
    vessel_p <- droplevels(vessel_pp[vessel_pp$BoatID == unique(vessel_pp$BoatID)[i],])
    
    
    
    #-----------------------  Rediscretized tracks according to ping rate ------------------------#
    
    Lon=vessel_p$Longitude;Lat=vessel_p$Latitude;Date=vessel_p$Time; ID=vessel_p$BoatID; burst=vessel_p$TripID
    t = median(vessel_p$dtime, na.rm=TRUE)
    data_traj = as.ltraj(as.data.frame(cbind(Lon,Lat)), Date, 
                         id = ID, burst=burst, slsp="remove", typeII = T)
    
    infolocs = cbind(vessel_p$BoatID,vessel_p$ping_rate_min)
    
    data_redis <- redisltraj(data_traj, u = t, type="time", burst=burst) # resample at xmin fix rate
    
    data_redis_df <- ld(data_redis)[,c(1:3,11,12)]
    names(data_redis_df) <- c('longitude','latitude','Date','ID','TripID')
    
    vessel_p <- data_redis_df
    #---------------------------------------------------------------------------------------------#
    
    
    #data formating and calculate angle/step length
    vessel_prep <- prepData(trackData = vessel_p[,c(1:2,4)], type = "LL",coordNames = c("longitude", "latitude"))
    
    #GM_edit: defining step and angle change using mixtools package
    
    message("Defining initial parameters for HMM using EM algorithm")
    
    # starting parameters
    #histD<-hist(vessel_prep$step,freq=T,breaks = 50,main="",
    #            ylab = "Frequency")
    
    tmp<-vessel_prep$step[!is.na(vessel_prep$step) & vessel_prep$step != 0]
    
    #cutsp<-quantile(tmp)
    cutsp<-seq(min(tmp),max(tmp),length.out=50)
    
    mxn<-makemultdata(tmp,cuts=cutsp)
    #sel<-multmixmodel.sel(mxn$x, comps = c(1:3))
    #k<-sel[2,"Winner"]
    
    res<-tryCatch(normalmixEM(mxn$x, k =3, arbmean = TRUE, arbvar = TRUE),
                  error=function(e) NULL)
    #plot(res,density=TRUE)
    if(!is.null(res)){
      res$mu2<-signif(sort(res$mu),2)
      res$sigma2<-signif(res$sigma[order(match(res$sigma,res$mu))],2)
    } else {
      message("Not convergence of EM algorthim. Initial parameter values fixed")
      
      res$mu2<-c(0.01,0.2,0.1)
      res$sigma2<-c(0.1,0.02,0.1)
    }
    
    
    if(length(which(vessel_prep$step==0))>0) { #check for zero inflation
      zeroInflation <- TRUE } else zeroInflation <- FALSE
    if (zeroInflation == TRUE) {
      stepPar0 <- c(res$mu2[2], res$mu2[1], res$mu2[3],  res$sigma2[2], res$sigma2[1], res$sigma2[3],  0.5, 0, 0)
    } else { 
      stepPar0 <- c(res$mu2[2], res$mu2[1], res$mu2[3],  res$sigma2[2], res$sigma2[1], res$sigma2[3]) 
    } 
    anglePar0 <- c(0,0,0,10,1,10)
    
    mod3 <- fitHMM(data = vessel_prep, nbStates = 3, stepPar0 = stepPar0, anglePar0 = anglePar0)
    list_mod[[i]] <- mod3
    names(list_mod)[i] <- paste0('HMM_pr',j,'_',i)  ### save an object with all models if visual checking necessary
    
    vessel_p$HMM3<- viterbi(mod3) #include model results in dataset

    if(length(unique(vessel_p$HMM3)) <= 2) {
      
      rm(res)
      
      message("Only 2 moving behaviours identified. Re-fitting HMM with 2 states.")
      
      cutsp<-quantile(tmp,probs = c(.25,.75))
      
      mxn<-makemultdata(tmp,cuts=cutsp)
      #sel<-multmixmodel.sel(mxn$x, comps = c(1:3))
      #k<-sel[2,"Winner"]
      
      res<-tryCatch(normalmixEM(mxn$x, k =2, arbmean = TRUE, arbvar = TRUE),
                    error=function(e) NULL)
      
      if(!is.null(res)){
        res$mu2<-signif(sort(res$mu),2)
        res$sigma2<-signif(res$sigma[order(match(res$sigma,res$mu))],2)
      } else {
        message("Not convergence of EM algorthim. Initial parameter values fixed")
        
        res$mu2<-c(0.01,0.2,)
        res$sigma2<-c(0.1,0.02,)
      }
      #plot(res,density=TRUE)
      
      vessel_p<-vessel_p[ , -which(names(vessel_p) %in% c("HMM3"))]
      
      if(length(which(vessel_prep$step==0))>0) { #check for zero inflation
        zeroInflation <- TRUE } else zeroInflation <- FALSE
        if (zeroInflation == TRUE) {
          stepPar0 <- c(res$mu2[1], res$mu2[2], res$sigma2[1], res$sigma2[2],  0.5, 0)
        } else { 
          stepPar0 <- c(res$mu2[1], res$mu2[2], res$sigma2[1], res$sigma2[2]) 
        } 
        anglePar0 <- c(0,0,1,10)
        
        mod3 <- fitHMM(data = vessel_prep, nbStates = 2, stepPar0 = stepPar0, anglePar0 = anglePar0)
        list_mod[[i]] <- mod3
        names(list_mod)[i] <- paste0('HMM_pr',j,'_',i)  ### save an object with all models if visual checking necessary
        
        vessel_p$HMM3<- viterbi(mod3)
    }
    
    ###________ Speed Filter
    #to reduce mode 2. All locations classified as 2 for which speed is > than mode 3 mean speed are reclassified as mode 3.
    #dtime
    vessel_d <- NULL
    #vessel$TripID <- NA
    for (k in 1:length(unique(vessel_p$TripID))) {
      d <- vessel_p[vessel_p$TripID == unique(vessel_p$TripID)[k],]
      d$dtime <- NA
      for (kk in 1:nrow(d)-1){
        d$dtime[kk] <- difftime(d$Date[kk+1],d$Date[kk],units="sec")
      }
      vessel_d <- rbind(vessel_d, d)
    }
    
    #speed
    vessel_d$step <- vessel_prep$step
    vessel_d$angle <- vessel_prep$angle
    vessel_d$speed <- (vessel_d$step*1000)/vessel_d$dtime
    
    
    #filter
    q3 = as.numeric(quantile(vessel_d$speed[vessel_d$HMM3 ==3], na.rm=TRUE)[2])
    #q11 = as.numeric(quantile(vessel_d$speed[vessel_d$HMM3 ==1], na.rm=TRUE)[2])
    #q12 = as.numeric(quantile(vessel_d$speed[vessel_d$HMM3 ==1], na.rm=TRUE)[4])
    vessel_d$HMM3_speedFilter <- ifelse(vessel_d$HMM3 == 2 & vessel_d$speed >= q3, 3, vessel_d$HMM3)
    #vessel_d$HMM3_speed_turn_Filter <- ifelse(vessel_d$HMM3_speedFilter == 2 & vessel_d$speed >= q11 & vessel_d$speed <= q12 & vessel_d$angle < pi/8, 1, vessel_d$HMM3_speedFilter)
    
    
    vessel_p <- vessel_d
    rm(vessel_d)
    
    ### PLOT ALL TRACKS FOR VESSEL i
    png(filename=paste0(path_plot,"/AllTrips_Vessel_",unique(vessel_pp$BoatID)[i],"_ping",n,"min.jpg"),width = 1920, height = 1200)
    colours <- ifelse(vessel_p$HMM3_speedFilter == 1, "red", ifelse(vessel_p$HMM3_speedFilter == 2, "orange", "grey40"))
    plot(vessel_p$longitude, vessel_p$latitude, asp=1, pch=20, cex = .8, type="p", xlab="x", ylab="y", 
         col=colours, main = paste0('Vessel ', unique(vessel_pp$BoatID)[i], 'all trips - ping rate ', n, 'min'),
         xlim=c(min(vessel_p$longitude),max(vessel_p$longitude)))
    legend("bottomleft", 
           legend = c("1 - Fishing", "2 - Mix Transit/Fishing", "3 - Transit"), 
           col = c("red","orange", "grey40"), 
           pch = c(20,20), 
           bty = "n", 
           pt.cex = 1, 
           cex = 0.9, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.01, 0.01))
    dev.off()
    
    ### PLOT INDIVIDUAL TRACKS
    
    for (k in 1:length(unique(vessel_p$TripID))){
      png(filename=paste0(path_plot,"/Trips_Vessel_",unique(vessel_p$TripID)[k],"_ping",n,"min.jpg"),width = 1280, height = 800)
      idt = unique(vessel_p$TripID)[k]
      colours <- ifelse(vessel_p$HMM3_speedFilter[vessel_p$TripID == idt] == 1, "red", ifelse(vessel_p$HMM3_speedFilter[vessel_p$TripID == idt] == 2, "orange", "grey40"))
      plot(vessel_p$longitude[vessel_p$TripID == idt], vessel_p$latitude[vessel_p$TripID == idt], asp=1, pch=20, cex = .8, type="p", xlab="x", ylab="y", col=colours, main = paste('Vessel ', unique(vessel_pp$BoatID)[i], 'trip #',idt,' - ping rate ', n, 'min'),
           xlim=c(min(vessel_p$longitude[vessel_p$TripID == idt]),max(vessel_p$longitude[vessel_p$TripID == idt])))
      segments(vessel_p$longitude[vessel_p$TripID == idt][-length(vessel_p$longitude[vessel_p$TripID == idt])],vessel_p$latitude[vessel_p$TripID == idt][-length(vessel_p$latitude[vessel_p$TripID == idt])],vessel_p$longitude[vessel_p$TripID == idt][-1L],vessel_p$latitude[vessel_p$TripID == idt][-1L],col=colours)
      points(vessel_p$longitude[vessel_p$TripID == idt], vessel_p$latitude[vessel_p$TripID == idt], asp=1, pch=20, cex = .8, type="p", xlab="x", ylab="y", col=colours)
      legend("bottomleft", 
             legend = c("1 - Fishing", "2 - Mix Transit/Fishing", "3 - Transit"), 
             col = c("red","orange", "grey40"), 
             pch = c(20,20), 
             bty = "n", 
             pt.cex = 1, 
             cex = 1, 
             text.col = "black", 
             horiz = F , 
             inset = c(0.01, 0.01))
      dev.off()
    }
    
    vessel_p$Activity_classification <- ifelse(vessel_p$HMM3_speedFilter==1, 'FISHING', ifelse(vessel_p$HMM3_speedFilter==2, 'MIX_TRANSIT', 'TRANSIT'))
    write.csv(vessel_p, paste0(path_res, "/VMS+HMMresults_Vessel_",unique(vessel_pp$BoatID)[i],"_ping",n,"min.csv"), row.names=FALSE)
    vessel_result_ping <- rbind(vessel_result_ping, vessel_p)
  }
  
  write.csv(vessel_result_ping, paste0(cwd,'/ping',n,'min_resultClassification.csv'), row.names = FALSE)
  #rm(i, j, vessel, vessel_p)
  
  save(list_mod, file = paste0(cwd,'/models_RawOutput_ping',n,'min.Rdata'))
  vessel_result_all <- rbind(vessel_result_all, vessel_result_ping)
}

write.csv(vessel_result_all, paste0(path_folder_result,'resultClassification.csv'), row.names = FALSE)

t2 <- Sys.time()

timing <- round(difftime(t2,t1,'mins'),1)

txt <- paste0('Task done! It took ',as.numeric(timing), ' minutes to process. A new file with the result of the classification model has been written in ', path_folder_result, ' You can ignore the following warnings.')

message(txt)
}

HMM_classification_function(vessel = dat,path_folder_result = path_folder_result )
