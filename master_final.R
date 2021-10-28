#' Varalika Jain
#' Project:
#' How common ravens (Corvus corax) exploit anthropogenic food sources (AFSs)
#' through time and space in a semi-transformed, alpine environment

#------------------------------------------------------------------------------
####----(i) LOAD LIBRARIES----####
#' For better understanding of the script, information from the associated 
#' guides or vignettes of the following packages have also been included
library(move)
library(lubridate)
library(tidyr)
library(recurse)
library(dplyr)
library(plyr)
library(sp)
library(geosphere)
library(sf)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(ii) DOWNLOAD DATA FROM MOVEBANK (MOVE OBJECT)----####
#' Permissions required to download data
#' Input movebank login details
login <- movebankLogin()

#' Source data from GPS tagged ravens from 1st July 2017 to 15th March 2020
ravens <- getMovebankData(study="Corvus corax, Common Raven - Eastern Alps", 
                          login=login,
                          timestamp_start = "20170701000000000", 
                          timestamp_end = "20200315000000000",
                          removeDuplicatedTimestamps=TRUE)

#' Assigning the correct time zone
timestamps(ravens) <- with_tz(timestamps(ravens), tz="Europe/Vienna")

#' Check the timezone
head(timestamps(ravens))

#' Check individuals
levels(ravens@trackId)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(iii) INDIVIDUAL OCCURRENCE DISTRIBUTION----####
#' For each individual-season-year combination, we estimated the occurrence distribution
#' as the 95% utilization distribution (UD) with dynamic Brownian
#' Bridge Movement Models (dBBMM) from the “move” package
#' The 95% UD results from the code below were manually entered into a spreadsheet 
#' containing demographic data available under the file name 'raven_demog.csv'

####----(1) Burst the dataset by season: season function----####
#' Split the move data by each individual (creates a large list)
indv <- split(ravens)

#' Since there are a large number of individuals, only select few at a time as
#' the dBBMM function can take a while.
#' Selecting one or a range of individuals from the list
select_indv <- indv[c(1:3)]

#' Function attributing dates to each season and year
season <- function(dates){
  ifelse(dates < "2017-09-23 00:00:00", "sum 2017",
         ifelse (dates >= "2017-09-23 00:00:00" & dates < "2017-12-21 00:00:00" , "aut 2017",
                 ifelse (dates >= "2017-12-21 00:00:00" & dates< "2018-03-21 00:00:00" , "win 2018",
                         ifelse (dates >= "2018-03-21 00:00:00" & dates < "2018-06-21 00:00:00" , "spr 2018",
                                 ifelse (dates >= "2018-06-21 00:00:00" & dates < "2018-09-23 00:00:00" , "sum 2018",
                                         ifelse (dates >= "2018-09-23 00:00:00" & dates < "2018-12-21 00:00:00" , "aut 2018",
                                                 ifelse (dates >= "2018-12-21 00:00:00" & dates < "2019-03-21 00:00:00" , "win 2019",
                                                         ifelse (dates >= "2019-03-21 00:00:00" & dates < "2019-06-21 00:00:00" , "spr 2019",
                                                                 ifelse (dates >= "2019-06-21 00:00:00" & dates < "2019-09-23 00:00:00" , "sum 2019",    
                                                                         ifelse (dates >= "2019-09-23 00:00:00" & dates < "2019-12-21 00:00:00" , "aut 2019",
                                                                                 ifelse (dates >= "2019-12-21 00:00:00" & dates < "2020-03-21 00:00:00" , "win 2020", "spr 2020")))))))))))
}

#' Setting up variables as lists to store results
ssn = indv_burst = indv_burst_prj = dBBv_burst = burst_split = list()

#' Splitting the object by season, and excluding time lags over 5hrs
for (i in 1:length(select_indv)){
  #' To burst the track, the length of this vector has to be one shorter than 
  #' the number of coordinates.
  ssn[[i]] <- season(timestamps(select_indv[[i]][-1]))
  #' The track is bursted by supplying a vector (ssn) with the specific 
  #' variable associated with each location.
  indv_burst[[i]] <- burst(x=select_indv[[i]], f=ssn[[i]])
  
  #' Remove the variance of the segments corresponding to the large time gaps
  indv_burst_prj[[i]] <- spTransform(indv_burst[[i]], center=TRUE)
  dBBv_burst[[i]] <- brownian.motion.variance.dyn(indv_burst_prj[[i]], 
                                                  location.error=20, 
                                                  window.size=31, 
                                                  margin=11)
  #' Segments corresponding to the large time gaps are set to FALSE
  dBBv_burst[[i]]@interest[timeLag(indv_burst_prj[[i]], "mins")>300]<- FALSE
  burst_split[[i]] <- split(dBBv_burst[[i]])}

####----(2) Calculate the UD using dBBMM----####
#' Setting up variables as lists to store results
dbb_burst_cor = ud_cor = vol_95 = area_95 = fixes = days = rapply(burst_split,length, how = "list")

for (j in seq_along(burst_split)){
  for (k in seq_along(burst_split[[j]])){
    #' Use the 'dBMvariance' object to calculate the dBBMM corrected for 
    #' long time lags
    dbb_burst_cor[[j]][[k]] <- brownian.bridge.dyn(burst_split[[j]][[k]], 
                                                   window_size = 31, margin = 13, 
                                                   raster=100,location.error=20, 
                                                   ext=1.5, time.step= 1)
    #' 95% UD
    ud_cor[[j]][[k]] <- getVolumeUD(dbb_burst_cor[[j]][[k]])
    vol_95[[j]][[k]] <- ud_cor[[j]][[k]] <=.95
    area_95[[j]][[k]] <- sum(values(vol_95[[j]][[k]]))
    #' Number of fixes and days for the bursted sections of the track
    fixes[[j]][[k]] <- n.locs(burst_split[[j]][[k]])
    days[[j]][[k]] <- as.Date(last(burst_split[[j]][[k]]@timestamps)) - 
      as.Date(first(burst_split[[j]][[k]]@timestamps))
  }
}
#' area_95 has the values of interest

####----(3) Export the contour plots from (2)----####
#' 95% UD
#' Setting up variables as lists to store results
cont_95 = prj_95 = rapply(vol_95,length, how = "list")
#' Convert into contour
for (l in seq_along(vol_95)){
  for (m in seq_along(vol_95[[l]]))
    cont_95[[l]][[m]]<- rasterToContour(vol_95[[l]][[m]], levels =.95)
}
#' Transforming projection
for (n in seq_along(cont_95)){
  for (o in seq_along(cont_95[[n]]))
    prj_95[[n]][[o]]<- spTransform(cont_95[[n]][[o]], CRSobj = "+proj=longlat +datum=WGS84 +nodefs")
}
#' Exporting the files
for (p in seq_along(unlist(prj_95))){
  writeOGR(unlist(prj_95)[[p]], paste0("1_contour_95_" ,p ,".kml"), 
           paste0("1_contour_95_" ,p ,".kml"), driver = "KML")}

####----(4) Troubleshooting: Calculate the UD----####
#' For some individual-season combinations, the same extent parameter cannot be used to
#' calculate the UD and make plots, and so indexing may be required.
index_indv <- burst_split[[1]][1]

#' Setting up variables as lists to store results
dbb_burst_cor = ud_cor = vol_95 = area_95 = fixes = 
  days = rapply(index_indv,length, how = "list")

for (j in seq_along(index_indv)){
  #' Use the 'dBMvariance' object to calculate the dBBMM corrected for 
  #' long time lags
  dbb_burst_cor[[j]] <- brownian.bridge.dyn(index_indv[[j]], 
                                            window_size = 31, margin = 13, 
                                            raster=100,location.error=20, 
                                            ext=2.5, time.step= 1)
  #' 95% UD
  ud_cor[[j]] <- getVolumeUD(dbb_burst_cor[[j]])
  vol_95[[j]] <- ud_cor[[j]] <=.95
  area_95[[j]] <- sum(values(vol_95[[j]]))
  #' Number of fixes and days for the bursted sections of the track
  fixes[[j]] <- n.locs(index_indv[[j]])
  days[[j]] <- as.Date(last(index_indv[[j]]@timestamps)) - 
    as.Date(first(index_indv[[j]]@timestamps))
}

####----(5) Export contour plots from (4)----####
#' 95% UD
#' Setting up variables as lists to store results
cont_95 = prj_95 = rapply(vol_95,length, how = "list")
#' Convert into contour
for (l in seq_along(vol_95)){
  cont_95[[l]]<- rasterToContour(vol_95[[l]], levels =.95)
}
#' Transforming projection
for (n in seq_along(cont_95)){
  prj_95[[n]]<- spTransform(cont_95[[n]], CRSobj = "+proj=longlat +datum=WGS84 +nodefs")
}
#' Exporting the files
for (p in seq_along(unlist(prj_95))){
  writeOGR(unlist(prj_95)[[p]], paste0("1_contour_95_" ,p ,".kml"), paste0("1_contour_95_" ,p ,".kml"), driver = "KML")}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(iv) CONVERT MOVE OBJECT TO DATAFRAME----####
ravens_df <- as(ravens, "data.frame")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(v) MODIFY DATAFRAME FOR REQUIRED INFO----####
####----(1) Creating columns: season & year & date----####
#' Check to see if timestamps are still in local time
head(ravens_df$timestamps)

#' Create a season-year column first as season
#' function for season in iii.1
ravens_df$season <- season(ravens_df$timestamps)

#' Splitting column into season and year
ravens_df <- separate(ravens_df, season, into = c("season", "year"), sep = "\\s")

#' Checking the columns are now split into season and year correctly
head(ravens_df)

#' Change to factor
class(ravens_df$season)
ravens_df$season <- as.factor(ravens_df$season)
class(ravens_df$year)
ravens_df$year <- as.factor(ravens_df$year)

# Create a column with the date 
ravens_df$date <- as.Date(ravens_df$timestamps)

####----(2) Selecting required columns----####
head(ravens_df)
ravens_filter <- ravens_df %>% dplyr::select(c("local_identifier", "location_lat", 
                                               "location_long", "timestamps", 
                                               "season", "year", "date"))
head(ravens_filter)

#' Local identifier as factor
class(ravens_filter$local_identifier)
ravens_filter$local_identifier <- as.factor(ravens_filter$local_identifier)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(vi) ESTIMATING THE BUFFER SIZE FOR AFSs----####
#' To measure the number of AFSs visited by individuals and their probability 
#' of being at an AFS, we set a buffer around each AFS, and considered 
#' individuals inside the buffers to be actively using the site.

####----(1) Convert raven data into UTM projection---####
head(ravens_filter)

ravens_coords <- ravens_filter[c("location_long", "location_lat", "timestamps", "local_identifier")]
coordinates(ravens_coords) <- c("location_long", "location_lat")
proj4string(ravens_coords) <- CRS("+proj=longlat +datum=WGS84")
#' The units for the radius are those of the (x,y) coordinates 
#' (e.g., meters in the case of a UTM projection)
ravens_coords_utm <- spTransform(ravens_coords, CRS("+proj=utm +zone=33 ellps=WGS84"))
ravens_coords_utm <- as.data.frame(ravens_coords_utm)
ravens_coords_utm <- ravens_coords_utm[c(3,4,1,2)] #long, lat, timestamp, id
head(ravens_coords_utm)

####----(2) Convert raven data (UTM projection) into an sf object----####
ravens_sf <- st_as_sf(x = ravens_coords_utm, 
                      coords = c("location_long", "location_lat"),
                      crs = "+proj=utm +zone=33 +datum=WGS84")

####----(3) Read in the AFS location data----####
#' AFS locations were identified using clusters of GPS fixes (from tracked ravens)
#' near obvious anthropogenic structures
afs <- read.csv('AFSs.csv') #locs - locations
head(afs)

####----(4) Select for columns: name, long, lat and type of AFS----####
afs_latlong <- afs[c(1,3,2,6)]
head(afs_latlong)

####----(5) Convert location data into UTM projection----####
afs_coords <- afs_latlong[c("FID", "Longitude", "Latitude", "Category")]
coordinates(afs_coords) <- c("Longitude", "Latitude")
proj4string(afs_coords) <- CRS("+proj=longlat +datum=WGS84")
afs_coords_utm <- spTransform(afs_coords, CRS("+proj=utm +zone=33 +datum=WGS84"))

####----(6) Convert location coordinates into dataframe----####
afs_coords_df <- as.data.frame(afs_coords_utm@coords)
head(afs_coords_df)

####----(7) Using recurse to determine buffer/ radius size----####
#' Calculate the number of revisits at AFSs at different radii
afs_revisits = list()
for (i in 1:150) {
  afs_revisits$revisits[i] = getRecursionsAtLocations(x = ravens_coords_utm, 
                                                      locations = afs_coords_df,
                                                      radius = i,
                                                      threshold = 0) #threshold in hours
}

afs_revisits <- as.data.frame(afs_revisits$revisits)

####----(8) Plotting the revisits according to radius size----####
radii <- c(1:150)
titles <- paste(afs_coords$FID)

#' Create a folder called revisits
folder <- "revisits"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}

#' Store the PDF graphs in that folder
par(mar = c(1,1,1,1))
par(mfrow = c(1,1))
for (j in 1:ncol(afs_revisits)){
  pdf(file = paste0("revisits/", main = titles[j], ".pdf"))
  plot(radii,afs_revisits[j,], main = titles[j],
       ylab = "Total number of revisits", xlab = "Radii")
  dev.off()
}

#' Site-specific buffers were  selected based on the best estimate 
#' of the distance above which revisitations no longer
#' increased with radius size, or at the maximum of 150m.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(vii) GPS FIX & BUFFER INTERSECTION ----####
#' An individual was considered to be using an AFS if its GPS fix intersected 
#' with a buffer

####----(1) Convert location data (UTM projection) into an sf object----####
afs_sf <- st_as_sf(x = afs_coords_utm, 
                   coords = c("long", "lat"),
                   crs = "+proj=utm +zone=33 +datum=WGS84")

####----(2) Filtering locations by type of AFS: Wildpark, Hut, Refuse & Farm----####
afs_wildparks <- afs_sf %>% filter(Category == c("Wildpark")) #buffer 80 for wildparks
afs_huts <- afs_sf %>% filter(Category == c("Hut")) #buffer 80 for huts
afs_dumps <- afs_sf %>% filter (Category == "Refuse") #buffer 100 for compost + dumps
afs_farms <- afs_sf %>% filter (Category == "Farm") #buffer 80 for farm

####----(3) Get the geom of location values----####
geom_wildparks = st_geometry(afs_wildparks)
geom_huts = st_geometry(afs_huts)
geom_dumps = st_geometry(afs_dumps)
geom_farms = st_geometry(afs_farms)

####----(4) Set buffers----####
#' For game parks and huts, optimal buffer radii had a median of 80m 
#' (range: 40 – 150m). Refuse sites had larger radii of median 100m 
#' (range: 50 – 150m)

#' In UTM, so dist is in meters
buf1 <- st_buffer(geom_wildparks, dist = 80)
buf2 <- st_buffer(geom_huts, dist = 80)
buf3<- st_buffer(geom_dumps, dist = 100)
buf4 <- st_buffer(geom_farms, dist = 80)

#' Combining the buffers into one list (but list number[[x]] does not correspond to FIDs!!!)
buf <- c(buf1, buf2, buf3, buf4)
#' As FIDs are lost in this not so neat process, find the centroids again.
buf_centre <- st_centroid(buf)

####----(5) Calculating the intersection between fixes and buffers----####
#' For each GPS fix, calculate whether or not it intersects with any of the bufs
#' #' If it does not intersect, in the intersection column, paste " ", else "AFS"
intersect_dat <- ravens_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, buf))
  , location = if_else(is.na(intersection), " ", paste0("AFS"))) 

#'By turning intersection into factor, we get the foraging sites 
#'But these correspond to buffer numbers. 
#'need to refer to buf_centre and rename the numbers according to FID

####----(6) Renaming intersection according to FID----####
intersect_dat$intersection <- as.character(intersect_dat$intersection)

#' Wildparks
intersect_dat$intersection <- revalue(intersect_dat$intersection, c("1"="F1", "2"="F35", 
                                                                    "3"="F39"))
#' Huts 
intersect_dat$intersection <- revalue(intersect_dat$intersection , c("4"="F21", "5"="F26", 
                                                                     "6"="F27", "7"="F28", "8"="F29", "9"="F32",
                                                                     "10"="F34", "11"="F36", "12"="F37", "13"="F38",
                                                                     "14"="F40", "15"="F41", "16"="F42", "17"="F43", 
                                                                     "18"="F45", "19"="F46", "20"="F47", "21"="F48"))
#' Refuse
intersect_dat$intersection <- revalue(intersect_dat$intersection, c("22"="F02", "23"="F03",
                                                                    "24"="F04", "25"="F05", "26"="F06", "27"="F08",
                                                                    "28"="F10", "29"="F11", "30"="F13", "31"="F14",
                                                                    "32"="F15", "33"="F16", "34"="F17", "35"="F18",
                                                                    "36"="F20", "37"="F22", "38"="F24", "39"="F25"))
#' Farms
intersect_dat$intersection <- revalue(intersect_dat$intersection, c("40"="F12", "41"="F19",
                                                                    "42"="F23", "43"="F30", "44"="F31"))


####----(7) Convert to dataframe ----####
intersect_df <- as.data.frame(intersect_dat)
head(intersect_df)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(viii) MERGE DATAFRAMES WITH DEMOGRAPHIC DATA----####
####----(1) First combine intersection data with raven move data----####
ravens_filter_intersect <- merge.data.frame(intersect_df, ravens_filter, 
                                            by = c("local_identifier", "timestamps"))
#' Check dataframe
head(ravens_filter_intersect)

####----(2) Read in the demography datafile, which includes the 95%UD----####
demog <- read.csv('raven_demog.csv') 
head(demog)

#' id name as factor
class(demog$id_name)
demog$id_name <- as.factor(demog$id_name)

####----(3) Checking to see if all the individual names match----####
levels(ravens_filter_intersect$local_identifier)
levels(demog$id_name)
setdiff(demog$id_name, ravens_filter$local_identifier)

#' Change the names in the lif_hist file as it's smaller & faster to do
demog$id_name <- revalue(demog$id_name, c("Hilary" = "Hillary", "Ilo" = "Io",
                                          "T.Rex" = "T-Rex", "Victor.Peter" = "Victor-Peter", 
                                          "Xeres" = "Xerxes"))

#' Check that there are no more differences
setdiff(demog$id_name, ravens_filter_intersect$local_identifier)

####----(4) Merging the dataframe with demography data----####
names(demog)
names(ravens_filter_intersect)
ravens_demog <- merge.data.frame(demog, ravens_filter_intersect, by.x = c("id_name","season", "year"),
                                 by.y = c('local_identifier', "season", "year"),
                                 all.y = TRUE, sort = TRUE)

head(ravens_demog)
write.csv(ravens_demog, "ravens_demog.csv")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(ix) CALCULATE THE DAILY DISPLACEMENT----####
#' We measured maximum daily displacement as the distance between an individual’s 
#' furthest GPS fix relative to its first GPS fix within a day

####----(1) Select the location from the first GPS fix in the day from dataset----####
first_location <- ravens_demog %>%
  group_by(id_name, date) %>%
  arrange(timestamps) %>%
  slice(1L)
head(first_location)

#' Rename latitude and longitude columns
names(first_location)[names(first_location) == "location_lat"] <- "first_lat"
names(first_location)[names(first_location) == "location_long"] <- "first_long"
names(first_location)

####----(2) Calculate the distance between each fix relative to the first fix----####
#'Select the columns needed
first_location <- first_location %>% dplyr::select(c("id_name", "date", 
                                                     "first_long", "first_lat"))

#' With the first GPS location in the day in unique columns, merge the dataframes
head(first_location)
head(ravens_demog)

daily_dist <- merge(ravens_demog,first_location, by=c("id_name", "date"), all = TRUE)
head(daily_dist)

#' Convert to spatial data to use distVincentyEllipsoid function & set projection
first_coords <- subset(daily_dist, select = c("id_name", "date", "timestamps",
                                              "first_long", "first_lat"))
daily_coords <- subset(daily_dist, select = c("id_name", "date", "timestamps",
                                              "location_long", "location_lat"))

#' set the projection
coordinates(first_coords) <-  c("first_long", "first_lat")
proj4string(first_coords) <- CRS("+proj=longlat +zone=33 +datum=WGS84")

coordinates(daily_coords) <-  c("location_long", "location_lat")
proj4string(daily_coords) <- CRS("+proj=longlat +zone=33 +datum=WGS84")

#' Difference in distance in meters stored in a vector
displacement <- distVincentyEllipsoid(first_coords, daily_coords)
head(displacement)

####----(3) Add daily displacement as a column in daily_dist----####
daily_dist$displacement <- displacement
head(daily_dist)

####----(4) Merge ravens demog and daily dist ----####
#' Select the columns
daily_dist_subset <- subset(daily_dist, select = c("id_name", "date", "timestamps", 
                                                   "location_long", "location_lat",
                                                   "first_long", "first_lat","displacement"))

head(daily_dist)
#' Check the arrangement of the dataframes 
identical(daily_dist_subset[["location_lat"]], ravens_demog[["location_lat"]])
identical(daily_dist_subset[["location_long"]], ravens_demog[["location_long"]])
identical(daily_dist_subset[["timestamps"]], ravens_demog[["timestamps"]])

#' Arrange the dataframes before merging
daily_dist_arrange <- arrange(daily_dist_subset, id_name, date, timestamps)
ravens_demog_arrange <- arrange(ravens_demog, id_name, date, timestamps)
identical(daily_dist_arrange[["location_lat"]], ravens_demog_arrange[["location_lat"]])
identical(daily_dist_arrange[["location_long"]], ravens_demog_arrange[["location_long"]])
identical(daily_dist_arrange[["timestamps"]], ravens_demog_arrange[["timestamps"]])



#' Merge the arranged dataframes
ravens_demog_dist <- merge(daily_dist_arrange, ravens_demog_arrange, 
                           by = c("id_name", "date", "timestamps",
                                  "location_long", "location_lat"))
names(ravens_demog_dist)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
####----(x) CREATE MASTER DATAFRAME FOR ANALYSIS----####
####----(1) Calculate the total number of fixes for each id, season & year ----####
fix_total <- ravens_demog_dist %>% group_by(id_name, season, year) %>%
  tally()
fix_total <- as.data.frame(fix_total)
head(fix_total)
names(fix_total)[names(fix_total) == "n"] <- "fix_total"

####----(2) Calculate the number of tracking days for each id, season & year----####
days <- ravens_demog_dist %>% group_by(id_name, season, year, date) %>% tally()
head(days)
days <- days[,-5]
tracking_days <- days %>% group_by(id_name, season, year) %>% tally()
head(tracking_days)
tracking_days <- as.data.frame(tracking_days)
names(tracking_days)[names(tracking_days) == "n"] <- "tracking_days"

####----(3) Merge total fixes and tracking days ----####
fixes_days <- merge.data.frame(fix_total, tracking_days, by = c("id_name", "season", "year"))
head(fixes_days)

####----(4) Calculate the ratio of fixes by days for each id,season & year----####
fixes_days$fixes_by_days <- (fixes_days$fix_total/fixes_days$tracking_days)

####----(5) Merge fixes_days with ravens_demog for master file----####
master <- merge(ravens_demog_dist, fixes_days, by = c("id_name", "season", "year"), all = TRUE)
head(master)

####----(6) Export master, (i.e. master file) ----####
#' some individuals don't have area_95/ it couldn't be computed so we want to
#' exclude those rows 
master_final <- master[!is.na(master$area_95),]
master_final <- droplevels(master_final)
levels(master_final$id_name)
names(master_final)

#' remove geometry column 
master_final <- master_final[,-26]

#' Export file
write.csv(master_final, "master_final.csv")
#------------------------------------------------------------------------------

####----END----####