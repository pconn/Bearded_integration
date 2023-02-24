## build and test bearded seal integrated model using spring data (with aerial surveys)
# and fall data (w/o aerial surveys)

#assume model run 2004-2018
#UDs available for all years
# aerial surveys for BOSS in 2012, 2013; CHESS for 2016
# CKMR for total abundance
# passive acoustic moorings in years, seasons where actual data available

# load grid that will be used for estimation
load('prediction_grid_ice_kriged.RData')

# clean up a little
list_rm =ls()
list_rm = list_rm[-which(list_rm %in% c("grid_sf","Ice"))]
rm(list=list_rm)
rm(list_rm)

yr_last = 2018 #will want to change if updating Ice year range

#### 1) generate underlying abundance intensity
library(sf)
#ice_2012 = Ice[,(2013-2004)*5 + 2]
#grid_sf$ice = ice_2012

#we'll make distribution conditional on depth, ice, distance from land,

depth_eff <- function(depth){
  -10*depth - 10*depth^2   #quadratic function maximized at -50 m
}

ice_eff <- function(ice){  #maximized at ice=0.7; weak effect
  2.8*ice - 2*ice^2
}

dist_land_eff <- function(dist_land){
  5*dist_land - 5*dist_land^2
}

grid_sf$depth = grid_sf$depth/100  #depth now in 100s of m's (numerical stability)
grid_sf$dist_land = grid_sf$dist_land/100 #now in 100s of km

easting_northing = st_coordinates(st_centroid(grid_sf))
grid_sf$easting = easting_northing[,1]
grid_sf$northing = easting_northing[,2]
grid_sf$northing = grid_sf$northing/2500000  #reduce extreme scale to prevent numerical issues
grid_sf$easting = grid_sf$easting/2500000  
northing_eff <- function(northing){ #make 63 degrees optimal
  -23*northing - 10 * northing^2
}

N_tot = 500000
Pi_s = Ice  #proportion of abundance in each grid cell (columns sum to one)
n_season = ncol(Ice)

Depth_eff = c(-10,-10)
Ice_eff = c(2.8,-2)
Land_eff = c(5, -5)
Northing_eff = c(-23, -10)

small = 0.00000001  #this will be added to Pi to try to prevent -inf loglik so probably want to add here too for debugging

#we'll set proportion land = 0 for all cells for simulation testing...will want to incorporate
# into final seal model though
grid_sf$land_cover = 0

for(iseason in 1:n_season){
  Pi_s[,iseason] = (1-grid_sf$land_cover)*exp(depth_eff(grid_sf$depth)+ice_eff(Ice[,iseason])+
                                                dist_land_eff(grid_sf$dist_land)+northing_eff(grid_sf$northing))+small
  Pi_s[,iseason]=Pi_s[,iseason]/sum(Pi_s[,iseason])
}
N_s = N_tot * Pi_s

Beta_s_true = c(Depth_eff,Land_eff,Northing_eff)
Beta_ice_true = Ice_eff

grid_plot = grid_sf
grid_plot$N_1 = N_s[,1]
library(ggplot2)
ggplot(grid_plot)+geom_sf(aes(fill=N_1),color=NA)


# generate acoustic data ... we'll load in actual mooring locations and try to 
# pattern data availability after our actual data!
library("RPostgreSQL")
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))
moorings_sf <- sf::st_read(con, 
                           query = "SELECT * FROM acoustics.geo_moorings", 
                           geometry_column = "geom")
Acc_detect <- dplyr::tbl(con, dplyr::sql("SELECT * FROM acoustics.tbl_detections"))
Acc_detect <- data.frame(Acc_detect)
Acc_detect = Acc_detect[which(lubridate::year(Acc_detect[,"detection_dt"])<=yr_last),]

Month = lubridate::month(Acc_detect[,"detection_dt"])
Year = lubridate::year(Acc_detect[,"detection_dt"])
Acc_detect$season = 1*(Month %in% c(1:3))+2*(Month %in% c(4:5))+3*(Month %in% c(6:7))+4*(Month %in% c(8:9))+5*(Month %in% c(10:12))
Acc_detect$tstep = (Year-2004)*5+Acc_detect$season
moorings_sf = st_transform(moorings_sf,st_crs(grid_sf))

#take out any moorings not within estimation grid
InGrid = st_within(moorings_sf,st_union(grid_sf),sparse=FALSE) 
moorings_sf = moorings_sf[which(InGrid==TRUE),]
Acc_detect = Acc_detect[which(Acc_detect$mooring_site_id %in% unique(moorings_sf$mooring_site_id_full)),]

#take out any moorings without acoustic detections
Unique_moorings_det = unique(Acc_detect$mooring_site_id)
moorings_sf = moorings_sf[which(moorings_sf$mooring_site_id_full %in% Unique_moorings_det),]

#combine data from same mooring, time step combo
Mooring_season_ID = paste(Acc_detect$mooring_site_id,Acc_detect$tstep)
Unique_MS = unique(Mooring_season_ID)
n_det = length(Unique_MS)
Acc_det_sum = data.frame(matrix(0,n_det,13))
colnames(Acc_det_sum)=c("mooring_site_id","calls_30s","effort_30s","calls_60s","effort_60s","calls_90s","effort_90s",
                        "calls_120s","effort_120s","calls_180s","effort_180s","season","tstep")  
for(ims in 1:n_det){
  Which_MS_ID=which(Mooring_season_ID==Unique_MS[ims])
  Acc_det_sum$mooring_site_id[ims]=Acc_detect[Which_MS_ID[1],"mooring_site_id"]
  Acc_det_sum$tstep[ims]=Acc_detect[Which_MS_ID[1],"tstep"]
  Acc_det_sum$season[ims]=Acc_detect[Which_MS_ID[1],"season"]
  Acc_det_sum[ims,c(2:11)]=colSums(Acc_detect[Which_MS_ID,8:17])
}

# put together a look up list for which mooring goes with which seasonal detection summary; attach latitude to detections
n_det = nrow(Acc_det_sum)
n_moor = nrow(moorings_sf)
Detect_to_mooring = rep(0,n_det)
Acc_det_sum$latitude = 0
for(imoor in 1:n_moor){
  cur_moor = moorings_sf$mooring_site_id_full[imoor]
  Which_cur = which(Acc_det_sum$mooring_site_id==cur_moor)
  Detect_to_mooring[Which_cur]=imoor
  Acc_det_sum[Which_cur,"latitude"] = as.numeric(moorings_sf[imoor,]$latitude)
}

# mapping from mooring location to analysis grid - integrate bivariate normal over grid cells that are 
# adjacent to the one including the mooring....we'll use the polyCube package for this
n_cells = nrow(grid_sf)
Mooring_to_grid = matrix(0,n_moor,n_cells)
library(polyCub)
intrfr <- function (R, sigma = 8000)  # the sigma here is designed to allow some (but very low) detection at 25km from mooring [max range reported by Cleator et al. 2011]
{
  (1 - exp(-R^2/2/sigma^2))/2/pi
}

#construct adjacency matrix - first identify the 'home' cell of each mooring
Home = rep(0,n_moor)
for(imoor in 1:n_moor) Home[imoor] = which(st_within(moorings_sf[imoor,],grid_sf,sparse=FALSE)==1)
Adj = st_is_within_distance(grid_sf[Home,],grid_sf,dist=1,sparse=FALSE)

for(imoor in 1:n_moor){
  center_pt = st_coordinates(moorings_sf[imoor,])
  colnames(center_pt)=c("x","y")
  Which_calc = which(Adj[imoor,]>0)  #cells within 25km of fixed mooring
  for(icell in 1:length(Which_calc)){
    Mooring_to_grid[imoor,Which_calc[icell]] = polyCub.iso(as(grid_sf[Which_calc[icell],],"Spatial"),intrfr=intrfr,center=center_pt)
  }
}

#each row of Area_acc gives the proportion of animals in each grid cell that are detectable by the associated mooring
Area_acc = matrix(0,n_det,n_cells)
for(idet in 1:n_det)Area_acc[idet,]=Mooring_to_grid[Detect_to_mooring[idet],]  

#attach ice to acoustic detections
Acc_det_sum$ice = 0
for(idet in 1:n_det){
  Acc_det_sum$ice[idet] = Ice[Home[Detect_to_mooring[idet]],Acc_det_sum$tstep[idet]]
}

#now we should be able to multiply [Mooring_to_grid %*% N_s[,iseason]] to get at expected number of detectable seals

# make acoustic calling rates vary by season, ice presence

#latitude runs from 56.9 to 72.6  we won't use this for simulations but will probably
#want to set up a season * latitude interaction for ultimate integrated model

# call rate = exponential hazard
Acc_det_sum$season = factor(Acc_det_sum$season,levels=as.character(1:5))
DM = model.matrix(~season+ice,data=Acc_det_sum)
Par = Beta_acc_true = matrix(c(-9,3,3,-1,-1,0),ncol=1)   #intercept, season2 eff ... season 5 eff, ice 
Acc_det_sum$rate = exp(DM %*% Par)  #log link on hazard rate
hist(1-exp(-Acc_det_sum$rate*60))  #this shows the implied probability of a single animal in detection range
#having a call detected in a 60 sec png file


save.image("bearded_sim_prep.RData")
