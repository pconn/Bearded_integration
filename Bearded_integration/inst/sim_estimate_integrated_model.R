## build and test bearded seal integrated model using spring data (with aerial surveys)
# and fall data (w/o aerial surveys)

#assume model run 2004-2018
#UDs available for all years
# aerial surveys for BOSS in 2012, 2013; CHESS for 2016
# CKMR for total abundance
# passive acoustic moorings in years, seasons where actual data available

# load grid that will be used for estimation
load('prediction_grid_ice_kriged.RData')
set.seed(2020)

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


# generate UD distributions

# we'll want these to sum to one and approximate true relative abundance pattern
n_cells = nrow(Pi_s)
Pi_UD = UD_mean_adj = Pi_s
W_st = Pi_s
w_max = 10000  # we'll want to have a TMB option for this to prevent really high weights for zero predictions
mean_UD = mean(Pi_UD)  #we'll standardize to mean to prevent numerical errors - need to account for this in objective function though!
for(iseason in 1:n_season){
  LP = rnorm(1,1,0.2)*depth_eff(grid_sf$depth)+
    rnorm(1,1,0.2)*ice_eff(Ice[,iseason])+
    rnorm(1,1,0.2)*dist_land_eff(grid_sf$dist_land)+
    rnorm(1,1,0.2)*northing_eff(grid_sf$northing)
  Exp_LP = exp(LP)
  Pi_UD[,iseason]=Exp_LP/sum(Exp_LP)
  UD_mean_adj[,iseason] = Pi_UD[,iseason]/mean_UD
  #assign CV = 20% - with minimum value f
  Var_UD = (Pi_UD[,iseason]*0.2/mean_UD)^2 #constant CV 0f 0.2 transformed to new scale
  Var_UD[which(Var_UD<small)]=small
  W_st[,iseason] = 1/Var_UD
}

ggplot(grid_plot)+geom_sf(aes(fill=UD_mean_adj[,1]),color=NA)


# generate aerial survey abundance maps (2012, 2013 in Bering; 2016 in Chukchi; spring only)

which_Bering = which(grid_sf$Bering_2013_smoothed>=0)
which_Chukchi = which(grid_sf$Chukchi_smoothed>0)

n_aerial = 100
frac_sampled = 0.02 # assume 2% of each cell sampled
Row = c(1:nrow(grid_sf))
sampled_Bering_2012 = sample(Row[which_Bering],n_aerial,replace=TRUE)
sampled_Bering_2013 = sample(Row[which_Bering],n_aerial,replace=TRUE)
sampled_Chukchi = sample(Row[which_Chukchi],n_aerial,replace=TRUE)

Counts_2012 = rpois(n_aerial,N_s[sampled_Bering_2012,(2012-2004)*5+2]*frac_sampled)
Counts_2013 = rpois(n_aerial,N_s[sampled_Bering_2013,(2013-2004)*5+2]*frac_sampled)
Counts_2016 = rpois(n_aerial,N_s[sampled_Chukchi,(2016-2004)*5+2]*frac_sampled)

#fit GAM models to counts
nugget = 0.0000001 #small value to add to VC diagonal to make invertible
# 2012 Bering
grid_sf$ice = Ice[,(2012-2004)*5+2]
Aerial_data = grid_sf[sampled_Bering_2012,]
Aerial_data$Count = Counts_2012
Aerial_data$off = log(frac_sampled)
Bering_data=grid_sf[which_Bering,]
library(mgcv)
aerial_gam =gam(Count~offset(off)+s(depth)+s(dist_land)+s(ice)+s(northing),family='poisson',data=Aerial_data)
Bering_data$off=0 #make inference to whole cell
Bering_data$aerial_N = predict(aerial_gam,newdata=Bering_data,type="response")
VC_beta = vcov(aerial_gam)
LP = predict(aerial_gam,newdata=Bering_data,type="lpmatrix")
#variance covariance matrix of predictions
VC_pred_log = LP %*% VC_beta %*% t(LP) # this is on log scale
Bering_2012_VC = diag(Bering_data$aerial_N) %*% VC_pred_log %*% diag(Bering_data$aerial_N) #real scale
Bering_2012_N = Bering_data$aerial_N
Bering_2012_logVC = VC_pred_log  #in case we want to do data integration on log scale
diag(Bering_2012_logVC)=diag(Bering_2012_logVC)+nugget
SE_2012 = predict(aerial_gam,newdata=Bering_data,type="response",se.fit=TRUE)$se.fit
W_2012 = 1/SE_2012^2

# grid_plot$se2012 = NA
# grid_plot$se2012[which_Bering]=sqrt(diag(Bering_2012_VC))
# ggplot(grid_plot)+geom_sf(aes(fill=se2012),color=NA)
# 
# grid_plot$sampled2012 = NA
# grid_plot$sampled2012[sampled_Bering_2012]=1
# ggplot(grid_plot)+geom_sf(aes(fill=sampled2012),color=NA)
# 
# grid_plot$log_se2012 = NA
# grid_plot$log_se2012[which_Bering]=sqrt(diag(VC_pred_log))
# ggplot(grid_plot)+geom_sf(aes(fill=log_se2012),color=NA)
# 
# rep1 = exp(rmvnorm(1,log(Bering_2012_N),Bering_2012_logVC))
# grid_plot$rep1 = NA
# grid_plot$rep1[which_Bering]=rep1
# ggplot(grid_plot)+geom_sf(aes(fill=rep1),color=NA)
# 
# grid_plot$Aerial_2012 = NA
# #grid_plot$Aerial_2012[which_Bering]=Bering_2012_N
# grid_plot$Aerial_2012[which_Bering]=exp(Data$log_2012_N)
# ggplot(grid_plot)+geom_sf(aes(fill=Aerial_2012),color=NA)
# 
# grid_plot$Est_2012 = NA
# grid_plot$Est_2012[which_Bering] = exp(Report$log_E_N_2012)
# ggplot(grid_plot)+geom_sf(aes(fill=Est_2012),color=NA)


#2013 Bering
grid_sf$ice = Ice[,(2013-2004)*5+2]
Aerial_data = grid_sf[sampled_Bering_2013,]
Aerial_data$Count = Counts_2013
Aerial_data$off = log(frac_sampled)
Bering_data=grid_sf[which_Bering,]
library(mgcv)
aerial_gam =gam(Count~offset(off)+s(depth)+s(dist_land)+s(ice)+s(northing),family='poisson',data=Aerial_data)
Bering_data$off=0 #make inference to whole cell
Bering_data$aerial_N = predict(aerial_gam,newdata=Bering_data,type="response")
VC_beta = vcov(aerial_gam)
LP = predict(aerial_gam,newdata=Bering_data,type="lpmatrix")
#variance covariance matrix of predictions
VC_pred_log = LP %*% VC_beta %*% t(LP) # this is on log scale
Bering_2013_VC = diag(Bering_data$aerial_N) %*% VC_pred_log %*% diag(Bering_data$aerial_N) #real scale
Bering_2013_N = Bering_data$aerial_N
Bering_2013_logVC = VC_pred_log  #in case we want to do data integration on log scale
diag(Bering_2013_logVC)=diag(Bering_2013_logVC)+nugget
SE_2013 = predict(aerial_gam,newdata=Bering_data,type="response",se.fit=TRUE)$se.fit
W_2013 = 1/SE_2013^2

grid_plot$Aerial_2013 = NA
grid_plot$Aerial_2013[which_Bering]=Bering_2013_N
ggplot(grid_plot)+geom_sf(aes(fill=Aerial_2013),color=NA)

# grid_plot$Est_2013 = NA
# grid_plot$Est_2013[which_Bering] = exp(Report$log_E_N_2013)
# ggplot(grid_plot)+geom_sf(aes(fill=Est_2013),color=NA)


#2016 Chukchi 
grid_sf$ice = Ice[,(2016-2004)*5+2]
Aerial_data = grid_sf[sampled_Chukchi,]
Aerial_data$Count = Counts_2016
Aerial_data$off = log(frac_sampled)
Chukchi_data=grid_sf[which_Chukchi,]
library(mgcv)
aerial_gam =gam(Count~offset(off)+s(depth)+s(dist_land)+s(ice)+s(northing),family='poisson',data=Aerial_data)
Chukchi_data$off=0 #make inference to whole cell
Chukchi_data$aerial_N = predict(aerial_gam,newdata=Chukchi_data,type="response")
VC_beta = vcov(aerial_gam)
LP = predict(aerial_gam,newdata=Chukchi_data,type="lpmatrix")
#variance covariance matrix of predictions
VC_pred_log = LP %*% VC_beta %*% t(LP) # this is on log scale
Chukchi_VC = diag(Chukchi_data$aerial_N) %*% VC_pred_log %*% diag(Chukchi_data$aerial_N) #real scale
Chukchi_N = Chukchi_data$aerial_N
Chukchi_logVC = VC_pred_log  #in case we want to do data integration on log scale
diag(Chukchi_logVC)=diag(Chukchi_logVC)+nugget
SE_2016 = predict(aerial_gam,newdata=Chukchi_data,type="response",se.fit=TRUE)$se.fit
W_2016 = 1/SE_2016^2


grid_plot$Aerial_2016 = NA
grid_plot$Aerial_2016[which_Chukchi]=Chukchi_N
ggplot(grid_plot)+geom_sf(aes(fill=Aerial_2016),color=NA)

# grid_plot$Est_2016 = NA
# grid_plot$Est_2016[which_Chukchi] = exp(Report$log_E_N_2016)
# ggplot(grid_plot)+geom_sf(aes(fill=Est_2016),color=NA)

#grid_sf$sampled = NA
#grid_sf$sampled[sampled]=1
#plot(grid_sf[,"sampled"])

#Ones = matrix(1,nrow=1,ncol=nrow(VC_pred_aerial))
#N_se_aerial = sqrt(Ones %*% VC_pred_aerial %*% t(Ones))  # calculate SE for aerial survey total abundance (just a check)
#right now CV is ~10%

# generate CKMR estimate 
N_ckmr = N_tot  #just set to true value for now
CKMR_cv = 0.4 
log_ckmr_N = log(N_ckmr)
log_ckmr_se = sqrt(log(CKMR_cv^2+1))


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

#now simulate acoustic detections
for(idet in 1:n_det){
  N_det = Mooring_to_grid[Detect_to_mooring[idet],] %*% N_s[,Acc_det_sum$tstep[idet]]
  Tot_rate = N_det * Acc_det_sum$rate[idet] 
  #30 minute png files
  if(Acc_det_sum$effort_30s[idet]>0){
    p = 1-exp(-30*Tot_rate)
    Acc_det_sum$calls_30s[idet]=rbinom(1,Acc_det_sum$effort_30s[idet],p)
  }
  #60 minute png files
  if(Acc_det_sum$effort_60s[idet]>0){
    p = 1-exp(-60*Tot_rate)
    Acc_det_sum$calls_60s[idet]=rbinom(1,Acc_det_sum$effort_60s[idet],p)
  }
  #90 minute png files
  if(Acc_det_sum$effort_90s[idet]>0){
    p = 1-exp(-90*Tot_rate)
    Acc_det_sum$calls_90s[idet]=rbinom(1,Acc_det_sum$effort_90s[idet],p)
  }  
  #120 minute png files
  if(Acc_det_sum$effort_120s[idet]>0){
    p = 1-exp(-120*Tot_rate)
    Acc_det_sum$calls_120s[idet]=rbinom(1,Acc_det_sum$effort_120s[idet],p)
  }  
  #180 minute png files
  if(Acc_det_sum$effort_180s[idet]>0){
    p = 1-exp(-180*Tot_rate)
    Acc_det_sum$calls_180s[idet]=rbinom(1,Acc_det_sum$effort_180s[idet],p)
  }  
}

save.image('sim_data.RData')

load('sim_data.RData')
#### Fit TMB models
library( TMB )
library(Matrix)
library(mgcv)

TmbFile1 = "C:/Users/paul.conn/git/Bearded_integration/Bearded_integration/src/fit_bearded_sim"
compile(paste0(TmbFile1,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile1) )

#design matrix for spatial effects (omitting Ice) - we'll make this GAM style

GAM_data = data.frame(Dummy=rep(1,nrow(grid_sf)),dist_land=grid_sf$dist_land,
                      depth=grid_sf$depth, easting = grid_sf$easting, northing=grid_sf$northing)
          
gam_setup = gam(Dummy ~ s(depth, bs = "cs",k=4) + s(dist_land, bs = "cs",k=4) + 
                  s(easting, bs = "cs",k=5) + s(northing,bs="cs",k=5),
                data = GAM_data,fit=FALSE)
S_depth = gam_setup$smooth[[1]]$S[[1]]
S_dist_land = gam_setup$smooth[[2]]$S[[1]]
S_easting = gam_setup$smooth[[3]]$S[[1]]
S_northing = gam_setup$smooth[[4]]$S[[1]]

S_list = list(S_depth,S_dist_land,S_easting,S_northing)

S_combined = .bdiag(S_list)         # join S's in sparse matrix
S_combined = as(S_combined,"generalMatrix")
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

#X_s = gam_setup$X[,-1]
X_s = matrix(0,n_cells,6)  #use fixed effects set-up for debugging 
Beta_s_true = c(Depth_eff,Land_eff,Northing_eff)
X_s[,1]=grid_sf$depth
X_s[,2]=grid_sf$depth^2
X_s[,3]=grid_sf$dist_land
X_s[,4]=grid_sf$dist_land^2
X_s[,5]=grid_sf$northing
X_s[,6]=grid_sf$northing^2

#For report, used for constructing plots----
depth = seq(min(GAM_data$depth),max(GAM_data$depth),by = (max(GAM_data$depth)-min(GAM_data$depth))/100)
dist_land = seq(min(GAM_data$dist_land),max(GAM_data$dist_land),by = (max(GAM_data$dist_land)-min(GAM_data$dist_land))/100)
easting = seq(min(GAM_data$easting),max(GAM_data$easting),by = (max(GAM_data$easting)-min(GAM_data$easting))/100)
northing = seq(min(GAM_data$northing),max(GAM_data$northing),by = (max(GAM_data$northing)-min(GAM_data$northing))/100)

depthReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(depth))
landReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(dist_land))
eastingReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(easting))
northingReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(northing))

designMatrixForReport = list(depthReport,landReport,eastingReport,northingReport)

#design matrices for Ice 
X_ice = array(0,dim=c(n_cells,2,n_season)) 
for(is in 1:n_season){
  cur_ice = data.frame("ice"=Ice[,is],"ice2"=Ice[,is]^2)
  X_ice[,,is]=model.matrix(~0+ice+ice2,data=cur_ice)
}

wt_ud = 1  #above 1 increases ud influence
wt_ckmr = 1  #above 1 decreases ckmr influence
wt_aerial = 1 #above 1 decreaes aerial influence
wt_acc=1  # above 1 increases acoustic influence
wt_spline =1

w_max = 1
W_st[which(W_st>w_max)]=w_max  #limit influence of zero UD predictions
W_2012[which(W_2012>w_max)]=w_max  #limit influence of zero aerial predictions
W_2013[which(W_2013>w_max)]=w_max
W_2016[which(W_2016>w_max)]=w_max


Data<- list("Acc_k"=as.matrix(Acc_det_sum[,c(2,4,6,8,10)]), 
            "Acc_n"=as.matrix(Acc_det_sum[,c(3,5,7,9,11)]),
            "Area_acc"=Area_acc,"X_acc"=DM,
            "Acc_tstep"=Acc_det_sum$tstep-1,
            "log_ckmr_N" = log(N_ckmr), "log_ckmr_se"=log_ckmr_se,
            "UD_mean_adj" = UD_mean_adj, "mean_UD" = mean_UD, "W_st" = W_st,
            "which_Chukchi"=which_Chukchi-1, "which_Bering"=which_Bering-1, #since C++ vector indexes start at 0
            "log_2012_N" = log(Bering_2012_N),
            "log_2012_VC_inv" = solve(Bering_2012_logVC), "log_2013_N" = log(Bering_2013_N),
            "log_2013_VC_inv" = solve(Bering_2013_logVC), "log_2016_N" = log(Chukchi_N),
            "log_2016_VC_inv" = solve(Chukchi_logVC), "log_2012_se" = sqrt(diag(Bering_2012_logVC)),
            "log_2013_se" = sqrt(diag(Bering_2013_logVC)),"log_2016_se" = sqrt(diag(Chukchi_logVC)),
            "N_2012" = Bering_2012_N, "N_2013"=Bering_2013_N, "N_2016"=Chukchi_N,
            "W_2012" = W_2012, "W_2013"=W_2013, "W_2016"=W_2016,
            "X_ice"=X_ice,"Land_cover"=grid_sf$land_cover, "X_s"=X_s,
            "designMatrixForReport"=.bdiag(designMatrixForReport), "S" = S_combined, "Sdims"=Sdims,
            "Wts" = c(wt_ud, wt_ckmr, wt_aerial, wt_acc, wt_spline),"n_tsteps"=n_season,
            "MVN" = 0, "Duty_lengths"=c(30,60,90,120,180),"est_acc"=1
)

# Params = list("log_N"=log(500000),"Beta_s" = rep(0,sum(Sdims)),"Beta_ice" = rep(0,ncol(X_ice[,,1])),
#               "Beta_acc"=rep(0,ncol(DM)),"log_lambda" = rep(0,length(Sdims))
# )

Params = list("log_N"=log(500000),"Beta_s" = rep(0,length(Beta_s_true)),"Beta_ice" = rep(0,ncol(X_ice[,,1])),
             "Beta_acc"=rep(0,ncol(DM)),"log_lambda" = rep(0,length(Sdims))
)

 # Params = list("log_N"=log(500000),"Beta_s" = Beta_s_true,"Beta_ice" = Beta_ice_true,
 #              "Beta_acc"=rep(0,ncol(DM)),"log_lambda" = rep(0,length(Sdims))
 # )

#Params = list("log_N"=log(500000),"Beta_s" = Beta_s_true,"Beta_ice" = Beta_ice_true,
#              "Beta_acc"=Beta_acc_true,"log_lambda" = rep(0,length(Sdims))
#)

Params$Beta_acc[1] = -9  #needs to start low to prevent numerical issues


#Random <- c("Beta_s","log_lambda")
Random <- NULL
#Data$est_acc=0
#Map <- list(Beta_acc=factor(rep(NA,ncol(DM))),log_lambda=factor(rep(NA,length(Sdims))))
Map <- list(log_lambda=factor(rep(NA,length(Sdims))))
#Map <- list(log_N=factor(NA),Beta_s=factor(rep(NA,length(Beta_s_true))),Beta_ice=factor(rep(NA,length(Beta_ice_true))),
#            log_lambda=factor(rep(NA,length(Sdims))))

Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_sim",silent=FALSE)
Obj$fn( Obj$par )
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
Report = Obj$report()
Opt$message
Opt2 = optim( par=Obj$par, fn=Obj$fn, gr=Obj$gr, method="BFGS",control=list(trace=1, maxit=10000))         #
Report2 = Obj$report()
# Params = list("log_N"=log(Report$N),"Beta_s" = Report$Beta_s,"Beta_ice" = Report$Beta_ice,
#               "Beta_acc"=Report$Beta_acc,"log_lambda" = rep(0,length(Sdims))
# )
# Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_sim",silent=FALSE)
# Obj$fn( Obj$par )
# Opt3 = optim( par=Obj$par, fn=Obj$fn, gr=Obj$gr, method="BFGS",control=list(trace=1, maxit=10000))         #
# Report3 = Obj$report()
SD=sdreport(Obj,bias.correct=FALSE)


n_iter = 4
Wts_history = matrix(0,n_iter+1,5)
Wts_history[1,]=Data$Wts
for(i in 1:n_iter){
  Data$Wts[1]=Data$Wts[1]*301425/Report$jnll_comp[1]
  Data$Wts[3]=Data$Wts[3]*5792/Report$jnll_comp[3]
  Data$Wts[4]=Data$Wts[4]*902/Report$chisq_acc
  Wts_history[i+1,]=Data$Wts
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_sim",silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
  Report = Obj$report()
}


SD2=sdreport(Obj,bias.correct=FALSE)

#matrix(c(-9,3,3,-1,-1,0),ncol=1)

library(ggplot2)
grid_plot$N1 = Report$Z_st[,1]
ggplot(grid_plot)+geom_sf(aes(fill=N1),color=NA)


