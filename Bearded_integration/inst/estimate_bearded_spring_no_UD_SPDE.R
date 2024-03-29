# estimate bearded seal abundance and spatial distribution using aerial survey data, CKMR, and 
# acoustics (omit UD for the moment).  Currently just spring, 2004-2018


load('prediction_grid_ice_kriged.RData')

# clean up a little
list_rm =ls()
list_rm = list_rm[-which(list_rm %in% c("grid_sf","Ice"))]
rm(list=list_rm)
rm(list_rm)

yr_last = 2021 #will want to change if updating Ice year range

library(sf)
grid_sf$depth = grid_sf$depth/100  #depth now in 100s of m's (numerical stability)
grid_sf$dist_land = grid_sf$dist_land/100 #now in 100s of km

easting_northing = st_coordinates(st_centroid(grid_sf))
grid_sf$easting = easting_northing[,1]
grid_sf$northing = easting_northing[,2]
grid_sf$northing = grid_sf$northing/2500000  #reduce extreme scale to prevent numerical issues
grid_sf$easting = grid_sf$easting/2500000  

which_Bering = which(grid_sf$Bering_2013_smoothed>=0)
which_Chukchi = which(grid_sf$Chukchi_smoothed>0)
which_Beaufort = which(grid_sf$JOBSS_smoothed>0)

#Ice = Ice
#n_season = ncol(Ice)
Seasons = c(1:ncol(Ice))
which_spring = which(Seasons%%5==2)
Ice = Ice[,which_spring]
n_season = ncol(Ice)  #now "n_season" is really the number of years, but keep consistency for when we put other seasons back in


#we'll set proportion land = 0 for all cells for simulation testing...will want to incorporate
# into final seal model though
grid_sf$land_cover = 0

N_ckmr =408651
N_ckmr_cv = 0.35
log_ckmr_se = sqrt(log(N_ckmr_cv^2+1))

# load acoustic data 
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
Acc_detect$season = 1*(Month %in% c(1:3))+2*(Month %in% c(4:5))+3*(Month %in% c(6:7))+4*(Month %in% c(8:9))+5*(Month %in% c(10:12))
#limit to spring
Acc_detect = Acc_detect[which(Acc_detect$season==2),]
#Acc_detect$tstep = (Year-2004)*5+Acc_detect$season
Year = lubridate::year(Acc_detect[,"detection_dt"])
Acc_detect$tstep=Year-2003  #2004 = time step 1
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

# make acoustic calling rates vary by latitude, ice presence
# call rate = exponential hazard
#Acc_det_sum$season = factor(Acc_det_sum$season,levels=as.character(1:5))
Acc_det_sum$ice2=Acc_det_sum$ice^2
Acc_DM = model.matrix(~latitude+ice+ice2,data=Acc_det_sum)

save.image('spring_no_UD_data.Rdata')
load('spring_no_UD_data.Rdata')

library( TMB )
library(Matrix)
library(mgcv)

TmbFile1 = "C:/Users/paul.conn/git/Bearded_integration/Bearded_integration/src/fit_bearded_real_noUD_SPDE"
compile(paste0(TmbFile1,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile1) )

#design matrix for spatial effects (omitting Ice) - we'll make this GAM style

GAM_data = data.frame(Dummy=rep(1,nrow(grid_sf)),dist_land=grid_sf$dist_land,
                      depth=grid_sf$depth, easting = grid_sf$easting, northing=grid_sf$northing)

gam_setup = gam(Dummy ~ s(depth, bs = "cs",k=4) + s(dist_land, bs = "cs",k=4) + 
                  s(easting, bs = "cs",k=6) + s(northing,bs="cs",k=6),
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
X_s = matrix(0,n_cells,4)  #use fixed effects set-up for debugging 
X_s[,1]=grid_sf$depth
X_s[,2]=grid_sf$depth^2
X_s[,3]=grid_sf$dist_land
X_s[,4]=grid_sf$dist_land^2

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

#set up SPDE basis using INLA
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
library(INLA)
Grid_locs = cbind(grid_sf$northing,grid_sf$easting) 
colnames(Grid_locs) = c("y","x")
mesh = inla.mesh.create( Grid_locs )
n_knots = mesh$n
Eta_index = mesh$idx$loc-1  #which spde REs to apply as random effects for each cell centroid
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]

#design matrices for Ice 
X_ice = array(0,dim=c(n_cells,3,n_season)) 
for(is in 1:n_season){
  cur_ice = data.frame("ice"=Ice[,is],"ice2"=Ice[,is]^2,"I_zero"=1*(Ice[,is]<0.01))
  X_ice[,,is]=model.matrix(~0+ice+ice2+I_zero,data=cur_ice)
}

wt_ud = 1  #above 1 increases ud influence
wt_ckmr = 1  #above 1 decreases ckmr influence
wt_aerial = 1 #above 1 decreaes aerial influence
wt_acc=0.0005  # above 1 increases acoustic influence; .0005 justified bassed on initial fits and chi-sq lack-of-fit adjustment
wt_spline =1
wt_spde=1

w_max = 1

low = -10
#high = abs(low)/4  #so low values will have SE's that make high counts imporobabe
log_2012_N = as.vector(log(grid_sf$Bering_2012_smoothed[which_Bering]))
log_2012_N[log_2012_N<low]=low #get rid of negative infinity
log_2013_N = as.vector(log(grid_sf$Bering_2013_smoothed[which_Bering]))
log_2013_N[log_2013_N<low]=low #get rid of negative infinity
log_2016_N = as.vector(log(grid_sf$Chukchi_smoothed[which_Chukchi]))
log_2016_N[log_2016_N<low]=low #get rid of negative infinity
log_2021_N = as.vector(log(grid_sf$JOBSS_smoothed[which_Beaufort]))
log_2021_N[log_2021_N<low]=low #get rid of negative infinity

log_2012_SE = sqrt(grid_sf$Bering_2012_smoothed_SE^2/grid_sf$Bering_2012_smoothed^2)[which_Bering]  #delta method
log_2013_SE = sqrt(grid_sf$Bering_2013_smoothed_SE^2/grid_sf$Bering_2013_smoothed^2)[which_Bering]  #delta method
log_2016_SE = sqrt(grid_sf$Chukchi_smoothed_SE^2/grid_sf$Chukchi_smoothed^2)[which_Chukchi]  #delta method
log_2021_SE = sqrt(grid_sf$JOBSS_smoothed_SE^2/grid_sf$JOBSS_smoothed^2)[which_Beaufort]  #delta method

#some SEs are infinity, since expectation in denom is zero
high = 2; 
log_2012_SE[log_2012_SE>high]=high
log_2013_SE[log_2013_SE>high]=high
log_2016_SE[log_2016_SE>high]=high
log_2021_SE[log_2021_SE>high]=high


#NOTE: I temporarily changed log_2012_VC_inv etc. so as not to have to repeat matrix inverse every sim
Data<- list("Acc_k"=as.matrix(Acc_det_sum[,c(2,4,6,8,10)]), 
            "Acc_n"=as.matrix(Acc_det_sum[,c(3,5,7,9,11)]),
            "Area_acc"=Area_acc,"X_acc"=Acc_DM,
            "Acc_tstep"=Acc_det_sum$tstep-1,
            "log_ckmr_N" = log(N_ckmr), "log_ckmr_se"=log_ckmr_se,
            #"UD_mean_adj" = UD_mean_adj, "mean_UD" = mean_UD, "W_st" = W_st,
            "which_Chukchi"=which_Chukchi-1, "which_Bering"=which_Bering-1, "which_Beaufort"=which_Beaufort-1,
            "log_2012_N" = log_2012_N, "log_2013_N"=log_2013_N,"log_2016_N"=log_2016_N,"log_2021_N"=log_2021_N,
            "log_2012_se" = log_2012_SE,"log_2013_se" = log_2013_SE,"log_2016_se" = log_2016_SE,"log_2021_se"=log_2021_SE,
            "N_2012" = exp(log_2012_N), "N_2013"=exp(log_2013_N), "N_2016"=exp(log_2016_N),"N_2021"=exp(log_2021_N),
            "X_ice"=X_ice,"Land_cover"=grid_sf$land_cover, "X_s"=X_s,
            "designMatrixForReport"=.bdiag(designMatrixForReport), "S" = S_combined, "Sdims"=Sdims,
            "Wts" = c(wt_ud, wt_ckmr, wt_aerial, wt_acc, wt_spline,wt_spde),"n_tsteps"=n_season,
            "Duty_lengths"=c(30,60,90,120,180),"est_acc"=1,
            "M0"=spde$M0,"M1"=spde$M1,"M2"=spde$M2,
            "Eta_index"=Eta_index,"options"=c(0,0),n_eta=n_knots,flag=1
)


Params = list("log_N"=log(500000),"Beta_s" = rep(0,ncol(X_s)),"Beta_ice" = rep(0,ncol(X_ice[,,1])),
              "Beta_acc"=rep(0,ncol(Acc_DM)),
              "log_lambda" = rep(0,length(Sdims)),"log_tau"=0, "log_kappa"=0,
               "Etainput_s"=rep(0,n_knots)
)


Params$Beta_acc[1] = -9  #needs to start low to prevent numerical issues

Random <- ("Etainput_s")
#Random <- c("Beta_s","log_lambda")
#Random <- NULL
Map <- list(log_lambda=factor(rep(NA,length(Sdims))))
#Map=NULL

Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_real_noUD_SPDE",silent=FALSE)
Obj$fn( Obj$par )
Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
Report = Obj$report()
Opt$message
Report1 = Report

save.image('bearded_SPDE_results.RData')

N_est[isim]=Report$N

Z_cat = c(Report$Z_score_2012,Report$Z_score_2013,Report$Z_score_2016)
prop_gt_90 = sum(abs(Z_cat)>1.64)/length(Z_cat)


# 
n_iter = 4
Wts_history = matrix(0,n_iter+1,5)
Wts_history[1,]=Data$Wts
for(i in 1:n_iter){
  Data$Wts[1]=Data$Wts[1]*301425/Report$jnll_comp[1]
  #Data$Wts[3]=Data$Wts[3]*5792/Report$jnll_comp[3]
  if(Data$Wts[3]<0.09 | Data$Wts[3]>0.11)Data$Wts[3] = Data$Wts[3]*(0.2+prop_gt_90)/0.3
  Data$Wts[4]=Data$Wts[4]*905/Report$chisq_acc
  Wts_history[i+1,]=Data$Wts
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_sim_all",silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
  Report = Obj$report()
  Z_cat = c(Report$Z_score_2012,Report$Z_score_2013,Report$Z_score_2016)
  prop_gt_90 = sum(abs(Z_cat)>1.64)/length(Z_cat)
  
}


#SD2=sdreport(Obj,bias.correct=FALSE)

#matrix(c(-9,3,3,-1,-1,0),ncol=1)

library(ggplot2)
png('Est_2012.png')
grid_sf$N = Report$Z_st[,9]
#grid_sf$N1[which(grid_sf$N1>1000)]=1000
ggplot(grid_sf)+geom_sf(aes(fill=N),color=NA)+theme(axis.text=element_blank())
dev.off()

png('Est_2018.png')
grid_sf$N = Report$Z_st[,15]
#grid_sf$N1[which(grid_sf$N1>1000)]=1000
ggplot(grid_sf)+geom_sf(aes(fill=N),color=NA)+theme(axis.text=element_blank())
dev.off()

grid_sf$z = NA
grid_sf$z[which_Bering]=Report$Z_score_2012
ggplot(grid_sf)+geom_sf(aes(fill=z),color=NA)
Which_high_z = which(Report$Z_score_2012>200)
Data$N_2012[Which_high_z]

ggplot(grid_sf)+geom_sf(aes(fill=Bering_2012),color=NA)

ggplot(grid_sf)+geom_sf(aes(fill=Bering_2013),color=NA)

#plot acoustic effects
Ice_conc = c(0:100)/100
Rate = exp(Report$Beta_acc[1]+Report$Beta_acc[3]*Ice_conc+Report$Beta_acc[4]*Ice_conc^2)
Prob60 = 1-exp(-60*Rate*100)
jpeg('Call_rate.jpeg')
par(mfrow=c(1,2))
plot(Ice_conc,Rate,type="l",lwd=2,ylab="Rate (per second, per capita)")
plot(Ice_conc,Prob60,type="l",lwd=2,ylab="Probability (>0 det in 60 secs for 100 seals")
dev.off()



save.image("bearded_SPDE_results_Mar2023.RData")
