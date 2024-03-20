# estimate bearded seal abundance and spatial distribution using aerial survey data, CKMR, and 
# acoustics (omit UD for the moment).  Currently just spring, 2004-2018


load('prediction_grid_ice_kriged_2season.RData')

# clean up a little
list_rm =ls()
list_rm = list_rm[-which(list_rm %in% c("grid_sf","Ice","Ice_mo"))]
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

Seasons = c(1:ncol(Ice))
n_season = ncol(Ice)  #now "n_season" is really the number of years * # of seasons


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
Year = lubridate::year(Acc_detect[,"detection_dt"])
Acc_detect$season = 1*(Month %in% c(1:5,12))+2*(Month %in% c(6:11))
Acc_detect$Month = Month #factor(as.character(Month),levels=as.character(c(1:12)))
Acc_detect$tstep = (Year-2004)*2+Acc_detect$season
Acc_detect$tstep[which(Month==12)]=Acc_detect$tstep[which(Month==12)]+2 #need to put december with next year's time step
Acc_detect$month_step = (Year-2004)*12+Acc_detect$Month

#Acc_detect = Acc_detect[-which(Month==12 & Year==2020),]  # don't need this - there are no records from Dec 2020 (they go through Oct 2019)
moorings_sf = st_transform(moorings_sf,st_crs(grid_sf))

#take out any moorings not within estimation grid
InGrid = st_within(moorings_sf,st_union(grid_sf),sparse=FALSE) 
moorings_sf = moorings_sf[which(InGrid==TRUE),]
Acc_detect = Acc_detect[which(Acc_detect$mooring_site_id %in% unique(moorings_sf$mooring_site_id_full)),]

#take out any moorings without processed recordings
Unique_moorings_det = unique(Acc_detect$mooring_site_id)
moorings_sf = moorings_sf[which(moorings_sf$mooring_site_id_full %in% Unique_moorings_det),]


#produce table of number of files per year by mooring
Acc_df = as.data.frame(Acc_detect)
Acc_df$year = lubridate::year(Acc_df$detection_dt)
library(doBy)
Acc_table = summaryBy(num_png_with_call~mooring_site_id+year,data=Acc_df,FUN=sum)

#combine data from same mooring, time step combo
# Mooring_season_ID = paste(Acc_detect$mooring_site_id,Acc_detect$tstep)
# Unique_MS = unique(Mooring_season_ID)
# n_det = length(Unique_MS)
# Acc_det_sum = data.frame(matrix(0,n_det,14))
Mooring_month_ID = paste(Acc_detect$mooring_site_id,Acc_detect$month_step)
Unique_MS = unique(Mooring_month_ID)
n_det = length(Unique_MS)
Acc_det_sum = data.frame(matrix(0,n_det,15))

colnames(Acc_det_sum)=c("mooring_site_id","calls_30s","effort_30s","calls_60s","effort_60s","calls_90s","effort_90s",
                        "calls_120s","effort_120s","calls_180s","effort_180s","season","tstep","month","month_step")  
for(ims in 1:n_det){
  Which_MS_ID=which(Mooring_month_ID==Unique_MS[ims])
  Acc_det_sum$mooring_site_id[ims]=Acc_detect[Which_MS_ID[1],"mooring_site_id"]
  Acc_det_sum$tstep[ims]=Acc_detect[Which_MS_ID[1],"tstep"]
  Acc_det_sum$season[ims]=Acc_detect[Which_MS_ID[1],"season"]
  Acc_det_sum$month[ims]=Acc_detect[Which_MS_ID[1],"Month"]
  Acc_det_sum$month_step[ims]=Acc_detect[Which_MS_ID[1],"month_step"]
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
  Acc_det_sum$ice[idet] = Ice_mo[Home[Detect_to_mooring[idet]],Acc_det_sum$tstep[idet]]
}

#now we should be able to multiply [Mooring_to_grid %*% N_s[,iseason]] to get at expected number of detectable seals

# make acoustic calling rates vary by latitude, ice presence
# call rate = exponential hazard
Acc_det_sum$season = factor(Acc_det_sum$season,levels=as.character(1:2))
Acc_det_sum$season2[which(Acc_det_sum$month %in% c(4:6))]
Acc_det_sum$ice2=Acc_det_sum$ice^2
Acc_det_sum$month=factor(as.character(Acc_det_sum$month),levels=as.character(c(1:12)))
Acc_det_sum$season2 = as.character(Acc_det_sum$month)
Acc_det_sum$season2[which(as.character(Acc_det_sum$month)%in%as.character(c(7:12)))]="open"
Acc_det_sum$season2 = factor(Acc_det_sum$season2,levels=unique(Acc_det_sum$season2))
Acc_DM = model.matrix(~season2+ice+ice2,data=Acc_det_sum)

load("logUD_output_2seasons.RData")

#format UD data
logUD = logUD_output$logUD
Pi_UD = exp(logUD_output$logUD)
Var_logUD = logUD_output$SE^2 
for(i in 1:n_season)Pi_UD[,i]=Pi_UD[,i]/sum(Pi_UD[,i])
UD_mean_adj = Pi_UD
W_st = Pi_UD
w_max = 10000  # we'll want to have a TMB option for this to prevent really high weights for zero predictions
mean_UD = mean(Pi_UD)  #we'll standardize to mean to prevent numerical errors - need to account for this in objective function though!
Index = c(1:n_cells)
for(iseason in 1:n_season){
  UD_mean_adj[,iseason] = Pi_UD[,iseason]/mean_UD
  Cur_exp = exp(logUD[,iseason])
  sumSq = (sum(Cur_exp))^2
  Grad_trans = matrix(0,n_cells,n_cells)
  for(is in 1:n_cells){
    Grad_trans[is,is]=Cur_exp[is]*sum(Cur_exp[-is])
    Off_diags = Index[-is]
    for(is2 in 1:(n_cells-1)){
      Grad_trans[is,Off_diags[is2]]= -Cur_exp[is]*Cur_exp[Off_diags[is2]]
    }
  }
  Grad_trans = Grad_trans /sumSq
  Var_UD = 1/mean_UD^2 * Grad_trans %*% diag(Var_logUD[,iseason]) %*% t(Grad_trans)
  W_st[,iseason] = 1/diag(Var_UD)
  if(sum(W_st[,iseason]>w_max)>0)W_st[which(W_st[,iseason]>w_max),iseason]=w_max
}

save.image('bearded_TMB_inputs_2season_Jan2024.RData')

load('bearded_TMB_inputs_2season_Jan2024.RData')

library( TMB )
library(Matrix)
library(mgcv)

TmbFile1 = "C:/Users/paul.conn/git/Bearded_integration/Bearded_integration/src/fit_bearded_betabinom_SPDErr"
compile(paste0(TmbFile1,".cpp"),"-O1 -g",DLLFLAGS="") 

dyn.load( dynlib(TmbFile1) )


#design matrix for spatial effects (omitting Ice) - we'll make this GAM style
grid_sf$depth[which(grid_sf$depth>0)]=0
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
X_s = matrix(0,n_cells,2)  #use fixed effects set-up for debugging 
X_s[,1]=grid_sf$depth
X_s[,2]=grid_sf$depth^2
#X_s[,3]=grid_sf$dist_land
#X_s[,4]=grid_sf$dist_land^2

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
# using a predictive process here; 200 knots.
library(INLA)
Grid_locs = cbind(grid_sf$northing,grid_sf$easting) 
knots <- stats::kmeans(x = Grid_locs, centers = 200)  # # knots will be more after triangulation
loc_centers <- knots$centers
mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
spde <- INLA::inla.spde2.matern(mesh)  
A <- INLA::inla.spde.make.A(mesh, loc = as.matrix(Grid_locs)) #prediction matrix

#design matrices for Ice - diff effects for each season
Ice[which(Ice>1.0)]=1.0
X_ice = array(0,dim=c(n_cells,6,n_season)) 
for(is in 1:n_season){
  cur_seas = factor(as.character((is %% 2)+1),levels=c("1","2")) #actual season!
  cur_ice = data.frame("ice"=Ice[,is],"ice2"=Ice[,is]^2,"I_zero"=1*(Ice[,is]<0.01),"season"=cur_seas)
  X_ice[,,is]=model.matrix(~0+ice:season+ice2:season+I_zero:season,data=cur_ice)
}

wt_ud = 1/13  #above 1 increases ud influence
wt_ckmr = 1 #1/10 #sqrt(1/4020)  #above 1 decreases ckmr influence
wt_aerial = 1 #above 1 decreaes aerial influence
wt_acc=1  # above 1 increases acoustic influence; .00005 justified based on initial fits and chi-sq lack-of-fit adjustment
wt_spline =1
wt_spde=1

# wt_ud = .07  #above 1 increases ud influence
# wt_ckmr = 1  #above 1 decreases ckmr influence
# wt_aerial = 67 #above 1 decreases aerial influence
# wt_acc=0.0005  # above 1 increases acoustic influence; .0005 justified bassed on initial fits and chi-sq lack-of-fit adjustment
# wt_spline =1
# wt_spde=1

low = -10
#high = abs(low)/4  #so low values will have SE's that make high counts improbable
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


load("BD_telem_m_y_table.RData")
Points_s_y = matrix(0,15,2)
Points_s_y[,1]=rowSums(Points_m_y[,1:5])
Points_s_y[,2]=rowSums(Points_m_y[,6:11])
Points_s_y[2:15,1]=Points_s_y[2:15,1]+Points_m_y[1:14,12] #adding in december to next year's season 1

Which_tsteps_UD = which(t(Points_s_y)>0)

#nb: 2 seasons/yr hardwired in Data and Params
Data<- list("Acc_k"=as.matrix(Acc_det_sum[,c(2,4,6,8,10)]), 
            "Acc_n"=as.matrix(Acc_det_sum[,c(3,5,7,9,11)]),
            "Area_acc"=Area_acc,"X_acc"=Acc_DM,
            "Acc_tstep"=Acc_det_sum$tstep-1,
            "log_ckmr_N" = log(N_ckmr), "log_ckmr_se"=log_ckmr_se,
            "UD_mean_adj" = UD_mean_adj, "mean_UD" = mean_UD, "which_tsteps_UD" = Which_tsteps_UD-1,"W_st" = W_st,
            "which_Chukchi"=which_Chukchi-1, "which_Bering"=which_Bering-1, "which_Beaufort"=which_Beaufort-1,
            "log_2012_N" = log_2012_N, "log_2013_N"=log_2013_N,"log_2016_N"=log_2016_N,"log_2021_N"=log_2021_N,
            "log_2012_se" = log_2012_SE,"log_2013_se" = log_2013_SE,"log_2016_se" = log_2016_SE,"log_2021_se"=log_2021_SE,
            "N_2012" = exp(log_2012_N), "N_2013"=exp(log_2013_N), "N_2016"=exp(log_2016_N),"N_2021"=exp(log_2021_N),
            "X_ice"=X_ice,"Land_cover"=grid_sf$land_cover, "X_s"=X_s,
            "designMatrixForReport"=.bdiag(designMatrixForReport), "S" = S_combined, "Sdims"=Sdims,
            "Wts" = c(wt_ud, wt_ckmr, wt_aerial, wt_acc, wt_spline,wt_spde),"n_tsteps"=n_season,
            "Duty_lengths"=c(30,60,90,120,180),"est_acc"=1,
            "M0"=spde$param.inla$M0,"M1"=spde$param.inla$M1,"M2"=spde$param.inla$M2,
            "A"=A,"options"=c(0,0),n_eta=ncol(A),n_seasons=2,flag=1
)


Params = list("log_N"=log(500000),"Beta_s" = rep(0,ncol(X_s)),"Beta_ice" = rep(0,ncol(X_ice[,,1])),
              "Beta_acc"=rep(0,ncol(Acc_DM)),
              "log_lambda" = rep(0,length(Sdims)),"log_tau"=0, "log_kappa"=0,
               "Etainput_s"=rep(0,2*ncol(A)),"logit_rho"=0
)

Params$Beta_acc[1] = -9  #needs to start low to prevent numerical issues

Random <- ("Etainput_s")
#Random <- c("Beta_s","log_lambda")
#Random <- NULL
Map <- list(log_lambda=factor(rep(NA,length(Sdims))),
            Etainput_s = factor(rep(NA,2*ncol(A))),
            log_tau = factor(NA),log_kappa = factor(NA)
            #Beta_acc=factor(rep(NA,ncol(Acc_DM)))
            #Beta_s = factor(rep(NA,2)),
            #Beta_ice = factor(rep(NA,dim(X_ice)[2]))
            )
#Map=NULL

Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_betabinom_SPDErr",silent=FALSE)
Obj$fn( Obj$par )
Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=100000, iter.max=100000))         #
Report = Obj$report()
Opt$message
Report0 = Report

Map <- list(log_lambda=factor(rep(NA,length(Sdims)))
            #Beta_acc=factor(rep(NA,ncol(Acc_DM)))
            #Beta_s = factor(rep(NA,2)),
            #Beta_ice = factor(rep(NA,dim(X_ice)[2]))
            #Etainput_s = factor(rep(NA,ncol(A))),
            #log_tau = factor(NA),log_kappa = factor(NA)
)

init_time <- Sys.time()

#Params$log_N=log(Report$N)
#Params$Beta_s = Report$Beta_s
#Params$Beta_ice = Report$Beta_ice  #for some reason, getting newton failure in Obj$fn(Obj$par) when starting at values from fixed effects model
Params$Beta_acc = Report$Beta_acc
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_betabinom_SPDErr",silent=FALSE)
Obj$fn( Obj$par )
Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=100000, iter.max=100000))         #
Report = Obj$report()
Opt$message
Report1 = Report


end_time <- Sys.time()
end_time - init_time

save.image('bearded_SPDErr_results_2season.RData')

Z_cat = c(Report$Z_score_2012,Report$Z_score_2013,Report$Z_score_2016)
prop_gt_90 = sum(abs(Z_cat)>1.64)/length(Z_cat)

# wt_ud = .6  #above 1 increases ud influence
# wt_ckmr = 1  #above 1 decreases ckmr influence
# wt_aerial = 5 #above 1 decreases aerial influence (the value of 5 is justified by examining variance ratios of Var(N) to Var(sum(Z)) )
# wt_acc=0.0005  # obsolete w/ beta-binom
# wt_spline =1
# wt_spde=1
# Data$Wts = c(wt_ud, wt_ckmr, wt_aerial, wt_acc, wt_spline,wt_spde)

Params$log_N = log(500000)
Params$Beta_s = Report$Beta_s
Params$Beta_ice = Report$Beta_ice
Params$Beta_acc = Report$Beta_acc
Params$log_lambda = log(Report$lambda)
Params$log_tau = Report$log_tau
Params$log_kappa = Report$log_kappa
Params$Etainput_s = Report$Etainput_s

# init_time <- Sys.time()
# Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_real_UD_SPDErr",silent=FALSE)
# Obj$fn( Obj$par )
# Obj <- normalize ( Obj , flag ="flag", value = 0)
# Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
# Upper = 50
# Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=100000, iter.max=100000))         #
# Report = Obj$report()
# Opt$message
# Report1 = Report
# end_time <- Sys.time()
# end_time - init_time
# 
# save.image('bearded_SPDErr_results2.RData')

Z_cat = c(Report$Z_score_2012,Report$Z_score_2013,Report$Z_score_2016)
prop_gt_90 = sum(abs(Z_cat)>1.64)/length(Z_cat)
52260/Report$jnll_comp[1]  #supposed UD weight adjustment
# 
n_iter = 4
N_history = rep(0,n_iter+1)
N_history[1]=Report$N

Wts_history = matrix(0,n_iter+1,6)
Wts_history[1,]=Data$Wts
Conv = rep(0,n_iter)
for(i in 1:n_iter){
  Params$log_N = Report$log_N
  Params$Beta_s = Report$Beta_s
  Params$Beta_ice = Report$Beta_ice
  Params$Beta_acc = Report$Beta_acc
  Params$log_lambda = c(0,0,0,0) #not used currently
  Params$log_tau = Report$log_tau
  Params$log_kappa = Report$log_kappa
  Params$Etainput_s = Report$Etainput_s
  Params$logit_rho = log(Report$rho / (1-Report$rho))
  Data$Wts[1]=Data$Wts[1]*52260/Report$jnll_comp[1]
  #Data$Wts[3]=Data$Wts[3]*5792/Report$jnll_comp[3]
  #if(Data$Wts[3]<0.09 | Data$Wts[3]>0.11)Data$Wts[3] = Data$Wts[3]*(0.2+prop_gt_90)/0.3
  #Data$Wts[4]=Data$Wts[4]*905/Report$chisq_acc
  Wts_history[i+1,]=Data$Wts
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="fit_bearded_betabinom_SPDErr",silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=10000, iter.max=10000))         #
  Report = Obj$report()
  N_history[i+1]=Report$N
  Z_cat = c(Report$Z_score_2012,Report$Z_score_2013,Report$Z_score_2016)
  prop_gt_90 = sum(abs(Z_cat)>1.64)/length(Z_cat)
  Conv[i]=1-Opt$convergence
  
}

save.image('temp_2season.RData')

#SD2=sdreport(Obj,bias.correct=FALSE)

#matrix(c(-9,3,3,-1,-1,0),ncol=1)

sum(Report$Z_st[which_Bering,17])  #number in Bering, winter 2012
sum(Report$Z_st[which_Bering,25])  #number in Bering, winter 2016
sum(Report$Z_st[which_Bering,29])  #number in Bering, winter 2018
sum(Report$Z_st[which_Bering,31])  #number in Bering, winter 2019

library(ggplot2)
png('Est_2012.png')
grid_sf$N = Report$Z_st[,17]
#grid_sf$N1[which(grid_sf$N1>1000)]=1000
ggplot(grid_sf)+geom_sf(aes(fill=N),color=NA)+theme(axis.text=element_blank())
dev.off()

png('Est_2018.png')
grid_sf$N = Report$Z_st[,30]
#winter_index = c(1:18)*2-1
#summer_index = winter_index+1
#grid_sf$N = rowMeans(Report$Z_st[,winter_index])
#grid_sf$N = rowMeans(Report$Z_st[,summer_index])
#grid_sf$N = rowMeans(Ice[,summer_index])
#grid_sf$N = Ice[,29]

#Eta_1 = as.matrix(Data$A) %*% Report$Etainput_s[1:208]
#Eta_2 = as.matrix(Data$A) %*% Report$Etainput_s[209:416]
#Eta_2[which(Eta_2<(-11))]=NA
#Eta_1[which(Eta_1<(-14))]=NA
#grid_sf$N=Eta_1
#grid_sf$N=exp(logUD[,18])
#grid_sf$N=grid_sf$Chukchi_smoothed
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

ggplot(grid_sf)+geom_sf(aes(fill=JOBSS_smoothed),color=NA)

#plot acoustic effects
Ice_conc = c(0:100)/100
Season = c("July-Dec","Jan","Feb","Mar","Apr","May","June")
Plot_df = expand.grid(Ice_conc=Ice_conc,Season=Season)
Plot_df$Rate = 0
for(i in 1:nrow(Plot_df)){
  which_season = which(Season==Plot_df$Season[i])
  beta_season = 0
  if(which_season>1)beta_season=Report$Beta_acc[which_season]
  Plot_df$Rate[i]=exp(Report$Beta_acc[1]+beta_season+Report$Beta_acc[8]*Plot_df$Ice[i]+Report$Beta_acc[9]*(Plot_df$Ice[i])^2)
}
Plot_df$Probability = 1-exp(-60*Plot_df$Rate*100)
ggplot(Plot_df)+geom_line(aes(x=Ice_conc,y=Probability,color=Season,group=Season),size=1.1)+
  xlab('Ice concentration (proportion)')+theme(axis.text=element_text(size=10),legend.text=element_text(size=10))

png('Call_rate.png',width=6,height=6,units="in",res=600)
ggplot(Plot_df)+geom_line(aes(x=Ice_conc,y=Probability,color=Season,group=Season),size=1.1)+
  xlab('Ice concentration (proportion)')+theme(axis.text=element_text(size=10),legend.text=element_text(size=10))
dev.off()


#plot depth effect
par(mfrow=c(1,1))
Depth = c(-80:0)/10
plot(Depth*100,Report$Beta_s[1]*Depth+Report$Beta_s[2]*Depth^2,main="Depth",xlim=c(-800,0),ylim=c(-10,0))

#plot ice effects
par(mfrow=c(1,2))
plot(Ice_conc,Report$Beta_ice[1]*Ice_conc+Report$Beta_ice[2]*Ice_conc^2,main='Winter')
plot(Ice_conc,Report$Beta_ice[3]*Ice_conc+Report$Beta_ice[4]*Ice_conc^2,main='Summer')

#how about a straight polynomial regression of Z_st on ice in winter?
Which_winter = c(1:18)*2-1
Z_winter = Report$Z_st[,Which_winter]
Z_winter = as.vector(Z_winter)
Ice_winter = Data$X_ice[,2,Which_winter]
Ice_winter = as.vector(Ice_winter)
Ice_winter2 = Ice_winter^2
crap = lm(log(Z_winter)~Ice_winter+Ice_winter2)
Coefs = coef(crap)
plot(Ice_conc,exp(Coefs[1]+Coefs[2]*Ice_conc+Coefs[3]*Ice_conc^2))

save.image("bearded_SPDE_results_Mar2023.RData")
