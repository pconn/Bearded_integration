library(sf)
library("RPostgreSQL")

#from J. London

#remotes::install_github("rstudio/pins")
#library(pins)
#board <- board_rsconnect(server = "connect.fisheries.noaa.gov")
#pin_read(board, "josh.london/bearded_seal_movement")

library(pins)
board <- board_url(c(
  bs_move_tble = 'https://connect.fisheries.noaa.gov/content/05b5b060-b0e1-451d-a0b9-077fc3576c06/',
  bs_pred_sf = 'https://connect.fisheries.noaa.gov/content/ed2f5b6d-4075-4285-96fe-2adbb4c89dc7/'))
#dat_tbl <- pin_read(board, "bs_move_tble")


BD_move <- readRDS("bearded_seal_movement.rds")

#access invidividual seal sfc objects using eg. BD_move[[4]][[1]]
library(tidyverse)
n_seal = length(BD_move[[1]])
BD_unlist = st_sf(BD_move[[4]][[1]])
for(iseal in 2:n_seal){
  BD_unlist <- rbind(BD_unlist,st_sf(BD_move[[4]][[iseal]]))
}
BD_unlist <- BD_unlist %>% mutate(year = lubridate::year(datetime)) %>% 
               mutate(month = lubridate::month(datetime))

BD_telem <- BD_unlist
save(BD_telem,file="BD_telem_NOAA.RData")

load("BD_telem_NOAA.RData")
#read in analysis grid, Ice data 
load('prediction_grid_ice_kriged.RData')

thin = 10  #use 1 out of every 10 hourly locations to start with...we can push this up later
n_rec = nrow(BD_telem)
BD_telem = BD_telem[which((c(1:n_rec)%%thin)==1),]

EN = st_coordinates(st_centroid(grid_sf))
grid_sf$easting=EN[,1]/1000000
grid_sf$northing=EN[,2]/1000000

Telem = vector("list",n_seasons)

Month = lubridate::month(BD_telem$datetime)
Year = lubridate::year(BD_telem$datetime)
BD_telem$ID = paste(BD_telem$deployid,Year)
get_which<-function(Row){
  which(Row==1)
}

N_ind = rep(0,n_seasons)
for(is in 1:n_seasons){
  if(is==1)Which = which(Month %in% c(12,1:3))
  if(is==2)Which = which(Month %in% c(4:5))
  if(is==3)Which = which(Month %in% c(6:7))
  if(is==4)Which = which(Month %in% c(8:9))
  if(is==5)Which = which(Month %in% c(10:11))
  Cur_telem = BD_telem[Which,]
  n_ind = length(unique(Cur_telem$ID))
  N_ind[is]=n_ind
  Telem[[is]]=vector("list",n_ind)
  IDs = unique(Cur_telem$ID)
  for(iind in 1:n_ind){
    Which_obs = which(Cur_telem$ID == IDs[iind])
    Telem_ind = Cur_telem[Which_obs,]
    #intersect telemetry locations with grid
    Inter = st_intersects(Telem_ind,grid_sf,sparse=FALSE)
    Which = unlist(apply(Inter,1,"get_which"))
    Telem[[is]][[iind]]=grid_sf[,c("easting","northing","depth","dist_land","land_cover")]
    Telem[[is]][[iind]]$Count = tabulate(Which,nbins=nrow(grid_sf))
    yr = lubridate::year(Telem_ind[1,]$datetime)-2003 #start in 2004
    Telem[[is]][[iind]]$sea_ice = Ice[,(yr-1)*n_seasons+is]
  }
}

#now group data by individual, including season as a categorical covariate
Cur_data = Telem[[1]][[1]]
n_ind = length(Telem[[1]])
Data_list = vector("list",n_ind)
Data_list[[1]]=Cur_data
Data_list[[1]]$season=1
Data_list[[1]]$season = factor(Data_list[[1]]$season, levels=c('1','2','3','4','5'))
Cur_data$ind=1
for(iind in 2:n_ind){
  Tmp = Telem[[1]][[iind]]
  Tmp$ind = iind
  Cur_data = rbind(Cur_data,Tmp)
  Data_list[[iind]]=Tmp
  Data_list[[iind]]$season=1
  Data_list[[iind]]$season = factor(Data_list[[iind]]$season, levels=c('1','2','3','4','5'))
}
Cur_data$season=1
Data = Cur_data
prev_ind = iind

for(is in 2:n_seasons){
  Cur_data = Telem[[is]][[1]]
  prev_ind=prev_ind+1
  Data_list[[prev_ind]]=Cur_data
  Data_list[[prev_ind]]$season=is
  Data_list[[prev_ind]]$season = factor(Data_list[[prev_ind]]$season, levels=c('1','2','3','4','5'))
  Cur_data$ind=prev_ind
  n_ind = length(Telem[[is]])
  for(iind in 2:n_ind){
    prev_ind=prev_ind+1
    Tmp = Telem[[is]][[iind]]
    Tmp$ind = prev_ind
    Cur_data = rbind(Cur_data,Tmp)
    Data_list[[prev_ind]]=Tmp
    Data_list[[prev_ind]]$season=is
    Data_list[[prev_ind]]$season = factor(Data_list[[prev_ind]]$season, levels=c('1','2','3','4','5'))
  }
  Cur_data$season=is
  Data = rbind(Data,Cur_data)
}
Data$season = as.factor(as.character(Data$season))
Data$ind = as.factor(Data$ind)

save(Data,file="telem_counts.RData")
save(Data_list,file="telem_counts_list.RData")
save.image("telem_analysis_workspace.RData")



load("telem_analysis_workspace.RData")

#analysis
library(TMB)
TmbFile = "C:/Users/paul.conn/git/Bearded_integration/Bearded_integration/src/fit_telem_counts_depthsm"
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 
dyn.load( dynlib(TmbFile) )


#set up SPDE basis using fmesher (these funcs were originally in INLA)
# Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
library(fmesher)
library(INLA)
library(sf)

set.seed(2020)

Grid_locs = st_coordinates(st_centroid(grid_sf))/1000
lat = st_coordinates(st_centroid(grid_sf))/1000
#mesh_dense = fm_mesh_2d_inla(Grid_locs)
#n_knots_dense = mesh_dense$n
#Mesh_index_dense = mesh_dense$idx$loc-1  #which spde REs to apply as random effects for each cell centroid
#spde_dense <- (inla.spde2.matern(mesh_dense, alpha=2)$param.inla)[c("M0","M1","M2")]


n_ind = length(Data_list)
Y_i = matrix(0,n_cells,n_ind)
for(i in 1:n_ind)Y_i[,i]=Data_list[[i]]$Count

Data$depth[which(Data$depth>0)]=0
Data$depth = Data$depth/100
Data$dist_land = Data$dist_land/100
Data$ice2 = Data$sea_ice^2
Data$depth2 = Data$depth^2
Data$dist_land2 = Data$dist_land^2
I_Bering = (grid_sf$cell<20000)
Data$I_Bering = rep(I_Bering,n_ind)
#X = model.matrix(~ sea_ice + ice2 + depth + depth2 + season:sea_ice + season:ice2 ,data=st_drop_geometry(Data))
X = model.matrix(~ sea_ice + ice2 + I_Bering + season + season:I_Bering + season:sea_ice + season:ice2 ,data=st_drop_geometry(Data))

Data$Dummy=1
gam_setup = mgcv::gam(Dummy ~ s(depth, bs = "cs",k=5),data = Data,fit=FALSE)
S_depth = gam_setup$smooth[[1]]$S[[1]]
depth = seq(min(Data$depth),max(Data$depth),by = 0.1)
depthReport = mgcv::PredictMat(gam_setup$smooth[[1]],data = data.frame(depth))
X_sm = mgcv::PredictMat(gam_setup$smooth[[1]],data = Data)


#format prediction data - looping over years, seasons
Data_pred = data.frame(matrix(0,n_yrs*n_seasons*n_cells,ncol(Data)))
colnames(Data_pred)=colnames(Data)
counter=1
for(iyr in 1:18){
  for(iseas in 1:5){
    Data_pred[counter:(counter+n_cells-1),]=Data[1:n_cells,]
    Data_pred[counter:(counter+n_cells-1),"sea_ice"]=Ice[,(iyr-1)*n_seasons+iseas]
    Data_pred[counter:(counter+n_cells-1),"ice2"]=Ice[,(iyr-1)*n_seasons+iseas]^2
    Data_pred[counter:(counter+n_cells-1),"season"]=iseas
    Data_pred[counter:(counter+n_cells-1),"I_Bering"]=I_Bering
    counter = counter + n_cells
  }
}
Data_pred$season = as.factor(Data_pred$season)
#X_pred = model.matrix(~ sea_ice + ice2 + depth + depth2 + season:sea_ice + season:ice2 ,data=Data_pred)
X_pred = model.matrix(~ sea_ice + ice2 + I_Bering + season + season:I_Bering + season:sea_ice + season:ice2 ,data=Data_pred)
X_sm_pred = mgcv::PredictMat(gam_setup$smooth[[1]],data = Data_pred)

 

#this block lifted from sdmTMB package
set.seed(12345)
knots <- stats::kmeans(x = Grid_locs, centers = 200)  #200 "center" locations (# knots will be more after triangulation)
loc_centers <- knots$centers
mesh <- INLA::inla.mesh.create(loc_centers, refine = TRUE)
spde <- INLA::inla.spde2.matern(mesh)  

A <- INLA::inla.spde.make.A(mesh, loc = as.matrix(Grid_locs)) #prediction matrix at data locs

n_knots = ncol(A)

Data_tmb <- list("Y_i"=Y_i, "X"=X,"X_pred"=X_pred, "N_samp"=colSums(Y_i),
                 "flag"=1, "n_ind"=n_ind, "n_s"=n_cells, "n_mesh"=n_knots, "n_yrs"=18,
                 "options"=c(0,0),"matern_pri"=c(0,0,0,0),"Season"=as.numeric(Data$season)-1,
                 "X_sm"=X_sm,"X_sm_pred"=X_sm_pred,"S"=as(S_depth,"dgTMatrix"),"depthReport"=depthReport)

Data_tmb$A = as(A,"dgTMatrix")
Data_tmb$M0 = spde$param.inla$M0
Data_tmb$M1 = spde$param.inla$M1
Data_tmb$M2 = spde$param.inla$M2
Data_tmb$Offset = log(colSums(Data_tmb$Y_i)/mean(colSums(Data_tmb$Y_i)))


Params <- list("Beta"=rep(0,ncol(X)), "log_tau_eta"=0, "log_kappa_eta"=0,
              "log_tau_xi" = 0, "log_kappa_xi" = 0, "Eta_s"=matrix(0,n_knots,n_seasons),
              "Xi_s"=matrix(0,n_knots,n_ind),"Beta_sm"=rep(0,nrow(S_depth)),"log_lambda"=0)      
Random <- c("Eta_s","Xi_s","Beta_sm","log_lambda")
Map <- list("Eta_s"=factor(matrix(NA,n_knots,n_seasons)),
            "Xi_s"=factor(matrix(NA,n_knots,n_ind)),
            "log_tau_eta"=factor(NA),"log_kappa_eta"=factor(NA),
            "log_tau_xi"=factor(NA),"log_kappa_xi"=factor(NA))

Data_tmb$options[0]=0  #use normalization trick if options[0]=1
Obj = MakeADFun( data=Data_tmb, parameters=Params, random=Random, map=Map, DLL="fit_telem_counts_depthsm",silent=FALSE)
Obj$fn( Obj$par )
#Obj <- normalize ( Obj , flag ="flag", value = 0)
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=4000, iter.max=4000))         #
Report = Obj$report()

#cat(paste0("time = ",Sys.time()-st_time))

SD=sdreport(Obj,bias.correct=FALSE,getReportCovariance=FALSE)


#output surfaces

logUD_output <- list(logUD = Report$log_UD_s, SE = matrix(SD$sd,n_cells,t_steps))
save(logUD_output,file="logUD_output.RData")

#make a plot with all 5 season and two very different years, 2012 and 2019

year = 2012-2004
library(ggplot2)
grid_sf_big = grid_sf
grid_sf_big$Season = "Winter"
grid_sf_big$Year = "2012"
grid_sf_big$RelAbund=exp(Report$log_UD_s[,year*5+1])/sum(exp(Report$log_UD_s[,year*5+1]))


Seasons = c("Winter","Spring","Early open","Open","Fall")
for(is in 2:5){
  grid_tmp = grid_sf
  grid_tmp$Season = Seasons[is]
  grid_tmp$Year = "2012"
  grid_tmp$RelAbund=exp(Report$log_UD_s[,year*5+is])/sum(exp(Report$log_UD_s[,year*5+is]))
  grid_sf_big = rbind(grid_sf_big,grid_tmp)
}
year = 2019-2004
for(is in 1:5){
  grid_tmp = grid_sf
  grid_tmp$Season = Seasons[is]
  grid_tmp$Year = "2019"
  grid_tmp$RelAbund=exp(Report$log_UD_s[,year*5+is])/sum(exp(Report$log_UD_s[,year*5+is]))
  grid_sf_big = rbind(grid_sf_big,grid_tmp)
}
grid_sf_big$Season = factor(grid_sf_big$Season, levels = Seasons)

library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))

Plot = ggplot(grid_sf_big)+geom_sf(aes(color=RelAbund,fill=RelAbund))+
  facet_grid(Season~Year)+
  theme(legend.position="none",axis.text.x = element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank())+
  scale_fill_gradientn(colours=myPalette(100))+
  scale_color_gradientn(colours=myPalette(100))

png("UD_plots.png",width=8,height=8,units="in",res=1200)
  Plot
dev.off()


