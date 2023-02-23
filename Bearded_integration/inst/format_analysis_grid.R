# Plot grid for bearded seal integration, make some data
# summary maps

###  Attach environmental covariates to "big" (US + Russia) CHESS grid

install_pkg <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Install libraries ----------------------------------------------
install_pkg("RPostgreSQL")
install_pkg("sf")
install_pkg("devtools")
install_pkg("dplyr")
install_pkg("sp")
install_pkg("rgeos")
install_pkg("automap")
install_pkg("rpostgis")
install_pkg("raster")
install_pkg("fasterize")
install_pkg("rgdal")

#function to drop the geometry of an st object (easier to work w/ the data.frame)
st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}

# Run code -------------------------------------------------------
# Extract data from DB ------------------------------------------------------------------
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))
grid_sf <- sf::st_read(con, 
                             query = "SELECT * FROM base.geo_analysis_grid_no_polarbear", 
                             geometry_column = "geom")

grid_sf <- grid_sf[-which(grid_sf$cell<5000),]#take out sea of Okhotsk


#calculate proportion area covered by land and remove cells with >99% of area on land
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK = st_transform(AK,st_crs(grid_sf))
Russia = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/russia_dcw.shp')
Russia = st_transform(Russia,st_crs(grid_sf))
Land=st_union(st_union(AK),st_union(Russia))
n_cells=nrow(grid_sf)
Land.area=rep(0,n_cells)
I.intersect = st_intersects(grid_sf,Land)
for(icell in 1:n_cells){
  if(length(I.intersect[[icell]])>0)Land.area[icell]=st_area(st_intersection(grid_sf[icell,],Land))
}
grid_sf$land_cover=Land.area/628381060
grid_sf=grid_sf[-which(grid_sf$land_cover>0.99),]
n_cells=nrow(grid_sf)

#attach depth
depth <- sf::st_read(con,
                     query = "SELECT *
                         FROM environ.geo_etopo1_sdr_bins", 
                     geometry_column = "geom")
depth = st_transform(depth,st_crs(grid_sf))
depth = depth[-which(depth$gridcode>0),]  #take out land
#take out features that don't intersect w/ chess grid
chess_union = st_union(grid_sf)
I.grid = st_intersects(depth,chess_union,sparse=FALSE)
depth=depth[I.grid>0,]
depth$mean = apply(cbind(depth$gridcode,depth$meandepth),1,'mean')  #using midpoint between min and max depth for each polygon
grid_sf$depth = NA
for(i in 1:nrow(grid_sf)){
  Pts = st_sample(grid_sf[i,],400,type='regular')  #lay down 400 regular points in each cell
  N.cells = colSums(st_intersects(Pts,depth,sparse=FALSE))
  grid_sf$depth[i]=sum(N.cells %*% depth$mean)/sum(N.cells)
} 
grid_sf$depth[which(is.na(grid_sf$depth))]=0  #in case some cells 400 points were all on land
#note depth truncated at 600m (deeper gets assigned 600m)

#distance to land
library(ptolemy) #from josh london's github
bbox = st_union(grid_sf)
Land_sf = extract_gshhg(
  data=bbox,
  resolution = "i",
  epsg = NULL,
  buffer = 5000,
  simplify = FALSE,
  warn = FALSE
)
grid_centroids = st_centroid(grid_sf)
Distances=st_distance(grid_centroids,Land_sf,byid=TRUE)
grid_sf$dist_land=apply(Distances,1,'min')/1000 #change to km



save.image("format_analysis_grid1.RData")  #the depth calcs take a bit of time....might just want to load workspace if needed
#load("format_analysis_grid1.RData")   #the above takes awhile; might want to load from here.


#now, attach abundance estimates!!!
load('c:/users/paul.conn/git/CHESS/Ests_CHESS_mean.RData')
load('c:/users/paul.conn/git/BOSSst/Ests_EBS_2012_mean.RData')
load('c:/users/paul.conn/git/BOSSrussia/Ests_WBS_2012_mean.RData')
load('c:/users/paul.conn/git/BOSSst/Ests_EBS_2013_mean.RData')
load('c:/users/paul.conn/git/BOSSrussia/Ests_WBS_2013_mean.RData')

#these estimates cover different spatial areas and different survey
#periods.  For simplicity, we'll use mean values (taken over day of survey)

#Chukchi
Est_CHESS$meta$grid$Bd = Est_CHESS$bd_mean
Est_CHESS$meta$grid$Bd_SE = sqrt(Est_CHESS$var_infl*diag(Est_CHESS$bd_vc))
# looks like cell #s are inconsistent among the CHESS grid and the bigger grid; resort to sf 
Within = st_within(st_centroid(Est_CHESS$meta$grid),grid_sf,sparse=FALSE)
Mapping_chukchi = rep(0,nrow(Est_CHESS$meta$grid))  #this object will be needed for integrated modeling too
for(icell in 1:nrow(Est_CHESS$meta$grid))Mapping_chukchi[icell]=which(Within[icell,]==1)
grid_sf$bd_Ch = grid_sf$bd_Ch_SE = NA
grid_sf$bd_Ch[Mapping_chukchi]=Est_CHESS$meta$grid$Bd
grid_sf$bd_Ch_SE[Mapping_chukchi] = Est_CHESS$meta$grid$Bd_SE 

#western Bering
Est_wBS_2012$meta$grid = st_as_sf(Est_wBS_2012$meta$grid)  
Est_wBS_2012$meta$grid$Bd = Est_wBS_2012$bd_mean
Est_wBS_2012$meta$grid$Bd_SE = sqrt(Est_wBS_2012$var_infl*(diag(Est_wBS_2012$bd_vc)+0.0001)) # a couple entries are slightly negative
Within = st_within(st_centroid(Est_wBS_2012$meta$grid),grid_sf,sparse=FALSE)
Mapping_wBS = rep(0,nrow(Est_wBS_2012$meta$grid))  #this object will be needed for integrated modeling too
for(icell in 1:nrow(Est_wBS_2012$meta$grid))Mapping_wBS[icell]=which(Within[icell,]==1)
grid_sf$bd_wBS_2012 = grid_sf$bd_wBS_SE_2012 = NA
grid_sf$bd_wBS_2012[Mapping_wBS]=Est_wBS_2012$meta$grid$Bd
grid_sf$bd_wBS_SE_2012[Mapping_wBS] = Est_wBS_2012$meta$grid$Bd_SE 

Est_wBS_2013$meta$grid = st_as_sf(Est_wBS_2013$meta$grid)  
Est_wBS_2013$meta$grid$Bd = Est_wBS_2013$bd_mean
Est_wBS_2013$meta$grid$Bd_SE = sqrt(Est_wBS_2013$var_infl*(diag(Est_wBS_2013$bd_vc)+0.0001)) # a couple entries are slightly negative
grid_sf$bd_wBS_2013 = grid_sf$bd_wBS_SE_2013 = NA
grid_sf$bd_wBS_2013[Mapping_wBS]=Est_wBS_2013$meta$grid$Bd
grid_sf$bd_wBS_SE_2013[Mapping_wBS] = Est_wBS_2013$meta$grid$Bd_SE 

#eastern Bering
Est_eBS_2012$meta$grid$Bd = Est_eBS_2012$bd_mean
Est_eBS_2012$meta$grid$Bd_SE = sqrt(Est_eBS_2012$var_infl*(diag(Est_eBS_2012$bd_vc)+0.0001)) # a couple entries are slightly negative
Within = st_within(st_centroid(Est_eBS_2012$meta$grid),grid_sf,sparse=FALSE)
Mapping_eBS = rep(0,nrow(Est_eBS_2012$meta$grid))  #this object will be needed for integrated modeling too
for(icell in 1:nrow(Est_eBS_2012$meta$grid))Mapping_eBS[icell]=which(Within[icell,]==1)
grid_sf$bd_eBS_2012 = grid_sf$bd_eBS_SE_2012 = NA
grid_sf$bd_eBS_2012[Mapping_eBS]=Est_eBS_2012$meta$grid$Bd
grid_sf$bd_eBS_SE_2012[Mapping_eBS] = Est_eBS_2012$meta$grid$Bd_SE 

Est_eBS_2013$meta$grid$Bd = Est_eBS_2013$bd_mean
Est_eBS_2013$meta$grid$Bd_SE = sqrt(Est_eBS_2013$var_infl*(diag(Est_eBS_2013$bd_vc)+0.0001)) # a couple entries are slightly negative
grid_sf$bd_eBS_2013 = grid_sf$bd_eBS_SE_2013 = NA
grid_sf$bd_eBS_2013[Mapping_eBS]=Est_eBS_2013$meta$grid$Bd
grid_sf$bd_eBS_SE_2013[Mapping_eBS] = Est_eBS_2013$meta$grid$Bd_SE 


#plot BOSS
library(ggplot2)
library(viridis)
library(cowplot)
grid_sf$Bering_2012 = grid_sf$bd_eBS_2012
grid_sf$Bering_SE_2012 = grid_sf$bd_eBS_SE_2012
grid_sf$Bering_2012[which(! is.na(grid_sf$bd_wBS_2012))]=grid_sf$bd_wBS_2012[which(! is.na(grid_sf$bd_wBS_2012))]
grid_sf$Bering_SE_2012[which(! is.na(grid_sf$bd_wBS_2012))]=grid_sf$bd_wBS_SE_2012[which(! is.na(grid_sf$bd_wBS_2012))]
plot_BOSS_2012 = ggplot(grid_sf)+geom_sf(aes(fill=Bering_2012,colour=Bering_2012))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) + 
  ggtitle("A. Bering 2012")

grid_sf$Bering_2013 = grid_sf$bd_eBS_2013
grid_sf$Bering_SE_2013 = grid_sf$bd_eBS_SE_2013
grid_sf$Bering_2013[which(! is.na(grid_sf$bd_wBS_2013))]=grid_sf$bd_wBS_2013[which(! is.na(grid_sf$bd_wBS_2013))]
grid_sf$Bering_SE_2013[which(! is.na(grid_sf$bd_wBS_2013))]=grid_sf$bd_wBS_SE_2013[which(! is.na(grid_sf$bd_wBS_2013))]

plot_BOSS_2013 = ggplot(grid_sf)+geom_sf(aes(fill=Bering_2013,colour=Bering_2013))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) +
  ggtitle("B. Bering 2013")

plot_Chukchi = ggplot(grid_sf)+geom_sf(aes(fill=bd_Ch,colour=bd_Ch))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) +
  ggtitle("C. Chukchi 2016")

png("Abundance_maps.png")
  plot_grid(plot_BOSS_2012,plot_BOSS_2013,plot_Chukchi,ncol=1)
dev.off()


### kernel code 
Centroids = st_centroid(grid_sf)
Dists = units::drop_units(st_distance(Centroids,Centroids))
Probs = matrix(0,n_cells,n_cells)
kern.sd = 20000 #in m; based in part on kernel home range analysis in eda_telemetry.R
cutoff = qnorm(0.9999,0,kern.sd)  #cutoff to limit # of normal calcs
for(icell in 1:n_cells){
  Which_calc = which(Dists[icell,]<cutoff)
  Probs[icell,Which_calc]=dnorm(Dists[icell,Which_calc],0,kern.sd)
  Probs[icell,Which_calc]=Probs[icell,Which_calc]/sum(Probs[icell,Which_calc])
}
#Probs holds transition matrix; now need to apply this to abundance values
#to smooth them into cells with no ice
grid_sf$Bering_2012_smoothed = grid_sf$Bering_2012
grid_sf$Bering_2012_smoothed[which(is.na(grid_sf$Bering_2012_smoothed))]=0
grid_sf$Bering_2012_smoothed = Probs %*% grid_sf$Bering_2012_smoothed
grid_sf$Bering_2012_smoothed[which(is.na(grid_sf$Bering_2012))]=NA #reset cells in different sea to NA
Mapping_Bering = which(! is.na(grid_sf$Bering_2012))
library(Matrix)
SE_B2012 = diag(Probs[Mapping_Bering,Mapping_Bering] %*% 
                  Diagonal(x=grid_sf$Bering_SE_2012[Mapping_Bering]) %*% 
                  t(Probs[Mapping_Bering,Mapping_Bering]))
grid_sf$Bering_2012_smoothed_SE = rep(NA,n_cells) 
grid_sf$Bering_2012_smoothed_SE[Mapping_Bering]=SE_B2012

grid_sf$Bering_2013_smoothed = grid_sf$Bering_2013
grid_sf$Bering_2013_smoothed[which(is.na(grid_sf$Bering_2013_smoothed))]=0
grid_sf$Bering_2013_smoothed = Probs %*% grid_sf$Bering_2013_smoothed
grid_sf$Bering_2013_smoothed[which(is.na(grid_sf$Bering_2013))]=NA #reset cells in different sea to NA
SE_B2013 = diag(Probs[Mapping_Bering,Mapping_Bering] %*% 
                  Diagonal(x=grid_sf$Bering_SE_2013[Mapping_Bering]) %*% 
                  t(Probs[Mapping_Bering,Mapping_Bering]))
grid_sf$Bering_2013_smoothed_SE = rep(NA,n_cells) 
grid_sf$Bering_2013_smoothed_SE[Mapping_Bering]=SE_B2013

Mapping_Chukchi = which(! is.na(grid_sf$bd_Ch))
grid_sf$Chukchi_smoothed = grid_sf$bd_Ch
grid_sf$Chukchi_smoothed[which(is.na(grid_sf$Chukchi_smoothed))]=0
grid_sf$Chukchi_smoothed = Probs %*% grid_sf$Chukchi_smoothed
grid_sf$Chukchi_smoothed[which(is.na(grid_sf$bd_Ch))]=NA #reset cells in different sea to NA
SE_Chukchi = diag(Probs[Mapping_Chukchi,Mapping_Chukchi] %*% 
                  Diagonal(x=grid_sf$bd_Ch_SE[Mapping_Chukchi]) %*% 
                  t(Probs[Mapping_Chukchi,Mapping_Chukchi]))
grid_sf$Chukchi_smoothed_SE = rep(NA,n_cells) 
grid_sf$Chukchi_smoothed_SE[Mapping_Chukchi]=SE_Chukchi

plot_BOSS_2012_smoothed = ggplot(grid_sf)+geom_sf(aes(fill=Bering_2012_smoothed,color=Bering_2012_smoothed))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) + 
  ggtitle("D. Bering 2012 - smoothed")

plot_BOSS_2013_smoothed =ggplot(grid_sf)+geom_sf(aes(fill=Bering_2013_smoothed,color=Bering_2013_smoothed))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) + 
  ggtitle("E. Bering 2013 - smoothed")

plot_Chukchi_smoothed =ggplot(grid_sf)+geom_sf(aes(fill=Chukchi_smoothed,color=Chukchi_smoothed))+
  scale_fill_viridis(name="Abundance",na.value="transparent",limits=c(0,3000)) + 
  scale_color_viridis(name="Abundance",na.value="lightgray",limits=c(0,3000)) + 
  ggtitle("F. Chukchi 2016 - smoothed")


png("Abundance_maps2.png")
plot_grid(plot_BOSS_2012,plot_BOSS_2012_smoothed,plot_BOSS_2013,
          plot_BOSS_2013_smoothed,plot_Chukchi,plot_Chukchi_smoothed,ncol=2)
dev.off()

save.image("format_analysis_grid2.RData")  #the depth calcs take a bit of time....might just want to load workspace if needed

load("format_analysis_grid2.RData")
library(ggplot2)
library(viridis)
library(cowplot)




#telemetry data currently 2004 - 2018; attach seasonal sea ice concentration values
start_yr = 2004
end_yr = 2018
n_yrs = end_yr-start_yr+1
season_start = c("01-01","04-01","06-01","08-01","10-01")
season_end = c("03-31","05-31","07-31","09-30","12-31")
n_seasons = length(season_end)
t_steps = n_seasons * n_yrs

qry <- "SELECT * FROM base.tbl_analysis_grid_cov_seaice WHERE cell>10000"
sea_ice = sf::st_read(con,query=qry)

Cell_IDs = grid_sf$cell

### average sea ice concentration values by cell, year, season 
start_date = end_date = rep('NA',t_steps)
Ice = matrix(NA,n_cells,t_steps)
for(iyr in 1:n_yrs){
  for(iseason in 1:n_seasons){
    date1 = paste0(iyr+start_yr,'-',season_start[iseason])
    date2 = paste0(iyr+start_yr,'-',season_end[iseason])
    Temp_ice = sea_ice[which(sea_ice$fdate>=date1 & sea_ice$fdate<=date2),]
    for(icell in 1:n_cells){
      Ice[icell,n_seasons*(iyr-1)+iseason] <- mean(Temp_ice[which(Temp_ice$cell==Cell_IDs[icell]),"rast_seaice"],na.rm=TRUE)
    }
  }
}

save.image(file='analysis_grid_ice_attached.RData')

#krige ice data
library(automap)
grid_sp <- as(st_centroid(grid_sf),'Spatial')
for(iyr in 1:n_yrs){
  for(iseason in 1:n_seasons){
    grid_sp$ice = Ice[,n_seasons*(iyr-1)+iseason]
    Which_problem = which(Ice[,n_seasons*(iyr-1)+iseason]>1.0 | grid_sf$land_cover>0.01)
    fit_points = grid_sp[-Which_problem,]
    pred_points = grid_sp[Which_problem,]
    krige_out=autoKrige(ice~1,input_data=fit_points,new_data=pred_points)$krige_output[["var1.pred"]] 
    if(sum(krige_out<0)>0){
      krige_out[which(krige_out<0)]=0
    }
    Ice[Which_problem,n_seasons*(iyr-1)+iseason]=krige_out
  }
}

save.image(file='prediction_grid_ice_kriged.RData')


#plot showing seasonal variability between years
Plots = vector("list",5)  #season on first level
Years = c(1,6,11,15)
for(is in 1:5){
  Plots[[is]]=vector("list",4) 
  for(iyr in 1:4){
    grid_sf$ice = Ice[,n_seasons*(Years[iyr]-1)+is]
    Plots[[is]][[iyr]]=ggplot(grid_sf)+geom_sf(aes(fill=ice,color=ice))+
      scale_fill_viridis(name="Ice",na.value="transparent") + 
      scale_color_viridis(name="Ice",na.value="lightgray") 
  }
}
#now, cowplot

#plot with 
png("Ice_maps.png")
plot_grid(Plots[[1]][[1]],Plots[[1]][[2]],Plots[[1]][[3]],Plots[[1]][[4]],
          Plots[[2]][[1]],Plots[[2]][[2]],Plots[[2]][[3]],Plots[[2]][[4]],
          ncol=4)
dev.off()


grid_sf$ice5 = Ice[,5]
ggplot(grid_sf)+geom_sf(aes(fill=ice5,color=ice5))+
  scale_fill_viridis(name="Ice",na.value="transparent",begin=0,end=1) + 
  scale_color_viridis(name="Ice",na.value="lightgray",begin=0,end=1) 


