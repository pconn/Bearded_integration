### EDA bearded telemetry data
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


#read in analysis grid
library(raster)
load('analysis_grid_ice_kriged.RData')


Land_cover = grid_sf[,"land_cover"]
Land_cover_sp = as(Land_cover, "Spatial")
Land_raster <- raster(Land_cover_sp, res = 25067.53) #create raster with same grid cell res
Land_raster <- rasterize(Land_cover_sp,Land_raster,field="land_cover")
Land_raster[is.na(Land_raster[])]=1 #treat missing cells as if on land

# the bedrock layer is 30 arcseconds
#NOAA National Centers for Environmental Information. 2022: ETOPO 2022 15 Arc-Second
#Global Relief Model. NOAA National Centers for Environmental Information.
#https://doi.org/10.25921/fd45-gt74 . Accessed [12-16-2022]
r1<-raster("./data/depth_etopo2022_bedrock_1.tiff") #W of dateline
r2<-raster("./data/depth_etopo2022_bedrock_2.tiff") #E of dateline
r <- merge(r1,r2)  #note this is really big because we're crossing the date line

depth_raster <- projectRaster(r, Land_raster) 

#distance to land
library(ptolemy) #from josh london's github
raster_sp = rasterToPolygons(Land_raster)
raster_sf = st_as_sf(raster_sp)
bbox = st_union(raster_sf)
Land_sf = extract_gshhg(
  data=bbox,
  resolution = "i",
  epsg = NULL,
  buffer = 5000,
  simplify = FALSE,
  warn = FALSE
)
raster_centroids = st_centroid(raster_sf)
Distances=st_distance(raster_centroids,Land_sf,byid=TRUE)
Dist_land=apply(Distances,1,'min')
dist_land_raster = depth_raster
dist_land_raster$layer=Dist_land

EN = st_coordinates(st_centroid(grid_sf))
grid_sf$easting=EN[,1] #/max(EN[,1])
grid_sf$northing=EN[,2] #/max(EN[,2])
Easting=grid_sf[,"easting"]
Northing=grid_sf[,"northing"]


Land_cover = grid_sf[,"land_cover"]
Land_cover_sp = as(Land_cover, "Spatial")
Land_raster <- raster(Land_cover_sp, res = 25067.53) #create raster with same grid cell res
Land_raster <- rasterize(Land_cover_sp,Land_raster,field="land_cover")
Land_raster[is.na(Land_raster[])]=1 #treat missing cells as if on land

Easting_sp = as(Easting, "Spatial")
Easting_raster <- raster(Easting_sp, res = 25067.53) #create raster with same grid cell res
Easting_raster <- rasterize(Easting_sp,Easting_raster,field="easting")
Easting_raster[is.na(Easting_raster[])]=0 
Northing_sp = as(Northing, "Spatial")
Northing_raster <- raster(Northing_sp, res = 25067.53) #create raster with same grid cell res
Northing_raster <- rasterize(Northing_sp,Northing_raster,field="northing")
Northing_raster[is.na(Northing_raster[])]=-2500000 


con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))


# now the hard part: we're going to need to come up with a rasterBrick for all 10 day 
# periods in the study.  Each 10 day period is going to need to be kriged first though!
qry <- "SELECT * FROM base.tbl_analysis_grid_cov_seaice WHERE cell>10000"
sea_ice = sf::st_read(con,query=qry)


start_yr = 2004
end_yr = 2018
n_yrs = end_yr-start_yr+1
season_start = c("01-01","01-11","01-21",
                 "02-01","02-11","02-21",
                 "03-01","03-11","03-21",
                 "04-01","04-11","04-21",
                 "05-01","05-11","05-21",
                 "06-01","06-11","06-21",
                 "07-01","07-11","07-21",
                 "08-01","08-11","08-21",
                 "09-01","09-11","09-21",
                 "10-01","10-11","10-21",
                 "11-01","11-11","11-21",
                 "12-01","12-11","12-21")
season_end = c("01-10","01-20","01-31",
               "02-10","02-20","02-28",
               "03-10","03-20","03-31",
               "04-10","04-20","04-30",
               "05-10","05-20","05-31",
               "06-10","06-20","06-30",
               "07-10","07-20","07-31",
               "08-10","08-20","08-31",
               "09-10","09-20","09-30",
               "10-10","10-20","10-31",
               "11-10","11-20","11-30",
               "12-10","12-20","12-31")
n_seasons = length(season_end)
t_steps = n_seasons * n_yrs


qry <- "SELECT * FROM base.tbl_analysis_grid_cov_seaice WHERE cell>10000"
sea_ice = sf::st_read(con,query=qry)

Cell_IDs = grid_sf$cell

### average sea ice concentration values by cell, year, 10 day interval
start_date = end_date = rep('NA',t_steps)
Ice = matrix(NA,n_cells,t_steps)
for(iyr in 1:n_yrs){
  for(iseason in 1:n_seasons){
    date1 = paste0(iyr+start_yr-1,'-',season_start[iseason])
    date2 = paste0(iyr+start_yr-1,'-',season_end[iseason])
    Temp_ice = sea_ice[which(sea_ice$fdate>=date1 & sea_ice$fdate<=date2),]
    for(icell in 1:n_cells){
      Ice[icell,n_seasons*(iyr-1)+iseason] <- mean(Temp_ice[which(Temp_ice$cell==Cell_IDs[icell]),"rast_seaice"],na.rm=TRUE)
    }
  }
}

save.image(file='create_movement_rasters_ice_attached.RData')


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

save.image(file='analysis_grid_ice_kriged.RData')

#convert kriged ice data into raster brick
n_layers = n_yrs*n_seasons
Raster_list = vector("list",n_layers)
for(i in 1:n_layers){
  grid_sf$ice=Ice[,i]
  grid_sf$ice[which(grid_sf$ice>1)]=1
  ice_sp = as(grid_sf[,"ice"],"Spatial")
  Raster_list[[i]]<-raster(ice_sp,res=25067.53)
  Raster_list[[i]] =rasterize(ice_sp,Raster_list[[i]],field="ice")
  Raster_list[[i]][is.na(Raster_list[[i]][])]=0
}
Ice_brick = brick(Raster_list)

# now, we need to link movement data with our covariate data; we'll do this
# for ice by giving the layers of the ice brick a "date" and associating
# telemetry data with one of these potential dates
Dates = Last = rep("",n_layers)
i=0
for(iyr in 1:n_yrs){
  for(iseason in 1:n_seasons){
    i=i+1
    Dates[i] = paste0(iyr+start_yr-1,'-',season_start[iseason])
    Last[i] = paste0(iyr+start_yr-1,'-',season_end[iseason])
  }
}
Ice_brick = setZ(Ice_brick,z=Dates,name="date")


BD_df = data.frame(BD_telem)
BD_df = BD_df[,c("deployid","datetime")]
Coords = st_coordinates(BD_telem)
BD_df$X = Coords[,"X"]
BD_df$Y = Coords[,"Y"]
BD_df$date = Dates[1]
BD_dates = as.Date(BD_df$datetime)
BD_dates[which(BD_dates=="2012-02-29")]="2012-02-28"

for(irec in 1:nrow(BD_df)){
  Match = (BD_dates[irec]>=Dates & BD_dates[irec]<=Last)
  BD_df$date[irec]=Dates[which(Match==1)]
}

#save.image('data_for_Brett.RData')
writeRaster(Ice_brick,"Ice_brick.grd",overwrite=T)
writeRaster(dist_land_raster,"dist_land_raster.grd",overwrite=T)
writeRaster(Northing_raster,"Northing_raster.grd",overwrite=T)
writeRaster(Easting_raster,"Easting_raster.grd",overwrite=T)
writeRaster(Land_raster,"Land_raster.grd",overwrite=T)
writeRaster(depth_raster,"depth_raster.grd",overwrite=T)


save(BD_df,file='bd_telem.RData')





