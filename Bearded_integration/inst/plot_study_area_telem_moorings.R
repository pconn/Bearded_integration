# plot study area, telemetry, acoustic mooring locations


library(ptolemy)
library(ggplot2)
library(sf)
alaska_bbox = st_bbox(ptolemy::alaska_bbox)

xmin = alaska_bbox$xmin - 2000000
xmax = alaska_bbox$xmax
ymin = alaska_bbox$ymin +100000
ymax = alaska_bbox$ymax + 300000

lon = c(xmin, xmax)
lat = c(ymin,ymax)
Poly_Coord_df = data.frame(lon, lat)

pol = st_polygon(
  list(
    cbind(
      Poly_Coord_df$lon[c(1,2,2,1,1)], 
      Poly_Coord_df$lat[c(1,1,2,2,1)])
  )
)
polc = st_sfc(pol, crs=st_crs(alaska_bbox))
land_sf = extract_gshhg(data = polc)


#input grid
load('prediction_grid_ice_kriged.RData')

#telemetry data
load("BD_telem_NOAA.RData")
BD_telem$Season="Winter"
BD_telem[which(BD_telem$month %in% c(6:11)),"Season"]="Summer"
BD_telem$Season = factor(BD_telem$Season,levels=c("Summer","Winter"))

# BD_telem[which(BD_telem$month %in% c(4,5)),"Season"]="Spring"
# BD_telem[which(BD_telem$month %in% c(6,7)),"Season"]="Early_open"
# BD_telem[which(BD_telem$month %in% c(8:9)),"Season"]="Open_water"
# BD_telem[which(BD_telem$month %in% c(10:11)),"Season"]="Fall"
# BD_telem$Season = factor(BD_telem$Season,levels=c("Spring","Early_open","Open_water","Fall","Winter"))



#Which_plot = c(1:10764)*10
#Plot_telem = BD_telem[Which_plot,]
Plot_telem = BD_telem

text_df = data.frame(x=c(-800000,1250000),y=c(-2600000,-2250000),Location=c("Russia","Alaska,USA") )

#acoustic mooring locations
yr_last = 2021 #will want to change if updating Ice year range

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
Acc_detect$season = 1*(Month %in% c(1:3,12))+2*(Month %in% c(4:5))+3*(Month %in% c(6:7))+4*(Month %in% c(8:9))+5*(Month %in% c(10:11))
#limit to spring
#Acc_detect = Acc_detect[which(Acc_detect$season==2),]
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


#one plot to rule them all - note the compass position
#entirely depends on teh size of the plot.  
my_plot = ggplot() + geom_sf(data=grid_sf,color='lightblue')+
  geom_sf(data=land_sf,fill='lightgray') +
  geom_sf(data=Plot_telem,aes(color=Season),size=0.1)+
  geom_sf(data=moorings_sf,color='darkred')+
  geom_text(data=text_df,aes(x=x,y=y,label=Location),size=3)+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggspatial::annotation_scale(location = "bl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book")+
  ggspatial::annotation_north_arrow(
    which_north = "grid",
    pad_x = unit(1.3, "in"), pad_y = unit(4.3, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),line_width=0.1,
      line_col = "grey20",
      text_family = "ArcherPro Book"))+
  xlim(c(-1600000,1540000))+
  ylim(c(-4000000,-1400000))+ #alaska_bbox$ymax))
  xlab('Easting')+ylab('Northing')

png("study_area_telem_moorings_2season.png",width=8,height=8,units="in",res=1200)
my_plot
dev.off()