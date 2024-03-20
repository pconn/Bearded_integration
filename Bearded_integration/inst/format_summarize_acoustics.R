# summarize/format acoustic detection data
load('prediction_grid_ice_kriged.RData')

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
Acc_detect$season = 1*(Month %in% c(1:3,12))+2*(Month %in% c(4:5))+3*(Month %in% c(6:7))+4*(Month %in% c(8:9))+5*(Month %in% c(10:11))
#limit to spring
#Acc_detect = Acc_detect[which(Acc_detect$season==2),]
#Acc_detect$tstep = (Year-2004)*5+Acc_detect$season
Year = lubridate::year(Acc_detect[,"detection_dt"])
Acc_detect$tstep=Year-2003  #2004 = time step 1
library(sf)
moorings_sf = st_transform(moorings_sf,st_crs(grid_sf))

#take out any moorings not within estimation grid
InGrid = st_within(moorings_sf,st_union(grid_sf),sparse=FALSE) 
moorings_sf = moorings_sf[which(InGrid==TRUE),]
Acc_detect = Acc_detect[which(Acc_detect$mooring_site_id %in% unique(moorings_sf$mooring_site_id_full)),]

#take out any moorings without acoustic detections
Unique_moorings_det = unique(Acc_detect$mooring_site_id)
moorings_sf = moorings_sf[which(moorings_sf$mooring_site_id_full %in% Unique_moorings_det),]

#produce plot of moorings with labels

#produce table of number of files per year by mooring
Acc_df = as.data.frame(Acc_detect)
Acc_df$year = lubridate::year(Acc_df$detection_dt)
library(doBy)
Acc_table = summaryBy(num_png_with_call~mooring_site_id+year,data=Acc_df,FUN=sum)
Acc_table = summaryBy(num_png_with_call~mooring_site_id+year,data=Acc_df,FUN=sum)

first_yr = min(Acc_df$year)
last_yr = max(Acc_df$year)
Unique_IDs = unique(Acc_table$mooring_site_id)


#produce plot of moorings
library(ptolemy)
alaska_bbox = st_bbox(ptolemy::alaska_bbox)

xmin = alaska_bbox$xmin + 1700000 #- 2000000
xmax = alaska_bbox$xmax -1200000
ymin = alaska_bbox$ymin +300000
ymax = alaska_bbox$ymax - 400000

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


library(ggplot2)
moorings_plot = ggplot()+geom_sf(data=land_sf)+geom_sf_text(data=moorings_sf,size=2,aes(label=mooring_site_id_full))+
  theme(axis.text=element_blank())+xlab('')+ylab('')

png('mooring_labels_plot.png',  width     = 6,
    height    = 6,
    units     = "in",
    res       = 1200,
    pointsize = 4)
  moorings_plot
dev.off()


Data_mat = matrix(0,length(Unique_IDs),last_yr-first_yr+1)
rownames(Data_mat)=Unique_IDs
colnames(Data_mat)=first_yr:last_yr
for(iid in 1:length(Unique_IDs)){
  for(iyr in first_yr:last_yr){
    Cur_data = Acc_df[which(Acc_df$mooring_site_id==Unique_IDs[iid] & Acc_df$year==iyr),]
    if(nrow(Cur_data)>0)Data_mat[iid,iyr-first_yr+1] = sum(Cur_data$num_png_with_effort)
  }
}

library(xtable)
xtable(Data_mat)


