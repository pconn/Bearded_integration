### EDA bearded telemetry data
library(sf)
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

#summary statistics for numbers of seals, records per year, season 
min_yr = min(BD_unlist$year)
max_yr = max(BD_unlist$year)
n_yrs = max_yr-min_yr+1
Month_year = matrix(0,n_yrs,12)
colnames(Month_year)=as.character(c(1:12))
rownames(Month_year)=as.character(c(min_yr:max_yr))

Points_m_y = Month_year
Seals_m_y = Month_year

for(irow in 1:nrow(BD_unlist)){
  Points_m_y[BD_unlist$year[irow]-min_yr+1,BD_unlist$month[irow]]=Points_m_y[BD_unlist$year[irow]-min_yr+1,BD_unlist$month[irow]]+1
}

for(iseal in 1:n_seal){
  BD_move[[4]][[iseal]] <- BD_move[[4]][[iseal]] %>% 
                         mutate(year = lubridate::year(datetime)) %>% 
                         mutate(month = lubridate::month(datetime)) %>%
                         mutate(my=paste(month,year))
  Unique_my = unique(BD_move[[4]][[iseal]]$my)
  for(imy in 1:length(Unique_my)){
    cur_pl = which(BD_move[[4]][[iseal]]$my==Unique_my[imy])[1]
    Seals_m_y[BD_move[[4]][[iseal]]$year[cur_pl]-min_yr+1,BD_move[[4]][[iseal]]$month[cur_pl]]=
      Seals_m_y[BD_move[[4]][[iseal]]$year[cur_pl]-min_yr+1,BD_move[[4]][[iseal]]$month[cur_pl]]+1
  }
}

#Records by month, year

Points_m_y

#Records by month

colSums(Points_m_y)

#Records by year

rowSums(Points_m_y)

#Seals by month, year

Seals_m_y

#Seal-years by month

colSums(Seals_m_y)

xtable::xtable(Points_m_y,digits=0)


xtable::xtable(Seals_m_y,digits=0)

unique_IDs = unique(BD_unlist$deployid)


# now look to see which of these are AFSC, NSB, ADF&G
library("RPostgreSQL")
library("sf")

# connect to DB 
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))
#query table structure of 'telem'
dbGetQuery(con,
           "SELECT table_name FROM information_schema.tables  
                   WHERE table_schema='telem'")

# let's try comparing capture spenos to the ones we have in Josh's dataset
nsb_capture <- sf::st_read(con,query="SELECT * FROM telem.nsb_capture")
nsb_capture$speno %in% unique_IDs
#hmm, no matches (?)  - potentially 2 more bearded seals?

adfg_capture <- sf::st_read(con,query="SELECT * FROM telem.adfg_capture")
adfg_capture$deployid %in% unique_IDs  #none here either; 24 more bearded seals?


#compute kernel SD from bearded seals in our dataset
Month = lubridate::month(BD_telem$datetime)
Telem_spring <- BD_telem[which(Month==4 | Month==5),]
library(adehabitatHR)
library(sp)
#convert to SpPointsDF (needed for adehabitatHR)
Telem_spdf = as(Telem_spring,"Spatial")
unique(paste(Telem_spdf$deployid,Telem_spdf$year))  #there aren't any IDs that span years
Telem_spdf$id=Telem_spdf$deployid
Telem_spdf = Telem_spdf[,"id"]
kernel.ref <- kernelUD(Telem_spdf, h = "href")  # href = the reference bandwidth
image(kernel.ref) # plot

n_indiv = length(unique(Telem_spdf$id))
kernel_SDs = Days = rep(0,n_indiv)
Unique_ids = unique(Telem_spdf$id)
for(i in 1:n_indiv){
  kernel_SDs[i]=kernel.ref[[i]]@h$h
  Cur_records = Telem_spring[which(Telem_spring$deployid==Unique_ids[i]),]
  Days[i] = max(lubridate::yday(Cur_records$datetime))-min(lubridate::yday(Cur_records$datetime))+1
}
lm(kernel_SDs~Days)
prediction = 7032 + 55.3*60

