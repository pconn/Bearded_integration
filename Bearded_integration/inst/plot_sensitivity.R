#plot sensitivity run differences

load('bearded_SPDErr_results_2season.RData')

N_df = data.frame("Data"=c("Acoustics","Aerial","CKMR","UD"),"Omit"=rep(0,4),"More"=rep(0,4),
                  "Omit2" = rep(0,4), "More2"=rep(0,4))
N_base = Report$N

grid_winter2012 = grid_summer2012 = grid_winter2019 = grid_summer2019 = grid_sf
grid_winter2012$Z_base = Report$Z_st[,17]
grid_summer2012$Z_base = Report$Z_st[,18]
grid_winter2019$Z_base = Report$Z_st[,31]
grid_summer2019$Z_base = Report$Z_st[,32]

Z_base = Report$Z_st


load('bearded_SPDErr_noAcoustics.RData')
N_df$Omit[1] = (Report$N-N_base)/N_base
N_df$Omit2[1] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$Omit_acoustics = Report$Z_st[,17]
grid_summer2012$Omit_acoustics = Report$Z_st[,18]


load('bearded_SPDErr_noAerial.RData')
N_df$Omit[2] = (Report$N-N_base)/N_base
N_df$Omit2[2] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$Omit_aerial = Report$Z_st[,17]
grid_summer2012$Omit_aerial = Report$Z_st[,18]

load('bearded_SPDErr_noCKMR.RData')
N_df$Omit[3] = (Report$N-N_base)/N_base
N_df$Omit2[3] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$Omit_CKMR = Report$Z_st[,17]
grid_summer2012$Omit_CKMR = Report$Z_st[,18]

N_df$Omit[4]=NA
N_df$Omit2[4] = NA
grid_winter2012$Omit_UD = NA
grid_summer2012$Omit_UD = NA

load('bearded_SPDErr_MoreAcoustics.RData')
N_df$More[1] = (Report$N-N_base)/N_base
N_df$More2[1] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$More_acoustics = Report$Z_st[,17]
grid_summer2012$More_acoustics = Report$Z_st[,18]

load('bearded_SPDErr_MoreAerial10.RData')
N_df$More[2] = (Report$N-N_base)/N_base
N_df$More2[2] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$More_aerial = Report$Z_st[,17]
grid_summer2012$More_aerial = Report$Z_st[,18]

load('bearded_SPDErr_MoreCKMR.RData')
N_df$More[3] = (Report$N-N_base)/N_base
N_df$More2[3] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$More_CKMR = Report$Z_st[,17]
grid_summer2012$More_CKMR = Report$Z_st[,18]

load('bearded_SPDErr_MoreUD.RData')
N_df$More[4] = (Report$N-N_base)/N_base
N_df$More2[4] = sqrt(sum((Report$Z_st - Z_base)^2))
grid_winter2012$More_UD = Report$Z_st[,17]
grid_summer2012$More_UD = Report$Z_st[,18]

library(xtable)
xtable(N_df,digits=2)

#plots

#1) base model abundance plots for 2012, 2019 winter/summer
library(ggplot2)
library(viridis)
library(reshape2)
library(scales)

n_cells = nrow(grid_sf)
grid_plot = rbind(grid_sf,grid_sf,grid_sf,grid_sf)
grid_plot$N = c(grid_winter2012$Z_base,grid_summer2012$Z_base,grid_winter2019$Z_base,grid_summer2019$Z_base)
grid_plot$Season = rep(c("Winter","Open water","Winter","Open water"),each=n_cells)
grid_plot$Year = rep(c("2012","2019"),each=n_cells*2)
ggplot(grid_plot)+geom_sf(aes(fill=N),colour=NA)+facet_grid(Season~Year)+scale_fill_viridis(name="Abundance")
png("Seal_N_plot.png",width=6,height=6,res=600,units="in")
ggplot(grid_plot)+geom_sf(aes(fill=N),colour=NA)+facet_grid(Season~Year)+scale_fill_viridis(name="Abundance")
dev.off()


# 2) differences between base and sensitivity runs

# 2a) winter

grid_plot=grid_sf
for(i in 2:8)grid_plot=rbind(grid_plot,grid_sf)
grid_plot$RelDiff = c(grid_winter2012$Omit_acoustics,grid_winter2012$Omit_aerial,
                      grid_winter2012$Omit_CKMR,rep(NA,n_cells),
                      grid_winter2012$More_acoustics,grid_winter2012$More_aerial,
                      grid_winter2012$More_CKMR,grid_winter2012$More_UD)
grid_plot$Change = rep(c("-","+"),each=n_cells*4)
grid_plot$Data = rep(rep(c("Acoustics","Aerial","CKMR","UD"),each=n_cells),2)
Z_base = rep(grid_winter2012$Z_base,8)
grid_plot$RelDiff=(grid_plot$RelDiff-Z_base)/Z_base
grid_plot$RelDiff[which(grid_plot$RelDiff>2)]=2
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())

png("Sens_winter_diff.png",width=6,height=8,res=600,units="in")
ggplot(grid_plot)+geom_sf(aes(fill=RelDiff),colour=NA)+facet_grid(Data~Change)+
  scale_fill_gradientn(colors=c("blue","white","orange"),
                       values=rescale(c(-1,0,2)),limits=c(-1,2),
                       breaks=c(-1,0,1,2),labels=c("-100%","0%","100%",">200%"))+
  tmp.theme
dev.off()

# 2b) open water

grid_plot=grid_sf
for(i in 2:8)grid_plot=rbind(grid_plot,grid_sf)
grid_plot$RelDiff = c(grid_summer2012$Omit_acoustics,grid_summer2012$Omit_aerial,
                      grid_summer2012$Omit_CKMR,rep(NA,n_cells),
                      grid_summer2012$More_acoustics,grid_summer2012$More_aerial,
                      grid_summer2012$More_CKMR,grid_summer2012$More_UD)
grid_plot$Change = rep(c("-","+"),each=n_cells*4)
grid_plot$Data = rep(rep(c("Acoustics","Aerial","CKMR","UD"),each=n_cells),2)
Z_base = rep(grid_summer2012$Z_base,8)
grid_plot$RelDiff=(grid_plot$RelDiff-Z_base)/Z_base
grid_plot$RelDiff[which(grid_plot$RelDiff>2)]=2
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())

png("Sens_summer_diff.png",width=6,height=8,res=600,units="in")
ggplot(grid_plot)+geom_sf(aes(fill=RelDiff),colour=NA)+facet_grid(Data~Change)+
  scale_fill_gradientn(colors=c("blue","white","orange"),
                       values=rescale(c(-1,0,2)),limits=c(-1,2),
                       breaks=c(-1,0,1,2),labels=c("-100%","0%","100%",">200%"))+
  tmp.theme
dev.off()


#goodness of fit: absolute differences in abundance between aerial surveys 
# and integrated predictions

grid_sf$Diff_2012 = grid_sf$Diff_2013 = grid_sf$Diff_2016 = grid_sf$Diff_2021 = NA
grid_sf$Diff_2013[which_Bering]=(Z_base[which_Bering,19]-exp(log_2013_N))/exp(log_2013_N)
grid_sf$Diff_2012[which_Bering]=(Z_base[which_Bering,17]-exp(log_2012_N))/exp(log_2012_N)
grid_sf$Diff_2016[which_Chukchi]=(Z_base[which_Chukchi,25]-exp(log_2016_N))/exp(log_2016_N)
grid_sf$Diff_2021[which_Beaufort]=(Z_base[which_Beaufort,35]-exp(log_2021_N))/exp(log_2021_N)
grid_plot=rbind(grid_sf[,1],grid_sf[,1],grid_sf[,1],grid_sf[,1])
grid_plot$Year = rep(c(2012,2013,2016,2021),each=n_cells)
grid_plot$Diff = c(grid_sf$Diff_2012,grid_sf$Diff_2013,grid_sf$Diff_2016,grid_sf$Diff_2021)
ggplot(grid_plot)+geom_sf(colour=NA,aes(fill=Diff))+facet_wrap(vars(Year))+
  tmp.theme+  scale_fill_gradientn(colors=c("blue","white","orange"),
                                   values=rescale(c(-2000,0,500)),limits=c(-2000,500),
                                   breaks=c(-2000,-1500,-1000,-500,0,500))


grid_sf$Diff_2012 = grid_sf$Diff_2013 = grid_sf$Diff_2016 = grid_sf$Diff_2021 = NA
grid_sf$Diff_2013[which_Bering]=Z_base[which_Bering,19]-exp(log_2013_N)
grid_sf$Diff_2012[which_Bering]=Z_base[which_Bering,17]-exp(log_2012_N)
grid_sf$Diff_2016[which_Chukchi]=Z_base[which_Chukchi,25]-exp(log_2016_N)
grid_sf$Diff_2021[which_Beaufort]=Z_base[which_Beaufort,35]-exp(log_2021_N)
grid_plot=rbind(grid_sf[,1],grid_sf[,1],grid_sf[,1],grid_sf[,1])
grid_plot$Year = rep(c(2012,2013,2016,2021),each=n_cells)
grid_plot$Diff = c(grid_sf$Diff_2012,grid_sf$Diff_2013,grid_sf$Diff_2016,grid_sf$Diff_2021)
ggplot(grid_plot)+geom_sf(colour=NA,aes(fill=Diff))+facet_wrap(vars(Year))+
  tmp.theme+  scale_fill_gradientn(colors=c("blue","white","orange"),
                                   values=rescale(c(-2000,0,500)),limits=c(-2000,500),
                                   breaks=c(-2000,-1500,-1000,-500,0,500))

