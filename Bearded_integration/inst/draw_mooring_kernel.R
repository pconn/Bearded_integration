### plot acoustic kernels for bearded seal integration methods white paper

library(sf)

box = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))))
grid = st_make_grid(box,cellsize=0.25)
point = st_centroid(box)
circle = st_buffer(point,dist=0.25)

library(ggplot2)
plot1 = ggplot()+geom_sf(data=grid)+
  geom_sf(data=circle,color="blue",alpha=0.3)+
  geom_sf(data=point,color="blue",size=2)+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )

#now plot kernel
pdf('hypo_grid_mooring.pdf')
 plot1
dev.off()

#now plot kernel
png('hypo_grid_mooring.png')
plot1
dev.off()

library(mvtnorm)

set.seed(0)
easting <- seq(0, 1, 0.01)
northing <- seq(0, 1, 0.01)
mean <- c(0.45, 0.45)
cov <- matrix(c(.0765^2, 0, 0, .0765^2), nrow=2)
f <- function(easting, northing) dmvnorm(cbind(easting, northing), mean, cov)
density <- outer(easting, northing, f)

density_df = data.frame(matrix(0,10000,3))
colnames(density_df)=c("easting","northing","density")
for(i in 1:100){
  for(j in 1:100){
    density_df$easting[(i-1)*100+j]=easting[i]
    density_df$northing[(i-1)*100+j]=northing[j]
    density_df$density[(i-1)*100+j]=density[i,j]
  }
}

density_df$density = density_df$density/25 
#countour plot
contour_plot = ggplot(density_df, aes(easting, northing)) + 
  geom_contour(aes(z = density, colour = after_stat(level)),bins=20)+
  xlim(c(.1,.9))+ylim(c(.1,.9))+theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank())+
  geom_vline(xintercept = c(0.125,0.375,0.625,0.875),col="darkgray",size=1.2)+
  geom_hline(yintercept = c(0.125,0.375,0.625,0.875),col="darkgray",size=1.2)+
  geom_point(aes(x=0.45,y=0.45),color='darkred',size=2)+
  scale_color_viridis_c(name="Density")+
  xlab('Easting')+ylab('Northing')
  
png('kernel_2d_contour.png',width=6,height=6,units="in",res=600)
 contour_plot
dev.off()

#create surface plot

res1 <- persp(data.lm,x~y, zlim=c(0,max(z)),contour=list(z="bottom",col="colors"),theta=-55,phi=25)    
xy <- matrix(c((-3-8)/5,-3,(3-8)/5,3),ncol=2,byrow = T)
lines(trans3d(xy[,2], xy[,1], 0, pmat = res1$`y ~ x`$transf), col = 3)

library(rsm)
png('kernel_2d.png')
persp(easting, northing, density, theta=-20, phi=20, col = 'blue',
      expand=0.8, ticktype='simple', contour=list(z="bottom",col="colors"))
dev.off()
