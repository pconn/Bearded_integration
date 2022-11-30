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

library(mvtnorm)

set.seed(0)
x1 <- seq(0, 1, 0.01)
x2 <- seq(0, 1, 0.01)
mean <- c(2/3, 1/3)
cov <- matrix(c(.0765^2, 0, 0, .0765^2), nrow=2)
f <- function(x1, x2) dmvnorm(cbind(x1, x2), mean, cov)
y <- outer(x1, x2, f)

#create surface plot
pdf('kernel_2d.pdf')
persp(x1, x2, y, theta=-20, phi=20, col = 'blue',
      expand=0.8, ticktype='simple')
dev.off()


