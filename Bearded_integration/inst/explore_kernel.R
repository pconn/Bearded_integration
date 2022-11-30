#### look at kernel inverse for ice --> sea

kern_sd = 2
n_cells = 225

n_x = sqrt(n_cells)

Ice_g=matrix(0,n_x,n_x)

for(i in 1:5){
  for(j in 1:5){
    Ice_g[i,j] = 1
  }
}
Pi = Ice_g /sum(Ice_g)
N = 100
N_g = N*Pi


XY = matrix(0,n_cells,2)
Dist = Psi = matrix(0,n_cells,n_cells)
counter=1
for(x in 1:n_x){
  for(y in 1:n_x){
     XY[counter,1]=x
     XY[counter,2]=y
     counter=counter+1
  }
}

for(icell1 in 1:n_cells){
  for(icell2 in 1:n_cells){
    Dist[icell1,icell2]=sqrt((XY[icell1,1]-XY[icell2,1])^2+(XY[icell1,2]-XY[icell2,2])^2) #Euclid norm
    Psi[icell1,icell2]=dnorm(Dist[icell1,icell2],0,kern_sd)
  }
  Psi[icell1,]=Psi[icell1,]/sum(Psi[icell1,])
}

N_ice = as.vector(N_g)
N_sea = N_ice %*% Psi


Plot_df = data.frame("x"=XY[,1],"y"=XY[,2],"N_ice"=N_ice,"N_sea"=as.numeric(N_sea))


library(ggplot2)
ggplot(Plot_df)+geom_raster(aes(x=x,y=y,fill=N_ice))
ggplot(Plot_df)+geom_raster(aes(x=x,y=y,fill=N_sea))

Psi = Psi+0.0001
for(icell in 1:n_cells)Psi[icell,]=Psi[icell,]/sum(Psi[icell,])

Psi_inv = solve(Psi)

Plot_df$Est = as.vector(as.vector(N_sea) %*% Psi_inv)

ggplot(Plot_df)+geom_raster(aes(x=x,y=y,fill=Est))

### take away : inverse kernels very sensitive to departures from expectation


Plot_df$Crap = as.numeric(Plot_df$Est %*% Psi)
ggplot(Plot_df)+geom_raster(aes(x=x,y=y,fill=Crap))



