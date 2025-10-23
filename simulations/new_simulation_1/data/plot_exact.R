
library(fdaPDE)
mesh_dir = "mesh/"
time_mesh = as.matrix(readMM("mesh/time_mesh.mtx"))*2e-2
n_times = length(time_mesh)
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)
coeff = as.matrix(readMM("mesh/exact.mtx"))
grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)


exact = matrix(0, nrow=length(grid$x), ncol=n_times)
for(t in 1:n_times){
  exact[,t] = eval.FEM(FEM(coeff[,t], FEMbasis), locations = grid)
}
range(coeff)
range(coeff[,1])
range(coeff[,2])
range(coeff[,n_times])

par = as.matrix(readMM("mesh/parabolic.mtx"))
max(abs(par - coeff))

parabolic = matrix(0, nrow=length(grid$x), ncol=n_times)
for(t in 1:n_times){
  parabolic[,t] = eval.FEM(FEM(par[,t], FEMbasis), locations = grid)
}

locs = as.matrix(readMM("100/0/locs.mtx"))
obs = as.matrix(readMM("100/0/obs.mtx"))
lims = range( c(range(exact), range(parabolic), range(obs))) 

{  
  pdf("exact.pdf")
  for(t in 1:n_times){
    data <- data.frame(x = grid$x, y = grid$y, z = exact[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = 10) +
      scale_fill_viridis_c(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

{  
  pdf("parabolic.pdf")
  for(t in 1:n_times){
    data <- data.frame(x = grid$x, y = grid$y, z = parabolic[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = 10) +
      scale_fill_viridis_c(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

difference = abs(parabolic-exact)
{  
  pdf("difference.pdf")
  for(t in 1:n_times){
    data <- data.frame(x = grid$x, y = grid$y, z = difference[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = 10) +
      scale_fill_viridis_c(limits=range(difference)) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

range(difference)

{ 
  n_breaks = 20
  
  mycolors<- function(x) {
    colors<-viridis( x + 1 )
    colors[1:x]
  }
  
  mycols <- mycolors(n_breaks + 2)
  vals   <- seq(lims[1], lims[2], length.out = length(mycols))
  
  linewidth=2
  locs = as.matrix(readMM("100/0/locs.mtx"))
  obs = as.matrix(readMM("100/0/obs.mtx"))
  pdf(paste0("100/pointwise_data.pdf"))
  for(t in 1:n_times){
    data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
    plt <- ggplot(data, aes(x, y, fill = z)) + 
      geom_rect(xmin = 0, ymin = 0, xmax = 1, ymax = 1,
                fill = "transparent",
                color = "black",
                linewidth = 1) +
      geom_point(shape = 21, size = 5, color = "transparent") + 
      scale_fill_viridis_c(limits = lims, oob = scales::squish) +
      coord_fixed() +
      theme_void() +
      theme(legend.position = "none")
    
    print(plt)
  }
  dev.off()
}


{
  lwd = 4
  cex.axis = 2.5
  idxs = c(1,2,3)
  colors = viridis(length(idxs))
  psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, locs)
  exact = psi%*%coeff
  system("mkdir -p imgs")
  pdf("imgs/trend.pdf", width=14)
  plot(times, obs[1,], type="l", lwd=1, col="gray", 
       ylim=c(-0.05, 1.15), frame.plot = F, xlab = "", ylab="", 
       cex.axis = cex.axis)
  for(n in 2:nrow(exact)){
    points(times, obs[n,], type="l", lwd=1, col="gray")
  }
  for(n in 1:length(idxs)){
    points(times, obs[idxs[n],], type="l", lwd=lwd, col=colors[idxs[n]])
  }
  
  plot(times, exact[1,], type="l", lwd=1, col="gray", 
       ylim=c(-0.05, 1.15), frame.plot = F, xlab = "", ylab="", 
       cex.axis = cex.axis)
  for(n in 2:nrow(exact)){
    points(times, exact[n,], type="l", lwd=1, col="gray")
  }
  for(n in 1:length(idxs)){
    points(times, exact[idxs[n],], type="l", lwd=lwd, col=colors[idxs[n]])
  }
  
  dev.off()
}

