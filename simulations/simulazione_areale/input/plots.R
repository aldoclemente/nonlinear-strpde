rm(list=ls())
library(fdaPDE)
library(ggplot2)
library(viridis)

mesh_dir = "mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)

n_bins = 10
exact_coeff = as.matrix(readMM("fisher_kpp.mtx"))
diff_coeff = as.matrix(readMM("diffusion.mtx"))
#adv_coeff = as.matrix(readMM("adv_diff.mtx"))
range(exact_coeff)
range(diff_coeff)
#range(adv_coeff)

grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)

n_times = ncol( exact_coeff )
exact = matrix(0, nrow=length(grid$x), ncol=n_times)
diff = exact
#adv = exact
psi_test = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, as.matrix(grid))
exact = psi_test%*%exact_coeff
diff = psi_test%*%diff_coeff
#adv = psi_test%*%adv_coeff
frame = 1:ncol(exact_coeff)
lims = c(0,1)
{
  pdf(paste0("fisher_kpp.pdf"))
  for(t in frame){
    data <- data.frame(x = grid$x, y = grid$y, z = exact[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = n_bins) +
      scale_fill_viridis_c(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}


{
  pdf(paste0("diffusion.pdf"))
  for(t in frame){
    data <- data.frame(x = grid$x, y = grid$y, z = diff[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = n_bins) +
      scale_fill_viridis_c(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

{
  pdf(paste0("adv_diff.pdf"))
  for(t in frame){
    data <- data.frame(x = grid$x, y = grid$y, z = adv[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = n_bins) +
      scale_fill_viridis_c(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

max(abs(exact_coeff- diff_coeff))

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}

n_breaks = 20

locs= as.matrix(readMM("1000/0/locs.mtx"))
obs = as.matrix(readMM("1000/0/obs.mtx"))
#obs = (obs - min(obs)) / (max(obs) - min(obs))
range(obs)
{ 
  mycols <- mycolors(n_breaks + 2)
  vals   <- seq(lims[1], lims[2], length.out = length(mycols))
  
  linewidth=2
  pdf("pointwise_data.pdf")
  for(t in frame){
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


plot(locs[,1], locs[,2])  
