rm(list=ls())
library(Matrix)
library(fdaPDE)

mesh_dir = "mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

dist = function(x,y){
  return( sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2) )
}


mesh = create.mesh.2D(nodes, triangles = (cells + 1)) 
plot(mesh, pch=".")
FEMbasis = create.FEM.basis(mesh)

h = abs(mesh$nodes[1,1] - mesh$nodes[2,1]) 
a=0
b=1
x = seq(a + h, b-h, by=2*h)
y = x

centers = expand.grid(x,y)
plot(mesh, pch=".")
points(centers, pch=16, col="blue")
incidence_matrix = matrix(0, nrow=nrow(centers), ncol=nrow(mesh$triangles))

plot(mesh,pch=".")
points(centers, col="blue", pch=16)
for(r in 1:nrow(centers)){
  for(e in 1:nrow(mesh$triangles)){
    point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))
    points(t(point), col="black")
    if( dist(point, t(centers[r,]) )  <=  h   ){
      incidence_matrix[r,e] = 1 
    }
  }
}

point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))

cols = 1:nrow(centers)
for( r in 1:nrow(incidence_matrix)){
  apply(mesh$triangles[which(incidence_matrix[r,] == 1),], 
        MARGIN=1, FUN= function(x){
          polygon(mesh$nodes[x,],col=cols[r])})
  
}

head(centers)
tail(centers)
# tmp_coords = mesh$nodes[mesh$triangles[which(incidence_matrix[1,]==1),],]
# tmp_coords2 = mesh$nodes[as.vector(mesh$triangles[which(incidence_matrix[1,]==1),]),]
# tmp_coords - tmp_coords2

n_dofs = nrow(nodes)
FEMbasis = create.FEM.basis(mesh)
Mass = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
one = as.matrix(rep(1, times=n_dofs))
D = as.matrix(rep((2*h)^2 , times=nrow(centers))) 

#sum( I[as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),]  ))

time_locs = as.matrix(readMM("mesh/time_mesh.mtx"))
#time_locs = as.matrix(time_locs[-1])
nsim=30

n_time_locs = length(time_locs)

set.seed(0)
exact = readMM("fisher_kpp.mtx")

for(i in 0:(nsim-1)){
sim_dir = paste0(i, "/")
if(!dir.exists(sim_dir)) dir.create(sim_dir)
int_exact = matrix(0, nrow=nrow(centers), ncol=n_time_locs)
for(r in 1:nrow(centers)){
    for(t in 1:n_time_locs){
        I = Mass %*% exact[,t]
        int_exact[r,t]= 1/D[r,]*sum(I[unique( as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),]))])
      }
    }
    
    observations = int_exact + rnorm(nrow(int_exact)*ncol(int_exact), sd = 0.05) # 0.04 circa ==  0.05*diff(range(abs(int_exact)))
    observations = as.matrix(observations)
    writeMM(Matrix(observations, sparse=T), paste0(sim_dir,"obs.mtx"))
    writeMM(Matrix(observations[,1], sparse=T), paste0(sim_dir,"obs.t0.mtx"))
}


# prova
1/D[r,]*sum(I[as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),])])

write.csv(format(incidence_matrix, digits=16), file = "incidence_matrix.csv")

library(Matrix)
incidence_matrix = read.csv("incidence_matrix.csv")
observations = readMM("0/obs.mtx")
library(ggplot2)
{
  poly_coords = list()
  for(r in 1:nrow(incidence_matrix)){
    triangles = mesh$triangles[which(incidence_matrix[r,] == 1),]
    #triangles = cbind(triangles, triangles[,1]) 
    poly_coords[[r]] = list()
    for(t in 1:nrow(triangles)){
      poly_coords[[r]][[t]] = mesh$nodes[triangles[t,], ]
      poly_coords[[r]][[t]] = rbind(poly_coords[[r]][[t]], 
                                    mesh$nodes[triangles[t,][1], ])
    }
  }
}
 
pdf("areal_data.pdf")
library(viridis)
for (t in seq_len(ncol(exact))) {
  
  ## Build one big data.frame with all polygons for this time t
  poly_list <- list()
  idx <- 1
  
  for (r in seq_len(nrow(incidence_matrix))) {
    for (i in seq_along(poly_coords[[r]])) {
      
      coords <- poly_coords[[r]][[i]]
      # coords is assumed to be a matrix with 2 columns (x, y)
      
      poly_list[[idx]] <- data.frame(
        x     = coords[, 1],
        y     = coords[, 2],
        z     = observations[r, t],
        group = paste(r, i, sep = "_")  # group polygons
      )
      idx <- idx + 1
    }
  }
  
  poly_df <- do.call(rbind, poly_list)
  
  plt_rect <- ggplot(poly_df, aes(x = x, y = y, group = group, fill = z)) +
    geom_polygon(color = NA) +
    scale_fill_viridis(
      option = "viridis",
      limits = range(observations),
      oob    = scales::squish
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    geom_rect(
      data = data.frame(xmin = 0, ymin = 0, xmax = 1, ymax = 1),
      aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
      fill = "transparent",
      color = "black",
      linewidth = 1,
      inherit.aes = FALSE
    )
  
  print(plt_rect)
}

dev.off()







