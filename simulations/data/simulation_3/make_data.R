rm(list=ls())
library(Matrix)
library(fdaPDE)
data_dir = "../simulation_3/"
mesh_dir = "../simulation_1/mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

dist = function(x,y){
  return( sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2) )
}

sds = c("0.00", "0.05", "0.10")

nodes = as.matrix(readMM(paste0(mesh_dir,"points.mtx")))
cells = as.matrix(readMM(paste0(mesh_dir, "cells.mtx")))
boundary = as.matrix(readMM(paste0(mesh_dir, "boundary.mtx")))
mesh = create.mesh.2D(nodes, triangles = (cells + 1)) 
plot(mesh)
FEMbasis = create.FEM.basis(mesh)

h = abs(mesh$nodes[1,1] - mesh$nodes[2,1]) 
a=0
b=1
x = seq(a + 2*h, b-2*h, by=3*h)
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

#L = sqrt(2)*0.625
point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))

for( r in 1:nrow(incidence_matrix)){
  apply(mesh$triangles[which(incidence_matrix[r,] == 1),], 
        MARGIN=1, FUN= function(x){
          polygon(mesh$nodes[x,],col="blue")})
  
}

# tmp_coords = mesh$nodes[mesh$triangles[which(incidence_matrix[1,]==1),],]
# tmp_coords2 = mesh$nodes[as.vector(mesh$triangles[which(incidence_matrix[1,]==1),]),]
# tmp_coords - tmp_coords2

n_dofs = nrow(nodes)
FEMbasis = create.FEM.basis(mesh)
Mass = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
one = as.matrix(rep(1, times=n_dofs))
D = as.matrix(rep(0, times=nrow(centers))) 
I =  Mass%*% one
for(r in 1:nrow(centers)){
  D[r,] = sum(I[unique( as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),]))])
}

time_locs = as.matrix(readMM("../simulation_2/time_mesh.mtx"))
time_locs = as.matrix(time_locs[-1])
nsim=30

n_time_locs = length(time_locs)

set.seed(0)
for(k in 1:length(sds)){
  sd = sds[k]
  sigma_dir = paste0(data_dir, "sigma_", sds[k], "/")
  for(i in 0:(nsim-1)){
    sim_dir = paste0(sigma_dir, i, "/")
    if(!dir.exists(sim_dir)) dir.create(sim_dir, recursive = T)
    
    exact_dir = paste0("../simulation_2/", "sigma_", sds[k], "/", i, "/")
    exact  = as.matrix(readMM(paste0(exact_dir, "exact.mtx")))
  
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
  }
}

write.csv(format(incidence_matrix, digits=16), file = "incidence_matrix.csv")
