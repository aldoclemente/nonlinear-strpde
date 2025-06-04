library(fdaPDE)
library(Matrix)
library(ggplot2)
library(viridis)
datadir = "data/"

time_locs = as.matrix(read.table("../simulation_0_dirichlet/data/time_locations.txt"))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("../simulation_0_dirichlet/data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("../simulation_0_dirichlet/data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("../simulation_0_dirichlet/data/mesh/elements.mtx"))
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 
FEMbasis = create.FEM.basis(mesh)
nsim = 30
IC = as.matrix( as.vector(read.table("../simulation_0_dirichlet/data/exact.txt")[,1]))
PDE_parameters = list(K=matrix(c(0.1,0,0,0.1), nrow=2,ncol=2), b = c(0,0), c=3.0)
#lambda = 1.
incidence_matrix = as.matrix(read.table("data/incidence_matrix.txt"))
lambda = 1. / (nrow(incidence_matrix) * 11)
PDE_parameters = list(K=matrix(c(0.1,0,0,0.1), nrow=2,ncol=2), b = c(0,0), c=3.0)

BC_indices = boundary %in% c(4)
BC_values = rep(0,times=)
plot(mesh, pch=".")
points(mesh$nodes[BC_indices,], pch=16)

dirichlet_bc = which(boundary==4)
BC= list(
  BC_indices = rep(dirichlet_bc, each=(n_time_locs-1)), 
  BC_values = rep(0, length=(length(dirichlet_bc)*(n_time_locs-1)))
)

for(i in 0:(nsim-1)){
  obs = as.matrix(read.table(paste0(datadir, "/", i, "/", "obs.txt")))
  
  invisible( capture.output(parabolic <- smooth.FEM.time(incidence_matrix = incidence_matrix, lambdaS = lambda, lambdaT = 1,
                                                         time_mesh = as.vector(time_locs), PDE_parameters = PDE_parameters,
                                                         observations = obs[,2:n_time_locs], FEMbasis = FEMbasis, 
                                                         IC = IC, BC = BC,
                                                         FLAG_PARABOLIC = T)))
  
  coeff = matrix(parabolic$fit.FEM$coeff, nrow=n_dofs, ncol=n_time_locs)
  write.table( coeff, 
               paste0(datadir, "/", i, "/", "parabolic.txt"),
               col.names = F, row.names = F)
  
}

f_parabolic = matrix(0,nrow=n_dofs, ncol=n_time_locs)
for(i in 0:(nsim-1)){
  f_parabolic <- f_parabolic + as.matrix(read.table(paste0(datadir, "/", i, "/", "parabolic.txt"), header = F))/nsim
}

write.table(f_parabolic, row.names = F, "data/parabolic.txt")


