library(fdaPDE)
library(Matrix)
library(ggplot2)
library(viridis)

datadir = "data/"

NLOCS = c(100, 250, 500, 1000)

time_locs = as.matrix(read.table(paste0(datadir, "time_locations.txt")))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("data/mesh/elements.mtx"))
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 
FEMbasis = create.FEM.basis(mesh)
nsim = 30
lambda = 1.
IC = as.matrix( as.vector(read.table(paste0(datadir, "exact.txt"))[,1]))
PDE_parameters = list(K=matrix(c(0.1,0,0,0.1), nrow=2,ncol=2), b = c(0,0), c=1.0)
for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
    
  for(i in 0:(nsim-1)){
    locs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "locs.txt")))
    obs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "obs.txt")))
    
    invisible( capture.output(parabolic <- smooth.FEM.time(locations=locs, lambdaS = lambda, lambdaT = 1,
                                time_mesh = as.vector(time_locs), PDE_parameters = PDE_parameters,
                              observations = obs, FEMbasis = FEMbasis,
                              FLAG_PARABOLIC = T)))
    coeff = matrix(parabolic$fit.FEM$coeff, nrow=n_dofs, ncol=n_time_locs)
    write.table( coeff, 
                 paste0(datadir, nlocs, "/", i, "/", "parabolic.txt"),
                 col.names = F, row.names = F)
  } 
}
