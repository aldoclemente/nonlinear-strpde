library(fdaPDE)
library(Matrix)
library(ggplot2)
library(viridis)
rm(list=ls())
source("../utils.R")
eval.locs = as.matrix(read.csv("data/eval_locs.csv")[,2:3])
datadir = "data/"

#NLOCS = c(100, 250, 500, 1000)
NLOCS = c(250)
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
PDE_parameters = list(K=matrix(c(0.1,0,0,0.1), nrow=2,ncol=2), b = c(0,0), c=1.0)
IC = as.matrix( as.vector(read.table(paste0(datadir, "exact.txt"))[,1]))

# exact 
exact = as.matrix( read.table(paste0(datadir, "exact.txt")))
test.evals = as.vector(eval.parabolic(exact, mesh, eval.locs, time_locs)[,2:n_time_locs])
#

results <- data.frame(rmse=vector(mode="numeric"), 
                      n_obs=vector(mode="integer"),
                      method=vector(mode="character"))
for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  errors <- data.frame(rmse=rep(0,times = nsim),
                       n_obs = rep(0L,times = nsim),
                       method = rep("",times = nsim))
  lambda = 1. / (nlocs * (n_time_locs))
  
  for(i in 0:(nsim-1)){
    locs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "locs.txt")))
    obs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "obs.txt")))
    
    invisible( capture.output(separable <- smooth.FEM.time(
      locations=locs, lambdaS = lambda, lambdaT = 1,
      time_mesh = as.vector(time_locs[2:n_time_locs]), 
      PDE_parameters = PDE_parameters,
      observations = obs[,2:n_time_locs], FEMbasis = FEMbasis,
      FLAG_PARABOLIC = F, FLAG_MASS = T)))
    coeff = eval.FEM.time(separable$fit.FEM.time, 
                          locations = dofs, time.instants = time_locs[2:n_time_locs])
    coeff = matrix(coeff, nrow=n_dofs, ncol=n_time_locs)
    write.table( coeff, 
                 paste0(datadir, nlocs, "/", i, "/", "separable.txt"),
                 col.names = F, row.names = F)
    evals = eval.FEM.time(separable$fit.FEM.time, 
                          locations = eval.locs, 
                          time.instants = time_locs[2:n_time_locs])
    errors$rmse[i+1] = rmse(evals, test.evals)
    errors$n_obs[i+1] = nlocs
    errors$method[i+1] = "STRPDE-SEP"
  }
  results = rbind(results, errors) 
}

write.table(results, paste0(datadir, "separable_rmse.txt") )
