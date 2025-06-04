library(fdaPDE)
library(Matrix)
library(ggplot2)
library(viridis)
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
nsim = 30 # 30
lambda = as.matrix(read.table(paste0(datadir, "lambda.txt"),
                    header = F))
IC = as.matrix( as.vector(read.table(paste0(datadir, "exact.txt"))[,1]))
PDE_parameters = list(K=matrix(c(0.1,0,0,0.1), nrow=2,ncol=2), b = c(0,0), c=3.0)

dirichlet_bc = which(boundary==4)
BC= list(
  BC_indices = rep(dirichlet_bc, each=(n_time_locs-1)), 
  BC_values = rep(0, length=(length(dirichlet_bc)*(n_time_locs-1)))
)

# exact 
exact = as.matrix( read.table(paste0(datadir, "exact.txt")))
test.evals = as.vector(eval.parabolic(exact, mesh, eval.locs, time_locs)[,2:n_time_locs])
#

results <- data.frame(rmse=vector(mode="numeric"), 
                      n_obs=vector(mode="integer"),
                      method=vector(mode="character"),
                      lambda = vector(mode="numeric"))
results_sep <- results
outdir = paste0(datadir, "output/")
if(!dir.exists(outdir)) dir.create(outdir)
# no BC neanche per il parabolic...
for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  errors <- data.frame(rmse=rep(0,times = nsim),
                       n_obs = rep(0L,times = nsim),
                       method = rep("",times = nsim),
                       lambda = rep(0, times=nsim))
  errors_sep <- errors
  #lambda = 1. / (nlocs * (n_time_locs))
  tmpdir = paste0(outdir, nlocs, "/")
  if(!dir.exists(tmpdir)) dir.create(tmpdir)
  for(l in 1:length(lambda)){
    lambdadir = paste0(tmpdir, l, "/")
    if(!dir.exists(lambdadir)) dir.create(lambdadir)
    for(i in 0:(nsim-1)){
      simdir = paste0(lambdadir, i, "/")
      if(!dir.exists(simdir)) dir.create(simdir)
      
      locs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "locs.txt")))
      obs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "obs.txt")))
      invisible( capture.output(parabolic <- 
                                  smooth.FEM.time(locations=locs, 
                                                  lambdaS = lambda[l]/((nlocs * (n_time_locs))),
                                                  lambdaT = 1,
                                                time_mesh = as.vector(time_locs), 
                                                PDE_parameters = PDE_parameters,
                                                observations = obs[,2:n_time_locs], 
                                                FEMbasis = FEMbasis, IC=IC,
                                                BC = BC,
                                                FLAG_PARABOLIC = T)))
    coeff = matrix(parabolic$fit.FEM$coeff, nrow=n_dofs, ncol=n_time_locs)
    write.table( coeff, 
                 paste0(simdir, "parabolic.txt"),
                 col.names = F, row.names = F)
    evals = as.vector(eval.parabolic(coeff, mesh, eval.locs, time_locs)[,2:n_time_locs])
    errors$rmse[i+1] = rmse(evals, test.evals)
    errors$n_obs[i+1] = nlocs
    errors$method[i+1] = "STRPDE-PAR"
    errors$lambda[i+1] = lambda[l]
    # ---
    invisible( capture.output(separable <- smooth.FEM.time(
      locations=locs, lambdaS = lambda[l]/((nlocs * (n_time_locs))), lambdaT = 1,
      time_mesh = as.vector(time_locs[2:n_time_locs]), 
      PDE_parameters = PDE_parameters,
      observations = obs[,2:n_time_locs], FEMbasis = FEMbasis,
      FLAG_PARABOLIC = F, FLAG_MASS = T)))
    coeff = eval.FEM.time(separable$fit.FEM.time, 
                          locations = dofs, time.instants = time_locs[2:n_time_locs])
    coeff = matrix(coeff, nrow=n_dofs, ncol=n_time_locs)
    write.table( coeff, 
                 paste0(simdir, "separable.txt"),
                 col.names = F, row.names = F)
    evals = eval.FEM.time(separable$fit.FEM.time, 
                          locations = eval.locs, 
                          time.instants = time_locs[2:n_time_locs])
    errors_sep$rmse[i+1] = rmse(evals, test.evals)
    errors_sep$n_obs[i+1] = nlocs
    errors_sep$method[i+1] = "STRPDE-SEP"
    errors_sep$lambda[i+1] = lambda[l]
  }
  results = rbind(results, errors)
  results_sep = rbind(results_sep, errors_sep)
  }
}

write.table(results, paste0(outdir, "parabolic_rmse.txt") )
write.table(results_sep, paste0(outdir, "separable_rmse.txt") )

# ---
library(dplyr)
rmse.par = results %>% group_by(lambda) %>% summarise( q1 = quantile(rmse, probs = 0.25),
                                                       q2 = median(rmse),
                                                       q3 = quantile(rmse, probs = 0.75))
rmse.sep = results_sep %>% group_by(lambda) %>% summarise( q1 = quantile(rmse, probs = 0.25),
                                                           q2 = median(rmse),
                                                           q3 = quantile(rmse, probs = 0.75))



myred = rgb( 139,0,0, max=255, alpha=125)
myblue = rgb( 65, 105, 225, max = 255, alpha = 125 )
mygreen = rgb( 0, 255, 0, max = 255, alpha = 125 )

xx = log(rmse.par$lambda/(n_time_locs*nlocs))
plot(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.par$q2,
     type="l", lwd=2, col="red4",
     xlab = expression(log(lambda)), ylab="")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.sep$q2, 
       type="l", lwd=2, col="green4")
polygon( x= c(rev( log(rmse.par$lambda/(n_time_locs*nlocs))),
              log(rmse.par$lambda/(n_time_locs*nlocs))),
         y= c( rev(rmse.par$q1), rmse.par$q3),
         col = myred, border = NA )
polygon( x= c(rev( log(rmse.par$lambda/(n_time_locs*nlocs))),
              log(rmse.par$lambda/(n_time_locs*nlocs))),
         y= c( rev(rmse.sep$q1), rmse.sep$q3),
         col = mygreen, border = NA )
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.par$q2,
     type="l", lwd=2, col="red4")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.sep$q2, 
       type="l", lwd=2, col="green4")

