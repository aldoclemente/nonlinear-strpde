library(Matrix)
library(mgcv)
library(fdaPDE)

args = commandArgs(trailingOnly = T)
LOCS = c("100", "250", "500", "1000", "5000", "10000")

data_dir = "input/"
mesh_dir = "input/mesh/"

nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)

sim_dir = paste0(data_dir, args[1], "/", args[2], "/")
locs = as.matrix(readMM(paste0(sim_dir, "locs.mtx")))

obs = as.matrix(readMM(paste0(sim_dir,"obs.mtx")))
time_locs = as.matrix(readMM(paste0(mesh_dir, "time_mesh.mtx")))

if(args[3] == "esatta"){
  out_dir = paste0("output-tps/", args[1], "/", args[2], "/")
  psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, locs)
  true = as.matrix(readMM(paste0(data_dir, "fisher_kpp.mtx")))[,1]
  obs[,1] = as.vector(psi%*%true)
} else if(args[3] == "stimata"){
  out_dir = paste0("output-tps-ic-stimata/", args[1], "/", args[2], "/")
}

if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
m = length(time_locs)
n = nrow(locs)

dat <- data.frame(
  y = as.vector(obs),
  p1 = rep(locs[,1], m),
  p2 = rep(locs[,2], m),
  time = rep(time_locs, each = n))

TPS <- gam( y ~ te(p1, p2, time, # k = c(k_tp_space, k_tp_time)
                   d = c(2, 1), bs = c("tp", "ps")),   # Note: with fix intercept 
            data = dat)


test_locs = as.matrix(readMM(paste0(mesh_dir, "test_locs.mtx")))
test_obs = as.matrix(readMM(paste0(mesh_dir,"test_obs.mtx")))

dat_test =  data.frame(p1 = rep(test_locs[,1],m),
                       p2 = rep(test_locs[,2],m),
                       time = rep(time_locs, each = nrow(test_locs)))

test_vals <- predict.gam(TPS,  newdata = dat_test) 

rmse = sqrt ( mean( (test_vals - as.vector(test_obs))^2 ) ) 
writeMM(Matrix(as.matrix(rmse), sparse = T), file = paste0(out_dir, "rmse_tps.mtx"))

dat_coeff = data.frame(p1 = rep(nodes[,1],m),
                       p2 = rep(nodes[,2],m),
                       time = rep(time_locs, each = nrow(nodes)))
coeff_tps = predict.gam(TPS, newdata = dat_coeff)
writeMM(Matrix( matrix(coeff_tps, nrow=nrow(nodes), ncol=m), sparse = T), 
        file = paste0(out_dir, "estimate_tps.mtx"))