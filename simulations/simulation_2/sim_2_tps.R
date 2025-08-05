library(Matrix)
library(mgcv)

args = commandArgs(trailingOnly = T)

sigma_lab = c("0.00", "0.05", "0.10")

sd = sigma_lab[ as.integer(args[1]) + 1] # coerenti con c++, numerazione parte da 0
sim = args[2] 

data_dir = "../data/simulation_2/"
mesh_dir = "../data/simulation_1/mesh/"

sigma_dir = paste0(data_dir, "sigma_", sigma_lab[ as.integer(args[1]) + 1], "/")

sim_dir = paste0(sigma_dir, sim, "/")

locs = as.matrix(readMM(paste0(sim_dir, "locs.mtx")))
obs = as.matrix(readMM(paste0(sim_dir,"obs.mtx")))
time_locs = as.matrix(readMM(paste0(data_dir, "time_mesh.mtx")))
time_locs = time_locs[2:nrow(time_locs),1]
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


test_locs = as.matrix(readMM(paste0(mesh_dir, "../test_locs.mtx")))
test_obs = as.matrix(readMM(paste0(sim_dir,"test_obs.mtx")))

dat_test =  data.frame(p1 = rep(test_locs[,1],m),
                       p2 = rep(test_locs[,2],m),
                       time = rep(time_locs, each = nrow(test_locs)))

test_vals <- predict.gam(TPS,  newdata = dat_test) 

rmse = sqrt ( mean( (test_vals - as.vector(test_obs))^2 ) ) 
writeMM(Matrix(as.matrix(rmse), sparse = T), file = paste0(sim_dir, "rmse_tps.mtx"))

nodes = as.matrix(readMM(paste0(mesh_dir, "points.mtx")))

dat_coeff = data.frame(p1 = rep(nodes[,1],m),
                       p2 = rep(nodes[,2],m),
                       time = rep(time_locs, each = nrow(nodes)))
coeff_tps = predict.gam(TPS, newdata = dat_coeff)
writeMM(Matrix( matrix(coeff_tps, nrow=nrow(nodes), ncol=m), sparse = T), 
        file = paste0(sim_dir, "estimate_tps.mtx"))
