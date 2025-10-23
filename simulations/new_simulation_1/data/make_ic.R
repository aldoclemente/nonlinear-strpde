library(Matrix)
library(fdaPDE)
library(ggplot2)
library(viridisLite)
source("../graphics.R")
data_dir = "data/simulation_2/"
mesh_dir = "mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
quad_nodes = readMM(paste0(mesh_dir, "quad_nodes.mtx"))
x11()
plot(quad_nodes[,1], quad_nodes[,2], pch=16)
coe <- function(x) {
  0.5 *  sin(5 * pi * x)  * exp(-x^2) + 1
}

mu <- function(x, y, t, tf) {
  a <- (9 / (5 * tf)) * t - 2
  
  abs(sin( 2 * pi * ( coe(y) * x * cos(a) - y * sin(a) ) ) *
        cos( 2 * pi * ( coe(y) * x * cos(a + pi/2) + coe(x) * y * sin( (pi/2) * a ) ) ))
}

times = as.matrix(readMM("mesh/time_mesh.mtx"))
n_times = length(times)
tf = 0.1
coeff = mu(rep(nodes[,1], n_times), rep(nodes[,2], n_times), rep(times, each=nrow(nodes)), tf=rep(tf,times=nrow(nodes)*length(time_mesh)))
#           rep(t, each=nrow(nodes)), tf = 1.)
coeff = matrix(coeff, nrow=nrow(nodes), ncol=n_times)

mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)
grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)

 
  vicini = matrix(0, nrow=length(grid$x), ncol=n_times)
  for(t in 1:n_times){
    vicini[,t] = eval.FEM(FEM(coeff[,t], FEMbasis), locations = grid)
  }

{  
  pdf("vicini.pdf")
  for(t in 1:n_times){
    data <- data.frame(x = grid$x, y = grid$y, z = vicini[,t])
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_raster(aes(fill = z)) +
      geom_contour(color = "black", bins = 10) +
      scale_fill_viridis_c(limits=c(0,1)) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

# ------------------------------------------------------------------------------
time_mesh = as.matrix(readMM("mesh/time_mesh.mtx"))
n_times = length(time_mesh)

forcing = matrix(nrow=nrow(quad_nodes), ncol=n_times)
tf = 0.1
for(t in 1:n_times){
  forcing[,t] = mu(quad_nodes[,1], quad_nodes[,2],time_mesh[t], tf)
}

writeMM(Matrix(forcing, sparse = T), "mesh/forcing_coeff.mtx")


IC = mvtnorm::dmvnorm(as.matrix(nodes), mean=c(0.5,0.5), sigma=0.025*diag(2))
IC = (IC - min(IC))/diff(range(IC))
range(IC)
{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = eval.FEM(FEM(IC, FEMbasis), grid) )
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = 10) +
    scale_fill_viridis_c(limits=c(0,1), oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}

writeMM(Matrix(IC, sparse = T), "mesh/IC.mtx")

grid <- expand.grid(
  x = seq(0.05, 0.95, length.out = 60),
  y = seq(0.05, 0.95, length.out = 60)
)

writeMM(Matrix(as.matrix(grid), sparse = T), "mesh/test_locs.mtx")

# estimate IC using fdaPDE
lambda
plot(mesh)

# writeMM(Matrix( matrix(coeff, nrow=nrow(nodes), ncol=length(t)), sparse = T), 
#         file = "ic.mtx")

#domain = triangulation(nodes, cells, boundary)
#IC = fe_function(domain = domain, type = "P1", coeff = coeff)
#plot.fem.base(list(IC))

mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)
grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)

{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = eval.FEM(FEM(coeff, FEMbasis), grid) )
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = 10) +
    scale_fill_viridis_c(limits=c(0,1), oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}


# ------------------------------------------------------------------------------

library(Matrix)
library(fdaPDE)
source("../graphics.R")
data_dir = "../data/simulation_2/"
mesh_dir = "../simulation_1/mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
plot(mesh)

std_devs = c("0.00", "0.05", "0.10")
lambda = 10^seq(from=-6,2,length.out=500)
FEMbasis = create.FEM.basis(mesh)

K = matrix(c(0.001, 0, 0, 0.001), nrow=2, ncol=2)
b = matrix(0, nrow=2, ncol=1)
c = 0.#-3.0
PDE_parameters = list(K = K, b = b, c = 0)
PDE_parameters2 = list(K = K, b = b, c = 1.)
set.seed(0)
for(sd in std_devs){
  for(i in 0:29){
    datadir = paste0("sigma_", sd, "/", i, "/")
    locs = as.matrix(readMM(paste0(datadir, "locs.mtx")))
    exact = readMM(paste0(datadir, "exact.mtx"))[,1]
    psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, locs)
    obs =  as.matrix(psi %*% exact + rnorm(nrow(locs), mean = 0, sd = 0.05))
    writeMM(Matrix(obs,sparse = T), file = paste0(datadir, "obs.t0.mtx"))
    sol = smooth.FEM(locations = locs, observations = obs, FEMbasis = FEMbasis, 
                     lambda.selection.criterion = "grid", 
                     lambda.selection.lossfunction = "GCV", PDE_parameters = PDE_parameters,
                     lambda = lambda, DOF.evaluation = "stochastic")
    sol2 = smooth.FEM(locations = locs, observations = obs, FEMbasis = FEMbasis, 
                     lambda.selection.criterion = "grid", 
                     lambda.selection.lossfunction = "GCV", PDE_parameters = PDE_parameters2,
                     lambda = lambda, DOF.evaluation = "stochastic")
    
    sol3 = smooth.FEM(locations = locs, observations = obs, FEMbasis = FEMbasis, 
                      lambda.selection.criterion = "grid", 
                      lambda.selection.lossfunction = "GCV",
                      lambda = lambda, DOF.evaluation = "stochastic")
    writeMM(obj = Matrix(sol$fit.FEM$coeff,sparse = T), file = paste0(datadir, "ic.mtx"))
    plot(sol$optimization$GCV_vector)
  }
}

x11()
plot(sol2$optimization$GCV_vector, type = "l", col="red")
points(sol3$optimization$GCV_vector, type="l", col="green3")
points(sol$optimization$GCV_vector, type="l")
estimates =list()
estimates$exactt0 = eval.FEM(FEM(exact, FEMbasis), locations = grid)
estimates$IC = eval.FEM(sol$fit.FEM, locations = grid)
estimates$IC2 = eval.FEM(sol2$fit.FEM, locations = grid)
estimates$IC3 = eval.FEM(sol3$fit.FEM, locations = grid)

lims = range(estimates$exactt0, estimates$IC, estimates$IC2, estimates$IC3 )
{  
  x11()
  data <- data.frame(x = grid$x, y = grid$y, z = estimates$exactt0)
  plt <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = n_bins) +
    scale_fill_viridis_c(limits=lims, oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt
}

{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = estimates$IC)
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = n_bins) +
    scale_fill_viridis_c(limits=lims, oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}

{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = estimates$IC)
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = n_bins) +
    scale_fill_viridis_c(limits=lims, oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}


{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = estimates$IC2)
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = n_bins) +
    scale_fill_viridis_c(limits=lims, oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}

{ x11()
  data <- data.frame(x = grid$x, y = grid$y, z = estimates$IC3)
  plt2 <-
    ggplot(aes(x = x, y = y, z = z), data = data) +
    geom_raster(aes(fill = z)) +
    geom_contour(color = "black", bins = n_bins) +
    scale_fill_viridis_c(limits=lims, oob = scales::squish) +
    coord_fixed() + theme_void() +
    theme(legend.position = "none")
  plt2
}

sqrt( mean( (estimates$exactt0 - estimates$IC)^2 ) )
sqrt( mean( (estimates$exactt0 - estimates$IC2)^2 ) )
sqrt( mean( (estimates$exactt0 - estimates$IC3)^2 ) )
