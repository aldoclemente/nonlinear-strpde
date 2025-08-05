library(Matrix)
library(fdaPDE2)
source("../graphics.R")
data_dir = "../data/simulation_2/"
mesh_dir = "../simulation_1/mesh/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))

coe <- function(x) {
  0.5 *  sin(5 * pi * x)  * exp(-x^2) + 1
}

mu <- function(x, y, t, tf) {
  a <- (9 / (5 * tf)) * t - 2
  
  abs(sin( 2 * pi * ( coe(y) * x * cos(a) - y * sin(a) ) ) *
    cos( 2 * pi * ( coe(y) * x * cos(a + pi/2) + coe(x) * y * sin( (pi/2) * a ) ) ))
}

t = 0.
coeff = mu(rep(nodes[,1], length(t)), rep(nodes[,2], length(t)), 
           rep(t, each=nrow(nodes)), tf = 1.)
writeMM(Matrix( matrix(coeff, nrow=nrow(nodes), ncol=length(t)), sparse = T), 
        file = "ic.mtx")

domain = triangulation(nodes, cells, boundary)
IC = fe_function(domain = domain, type = "P1", coeff = coeff)
plot.fem.base(list(IC))
