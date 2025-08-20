setwd("simulations/data/simulation_1")

library(Matrix)
library(fdaPDE2)
source("../graphics.R")
nodes = readMM("mesh/points.mtx")
cells = readMM("mesh/cells.mtx")
boundary = readMM("mesh/boundary.mtx")

domain = triangulation(nodes, cells, boundary)
plot(domain)

tps_coeff = as.matrix(readMM("sigma_0.00/0/estimate_tps.mtx"))
exact_coeff = as.matrix(readMM("sigma_0.00/0/exact.mtx"))
nonlin_coeff = as.matrix(readMM("sigma_0.00/0/estimate_nonlinear.mtx"))
par_coeff = as.matrix(readMM("sigma_0.00/0/estimate_parabolic.mtx"))
exact_coeff = exact_coeff[,-1]

nonlin_coeff = matrix(nonlin_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
par_coeff = matrix(par_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))

n_times = ncol( tps_coeff )

tps = list()
exact = list()
nonlin = list()
parabolic = list()
for(t in 1:n_times){
  tps[[t]] = fe_function(domain, type="P1", tps_coeff[,t])
  exact[[t]] = fe_function(domain, type="P1", exact_coeff[,t])
  nonlin[[t]] = fe_function(domain, type="P1", nonlin_coeff[,t])
  #parabolic[[t]] = fe_function(domain, type="P1", par_coeff[,t])
}

FE_list = c(exact, nonlin, tps) # !!!

plot.fem.base(FElist = FE_list, filename = "out.pdf")


# ---

coeff = as.matrix(readMM("simulation_1/sigma_0.00/0/fs.test.mtx"))
fs_test = list()
for(t in 1:ncol(coeff)) fs_test[[t]] =fe_function(domain, type="P1", coeff=coeff[,t])

plot.fem.base(fs_test, filename = "vicini.pdf")

