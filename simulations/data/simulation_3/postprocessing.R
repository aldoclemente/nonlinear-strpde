library(Matrix)
library(fdaPDE2)
rm(list=ls())
source("../graphics.R")
nodes = readMM("../simulation_1/mesh/points.mtx")
cells = readMM("../simulation_1/mesh/cells.mtx")
boundary = readMM("../simulation_1/mesh/boundary.mtx")

domain = triangulation(nodes, cells, boundary)
plot(domain)


parabolic_coeff = as.matrix(readMM("sigma_0.00/0/estimate_parabolic.mtx"))
nonlinear_coeff = as.matrix(readMM("sigma_0.00/0/estimate_nonlinear.mtx"))
exact_coeff = as.matrix(readMM("../simulation_2/sigma_0.00/0/exact.mtx"))
exact_coeff = exact_coeff[,-1]

range(exact_coeff)
range(parabolic_coeff)
range(nonlinear_coeff)

n_times = ncol( exact_coeff )
parabolic_coeff = matrix(parabolic_coeff, nrow = nrow(nodes), ncol=n_times)
nonlinear_coeff = matrix(nonlinear_coeff, nrow = nrow(nodes), ncol=n_times)
#n_times = ncol( tps_coeff )

exact = list()
parabolic = list()
nonlinear = list()
for(t in 1:n_times){
  exact[[t]] = fe_function(domain, type="P1", exact_coeff[,t])
  parabolic[[t]] = fe_function(domain, type="P1", parabolic_coeff[,t])
  nonlinear[[t]] = fe_function(domain, type="P1", nonlinear_coeff[,t])
}

FE_list = c(exact, nonlinear, parabolic) 

plot.fem.base(FElist = FE_list, filename = "out.pdf")
