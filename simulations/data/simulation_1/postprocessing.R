# boxplots ---------------------------------------------------------------------
setwd("/home/user/simulations/data/simulation_1")

library(Matrix)
rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"),
                  sd = vector(mode="character"))

nsim = 30
sds = c("0.00", "0.05", "0.10")

for(i in 1:length(sds)){
  sd = sds[i]
  sigmadir = paste0("sigma_",sd,"/")
  
  for(j in 0:(nsim-1)){
    resdir = paste0(sigmadir, j, "/")
    err_nonlin = readMM(paste0(resdir, "rmse_nonlinear.mtx"))[1,1]
    err_par = readMM(paste0(resdir, "rmse_parabolic.mtx"))[1,1]
    err_tps = readMM(paste0(resdir, "rmse_tps.mtx"))[1,1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="nonlinear", sd=sd))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "parabolic", sd=sd))
    rmse = rbind(rmse, data.frame(rmse= err_tps, method= "TPS", sd=sd))
  }
}


myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

if(!dir.exists("imgs/")) dir.create("imgs/")

pdf(paste0("imgs/rmse-boxplots.pdf"))
par(xpd=TRUE)
boxplot(rmse$rmse~ rmse$method + rmse$sd,
        col=c(myblue,myred,mygreen), xaxt="n", yaxt="n", ylab="", xlab="")
axis(1, at = seq(2, 8, by = 3), 
     labels = sds, cex.axis=2.5)
axis(2, cex.axis=2.5)
mtext(expression(sigma), side=1, line=3,cex=2.5)
dev.off()

pdf("imgs/rmse-legend.pdf", width = 14, height = 2.5)
plot.new()
legend(x="center",fill = c(myblue, myred, mygreen),
       horiz = T,
       legend = c("Nonlinear", "Parabolic", "TPS"),
       bty = "n", cex = 3)
dev.off()

# estimates --------------------------------------------------------------------
setwd("/home/user/simulations/data/simulation_1")
library(Matrix)
library(fdaPDE2)
source("../graphics.R")
nodes = readMM("mesh/points.mtx")
cells = readMM("mesh/cells.mtx")
boundary = readMM("mesh/boundary.mtx")

domain = triangulation(nodes, cells, boundary)
#plot(domain)

sd = "0.00"
nsim=30
exact_coeff = as.matrix(readMM("sigma_0.00/0/exact.mtx"))
exact_coeff = exact_coeff[,-1]
tps_coeff = matrix(0, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
nonlin_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
par_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
for(j in 0:(nsim-1)){
  tps_coeff = tps_coeff + as.matrix(readMM(paste0("sigma_", sd, "/", j,"/", "estimate_tps.mtx"))) / nsim
  nonlin_coeff = nonlin_coeff + as.matrix(readMM(paste0("sigma_", sd, "/", j,"/", "estimate_nonlinear.mtx"))) / nsim
  par_coeff = par_coeff + as.matrix(readMM(paste0("sigma_", sd, "/", j,"/", "estimate_parabolic.mtx"))) / nsim
}

nonlin_coeff = matrix(nonlin_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
par_coeff = matrix(par_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))

n_times = ncol( exact_coeff )

tps = list()
exact = list()
nonlin = list()
parabolic = list()
for(t in 1:n_times){
  exact[[t]] = fe_function(domain, type="P1", exact_coeff[,t])
  nonlin[[t]] = fe_function(domain, type="P1", nonlin_coeff[,t])
  parabolic[[t]] = fe_function(domain, type="P1", par_coeff[,t])
  tps[[t]] = fe_function(domain, type="P1", tps_coeff[,t])
}

FE_list = c(exact, nonlin, parabolic, tps) 

data_list = list()
for(t in 1:n_times){
  data_list$locations[[t]] = as.matrix(readMM("sigma_0.00/0/locs.mtx"))
  data_list$data[[t]] = as.matrix(readMM("sigma_0.00/0/obs.mtx"))[,t]
}

if(!dir.exists("imgs/")) dir.create("imgs/")
plot.fem.base(FElist = FE_list, filename = "imgs/estimates.pdf",
              data_list = data_list, filename_data = "imgs/data.pdf")

