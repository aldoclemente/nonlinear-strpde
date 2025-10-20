# boxplots ---------------------------------------------------------------------
setwd("/home/user/simulations/data/simulation_3")

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
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="nonlinear", sd=sd))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "parabolic", sd=sd))
  }
}


myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

if(!dir.exists("imgs/")) dir.create("imgs/")

pdf(paste0("imgs/rmse-boxplots.pdf"))
par(xpd=TRUE)
boxplot(rmse$rmse~ rmse$method + rmse$sd,
        col=c(myblue,myred), xaxt="n", yaxt="n", ylab="", xlab="")
axis(1, at = seq(2, 8, by = 3), 
     labels = sds, cex.axis=2.5)
axis(2, cex.axis=2.5)
mtext(expression(sigma), side=1, line=3,cex=2.5)
dev.off()

# estimates --------------------------------------------------------------------
setwd("/home/user/simulations/data/simulation_3")
library(Matrix)
rm(list=ls())
source("../graphics.R")
nodes = readMM("../simulation_1/mesh/points.mtx")
cells = readMM("../simulation_1/mesh/cells.mtx")
boundary = readMM("../simulation_1/mesh/boundary.mtx")

domain = create.mesh.2D(nodes, triangles=cells+1)
FEMbasis = create.FEM.basis(domain)
plot(domain)

sd = "0.00"
nsim=30
exact_coeff = as.matrix(readMM("../simulation_2/sigma_0.00/0/exact.mtx"))
exact_coeff = exact_coeff[,-1]
nonlin_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
par_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
for(j in 0:(nsim-1)){
  nonlin_coeff = nonlin_coeff + as.matrix(readMM(paste0("sigma_", sd, "/", j,"/", "estimate_nonlinear.mtx"))) / nsim
  par_coeff = par_coeff + as.matrix(readMM(paste0("sigma_", sd, "/", j,"/", "estimate_parabolic.mtx"))) / nsim
}

nonlin_coeff = matrix(nonlin_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
par_coeff = matrix(par_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))

n_times = ncol( exact_coeff )

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}

n_breaks = 20

data_list = list()
data_list$data = as.matrix(readMM("sigma_0.00/0/obs.mtx"))

# 0.2, 0.4, 0.6, 0.8, 1.0
frame = c(2,4,6,8,10)
obs_range = range(data_list$data[,frame])

lims = c( min(exact_coeff[,frame], nonlin_coeff[,frame], par_coeff[,frame], obs_range), 
          max(exact_coeff[,frame], nonlin_coeff[,frame], par_coeff[,frame], obs_range) )

mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 

dofs = nodes
outdir = "imgs/"
n_time_locs = ncol(exact_coeff)
{
  pdf(paste0(paste0(outdir,"nonlinear.pdf")))
  for(t in frame){
    #coef = eval.FEM(FEM(as.matrix(f[,t]), FEMbasis), locations=dofs.ref)
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(nonlin_coeff[,t]))
    #data <- data.frame(x = round(dofs.ref[,1], 10), y = round(dofs.ref[,2], 10), z = as.matrix(coef))
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_contour_filled(breaks = mybreaks) +
      scale_color_viridis(limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

{
  pdf(paste0(outdir, "exact.pdf"))
  for(t in frame){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(exact_coeff[,t]))    
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_contour_filled(breaks = mybreaks) +
      scale_color_viridis(limits=lims) + 
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

{
  pdf(paste0(paste0(outdir,"parabolic.pdf")))
  for(t in frame){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(par_coeff[,t]))
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_contour_filled(breaks = mybreaks) +
      scale_color_viridis(limits=lims) + 
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

incidence_matrix = as.matrix(read.csv("incidence_matrix.csv")[,-1])
poly_coords = list()
for(r in 1:nrow(incidence_matrix)){
  triangles = domain$triangles[which(incidence_matrix[r,] == 1),]
  #triangles = cbind(triangles, triangles[,1]) 
  poly_coords[[r]] = list()
  for(t in 1:nrow(triangles)){
    poly_coords[[r]][[t]] = domain$nodes[triangles[t,], ]
    poly_coords[[r]][[t]] = rbind(poly_coords[[r]][[t]], 
                                  domain$nodes[triangles[t,][1], ])
  }
} 


{
  pdf(paste0(outdir,"areal_data.pdf"))
  for(t in frame){
    plt <- ggplot()
    for(r in 1:nrow(incidence_matrix)){
      for(i in 1:length(poly_coords[[r]])){
        plt <- plt +
          geom_polygon(data=data.frame(x=poly_coords[[r]][[i]][,1], 
                                       y=poly_coords[[r]][[i]][,2],
                                       z=rep(data_list$data[r,t], times=4)), 
                       aes(x = x, y = y, fill = z), color=NA) 
      }
    }
    plt <- plt  +
      #scale_fill_gradientn(colors =mycolors(n_breaks+2), limits=lims) +
      scale_fill_viridis(option = "viridis", limits=lims) + 
      #scale_color_gradientn(colors =mycolors(n_breaks+2), limits=lims) +
      coord_fixed() + theme_void() +  
      theme(legend.position = "none")
    
    plt_rect <- plt + geom_rect(data = data.frame(xmin =0, ymin = 0, xmax = 1,ymax =1), 
                                aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                                fill = "transparent", 
                                color = "black",
                                linewidth = 1)
    print(plt_rect)
  }
  dev.off()
}

plt_rect
