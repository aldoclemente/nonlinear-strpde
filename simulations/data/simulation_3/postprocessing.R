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
library(fdaPDE2)
rm(list=ls())
source("../graphics.R")
nodes = readMM("../simulation_1/mesh/points.mtx")
cells = readMM("../simulation_1/mesh/cells.mtx")
boundary = readMM("../simulation_1/mesh/boundary.mtx")

domain = triangulation(nodes, cells, boundary)

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

exact = list()
nonlin = list()
parabolic = list()
for(t in 1:n_times){
  exact[[t]] = fe_function(domain, type="P1", exact_coeff[,t])
  nonlin[[t]] = fe_function(domain, type="P1", nonlin_coeff[,t])
  parabolic[[t]] = fe_function(domain, type="P1", par_coeff[,t])
}

FElist = c(exact, nonlin, parabolic) 

#
incidence_matrix = read.csv("incidence_matrix.csv")[,-1]
areal_list = list()
areal_list$data = list()
areal_list$incidence_matrix = list()

for(t in 1:n_times){
  areal_list$incidence_matrix[[t]] = incidence_matrix
  areal_list$data[[t]] = as.matrix(readMM("sigma_0.00/0/obs.mtx"))[,t]
}

plot.fem.base <- function(FElist, data_list = NULL, raster_list =NULL,
                          filename=NULL, filename_data=NULL, filename_raster=NULL, 
                          areal_list = NULL, filename_areal=NULL,
                          palette=viridis){
  
  lims = c(1e10, -1e10)
  for(i in 1:length(FElist)){
    lims = c(min(c(FElist[[i]]$coeff, lims[1]), na.rm = T),
             max(c(FElist[[i]]$coeff, lims[2]), na.rm = T))
  }
  
  if(!is.null(data_list)){
    for(i in 1:length(data_list$data)){
      lims = c(min(c(as.vector(data_list$data[[i]]), lims[1]), na.rm = T),
               max(c(as.vector(data_list$data[[i]]), lims[2]), na.rm = T))
    }
  }
  
  if(!is.null(areal_list)){
    for(i in 1:length(areal_list$data)){
    lims = c(min(c(as.vector(areal_list$data[[i]]), lims[1]), na.rm = TRUE),
             max(c(as.vector(areal_list$data[[i]]), lims[2]), na.rm = TRUE))
    }
  }
  
  palette_ <- palette(100)
  if(!is.null(filename)) pdf(file=filename)
  for(i in 1:length(FElist)){
    plot(FElist[[i]], palette = palette, lims = lims,
         frame.plot = FALSE, xaxt = "n", yaxt = "n")
    fields::image.plot(legend.only = TRUE, zlim = range(lims, na.rm = TRUE),
                       col = palette_, legend.args = list(side = 4),
                       legend.lab = "", legend.mar = 4)
  }
  if(!is.null(filename)) dev.off()
  
  if(!is.null(data_list)){
    n_col = 100
    breaks <- seq(lims[1], lims[2], length.out = n_col + 1)
    mesh = FElist[[1]]$geometry
    nodes = mesh$nodes
    edges = mesh$edges
    bd_edges = edges[mesh$boundary_edges == 1,]
    if(!is.null(filename_data)) pdf(file=filename_data)
    for(i in 1:length(data_list$data)){
      locations = data_list$locations[[i]]
      vals = data_list$data[[i]]
      
      col_idx <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      col <- palette_[col_idx]
      plot(locations[, 1], locations[, 2], xlab = "", ylab = "", pch = 16,
           col = col, asp = 1, cex = 1.5, frame.plot=F, xaxt="n", yaxt="n")
      segments(nodes[bd_edges[,1],1], nodes[bd_edges[,1],2],
               nodes[bd_edges[,2],1], nodes[bd_edges[,2],2], asp=1)
    }
    plot.new()
    fields::image.plot(legend.only = TRUE, zlim = range(lims, na.rm = TRUE),
                       col = palette_, legend.args = list(side = 4),
                       legend.lab = "", legend.mar = 4)
    if(!is.null(filename_data)) dev.off()
  }
  
  if(!is.null(areal_list)){
    n_col = 100
    breaks <- seq(lims[1], lims[2], length.out = n_col + 1)
    mesh = FElist[[1]]$geometry
    nodes = mesh$nodes
    edges = mesh$edges
    bd_edges = edges[mesh$boundary_edges == 1,]
    if(!is.null(filename_areal)) pdf(file=filename_areal)
    for(i in 1:length(areal_list$data)){
      incidence_matrix = areal_list$incidence_matrix[[i]]
      vals = areal_list$data[[i]]
      
      col_idx <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      col <- palette_[col_idx]
      
      plot.new()
      for( r in 1:nrow(incidence_matrix)){
        apply(mesh$cells[which(incidence_matrix[r,] == 1),], 
              MARGIN=1, FUN= function(x){
                polygon(mesh$nodes[x,],col=col[r], , border=NA)})
        }
      segments(nodes[bd_edges[,1],1], nodes[bd_edges[,1],2],
               nodes[bd_edges[,2],1], nodes[bd_edges[,2],2], asp=1)
    }
    plot.new()
    fields::image.plot(legend.only = TRUE, zlim = range(lims, na.rm = TRUE),
                       col = palette_, legend.args = list(side = 4),
                       legend.lab = "", legend.mar = 4)
    if(!is.null(filename_areal)) dev.off()
  }
  
}

if(!dir.exists("imgs/")) dir.create("imgs/")
plot.fem.base(FElist = FE_list, filename = "imgs/estimates.pdf",
              areal_list = areal_list, filename_areal = "imgs/data.pdf")




