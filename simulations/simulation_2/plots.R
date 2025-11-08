
args = commandArgs(trailingOnly = T)
if(length(args) == 0) stop("Pass 'esatta' or 'stimata' as argument")

library(fdaPDE)
library(ggplot2)
library(viridis)
library(dplyr)

mesh_dir = "input/mesh/"
input_dir = "input/"
nodes = readMM(paste0(mesh_dir,"points.mtx"))
cells = readMM(paste0(mesh_dir, "cells.mtx"))
boundary = readMM(paste0(mesh_dir, "boundary.mtx"))
mesh = create.mesh.2D(nodes=nodes, triangles = cells + 1)
FEMbasis = create.FEM.basis(mesh)


rmse = data.frame(rmse=vector(mode="numeric"),
                  method=vector(mode="character"),
                  n_obs = vector(mode="character"))

locs = c("100", "250", "500","1000", "5000", "10000")

reactions = matrix(nrow = 30, ncol=length(locs))
#err_nonlin = reactions
#rmse_para_kFold = reactions
#rmse_para_gcv = reactions

for(j in 1:length(locs)){
  for(i in 1:30){
    if(args[1] == "esatta"){
      outdir = paste0("output/", locs[j], "/", i-1, "/")
      outdir_tps = paste0("output-tps/",locs[j], "/", i-1, "/")
    }else{
      outdir = paste0("output-ic-stimata/", locs[j], "/", i-1, "/")
      outdir_tps = paste0("output-tps-ic-stimata/",locs[j], "/", i-1, "/")
    }
    reactions[i,j] = readMM(paste0(outdir, "cv_optim.mtx"))[2]
    err_nonlin = readMM(paste0(outdir, "rmse_iterative.mtx"))[1]
    err_par = readMM(paste0(outdir, "rmse_diffusion.mtx"))[1]
    err_tps = readMM(paste0(outdir_tps, "rmse_tps.mtx"))[1]
    
    rmse = rbind(rmse, data.frame(rmse=err_nonlin, method="nonlinear", n_obs=locs[j]))
    rmse = rbind(rmse, data.frame(rmse= err_par, method= "parabolic", n_obs=locs[j]))
    rmse = rbind(rmse, data.frame(rmse= err_tps, method= "TPS", n_obs=locs[j]))
  }
}

head(rmse)

myblue = rgb(23,55,94, maxColorValue = 255)
myred = colorspace::lighten(myblue, amount = 0.35)
mygreen = colorspace::lighten(myblue, amount = 0.9)

for(j in 1:length(locs)){
  
  if(args[1] == "esatta"){
    outdir = paste0("output/", locs[j], "/")
  }else if(args[1] == "stimata"){
    outdir = paste0("output-ic-stimata/", locs[j], "/")
  }
  
  rmse2 = rmse %>% filter(n_obs == locs[j])
  
  pdf(paste0(outdir, "rmse-boxplots.pdf"))
  par(xpd=TRUE)
  boxplot(rmse2$rmse~ rmse2$method,
          col=c(myblue,myred, mygreen), xaxt="n", yaxt="n", ylab="", xlab="")
  axis(2, cex.axis=2.5)
  dev.off()
}

pdf(ifelse(args[1] == "esatta", "output/rmse-boxplots.pdf", "output-ic-stimata/rmse-boxplots.pdf"))
par(xpd=TRUE)
boxplot(rmse$rmse~ rmse$method + rmse$n_obs,
        col=c(myblue,myred,mygreen), xaxt="n", yaxt="n", ylab="", xlab="")
axis(1, at = seq(2, 3*length(locs), by = 3),  
     labels = locs, cex.axis=2.5, tick = FALSE)
axis(2, cex.axis=2.5)
dev.off()

pdf(ifelse(args[1] == "esatta", "output/rmse-legend.pdf", "output-ic-stimata/rmse-legend.pdf"))
plot.new()
legend(x="center",fill = c(myblue, myred, mygreen),
       horiz = T,
       legend = c("Nonlinear", "Parabolic", "TPS"),
       bty = "n", cex = 3)
dev.off()


n_bins = 10
grid <- expand.grid(
  x = seq(0, 1, length.out = 250),
  y = seq(0, 1, length.out = 250)
)
Psi = fdaPDE:::CPP_get.psi.Matrix(FEMbasis,as.matrix(grid))

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}

n_breaks = 20
n_bins = 10 
for(j in 1:length(locs)){
  if(args[1]=="esatta"){
    imgdir = paste0("output/", locs[j], "/")
  }else{
    imgdir = paste0("output-ic-stimata/", locs[j], "/") 
  }
  exact_coeff = as.matrix(readMM(paste0(input_dir, "fisher_kpp.mtx")))
  tps_coeff = matrix(0, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
  nonlin_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
  par_coeff = matrix(0, nrow = nrow(exact_coeff)*ncol(exact_coeff), ncol=1)
  
  for(i in 1:30){
    if(args[1] == "esatta"){
      outdir = paste0("output/", locs[j], "/", i-1, "/")
      outdir_tps = paste0("output-tps/",locs[j], "/", i-1, "/")
    }else{
      outdir = paste0("output-ic-stimata/", locs[j], "/", i-1, "/")
      outdir_tps = paste0("output-tps-ic-stimata/",locs[j], "/", i-1, "/")
    }
    
    tps_coeff = tps_coeff + as.matrix(readMM(paste0(outdir_tps, "estimate_tps.mtx"))) / 30
    nonlin_coeff = nonlin_coeff + as.matrix(readMM(paste0(outdir, "estimate_iterative.mtx"))) / 30
    par_coeff = par_coeff + as.matrix(readMM(paste0(outdir, "estimate_diffusion.mtx"))) / 30
  }
  
  nonlin_coeff = matrix(nonlin_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
  par_coeff = matrix(par_coeff, nrow=nrow(exact_coeff), ncol=ncol(exact_coeff))
  
  n_times = ncol( exact_coeff )
  estimates = list()
  estimates$exact = matrix(0, nrow=length(grid$x), ncol=n_times)
  estimates$nonlin = estimates$par = estimates$tps = estimates$exact
  
  estimates$exact = Psi %*% exact_coeff
  estimates$nonlin = Psi %*% nonlin_coeff
  estimates$par = Psi %*% par_coeff
  estimates$tps = Psi %*% tps_coeff

  data_list = list()
  data_list$locations = as.matrix(readMM(paste0(input_dir, "250/0/locs.mtx")))
  data_list$data = as.matrix(readMM(paste0(input_dir, "250/0/obs.mtx")))
  
  obs_range = range(data_list$data)
  
  # NON fittato su t0 = 0
  # lims = c( min(estimates$exact, estimates$nonlin, estimates$par, estimates$tps),#, obs_range), 
  #           max(estimates$exact, estimates$nonlin, estimates$par, estimates$tps)) #, obs_range))
  # 
  lims = c(0,1)
  mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 
  
  dofs = nodes
  n_time_locs = ncol(exact_coeff)
  { 
    pdf(paste0(imgdir,"nonlinear.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = grid$x, y = grid$y, z = estimates$nonlin[,t])
      plt <-
        ggplot(aes(x = x, y = y, z = z), data = data) +
        geom_raster(aes(fill = z)) +
        geom_contour(color = "black", bins = n_bins) +
        scale_fill_viridis_c(limits=lims) +
        coord_fixed() + theme_void() +
        theme(legend.position = "none")
      print(plt)
    }
    dev.off()
  }
  
  {
    pdf(paste0(imgdir,"exact.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = grid$x, y = grid$y, z = estimates$exact[,t])
      plt <-
        ggplot(aes(x = x, y = y, z = z), data = data) +
        geom_raster(aes(fill = z)) +
        geom_contour(color = "black", bins = n_bins) +
        scale_fill_viridis_c(limits=lims) +
        coord_fixed() + theme_void() +
        theme(legend.position = "none")
      print(plt)
    }
    dev.off()
  }
  
  {
    pdf(paste0(imgdir,"parabolic.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = grid$x, y = grid$y, z = estimates$par[,t])
      plt <-
        ggplot(aes(x = x, y = y, z = z), data = data) +
        geom_raster(aes(fill = z)) +
        geom_contour(color = "black", bins = n_bins) +
        scale_fill_viridis_c(limits=lims) +
        coord_fixed() + theme_void() +
        theme(legend.position = "none")
      print(plt)
    }
    dev.off()
  }
  
  {
    pdf(paste0(imgdir,"tps.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = grid$x, y = grid$y, z = estimates$tps[,t])
      plt <-
        ggplot(aes(x = x, y = y, z = z), data = data) +
        geom_raster(aes(fill = z)) +
        geom_contour(color = "black", bins = n_bins) +
        scale_fill_viridis_c(limits=lims) +
        coord_fixed() + theme_void() +
        theme(legend.position = "none")
      print(plt)
    }
    dev.off()
  }
  
  
  { 
    mycols <- mycolors(n_breaks + 2)
    vals   <- seq(lims[1], lims[2], length.out = length(mycols))
    
    linewidth=2
    pdf(paste0(imgdir,"pointwise_data.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = data_list$locations[,1], y = data_list$locations[,2], 
                         z = as.matrix(data_list$data[,t]))    
      plt <- ggplot(data, aes(x, y, fill = z)) + 
        geom_rect(xmin = 0, ymin = 0, xmax = 1, ymax = 1,
                  fill = "transparent",
                  color = "black",
                  linewidth = 1) +
        geom_point(shape = 21, size = 5, color = "transparent") + 
        scale_fill_viridis_c(limits = lims, oob = scales::squish) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "none")
      
      print(plt)
    }
    dev.off()
  }
}



# for(n_locs in locs){
# CV_ERRORS = matrix(0,nrow=85,ncol=1)
#     for(sim in 0:29){
#     outdir = paste0(n_locs, "/", sim, "/")
#     cv_errors = as.matrix(readMM(paste0(outdir, "cv_errors.mtx")))
#     cv_errors
#     cv_grids = as.matrix(readMM(paste0(outdir, "cv_grids.mtx")))
#     lambdas = unique(cv_grids[,1])
#     reactions = unique(cv_grids[,2])
#     
#     locs = as.matrix(readMM(paste0("../input/",n_locs,"/",sim, "/", "locs.mtx")))
#     obs = as.matrix(readMM(paste0("../input/",n_locs,"/",sim, "/", "obs.mtx")))
#     n_times = ncol(obs)
#     obs = obs[,2:n_times]
#     
#     lambda_para = lambdas[ which.min(cv_errors[1:length(lambdas),]) ]
#     sol_para = smooth.FEM.time(locations = locs,  time_mesh = ,observations=obs,FEMbasis = FEMbasis, lambda = lambda_para, 
#                           FLAG_PARABOLIC =TRUE)
#     
#     
#     CV_ERRORS = CV_ERRORS + cv_errors
#     
#   }
#   CV_ERRORS = CV_ERRORS / 30
#   
#   res = c()
#   for(i in 1:length(reactions)){
#     res = c(res, rep(i, each=length(lambdas)))
#   }
#   
#   df <- data.frame(x = rep(1:length(lambdas), times=length(reactions)),  # column index
#                    y = res,   
#                    cv_error = CV_ERRORS)
#   
#   # Plot as contiguous tiles
#   pdf(paste0(n_locs, "/", "mean_cv_errors.pdf"), width=15)
#   print(
#   ggplot(df, aes(x, y, fill = cv_error)) +
#     geom_tile() +
#     geom_text(aes(label = sprintf("%.4f", cv_error)), color = "black", size = 5) +
#     scale_fill_gradientn(colours = hcl.colors(100, "YlOrRd", rev = TRUE)) +
#     scale_x_continuous(breaks = 1:ncol(cv_matrix), 
#                        labels = as.character( seq(-5,3,by=0.5)) ) +
#     scale_y_continuous(breaks = 1:nrow(cv_matrix), 
#                        labels = as.character(reactions)) +
#     coord_fixed() +
#     labs(x = expression(log[10](lambda)), y = expression(r), fill = "") +
#     theme_minimal(base_size = 21))
#   dev.off()
# }
