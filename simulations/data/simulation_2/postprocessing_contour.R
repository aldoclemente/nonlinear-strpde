# boxplots ---------------------------------------------------------------------
setwd("/home/user/simulations/data/simulation_2")

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

pdf("imgs/rmse-legend.pdf", width = 14, height = 3)
plot.new()
legend(x="center",fill = c(myblue, myred, mygreen),
       horiz = T,
       legend = c("Nonlinear", "Parabolic", "TPS"),
       bty = "n", cex = 3)
dev.off()

# estimates --------------------------------------------------------------------
library(Matrix)
library(fdaPDE)
source("../graphics.R")
nodes = readMM("../simulation_1/mesh/points.mtx")
cells = readMM("../simulation_1/mesh/cells.mtx")
boundary = readMM("../simulation_1/mesh/boundary.mtx")

domain = create.mesh.2D(nodes, triangles=cells+1)
FEMbasis = create.FEM.basis(domain)
plot(domain)
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

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}

n_breaks = 20

data_list = list()
data_list$locations = as.matrix(readMM("sigma_0.00/0/locs.mtx"))
data_list$data = as.matrix(readMM("sigma_0.00/0/obs.mtx"))

# 0.2, 0.4, 0.6, 0.8, 1.0
frame = c(2,4,6,8,10)
obs_range = range(data_list$data[,frame])

# NON fittato su t0 = 0
lims = c( min(exact_coeff[,frame], nonlin_coeff[,frame], par_coeff[,frame], tps_coeff[,frame], obs_range), 
          max(exact_coeff[,frame], nonlin_coeff[,frame], par_coeff[,frame], tps_coeff[,frame], obs_range) )

mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 

# {
#   linewidth=2
#   locs =  as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/locs.txt")))
#   obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt")))
#   pdf(paste0(datadir, nlocs,"/data.pdf"))
#   for(t in 1:n_time_locs){
#     data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
#     plt <- ggplot() +
#       geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
#       geom_segment(data=data.frame(x=0, xend=1, y=0,yend=0),
#                    aes(x=x, xend=xend, y=y, yend=yend), color="black", 
#                    linewidth=linewidth) +
#       geom_segment(data=data.frame(x=0, xend=0, y=0,yend=1),
#                    aes(x=x, xend=xend, y=y, yend=yend), color="red",
#                    linewidth=linewidth) +
#       geom_segment(data=data.frame(x=0, xend=1, y=1,yend=1),
#                    aes(x=x, xend=xend, y=y, yend=yend), color="black",
#                    linewidth=linewidth) +
#       geom_segment(data=data.frame(x=1, xend=1, y=0,yend=1),
#                    aes(x=x, xend=xend, y=y, yend=yend), color="black",
#                    linewidth=linewidth) +
#       scale_color_gradientn(colors =mycolors(n_breaks + 2), limits=lims) +
#       coord_fixed() + theme_void() +
#       theme(legend.position = "none")
#     print(plt)
#   }
#   dev.off()
# }
# 
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
# scale_fill_manual(
#   aesthetics = "fill",
#   values = mycolors(n_breaks + 2), limits=lims
# )
# 

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

# scale_fill_manual(
#   aesthetics = "fill",
#   values = mycolors(n_breaks + 2), limits=lims
# )
# 
# 
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

{
  pdf(paste0(paste0(outdir,"tps.pdf")))
  for(t in frame){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(tps_coeff[,t]))
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
  linewidth=2
  locs = as.matrix(readMM("sigma_0.00/0/locs.mtx"))
  obs = as.matrix(readMM("sigma_0.00/0/obs.mtx"))
  pdf(paste0(outdir,"pointwise_data.pdf"))
  for(t in frame){
    data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
    plt <- ggplot() +
      geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
      geom_rect(data = data.frame(xmin = 0,ymin = 0,xmax = +1,ymax = +1), 
                aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                fill = "transparent", 
                color = "black",
                linewidth = 1) +
      scale_color_gradientn(colors =mycolors(n_breaks+2), limits=lims, 
                            values = seq(from=lims[1], to=lims[2], length.out=n_breaks)) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}
