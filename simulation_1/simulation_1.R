library(Matrix)
library(ggplot2)
library(viridis)
rm(list=ls())
NLOCS = c(100, 250, 500, 1000) #c(1000) # c(250)

datadir = "data/"

set.seed(314156)
N = 1000
eval.locs = matrix(0,nrow=N,ncol=2)
eval.locs[,1] = runif(1000)
eval.locs[,2] = runif(1000)
write.csv(eval.locs, paste0(datadir, "eval_locs.csv"))

time_locs = as.matrix(read.table(paste0(datadir, "time_locations.txt")))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("data/mesh/elements.mtx"))

source("simulation_1_parabolic.R")
source("simulation_1_separable.R")
source("simulation_1_soap.R")

library(fdaPDE)
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 

nsim = 30

f_exact =  as.matrix(read.table(paste0(datadir, "exact.txt")))

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}


for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  f <- matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f <- f + matrix(as.matrix(readMM(paste0(datadir, nlocs, "/", i, "/", "y.mtx"))), nrow=n_dofs, ncol= n_time_locs)/nsim
  }
  
  n_breaks <- 50
  # qua andrebbe considerata anche la "exact" per fare grafici coerenti (?)
  f_separable = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f_separable <- f_separable + as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "separable.txt"), header = F))/nsim
  }
  
  f_parabolic = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f_parabolic <- f_parabolic + as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "parabolic.txt"), header = F))/nsim
  }
  #obs_range = c(1e10,-1e10)
  obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt"))) # uso solo la prima...
  obs_range = range(obs)
  
  # NON fittato su t0 = 0
  f_tps = matrix(0,nrow=n_dofs, ncol=n_time_locs-1)
  for(i in 0:(nsim-1)){
    f_tps <- f_tps + as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "tps.txt"), header = F))/nsim
  }
  
  lims = c( min(f, f_exact, f_separable,obs_range, f_tps), 
            max(f, f_exact, f_separable,obs_range, f_tps) )
  mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 
  
  {
  linewidth=2
  locs =  as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/locs.txt")))
  obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt")))
  pdf(paste0(datadir, nlocs,"/data.pdf"))
  for(t in 1:n_time_locs){
    data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
    plt <- ggplot() +
      geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
       geom_segment(data=data.frame(x=0, xend=1, y=0,yend=0),
                    aes(x=x, xend=xend, y=y, yend=yend), color="black", 
                    linewidth=linewidth) +
      geom_segment(data=data.frame(x=0, xend=0, y=0,yend=1),
                   aes(x=x, xend=xend, y=y, yend=yend), color="red",
                   linewidth=linewidth) +
       geom_segment(data=data.frame(x=0, xend=1, y=1,yend=1),
                    aes(x=x, xend=xend, y=y, yend=yend), color="black",
                    linewidth=linewidth) +
       geom_segment(data=data.frame(x=1, xend=1, y=0,yend=1),
                    aes(x=x, xend=xend, y=y, yend=yend), color="black",
                    linewidth=linewidth) +
      scale_color_gradientn(colors =mycolors(n_breaks + 2), limits=lims) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
  }

{
  pdf(paste0(paste0(datadir, nlocs, "/","y.pdf")))
  for(t in 1:n_time_locs){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f[,t]))    
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
  pdf(paste0(datadir, nlocs, "/", "y_exact.pdf"))
  for(t in 1:n_time_locs){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_exact[,t]))    
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
  pdf(paste0(paste0(datadir, nlocs, "/","y_separable.pdf")))
  for(t in 1:n_time_locs){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_separable[,t]))
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
    pdf(paste0(paste0(datadir, nlocs, "/","y_parabolic.pdf")))
    for(t in 1:n_time_locs){
      data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_parabolic[,t]))
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
    pdf(paste0(paste0(datadir, nlocs, "/","y_tps.pdf")))
    for(t in 1:(n_time_locs-1)){
      data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_tps[,t]))    
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
  
  
  
}

## boxplots

source("../utils.R")
eval.locs = as.matrix(read.csv("data/eval_locs.csv")[,2:3])

exact = as.matrix( read.table(paste0(datadir, "exact.txt")))
test.evals = as.vector(eval.parabolic(exact, mesh, eval.locs, time_locs)[,2:n_time_locs])

results <- data.frame(rmse=vector(mode="numeric"), 
              n_obs=vector(mode="integer"),
              method=vector(mode="character"))

for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  errors <- data.frame(rmse=rep(0,times = nsim),
                       n_obs = rep(0L,times = nsim),
                       method = rep("",times = nsim))
  for(i in 0:(nsim-1)){
  f <- matrix(as.matrix(readMM(paste0(datadir, nlocs, "/", i, "/", "y.mtx"))), 
              nrow=n_dofs, ncol= n_time_locs)
  
 
  evals = as.vector(eval.parabolic(f, mesh, eval.locs, time_locs)[,2:n_time_locs])
  errors$rmse[i+1] = rmse(evals, test.evals)
  errors$n_obs[i+1] = nlocs
  errors$method[i+1] = "STRPDE-NL"
  }
results = rbind(results, errors) 
}

write.table(results, paste0(datadir, "nonlinear_rmse.txt") )
#

nonlinear <- read.table(paste0(datadir, "nonlinear_rmse.txt"))[c("rmse", "n_obs", "method")]
parabolic <- read.table(paste0(datadir, "parabolic_rmse.txt"))[c("rmse", "n_obs", "method")]
separable <- read.table(paste0(datadir, "separable_rmse.txt"))[c("rmse", "n_obs", "method")]
tps <- read.table(paste0(datadir, "tps_rmse.txt"))[c("rmse", "n_obs", "method")]

rmse <- rbind(nonlinear, parabolic, separable, tps)
rmse$n_obs <- as.factor(rmse$n_obs)
rmse$method <- factor(rmse$method, levels = c("STRPDE-NL", "STRPDE-PAR", 
                                              "STRPDE-SEP", "TPS"))

table(rmse$method)

plt <- boxplot.rmse(rmse) + MyTheme

{
  tmp <- plt + 
    labs(title="RMSE", x="observations") +
    theme(legend.position = "bottom")  
  
  leg <- cowplot::get_plot_component(tmp, 'guide-box', return_all = TRUE)
  
  pdf(paste0(datadir, "RMSE-legend.pdf"), width = 16, height = 1)
  grid.draw(leg[[3]])
  dev.off()
}
{
  tmp <- plt + 
    labs(title="RMSE", x="observations") +
    theme(legend.position = "none")
  ggsave(paste0(datadir, "RMSE-no-legend.pdf"), plot=tmp, width = 8, height = 7)
}
