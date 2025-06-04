library(Matrix)
library(ggplot2)
library(viridis)
rm(list=ls())
#NLOCS = c(100, 250, 500, 1000) #c(1000) # c(250)

NLOCS = c(250)

datadir = "data/"

set.seed(314156)
N = 1000
eval.locs = matrix(0,nrow=N,ncol=2)
eval.locs[,1] = -2.5 + 5*runif(1000)
eval.locs[,2] = -2.5 + 5*runif(1000)
write.csv(eval.locs, paste0(datadir, "eval_locs.csv"))
lambda = 10^seq(0,5, by=0.25) # -2,0 nope, -2,3
write.table(lambda, paste0(datadir, "lambda.txt"), 
            row.names = F, col.names = F)
time_locs = as.matrix(read.table(paste0(datadir, "time_locations.txt")))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("data/mesh/elements.mtx"))

#source("parabolic.R")
#source("separable.R")
#source("soap.R")

library(fdaPDE)
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 

nsim = 30

f_exact =  as.matrix(read.table(paste0(datadir, "exact.txt")))

## nonlinear rmse ---

source("../utils.R")
eval.locs = as.matrix(read.csv("data/eval_locs.csv")[,2:3])

exact = as.matrix( read.table(paste0(datadir, "exact.txt")))
test.evals = as.vector(eval.parabolic(exact, mesh, eval.locs, time_locs[2:n_time_locs,]))

results <- data.frame(rmse=vector(mode="numeric"), 
                      n_obs=vector(mode="integer"),
                      method=vector(mode="character"),
                      lambda = vector(mode="numeric"))

outdir = "data/output/250/"
for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  errors <- data.frame(rmse=rep(0,times = nsim),
                       n_obs = rep(0L,times = nsim),
                       method = rep("",times = nsim),
                       lambda = rep(0, times=nsim))
  for(l in 1:length(lambda)){
  for(i in 0:(nsim-1)){
    f <- matrix(as.matrix(readMM(paste0(outdir, l,"/",i, "/", "nonlinear.mtx"))), 
                nrow=n_dofs, ncol= n_time_locs)
    
    
    evals = as.vector(eval.parabolic(f, mesh, eval.locs, time_locs)[,2:n_time_locs])
    errors$rmse[i+1] = rmse(evals, test.evals)
    errors$n_obs[i+1] = nlocs
    errors$method[i+1] = "STRPDE-NL"
    errors$lambda[i+1] = lambda[l]
  }
  results = rbind(results, errors) 
  }
}

write.table(results, paste0(outdir, "nonlinear_rmse.txt") )
#

# nonlinear <- read.table(paste0(datadir, "nonlinear_rmse.txt"))[c("rmse", "n_obs", "method")]
# parabolic <- read.table(paste0(datadir, "parabolic_rmse.txt"))[c("rmse", "n_obs", "method")]
# separable <- read.table(paste0(datadir, "separable_rmse.txt"))[c("rmse", "n_obs", "method")]
#tps <- read.table(paste0(datadir, "tps_rmse.txt"))[c("rmse", "n_obs", "method")]

# rmse <- rbind(nonlinear, parabolic) #, separable)#, tps)
# rmse$n_obs <- as.factor(rmse$n_obs)
# rmse$method <- factor(rmse$method, levels = c("STRPDE-NL", "STRPDE-PAR", 
#                                               "STRPDE-SEP")) #, "TPS"))

# table(rmse$method)
# 
# pdf(paste0(paste0(datadir, nlocs, "/", "rmse.pdf")))
# plt <- boxplot.rmse(rmse, easy = T) + MyTheme + theme(legend.position = "none")
# print(plt)
# dev.off()
# noTPS = rmse[rmse$method!= "TPS",]
# boxplot(noTPS$rmse~ noTPS$method)
# 
# noSEP = noTPS[noTPS$method!= "STRPDE-SEP",]
# boxplot(noSEP$rmse ~ noSEP$method)
# {
#   tmp <- plt + 
#     labs(title="RMSE", x="observations") +
#     theme(legend.position = "bottom")  
#   
#   leg <- cowplot::get_plot_component(tmp, 'guide-box', return_all = TRUE)
#   
#   pdf(paste0(datadir, "RMSE-legend.pdf"), width = 16, height = 1)
#   grid.draw(leg[[3]])
#   dev.off()
# }
# {
#   tmp <- plt + 
#     labs(title="RMSE", x="observations") +
#     theme(legend.position = "none")
#   ggsave(paste0(datadir, "RMSE-no-legend.pdf"), plot=tmp, width = 8, height = 7)
# }


# ---

parabolic = read.table(paste0(outdir, "parabolic_rmse.txt") )
separable = read.table(paste0(outdir, "separable_rmse.txt") )
nonlinear = read.table(paste0(outdir, "nonlinear_rmse.txt") )

library(dplyr)
rmse.par = parabolic %>% group_by(lambda) %>% summarise( q1 = quantile(rmse, probs = 0.25),
                                                       q2 = median(rmse),
                                                       q3 = quantile(rmse, probs = 0.75))
rmse.sep = separable %>% group_by(lambda) %>% summarise( q1 = quantile(rmse, probs = 0.25),
                                                           q2 = median(rmse),
                                                           q3 = quantile(rmse, probs = 0.75))
rmse.nl = nonlinear %>% group_by(lambda) %>% summarise( q1 = quantile(rmse, probs = 0.25),
                                                         q2 = median(rmse),
                                                         q3 = quantile(rmse, probs = 0.75))


myred = rgb( 139,0,0, max=255, alpha=125)
myblue = rgb( 65, 105, 225, max = 255, alpha = 125 )
mygreen = rgb( 0, 255, 0, max = 255, alpha = 125 )

xx = log(rmse.par$lambda/(n_time_locs*nlocs))
{
pdf(paste0(outdir, "rmse.pdf"))
plot(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.par$q2,
     type="l", lwd=3, col="red4",
     xlab = expression(log(lambda)), ylab="")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.sep$q2, 
       type="l", lwd=3, col="green4")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.nl$q2, 
       type="l", lwd=3, col="blue4")
polygon( x= c(rev( log(rmse.par$lambda/(n_time_locs*nlocs))),
              log(rmse.par$lambda/(n_time_locs*nlocs))),
         y= c( rev(rmse.par$q1), rmse.par$q3),
         col = myred, border = NA )
polygon( x= c(rev( log(rmse.par$lambda/(n_time_locs*nlocs))),
              log(rmse.par$lambda/(n_time_locs*nlocs))),
         y= c( rev(rmse.nl$q1), rmse.nl$q3),
         col = myblue, border = NA )
polygon( x= c(rev( log(rmse.par$lambda/(n_time_locs*nlocs))),
              log(rmse.par$lambda/(n_time_locs*nlocs))),
         y= c( rev(rmse.sep$q1), rmse.sep$q3),
         col = mygreen, border = NA )
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.par$q2,
       type="l", lwd=2, col="red4")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.nl$q2, 
       type="l", lwd=3, col="blue4")
points(log(rmse.par$lambda/(n_time_locs*nlocs)), rmse.sep$q2, 
       type="l", lwd=2, col="green4")
legend(-7.5, y=0.35, legend=c("NL", "PAR","SEP"), 
       col = c("blue4", "red4", "green4"), lty=1, cex=2,
       lwd=4,
       box.lty=0)
dev.off()
}
# --- plots of the estimates ---

par.idx = which.min(rmse.par$q2)
sep.idx = which.min(rmse.sep$q2)
nl.idx = which.min(rmse.nl$q2)

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}

for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  f <- matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f <- f + matrix(as.matrix(readMM(paste0(outdir, nl.idx, "/", i, "/", "nonlinear.mtx"))), 
                    nrow=n_dofs, ncol= n_time_locs)/nsim
  }
  
  n_breaks <- 50
  # qua andrebbe considerata anche la "exact" per fare grafici coerenti (?)
  f_separable = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f_separable <- f_separable + as.matrix(read.table(paste0(outdir, sep.idx, "/", i, "/", "separable.txt"), header = F))/nsim
  }
  
  f_parabolic = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f_parabolic <- f_parabolic + as.matrix(read.table(paste0(outdir, "/", par.idx, "/", i, "/", "parabolic.txt"), header = F))/nsim
  }
  #obs_range = c(1e10,-1e10)
  obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt"))) # uso solo la prima...
  obs_range = range(obs)
  
  # NON fittato su t0 = 0
  # f_tps = matrix(0,nrow=n_dofs, ncol=n_time_locs-1)
  # for(i in 0:(nsim-1)){
  #   f_tps <- f_tps + as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "tps.txt"), header = F))/nsim
  # }
  # 
  lims = c( min(f, f_exact, f_separable,obs_range),  # f_tps), 
            max(f, f_exact, f_separable,obs_range)) #f_tps) )
  mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 
  
  {
    linewidth=2
    locs =  as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/locs.txt")))
    obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt")))
    pdf(paste0(outdir,"/data.pdf"))
    for(t in 1:n_time_locs){
      data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
      plt <- ggplot() +
        geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
        geom_segment(data=data.frame(x=min(mesh$nodes[,1]), xend=max(mesh$nodes[,1]),
                                     y=max(mesh$nodes[,2]), yend=max(mesh$nodes[,2])),
                     aes(x=x, xend=xend, y=y, yend=yend), color="black", 
                     linewidth=linewidth) +
        geom_segment(data=data.frame(x=min(mesh$nodes[,1]), xend=min(mesh$nodes[,1]),
                                     y=min(mesh$nodes[,2]),yend=max(mesh$nodes[,2])),
                     aes(x=x, xend=xend, y=y, yend=yend), color="red",
                     linewidth=linewidth) +
        geom_segment(data=data.frame(x=min(mesh$nodes[,1]), xend=max(mesh$nodes[,1]), 
                                     y=min(mesh$nodes[,2]),yend=min(mesh$nodes[,2])),
                     aes(x=x, xend=xend, y=y, yend=yend), color="black",
                     linewidth=linewidth) +
        geom_segment(data=data.frame(x=max(mesh$nodes[,1]), xend=max(mesh$nodes[,1]),
                                     y=min(mesh$nodes[,2]),yend=max(mesh$nodes[,2])),
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
    pdf(paste0(paste0(outdir, "/","nonlinear.pdf")))
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
    pdf(paste0(outdir, "/", "exact.pdf"))
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
    pdf(paste0(paste0(outdir, "/","separable.pdf")))
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
    pdf(paste0(paste0(outdir, "/","parabolic.pdf")))
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
  
  # {
  #   pdf(paste0(paste0(datadir, nlocs, "/","y_tps.pdf")))
  #   for(t in 1:(n_time_locs-1)){
  #     data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_tps[,t]))    
  #     plt <-
  #       ggplot(aes(x = x, y = y, z = z), data = data) +
  #       geom_contour_filled(breaks = mybreaks) +
  #       scale_color_viridis(limits=lims) +
  #       coord_fixed() + theme_void() +
  #       theme(legend.position = "none")
  #     print(plt)
  #   }
  #   dev.off()
  # }
  # 
  
  
}


