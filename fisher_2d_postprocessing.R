#setwd("~/Desktop/nonlinear-strpde/")

rm(list=ls())
source("~/Desktop/readMesh/read.mesh.R")
library(fdaPDE)
nlocs = c(250) #c(1000) # c(250)

datadir = "data/"

Nt = as.integer(read.table( paste0(datadir,"Nt.txt")))

filename <- paste0(datadir, "square.mesh")
domain <- read.mesh(filename)
mesh <- create.mesh.2D(nodes=domain$nodes, triangles = domain$elements)
fembasis = create.FEM.basis(mesh)
lambda_aldo = 1.
#lambda_aldo = 0.001
#lambda_aldo = 0.000363636
#lambda_aldo = 0.000190476
#lambda_aldo = 9.09091e-05
nsim = 1
alpha = 3.0
time_mesh = seq(0,1, length.out = Nt)
lambda = 1.
#fisher = matrix(nrow=fembasis$nbasis, ncol=Nt)
# for(i in 1:Nt){
#   file = paste0("freefem_data/fisher_",i-1,".txt")
#   fisher[,i] = as.matrix( read.table(file) )
# }

fisher = read.table( paste0(datadir,"exact.txt"))

for(n in 1:length(nlocs)){
  resdir = paste0(datadir, nlocs[n],"/")
  parabolic = matrix(0, nrow=fembasis$nbasis*Nt, ncol=nsim)
  misfit = parabolic
  for(sim in 1:nsim){
    simdir = paste0( resdir, sim-1, "/")
    locs = as.matrix(read.table(paste0(simdir, "locs.txt")))
    obs = as.matrix(read.table(paste0(simdir, "obs.txt")))
    #obs = matrix(0, nrow=nrow(locs), ncol=Nt)
    # for(t in 1:Nt){
    #   obs[,t] = as.matrix(read.table(paste0(datadir, "obs_",t-1,".txt")))
    # }
    output = smooth.FEM.time(observations = obs[,2:Nt], # ok la prima deve saltare
                             time_mesh = time_mesh, lambdaS = lambda, lambdaT = 1,
                             locations = locs, FEMbasis = fembasis, IC = fisher[,1],
                             PDE_parameters = list(K=rbind(c(0.1,0), 
                                                           c(0., 0.1)), 
                                                   b=c(0,0), 
                                                   c=-alpha, u = NULL), 
                             FLAG_PARABOLIC = TRUE)
    parabolic[,sim] = output$fit.FEM.time$coeff
    misfit[,sim] = output$PDEmisfit.FEM.time$coeff
  }
  write.table(parabolic, paste0(resdir, "parabolic.txt"))
  misfit = as.matrix(misfit[nrow(misfit):1,])
  write.table(misfit, paste0(resdir, "misfit.txt"), row.names = F, col.names = F)
}

# ---
#write.csv(locs, "freefem_data/datacpp/locs.csv")
library(viridis)
nlevels = 20
color.palette = viridis

for(n in 1:length(nlocs)){
  aldo = matrix(0, nrow=fembasis$nbasis*Nt, ncol=nsim)
  resdir = paste0(datadir, nlocs[n],"/")
  for(sim in 1:nsim){
    simdir = paste0( resdir, sim-1, "/", "lambda_", lambda_aldo, "/")
    #tmp = matrix(0, nrow=fembasis$nbasis, ncol=Nt)
    tmp = as.matrix(read.table(paste0(simdir, "iterative.txt")))
    aldo[,sim] = as.vector(tmp)
  }
  write.table(aldo, paste0(resdir, "iterative_lambda_", lambda_aldo, ".txt"))
}

range(aldo) 
range(parabolic)
errors = matrix(nrow=0, ncol=3)
colnames(errors) <- c("err", "n_obs", "method") 
#errors <- as.data.frame(errors)
for( n in 1:length(nlocs)){
  resdir = paste0(datadir, nlocs[n],"/")
  aldo = read.table(paste0(resdir, "iterative_", lambda_aldo, ".txt"))
  parabolic = read.table(paste0(resdir, "parabolic.txt"))
  for(sim in 1:nsim){
    errors = rbind(errors, c(
      sqrt( mean( (aldo[,sim]-as.vector(fisher))^2 ) ),
      nlocs[n], "STR-PDE-NL"))
    errors = rbind(errors, c(
      sqrt( mean( (parabolic[,sim]-as.vector(fisher))^2 ) ),
      nlocs[n], "STR-PDE"))
  }
  aldo = rowMeans(aldo); 
  parabolic = rowMeans(parabolic)
  aldo = matrix(aldo, nrow=fembasis$nbasis, ncol=Nt)
  parabolic = matrix(parabolic, nrow=fembasis$nbasis, ncol=Nt)
  
  locs = read.table(paste0(resdir, "0/locs.txt"))
  obs = read.table(paste0(resdir, "0/obs.txt"))

  lims = range(c(range(fisher), range(aldo), range(parabolic)), range(obs))
  #lims = range(c(range(fisher), range(aldo), range(obs)))
  
  pdf(paste0(resdir, "aldo_lambda_", lambda_aldo, ".pdf"), family="serif", w=8,h=8)
  for(t in 1:Nt){
    filled.contour(z = matrix(aldo[,t], 
                              nrow=sqrt(fembasis$nbasis), 
                              ncol=sqrt(fembasis$nbasis)),nlevels = nlevels,
                   color.palette = color.palette, zlim = lims,axes=F)
  }
  dev.off()
  
  
  pdf(paste0(resdir, "exact.pdf"), family="serif", w=8,h=8)
  for(t in 1:Nt){
    filled.contour(x = seq(-2.5, 2.5, length.out=sqrt(fembasis$nbasis)),
                   y = seq(-2.5, 2.5, length.out=sqrt(fembasis$nbasis)),
                   z = matrix(fisher[,t], 
                              nrow=sqrt(fembasis$nbasis), 
                              ncol=sqrt(fembasis$nbasis)),nlevels = nlevels,
                   color.palette = color.palette, zlim = lims,axes=F)
  }
  dev.off()
  
  pdf(paste0(resdir, "parabolic.pdf"), family="serif", w=8,h=8)
  for(t in 1:Nt){
    filled.contour(z = matrix(parabolic[,t],
                              nrow=sqrt(fembasis$nbasis),
                              ncol=sqrt(fembasis$nbasis)),nlevels = nlevels,
                   color.palette = color.palette, zlim = lims,axes=F)
  }
  dev.off()
  # 
  ncol=nlevels
  col = color.palette(ncol)
  #locs = read.table(paste0(resdir, "0/locs.txt"))
  
  {
    pdf(paste0(resdir, "data_0.pdf"),family="serif", w=8, h=8)
    for(t in 1:Nt){
      filled.contour(x = seq(-2.5, 2.5, length.out=sqrt(fembasis$nbasis)),
                     y = seq(-2.5, 2.5, length.out=sqrt(fembasis$nbasis)),
                     z = matrix(0, 
                                nrow=sqrt(fembasis$nbasis), 
                                ncol=sqrt(fembasis$nbasis)), 
                     axes = F, ann=F, col="white", 
                     plot.axes={points(locs, pch=16, cex=3, 
                                       col=col[round((ncol-1)*(obs[,t]-min(lims))/diff(lims))+1])
                       rect(xleft = -2.5, xright = 2.5, ybottom = -2.5, ytop = 2.5,lwd = 4)})
    }
    dev.off()
  }
  
}

# boxplots ---
# errors = as.data.frame(errors)
# errors$err <- as.numeric(errors$err)
# errors$n_obs <- factor(errors$n_obs, levels =  c("100", "250", "500", "1000"))
# errors$method <- factor(errors$method, levels = c("STR-PDE-NL", "STR-PDE"))
# head(errors)
# 
# n_obs_ = as.numeric(levels(errors[["n_obs"]]))
# methods_ = levels(errors[["method"]])
# at_ <- c()
# for(i in 1:length(n_obs_)){
#   at_ <-  c(at_, ((i-1)*(1+length(methods_)) + (1:length(methods_))))
# }
# 
# fill_col = viridis::viridis((length(levels(errors[["method"]]))+1), begin=0.25, end=0.95)
# fill_col = fill_col[1:length(methods_)]
# facs = names(Filter(is.factor, errors))
# which( ! names(errors) %in% facs )
# doplot = names(errors)[which( ! names(errors) %in% facs )]
# 
# if(length(methods_)%%2 != 0){
#   at_label = seq(ceiling(length(methods)/2), at_[length(at_)], by=(length(methods_) + 1))
# }else{
#   at_label = seq(length(methods_)/2, at_[length(at_)], by=(length(methods_) + 1))
# }
# 
# {
#   pdf("freefem_data/errors.pdf", family = "serif", width = 7, height = 7)
#   for(i in doplot){
#     boxplot(errors[[i]] ~  errors[["method"]] + as.numeric(errors[["n_obs"]]),
#             ylab="", xlab="observations", at = at_, xaxt="n",
#             #ylim=c(min(errors[[i]]), max(errors[[i]])*(1.1)),
#             col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
#             main = i, ylim=c(0,0.185))
#     axis(side = 1, at = at_label, labels = n_obs_, cex.lab = 2, cex.axis = 2)
#     legend("topright",legend=methods_, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
#            bty="n")
#     
#   }
#   dev.off()
# }
# 
