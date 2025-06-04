library(mgcv)
library(Matrix)
library(fdaPDE)
#rm(list=ls())

source("../utils.R")
eval.locs = as.matrix(read.csv("data/eval_locs.csv")[,2:3])

datadir = "data/"

#NLOCS = c(100, 250, 500, 1000)
NLOCS = c(250)
time_locs = as.matrix(read.table(paste0(datadir, "time_locations.txt")))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)

boundary = as.matrix(readMM("data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("data/mesh/elements.mtx"))

bd1 = dofs[boundary==1, ]
bd1 = bd1[nrow(bd1):1,]
bd2 = dofs[boundary==2, ]
bd2 = bd2[nrow(bd2):1,]
bd3 = dofs[boundary==3, ]
#bd3 = bd3[nrow(bd3):1, ]
bd4 = dofs[boundary==4, ]
#bd4 = bd4[nrow(bd4):1, ]

bnd = list()
bnd[[1]] = list()
#bd[[1]][[1]] = c(bd1[,1], bd2[,1], bd3[,1], bd4[,1])
bnd[[1]][[1]] = c(bd4[,1], bd3[,1], bd2[,1], bd1[,1])

#bd[[1]][[2]] = c(bd1[,2], bd2[,2], bd3[,2], bd4[,2])
bnd[[1]][[2]] = c(bd4[,2], bd3[,2], bd2[,2], bd1[,2])

names(bnd[[1]]) = c("x0", "x1")
bnd[[1]]$x0 = c(bnd[[1]]$x0, bnd[[1]]$x0[1])
bnd[[1]]$x1 = c(bnd[[1]]$x1, bnd[[1]]$x1[1])

bnd

mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 
FEMbasis = create.FEM.basis(mesh)
nsim = 30

# exact 
exact = as.matrix( read.table(paste0(datadir, "exact.txt")))
test.evals = as.vector(eval.parabolic(exact, mesh, eval.locs, time_locs)[,2:n_time_locs])
#

eval.locs <- data.frame(x0=eval.locs[,1], 
                        x1=eval.locs[,2],
                        t = rep(time_locs[2:n_time_locs], each=nrow(eval.locs)))


k_space <- 50 #30
k_time <- 8       
j_knots = 8      
eps_tp=0.1  # come arnone-vicini 
knots <- data.frame(x0=rep(seq(0+2*eps_tp,1-2*eps_tp,length.out = j_knots),j_knots),
                    x1=rep(seq(0+2*eps_tp,1-2*eps_tp,length.out = j_knots), each=j_knots))

marks = mesh$nodesmarkers == 0 
# dat.plot <- data.frame(
#   x0 = rep(mesh$nodes[marks,1], times=n_time_locs-1),
#   x1 = rep(mesh$nodes[marks,2], times=n_time_locs-1),
#   t = rep(time_locs[2:n_time_locs], each=sum(marks)))

dat.plot <- data.frame(
  x0 = rep(mesh$nodes[,1], times=n_time_locs-1),
  x1 = rep(mesh$nodes[,2], times=n_time_locs-1),
  t = rep(time_locs[2:n_time_locs], each=n_dofs))

results <- data.frame(rmse=vector(mode="numeric"), 
                      n_obs=vector(mode="integer"),
                      method=vector(mode="character"))

for(k in 1:length(NLOCS)){
  nlocs = NLOCS[k]
  
  errors <- data.frame(rmse=rep(0,times = nsim),
                       n_obs = rep(0L,times = nsim),
                       methd = rep("",times = nsim))
  
  
  for(i in 0:(nsim-1)){
    locs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "locs.txt")))
    obs = as.matrix(read.table(paste0(datadir, nlocs, "/", i, "/", "obs.txt")))
    
    dat <- data.frame(y = as.vector(obs[,2:n_time_locs]),
                      x0 = rep(locs[,1], n_time_locs-1),
                      x1 = rep(locs[,2], n_time_locs-1),
                      t = rep(time_locs[2:n_time_locs], each = nlocs))
    # 
    # dat <- data.frame(y = as.vector(obs[,1:n_time_locs]),
    #                   x0 = rep(locs[,1], n_time_locs),
    #                   x1 = rep(locs[,2], n_time_locs),
    #                   t = rep(time_locs[1:n_time_locs], 
    #                           each = nlocs))
    
    tps <- gam(y ~ te(x0, x1, t, k = c(k_space, k_time),
                      d = c(2, 1), bs = c("tp", "ps"), xt=list(bnd=bnd)), 
               knots=knots, data = dat)
    
    # soap <- gam(y ~ te(x0, x1, t, k = c(k_space, k_time),
    #                   d = c(2, 1), bs = c("sf", "ps"), xt=list(bnd=bnd)),
    #            knots=knots, data = dat)
    
    # 
    # soap.sw <- gam(y ~ te(x0, x1, t, k = c(k_space, k_time),
    #                    d = c(2, 1), bs = c("sw", "ps"), xt=list(bnd=bnd)), 
    #             knots=knots, data = dat)
    
    #invisible( capture.output(separable <-  mgcv ))
    #coeff = eval.FEM.time(separable$fit.FEM.time, 
    #                      locations = dofs, time.instants = time_locs)
    #gam.check(soap)
    #gam.check(tps)
    #gam.check(soap.sw)
    coeff = predict(tps, newdata = dat.plot)
    #coeff = matrix(coeff, nrow=sum(marks), ncol=n_time_locs-1)
    coeff = matrix(coeff, nrow=n_dofs, ncol=n_time_locs-1)
    
    # coeff = predict(soap, newdata = dat.plot,
    #         type = "terms", terms = "te(x0,x1,t)") + coef(soap)["(Intercept)"]
    # 
    write.table( coeff, 
                 paste0(datadir, nlocs, "/", i, "/", "tps.txt"),
                 col.names = F, row.names = F)
    evals = predict(tps, newdata = eval.locs)
    errors$rmse[i+1] = rmse(evals, test.evals)
    errors$n_obs[i+1] = nlocs
    errors$method[i+1] = "TPS"
  }
  results = rbind(results, errors)
  
}

write.table(results, paste0(datadir, "tps_rmse.txt") )


# 
# mycolors<- function(x) {
#   colors<-viridis( x + 1 )
#   colors[1:x]
# }
# 
# 
# library(ggplot2)
# library(viridis)
# n_breaks <- 50
# n_time_locs = 11
# breaks <- seq(min(coeff,na.rm = T),max(coeff, na.rm = T), length.out = n_breaks) 
# {
#   pdf(paste0("/home/aldoclemente/Desktop/nonlinear-strpde/tps.pdf"))
#   for(t in 1:(n_time_locs-1)){
#     data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), 
#                        z = as.matrix(coeff[,t]))    
#     plt <-
#       ggplot(aes(x = x, y = y, z = z), data = data) +
#       geom_contour_filled(breaks = breaks) +
#       coord_fixed() + theme_void() +
#       theme(legend.position = "none")
#     print(plt)
#   }
#   dev.off()
# }
# 
# 
# {
#   pdf(paste0("/home/aldoclemente/Desktop/nonlinear-strpde/soap.pdf"))
#   for(t in 1:(n_time_locs-1)){
#     data <- data.frame(x = round(dofs[marks,1], 10), y = round(dofs[marks,2], 10), 
#                        z = as.matrix(coeff[,t]))    
#     plt <-
#       ggplot(aes(x = x, y = y, z = z), data = data) +
#       geom_contour_filled(breaks = breaks) +
#       coord_fixed() + theme_void() +
#       theme(legend.position = "none")
#     print(plt)
#   }
#   dev.off()
# }
# 
# 
