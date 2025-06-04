library(Matrix)
library(ggplot2)
library(viridis)
rm(list=ls())
datadir = "data/"

time_locs = as.matrix(read.table(paste0("../simulation_1/data/time_locations.txt")))
n_time_locs = length(time_locs)
dofs <- as.matrix(readMM("../simulation_1/data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("../simulation_1/data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("../simulation_1/data/mesh/elements.mtx"))

library(fdaPDE)
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 

nsim = 30

f_exact =  as.matrix(read.table(paste0("../simulation_1/data/exact.txt")))

mycolors<- function(x) {
  colors<-viridis( x + 1 )
  colors[1:x]
}


f <- matrix(0,nrow=n_dofs, ncol=n_time_locs)
for(i in 0:(nsim-1)){
  f <- f + matrix(as.matrix(readMM(paste0(datadir, i, "/", "y.mtx"))), nrow=n_dofs, ncol= n_time_locs)/nsim
}
  
n_breaks <- 50
# qua andrebbe considerata anche la "exact" per fare grafici coerenti (?)
  # f_separable = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  # for(i in 0:(nsim-1)){
  #   f_separable <- f_separable + as.matrix(read.table(paste0(datadir, "/", i, "/", "separable.txt"), header = F))/nsim
  # }
  
  f_parabolic = matrix(0,nrow=n_dofs, ncol=n_time_locs)
  for(i in 0:(nsim-1)){
    f_parabolic <- f_parabolic + as.matrix(read.table(paste0(datadir, "/", i, "/", "parabolic.txt"), header = F))/nsim
  }
  #obs_range = c(1e10,-1e10)
  obs = as.matrix(read.table(paste0(datadir, "/", 0,"/obs.txt"))) # uso solo la prima...
  obs_range = range(obs)
  f_separable = 0
  lims = c( min(f, f_exact, f_parabolic, obs_range), 
            max(f, f_exact, f_parabolic, obs_range) )
  mybreaks <- seq( lims[1], lims[2], length.out = n_breaks) 
  
  # {
  #   locs =  as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/locs.txt")))
  #   obs = as.matrix(read.table(paste0(datadir, "/",nlocs, "/", 0,"/obs.txt")))
  #   pdf(paste0(datadir, nlocs,"/data.pdf"))
  #   for(t in 1:n_time_locs){
  #     data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(obs[,t]))    
  #     plt <- ggplot() +
  #       geom_rect(data = data.frame(xmin = 0,ymin = 0,xmax = 1,ymax = 1), 
  #                 aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
  #                 fill = "white", 
  #                 color = "black",
  #                 linewidth = 1) +
  #       geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
  #       scale_color_gradientn(colors =mycolors(n_breaks + 2), limits=lims) +
  #       coord_fixed() + theme_void() +
  #       theme(legend.position = "none")
  #     print(plt)
  #   }
  #   dev.off()
  # }
  # 
  
  {
    pdf(paste0(paste0(datadir,"y.pdf")))
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

  {
    pdf(paste0(datadir, "y_exact.pdf"))
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
  
  # {
  #   pdf(paste0(paste0(datadir,"y_separable.pdf")))
  #   for(t in 1:n_time_locs){
  #     data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), z = as.matrix(f_separable[,t]))
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
  {
    pdf(paste0(paste0(datadir,"y_parabolic.pdf")))
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
  
 
  incidence_matrix = as.matrix(read.table("data/incidence_matrix.txt", header = F))
  poly_coords = list()
  for(r in 1:nrow(incidence_matrix)){
    triangles = mesh$triangles[which(incidence_matrix[r,] == 1),]
    #triangles = cbind(triangles, triangles[,1]) 
    poly_coords[[r]] = list()
    for(t in 1:nrow(triangles)){
      poly_coords[[r]][[t]] = mesh$nodes[triangles[t,], ]
      poly_coords[[r]][[t]] = rbind(poly_coords[[r]][[t]], 
                                    mesh$nodes[triangles[t,][1], ])
    }
    } 
{
  pdf(paste0(datadir,"data.pdf"))
  for(t in 1:n_time_locs){
  plt <- ggplot()
  for(r in 1:nrow(incidence_matrix)){
    for(i in 1:length(poly_coords[[r]])){
      plt <- plt +
        geom_polygon(data=data.frame(x=poly_coords[[r]][[i]][,1], 
                                     y=poly_coords[[r]][[i]][,2],
                                     z=rep(obs[r,t], times=4)), 
                 aes(x = x, y = y, fill = z), color=NA)  
    }
  }
  plt <- plt  +
    scale_fill_viridis(option = "viridis", limits=lims) + 
    coord_fixed() + theme_void() +  
    theme(legend.position = "none")
  
  plt_rect <- plt + geom_rect(data = data.frame(xmin = 0,ymin = 0,xmax = 1,ymax = 1), 
                              aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                              fill = "transparent", 
                              color = "black",
                              linewidth = 1)
  print(plt_rect)
  }
  dev.off()
}
  