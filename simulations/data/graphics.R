library(ggplot2)
library(viridis)
plot.plottable_function <- function(x, palette = NULL, lims = NULL, ...) {
  n_col <- 100
  if (is.null(palette)) {
    palette_ <- colorRampPalette(colors = c("lightyellow", "darkred"))(n_col)
  } else if(is.function(palette)) {
    palette_ <- palette(n_col)
  } else { # user-defined palette
    palette_ <- palette_
  }
  
  nodes <- x$geometry$nodes
  x_grid <- seq(min(nodes[, 1]), max(nodes[, 1]), length.out = 250)
  y_grid <- seq(min(nodes[, 2]), max(nodes[, 2]), length.out = 250)
  xy_grid <- expand.grid(x_grid, y_grid)
  ## evaluate fe_function at fine grid
  vals <- x$eval(xy_grid)
  
  ## use provided limits or compute from the data
  if (is.null(lims)) {
    lims <- range(vals, na.rm = TRUE)
  }
  
  ## cut the values into bins according to lims
  breaks <- seq(lims[1], lims[2], length.out = n_col + 1)
  col_idx <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  col <- palette_[col_idx]
  #col <- palette_[as.numeric(cut(vals, breaks = n_col))]
  plot(
    xy_grid[, 1],
    xy_grid[, 2],
    xlab = "",
    ylab = "",
    pch = 15,
    col = col,
    asp = 1,
    cex = .6, ...
  )
}

# FElist -> list of fe_function 
# data_list -> data_list = list(data = #n_locs x #..., locations = #n_locs x ...)
plot.fem.base <- function(FElist, data_list = NULL, raster_list =NULL,
                          filename=NULL, filename_data=NULL, filename_raster=NULL, palette=viridis){
  
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
  
  if(!is.null(raster_list)){
    for(i in 1:length(raster_list)){
      lims = c(min(c( range(values(raster_list[[i]]), na.rm = T), lims[1]) , na.rm = T),
               max(c( range(values(raster_list[[i]]), na.rm = T), lims[2]), na.rm = T))
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
  
  if(!is.null(raster_list)){
    n_col = 100
    breaks <- seq(lims[1], lims[2], length.out = n_col + 1)
    mesh = FElist[[1]]$geometry
    nodes = mesh$nodes
    edges = mesh$edges
    bd_edges = edges[mesh$boundary_edges == 1,]
    if(!is.null(filename_raster)) pdf(file=filename_raster)
    for(i in 1:length(raster_list)){
      
      vals = as.data.frame(raster_list[[i]])
      vals = vals[!is.na(vals),1]
      col_idx <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      col <- palette_[col_idx]
      plot(raster_list[[i]], xlab = "", ylab = "", asp = 1,  col=palette(n_col), box=FALSE, legend=FALSE,
           frame.plot=F, xaxt="n", yaxt="n", axes=F)
      segments(nodes[bd_edges[,1],1], nodes[bd_edges[,1],2],
               nodes[bd_edges[,2],1], nodes[bd_edges[,2],2], asp=1)
    }
    if(!is.null(filename_raster)) dev.off()
  }
}


plot.fem <- function(coeff, dofs, 
                     filename=NULL){
  
  mycolors<- function(x) {
    colors<-viridis( x + 1 )
    colors[1:x]
  }
  
  if(!is.null(filename)) pdf(file=filename)
  for(t in 1:n_time_locs){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), 
                       z = as.matrix(f[,t]))    
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_contour_filled(breaks = breaks) +
      scale_color_viridis_d(limits=c(min(coeff),max(coeff))) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  if(!is.null(filename_raster))dev.off()
  
  
}
