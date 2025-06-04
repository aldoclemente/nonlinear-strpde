library(viridis)
library(colorspace)
library(ggplot2)
library(gridExtra)
library(grid)
title.size <- 26
MyTheme <- theme(
  axis.text = element_text(size=title.size-5),
  axis.title = element_text(size=title.size),
  title = element_text(size=title.size),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size=title.size-5),
  legend.key.size = unit(1,"cm"),
  legend.key.height = unit(1,"cm"),
  legend.title = element_blank(),
  legend.background = element_rect(fill="white", color="white",
                                   linewidth =c(1,0.5))
)

# .... sqrt ( ... ) ....
rmse <- function(x,y){
  return ( sqrt( mean( (x-y)^2 ) ) )
}

eval.parabolic = function(coeff, mesh, locs, time_mesh){
  res = matrix(0, nrow=nrow(locs), ncol=length(time_mesh))
  FEMbasis = create.FEM.basis(mesh)
  for(t in 1:length(time_mesh)){
    res[,t] = eval.FEM(FEM(coeff = coeff[,t], FEMbasis), locations = locs)
  }
  return(res)
}

boxplot.rmse = function(rmse, easy = F){
    begin=0.25
    end=0.95
    border_col = darken(viridis(length(unique(rmse$method)), begin=begin,end=end), amount=0.25)
    fill_col = viridis(length(unique(rmse$method)), begin=begin, end=end)
    BORDER = c()
    FILL = c()
    for(i in 1:length(unique(rmse$method))){
      FILL = append(FILL, fill_col[i])
      BORDER = append(BORDER, border_col[i])
    }
    ggFILL <-scale_fill_manual(values = FILL)
    ggBORDER <- scale_color_manual(values= BORDER)
    
    if ( !easy) {
    p<-ggplot(rmse)+
      geom_boxplot(aes(x=n_obs,
                       y=rmse, group=interaction(method,n_obs),
                       fill=method, color = method))+
      scale_x_discrete(limits=as.character(unique(rmse$n_obs)))+
      labs(x="", y="") +
      theme(
        axis.ticks.x = element_blank(),
        legend.position = c(0.95,0.95),
        legend.background = element_rect(fill="white", color="black",
                                         linewidth =c(1,0.5)),
        legend.title = element_blank())
    
    p <- p + ggFILL + ggBORDER
    }else{
      p<-ggplot(rmse)+
        geom_boxplot(aes(x=method,
                         y=rmse,
                         fill=method, color = method))+
        #scale_x_discrete(limits=as.character(unique(rmse$n_obs)))+
        labs(x="", y="") +
        theme(
          axis.ticks.x = element_blank(),
          legend.position = c(0.95,0.95),
          legend.background = element_rect(fill="white", color="black",
                                           linewidth =c(1,0.5)),
          legend.title = element_blank())
      
      p <- p + ggFILL + ggBORDER
    }
    p
}
  
plot.fem <- function(coeff, dofs, 
                     filename="/home/aldoclemente/Desktop/nonlinear-strpde/fem.pdf"){
  
  mycolors<- function(x) {
    colors<-viridis( x + 1 )
    colors[1:x]
  }
  
  n_breaks <- 50
  n_time_locs = ncol(coeff)
  #breaks <- seq(0,1, length.out = n_breaks)
  breaks <- seq(min(coeff),max(coeff), length.out = n_breaks)
  f <- matrix(coeff, nrow=nrow(coeff), ncol=n_time_locs)
  pdf(file=filename)
  for(t in 1:n_time_locs){
    data <- data.frame(x = round(dofs[,1], 10), y = round(dofs[,2], 10), 
                       z = as.matrix(f[,t]))    
    plt <-
      ggplot(aes(x = x, y = y, z = z), data = data) +
      geom_contour_filled(breaks = breaks) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}

plot.data <- function (coeff, locs,
              filename= "/home/aldoclemente/Desktop/nonlinear-strpde/data.pdf"){
  mycolors<- function(x) {
    colors<-viridis( x + 1 )
    colors[1:x]
  }
  
  n_breaks <- 50
  n_time_locs = ncol(coeff)
  breaks <- seq(min(coeff),max(coeff), length.out = n_breaks)
  pdf(file=filename)
  for(t in 1:n_time_locs){
    data <- data.frame(x = locs[,1], y = locs[,2], z = as.matrix(coeff[,t]))    
    plt <- ggplot() +
      geom_point(aes(x=x, y=y, color=z), data=data, size = 5) +
      scale_color_gradientn(colors =mycolors(n_breaks + 2)) +
      coord_fixed() + theme_void() +
      theme(legend.position = "none")
    print(plt)
  }
  dev.off()
}


