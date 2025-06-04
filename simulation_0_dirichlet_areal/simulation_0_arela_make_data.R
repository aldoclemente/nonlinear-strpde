
library(fdaPDE)
dist = function(x,y){
  return( sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2) )
}
dofs <- as.matrix(readMM("../simulation_0_dirichlet/data/mesh/nodes.mtx"))
n_dofs = nrow(dofs)
boundary = as.matrix(readMM("../simulation_0_dirichlet/data/mesh/boundary.mtx"))
elements =  as.matrix(readMM("../simulation_0_dirichlet/data/mesh/elements.mtx"))
mesh = create.mesh.2D(dofs, triangles = (elements + 1)) 
FEMbasis = create.FEM.basis(mesh)

#x = seq(0.125,0.875, by=0.25)
#y = x 
h = abs(mesh$nodes[1,1] - mesh$nodes[2,1]) 
x = seq(-2.5 + 2*h, 2.5-2*h, by=3*h)
y = x

centers = expand.grid(x,y)
plot(mesh, pch=".")
points(centers, pch=16, col="red")
incidence_matrix = matrix(0, nrow=nrow(centers), ncol=nrow(mesh$triangles))

plot(mesh,pch=".")
points(centers, col="red", pch=16)
for(r in 1:nrow(centers)){
  for(e in 1:nrow(mesh$triangles)){
    point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))
    points(t(point), col="black")
    if( dist(point, t(centers[r,]) )  <=  h   ){
      incidence_matrix[r,e] = 1 
    }
  }
}

#L = sqrt(2)*0.625
point = as.matrix(apply(mesh$nodes[mesh$triangles[e,],],FUN= mean, MARGIN=2))

for( r in 1:nrow(incidence_matrix)){
  apply(mesh$triangles[which(incidence_matrix[r,] == 1),], 
        MARGIN=1, FUN= function(x){
          polygon(mesh$nodes[x,],col="blue")})
  
}
# prova 
tmp_coords = mesh$nodes[mesh$triangles[which(incidence_matrix[1,]==1),],]
tmp_coords2 = mesh$nodes[as.vector(mesh$triangles[which(incidence_matrix[1,]==1),]),]
tmp_coords - tmp_coords2


FEMbasis = create.FEM.basis(mesh)
Mass = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
one = as.matrix(rep(1, times=n_dofs))
D = as.matrix(rep(0, times=nrow(centers))) 
I =  Mass%*% one
for(r in 1:nrow(centers)){
  D[r,] = sum(I[unique( as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),]))])
}

exact  = as.matrix(read.table("../simulation_0_dirichlet/data/exact.txt", header=F))
time_locs = as.matrix(read.table("../simulation_0_dirichlet/data/time_locations.txt"))
n_time_locs = length(time_locs)
int_exact = matrix(0, nrow=nrow(centers), ncol=n_time_locs)
for(r in 1:nrow(centers)){
  for(t in 1:n_time_locs){
    I = Mass %*% exact[,t]
    int_exact[r,t]= 1/D[r,]*sum(I[unique( as.vector(mesh$triangles[which(incidence_matrix[r,] == 1),]))])
  }
}
set.seed(0)
nsim = 30
for(i in 0:(nsim-1)){
  observations = int_exact + rnorm(nrow(int_exact)*ncol(int_exact), sd = 0.05) # 0.04 circa ==  0.05*diff(range(abs(int_exact)))
  datadir = paste0("data/",i, "/")
  if(! dir.exists(datadir) )  dir.create(datadir) 
  
  write.table(observations, paste0(datadir, "obs.txt"),
              col.names = F, row.names = F)
}

write.table(incidence_matrix,file = "data/incidence_matrix.txt",
            row.names = F, col.names = F)

write.table(D,file = "data/area.txt",
            row.names = F, col.names = F)
