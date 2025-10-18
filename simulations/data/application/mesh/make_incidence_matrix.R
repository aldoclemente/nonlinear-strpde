library(fdaPDE)
rm(list=ls())
nodes = as.matrix(read.table("nodes.txt", header = F))
cells = as.matrix(read.table("cells.txt", header = F))
boundary_idx = as.matrix(read.table("boundary.txt", header = F))
regions = as.matrix(read.table("regions.txt", header = F))

#  #(1001:1035, 2001:2035) x #cells, see $FREESURFER/FreeSurferColorLUT.txt  

ctx = c(1001:1003, 1005:1035, 2001:2003, 2005:2035)
#corpocallosum -> 0, makes sense (1004 & 2004)
incidence_matrix = Matrix(nrow=length(ctx),ncol=nrow(cells), sparse = T)

for(i in 1:length(ctx)){
  reg = which(regions == ctx[i])
  for(j in reg){
    incidence_matrix[i, j] = 1 #as.integer( regions == ctx[i] ) 
  }
}

writeMM(incidence_matrix, "incidence_mat.mtx")

rowSums(incidence_matrix)
#write.csv(incidence_matrix, file = "incidence_matrix.csv")
write.csv(format(nodes, digits=16), file="nodes.csv")

boundary = matrix(0, nrow = nrow(nodes), ncol=1)
boundary[boundary_idx] = 1
write.csv(boundary, file="boundary.csv")
write.csv(cells, file="cells.csv")

# plot (?)
incidence_matrix = as.matrix(read.csv("incidence_matrix.csv")[,-1])
library(Matrix)
# BUG writeMM ?!
writeMM(Matrix(incidence_matrix+1, sparse=T), "incidence_matrix.mtx")


mtx_tmp = matrix(c(1,0,0,1), nrow=2,ncol=2, byrow = T)
#storage.mode(mtx_tmp) <- "integer"
writeMM(Matrix(mtx_tmp, sparse=T), "prova.mtx")


writeMM(Matrix(as.numeric(incidence_matrix, sparse=T), "incidence_matrix.mtx")
