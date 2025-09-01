library(fdaPDE)

nodes = as.matrix(read.table("nodes.txt", header = F))
cells = as.matrix(read.table("cells.txt", header = F))
boundary_idx = as.matrix(read.table("boundary.txt", header = F))
regions = as.matrix(read.table("regions.txt", header = F))

#  #(1001:1035, 2001:2035) x #cells, see $FREESURFER/FreeSurferColorLUT.txt  

ctx = c(1001:1003, 1005:1035, 2001:2003, 2005:2035)
#corpocallosum -> 0, makes sense (1004 & 2004)
incidence_matrix = matrix(nrow=length(ctx),ncol=nrow(cells))
for(i in 1:length(ctx)){
  incidence_matrix[i,] = as.integer( regions == ctx[i] ) 
}
rowSums(incidence_matrix)
write.csv(incidence_matrix, file = "incidence_matrix.csv")
write.csv(format(nodes, digits=16), file="nodes.csv")

boundary = matrix(0, nrow = nrow(nodes), ncol=1)
boundary[boundary_idx] = 1
write.csv(boundary, file="boundary.csv")
write.csv(cells, file="cells.csv")
