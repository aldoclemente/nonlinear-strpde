args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop(" .mesh files must be supplied (input file).n", call.=FALSE)
} 

library(fdaPDE)
source("~/Desktop/readMesh/read.mesh.R")
datadir = "data/mesh/"
filename <- paste0(datadir, args[1])
domain <- read.mesh(filename)
writeMM(Matrix(domain$nodes, sparse = T), paste0(datadir,"nodes.mtx"))
writeMM(Matrix(domain$elements-1, sparse = T), paste0(datadir,"elements.mtx"))
writeMM(Matrix(domain$nodesmarkers, sparse = T), paste0(datadir,"boundary.mtx"))
