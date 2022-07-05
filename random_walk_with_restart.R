#Random Walk with Restart
BiocManager::install("RandomWalkRestartMH")
library(RandomWalkRestartMH)
library(igraph)

data(Pathway_Network) # We load the Pathway Network

## We create a 2-layers Multiplex object
PPI_PATH_Multiplex <- 
  create.multiplex(list(PPI=PPI_Network,PATH=Pathway_Network))
PPI_PATH_Multiplex

AdjMatrix_PPI_PATH <- compute.adjacency.matrix(PPI_PATH_Multiplex)
AdjMatrixNorm_PPI_PATH <- normalize.multiplex.adjacency(AdjMatrix_PPI_PATH)


bp_go <- godata('org.Dr.eg.db', ont="BP", keytype = "SYMBOL")
mf_go <- godata('org.Dr.eg.db', ont="MF", keytype = "SYMBOL")
cc_go <- godata('org.Dr.eg.db', ont="CC", keytype = "SYMBOL")

gene_name <- read.delim("genes_in_ld_blocks_20220623.txt")
gene_name<- gene_name$X0
typeof(c(1,2))
method <- "Resnik"
comb <- "BMA"

bp_matrix <- mgeneSim(genes = gene_name, semData = bp_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)
mf_matrix <- mgeneSim(genes = gene_name, semData = mf_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)
cc_matrix <- mgeneSim(genes = gene_name, semData = cc_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)

library(Matrix)
test_bp <- as(bp_matrix, "dgCMatrix")

bp_normal <- normalize.multiplex.adjacency(test_bp)

#make an igraph object from dataframe

g <- graph.adjacency(bp_matrix, mode="undirected", weighted=TRUE) # For undirected networks

SeedGene <- c("nrp1b")
multiplex_bp <- create.multiplex(list(PPI=g))
multiplex_bp
## We launch the algorithm with the default parameters (See details on manual)
data(PPI_Network) # We load the PPI_Network

## We create a Multiplex object composed of 1 layer (It's a Monoplex Network) 
## and we display how it looks like
PPI_MultiplexObject <- create.multiplex(list(PPI=PPI_Network))
PPI_MultiplexObject
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
adjm_bp <- compute.adjacency.matrix(multiplex_bp)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
adjm_nor_bp <- normalize.multiplex.adjacency(adjm_bp)
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                 PPI_MultiplexObject,SeedGene)
bp_normal
SeedGene <- c("PIK3R1")
RWR_Results <- Random.Walk.Restart.Multiplex(adjm_nor_bp, multiplex_bp, SeedGene, r = 0.3)



is.matrix(adjm)
# We display the results
RWR_PPI_Results

PPI_Network

#https://alberto-valdeolivas.github.io/RandomWalkRestartMH/articles/RandomWalkRestartMH.html#random-walk-with-restart-on-a-heterogeneous-network-1





