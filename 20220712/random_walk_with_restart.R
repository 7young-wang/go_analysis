#Random Walk with Restart
BiocManager::install("RandomWalkRestartMH")
library(RandomWalkRestartMH)
library(igraph)
library(GOSemSim)

data(Pathway_Network) # We load the Pathway Network

## We create a 2-layers Multiplex object
PPI_PATH_Multiplex <- 
  create.multiplex(list(PPI=PPI_Network,PATH=Pathway_Network))
PPI_PATH_Multiplex

AdjMatrix_PPI_PATH <- compute.adjacency.matrix(PPI_PATH_Multiplex)
AdjMatrixNorm_PPI_PATH <- normalize.multiplex.adjacency(AdjMatrix_PPI_PATH)


bp_go <- godata('org.Dr.eg.db', ont="BP", keytype = "SYMBOL") ##download BP GO for all zebrafish genes 
mf_go <- godata('org.Dr.eg.db', ont="MF", keytype = "SYMBOL") ##download MF GO for all zebrafish genes
cc_go <- godata('org.Dr.eg.db', ont="CC", keytype = "SYMBOL") ##download CC GO for all zebrafish genes

gene_name <- read.delim("genes_in_ld_blocks.txt")
gene_name<- gene_name$genes
method <- "Resnik"
comb <- "BMA"

# make a matrix for each GO type 
bp_matrix <- mgeneSim(genes = gene_name, semData = bp_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)
mf_matrix <- mgeneSim(genes = gene_name, semData = mf_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)
cc_matrix <- mgeneSim(genes = gene_name, semData = cc_go, measure = method, drop = NULL, combine = comb, verbose = TRUE)

heatmap(bp_matrix, Rowv = NA, Colv = NA)
library(Matrix)
test_bp <- as(bp_matrix, "dgCMatrix")

bp_normal <- normalize.multiplex.adjacency(test_bp)

#make an igraph object from dataframe

g <- graph.adjacency(bp_matrix, mode="undirected", weighted=TRUE) # For undirected networks
plot(g, vertex.label=NA, vertex.size=3, edge.width=E(g)$width)

SeedGene <- c("htr1b")
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
RWR_Results <- Random.Walk.Restart.Multiplex(adjm_nor_bp, multiplex_bp,"htr1b" , r = 0.3)
TopResults_PPI <-
  create.multiplexNetwork.topResults(RWR_Results,multiplex_bp,
                                     k=5)
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI, vertex.label.color="black",vertex.frame.color="#ffffff",
     vertex.size= 20, edge.curved=.2,
     vertex.color = ifelse(igraph::V(TopResults_PPI)$name == "htr1b","yellow",
                           "#00CCFF"), edge.color="blue",edge.width=0.8)





is.matrix(adjm)
# We display the results
RWR_PPI_Results

PPI_Network

#https://alberto-valdeolivas.github.io/RandomWalkRestartMH/articles/RandomWalkRestartMH.html#random-walk-with-restart-on-a-heterogeneous-network-1


#draw a graph
target_matrix <- target_matrix[genes, genes]
target_matrix <- as.matrix(target_matrix)
selected_genes <- read.csv("BP_Resnik_BMA_final_result.csv")$gene
need <- NULL
all_genes <- sample(read.delim("genes_in_ld_blocks.txt")$genes, size = 10, replace = F)
need <- append(all_genes, selected_genes)
for (i in selected_genes) {
  RWR_Results <- Random.Walk.Restart.Multiplex(adjm_nor_bp, multiplex_bp, i, r = 0.3)
  need<- append(need, append(i, head(RWR_Results$RWRM_Results$NodeNames, 5)))
}
need <- unique(need)
need <- intersect(colnames(store_matrix), need)
store_matrix <- target_matrix
target_matrix <- store_matrix[need, need]
target_matrix <- as.matrix(target_matrix)
g <- graph.adjacency(target_matrix, mode="undirected", weighted=TRUE) # For undirected networks
plot(g, vertex.label=NA, vertex.size=3, edge.width=E(g)$width)
g.copy <- delete.edges(g, which(E(g)$weight == 0.3))
plot(g.copy, vertex.label.color="black",vertex.frame.color="#ffffff",
     vertex.size= 15, edge.curved=.2, edge.color="blue",edge.width=percentile_weight,
     vertex.color = V(g)$color)
V(g)$color <- "yellow"
V(g)[selected_genes]$color<-"red"

log(E(g.copy)$weight, 0.1)
new_weight <- ((E(g.copy)$weight*100000))^2/10/5
new_weight
E(g.copy)$weight
df <- data.frame(weight = E(g.copy)$weight)
df$percentile <- ecdf(df$weight)(df$weight)
df$edge <- percentile_weight
keep_percentile <- df$percentile
percentile_weight <- NULL
for (i in keep_percentile) {
  if (i < 0.25) {
    percentile_weight <- append(percentile_weight, 0.01)
  } else if (i < 0.5) {
    percentile_weight <- append(percentile_weight, 0.1)
  } else if (i < 0.75) {
    percentile_weight <- append(percentile_weight, 0.5)
  } else if (i < 0.95){
    percentile_weight <- append(percentile_weight, 1)
  } else {
    percentile_weight <- append(percentile_weight, 5)
  }
}
percentile_weight <- round(keep_percentile, 1) * 2
V(g)$color


#Color scaling function
c_scale <- colorRamp(c('red','yellow','cyan','blue'))

#Applying the color scale to edge weights.
#rgb method is to convert colors to a character vector.
E(g)$color = apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )

#plot using igraph
plot.igraph(g)

