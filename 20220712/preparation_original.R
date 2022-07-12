library(GOSemSim)
library(RandomWalkRestartMH)
library(igraph)
args<-commandArgs(TRUE)

#for i in {1..100}; do RScript Gibbs.R BP ./BP $i; done

preparation <- function(args) {
  #sample command BP Resnik BMA genes to make matrix(e.g. bp_final.csv, i.e. genes you wanna work with)
  ontology <- args[1]
  ic_method <- args[2]
  comb_method <- args[3]
  origin <- args[4]
  
  #when only one gene is found in one ld block in the resulting list(list of genes with annotation), exclude it
  exclude_single_gene <- function(ld_information) {
    count <- 0
    unique_snp <- unique(ld_information$peak_snp)
    for (i in unique_snp) {
      if (length(ld_information[ld_information$peak_snp == i,]$gene) == 1) {
        count <- count + 1
        ld_information <- subset(ld_information, peak_snp != i)
      }
    }
    cat (paste("excluded", count, "LD blocks \n"))
    ld_information
  }
  
  make_rwr_matrix <- function(original_matrix) {
    gt <- graph.adjacency(original_matrix, mode="undirected", weighted=TRUE)
    multiplex <- create.multiplex(list(PPI=gt))
    adjm <- compute.adjacency.matrix(multiplex)
    adjm_nor <- normalize.multiplex.adjacency(adjm)
    result <- rbind(Random.Walk.Restart.Multiplex(adjm_nor, multiplex, colnames(original_matrix)[1], r = 0.3)$RWRM_Results, c(colnames(original_matrix)[1],0.3))
    result[,-1] <- as.numeric(result[,-1])
    colnames(result) <- c("NodeNames", colnames(original_matrix)[1])
    for (i in 2:length(colnames(original_matrix))) {
      SeedGene <- colnames(original_matrix)[i]
      new_column <- rbind(Random.Walk.Restart.Multiplex(adjm_nor, multiplex, SeedGene, r = 0.3)$RWRM_Results, c(SeedGene, (0.3)))
      new_column[,-1] <- as.numeric(new_column[,-1])
      result <- merge(result, new_column, by = "NodeNames")
      colnames(result) <- c(head(colnames(result), -1), SeedGene)
    }
    return (result)
  }
  
  #generate the semetic matrix
  #ld fine must be in the current directory
  generate_matrix <- function(ldfile, ont1, method, comb) {
    if(file.exists(paste(ont1, "_GO.RData", sep = ""))) {
      load(paste(ont1, "_GO.RData", sep = ""))
    } else {
      hsgo <- godata("org.Dr.eg.db", keytype = "SYMBOL", ont = ont1)
      save(hsgo, file = paste(ont1, "_GO.RData", sep = ""))
    }
    if (ldfile == "all") {
      gene_name <- hsgo@keys
    } else {
      gene_name <- read.delim(ldfile)$genes #genes is the column name for txt file
    }
    #gene_name <- read.csv(ldfile)$gene
    target_matrix <- mgeneSim(genes = gene_name, semData = hsgo, measure = method, drop = NULL, combine = comb, verbose = TRUE)
    target_matrix <- make_rwr_matrix(target_matrix)
    rownames(target_matrix) <- target_matrix[,1]
    target_matrix <- target_matrix[,-1]
    target_matrix <- target_matrix[ , order(names(target_matrix))]
    if(!file.exists(ont1)) dir.create(ont1,showWarnings=T,recursive=T)
    if(!file.exists("./supporting_files")) dir.create("./supporting_files",showWarnings=T,recursive=T)
    #filename <- paste("./supporting_files/target_matrix", ont1, sep = "_")
    filename <<- paste(ont1, method, comb, "matrix", sep = "_")
    save(target_matrix, file = paste("./supporting_files/", filename, ".RData", sep = "")) 
    cat("generate_matrix was runned \n")
  }
  
  generate_matrix(origin, ontology, ic_method, comb_method)
  cat(paste("This is running", args[2], "method with", args[3], "combination method \n"))
  

  convert_all_to_target <- function(big, target_list) {
    load(big)
    genes <- read.delim(target_list)$genes
    genes <- intersect(genes, colnames(target_matrix))
    target_matrix <- target_matrix[genes, genes]
    if(!file.exists("./supporting_files")) dir.create("./supporting_files",showWarnings=T,recursive=T)
    save(target_matrix, file = paste("./supporting_files/", filename, ".RData", sep = ""))
  }
  filename <- paste(ontology, ic_method, comb_method, "matrix", sep = "_")
  #convert_all_to_target(paste(filename, "all.RData", sep = ""), origin)
  
  #generate a list of genes that appear in the matrix
  generate_gene_list <- function(matrixFile, original, ont) {
    #ont is BP, CC or MF
    load(matrixFile)
    gene_list <- as.data.frame(colnames(target_matrix))
    colnames(gene_list) <- c("gene")
    origin1 <- read.delim(original)
    gene_list <- merge(gene_list, origin1, by.x = "gene", by.y = "genes")
    gene_list <- gene_list[ , !names(gene_list) %in% c("trait")]
    gene_list <- exclude_single_gene(gene_list)
    if(!file.exists(ont)) dir.create(ont,showWarnings=T,recursive=T)
    write.table(gene_list, paste("./", ont, "/rank_table_", ont, sep = "") ,quote=F,row.names=F,sep="\t")
  }
  
  if (origin != "all") {
    generate_gene_list(paste("./supporting_files/", filename, ".RData", sep = ""), origin, ontology)
  }
  
  convert_to_csv <- function(result, name) {
    temp <- read.delim(result)
    write.csv(temp, paste(name, ".csv", sep = ""))
  }
  
  
  
  
}

