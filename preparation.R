library(GOSemSim)
args<-commandArgs(TRUE)
#sample command: just do BP is enough; we only use args1

#for i in {1..100}; do RScript Gibbs.R BP ./BP $i; done

preparation <- function(args) {
  ontology <- args[1]
  #sample command: CC ./CC $i (which ontology to use, where the files will be stored, the name of each file)
  
  #function used to normalize the semsim matrix
  normalize <- function(x) {
    sum_column <- sum(x)
    return (x/sum_column)
  }
  
  exclude_single_gene <- function(ld_information) {
    count <- 0
    unique_snp <- unique(ld_information$peak_snp)
    for (i in unique_snp) {
      if (length(ld_information[ld_information$peak_snp == i,]$gene) == 1) {
        count <- count + 1
        ld_information <- subset(ld_information, peak_snp != i)
      }
    }
    cat (paste("excluded", count, "LD blocks"))
    ld_information
  }
  
  #generate the semetic matrix
  #ld fine must be in the current directory and is named "genes_in_ld_blocks.txt"
  generate_matrix <- function(ldfile, ont1) {
    hsgo <- godata("org.Dr.eg.db", keytype = "SYMBOL", ont = ont1) #preapre the ontology map; takes time
    gene_name <- read.delim(ldfile)$genes
    target_matrix <- mgeneSim(genes = gene_name, semData = hsgo, measure = "Resnik", drop = NULL, combine = "BMA", verbose = TRUE)
    target_matrix <- apply(target_matrix, 2, normalize)
    if(!file.exists("./supporting_files")) dir.create("./supporting_files",showWarnings=T,recursive=T)
    filename <- paste("./supporting_files/target_matrix", ont1, sep = "_")
    save(target_matrix, file = paste(filename, ".RData", sep = "")) 
    cat("generate matrix was runned")
  }
  generate_matrix("genes_in_ld_blocks.txt", ontology)
  
  #generate a list of genes that appear in the matrix
  generate_gene_list <- function(matrixFile, original, ont) {
    #ont is BP, CC or MF
    load(matrixFile)
    gene_list <- as.data.frame(colnames(target_matrix))
    colnames(gene_list) <- c("gene")
    origin <- read.delim(original)
    gene_list <- merge(gene_list, origin, by.x = "gene", by.y = "genes")
    gene_list <- gene_list[ , !names(gene_list) %in% c("trait")]
    gene_list <- exclude_single_gene(gene_list)
    if(!file.exists(ont)) dir.create(ont,showWarnings=T,recursive=T)
    write.table(gene_list, paste("./", ont, "/rank_table_", ont, sep = "") ,quote=F,row.names=F,sep="\t")
  }
  
  generate_gene_list(paste("./supporting_files/target_matrix_", ontology, ".RData", sep = ""), "genes_in_ld_blocks.txt", ontology)
  
  convert_to_csv <- function(result, name) {
    temp <- read.delim(result)
    write.csv(temp, paste(name, ".csv", sep = ""))
  }
}

preparation(args)