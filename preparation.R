library(GOSemSim)
args<-commandArgs(TRUE)

#for i in {1..100}; do RScript Gibbs.R BP ./BP $i; done

preparation <- function(args) {
  #sample command BP Resnik BMA filename(for example lin_have_a_try) genes to make matrix(e.g. bp_final.csv, i.e. genes you wanna work with)
  ontology <- args[1]
  ic_method <- args[2]
  comb_method <- args[3]
  origin <- args[4]
  
  #function used to normalize the semsim matrix
  normalize <- function(x) {
    sum_column <- sum(x)
    return (x/sum_column)
  }
  
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
  
  #generate the semetic matrix
  #ld fine must be in the current directory
  generate_matrix <- function(ldfile, ont1, method, comb) {
    load("hsgo.RData")
    gene_name <- read.delim(ldfile)$genes #genes is the column name for txt file
    #gene_name <- read.csv(ldfile)$gene
    target_matrix <- mgeneSim(genes = gene_name, semData = hsgo, measure = method, drop = NULL, combine = comb, verbose = TRUE)
    target_matrix <- apply(target_matrix, 2, normalize)
    if(!file.exists(ont1)) dir.create(ont1,showWarnings=T,recursive=T)
    if(!file.exists("./supporting_files")) dir.create("./supporting_files",showWarnings=T,recursive=T)
    #filename <- paste("./supporting_files/target_matrix", ont1, sep = "_")
    filename <<- paste(method, comb, "matrix", sep = "_")
    save(target_matrix, file = paste("./supporting_files/", filename, ".RData", sep = "")) 
    cat("generate_matrix was runned \n")
  }
  
  generate_matrix(origin, ontology, ic_method, comb_method)
  cat(paste("This is running", args[2], "method with", args[3], "combination method \n"))
  
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
    
    write.table(gene_list, paste("./", ont, "/rank_table_", ont, sep = "") ,quote=F,row.names=F,sep="\t")
  }
  
  generate_gene_list(paste("./supporting_files/", filename, ".RData", sep = ""), origin, ontology)
  
  convert_to_csv <- function(result, name) {
    temp <- read.delim(result)
    write.csv(temp, paste(name, ".csv", sep = ""))
  }
}

preparation(args)