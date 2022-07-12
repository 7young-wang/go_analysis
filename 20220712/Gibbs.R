 
 args<-commandArgs(TRUE)

 gibbs<-function(args)
 {
  if(length(args)==0)
   stop("No arguments!\nFor help, please type \"Rscript Gibbs.R -h\" or \"Rscript Gibbs.R --help\".\n")

  if(length(args)<2 & args[1]!="-h" & args[1]!="-help")
   stop("No enough arguments!\nFor help, please type \"Rscript Gibbs.R -h\" or \"Rscript Gibbs.R --help\".\n")

  if(args[1]=="-h" | args[1]=="--help")
  {
   cat("\nUsage: Rscript Gibbs.R SNP_file flank res_path res_pref\nPlease desinage at least the first two parameters!\n") 
   cat("\nExplanation of parameters:\n")
   cat("  --SNP_file   File name of GWAS SNPs. The file must include three columns: SNP, Chr, and Pos_hg19.\n")
   cat("  --flank      An integer indicating the flanking region (basepair) of a SNP when identifying candidate genes.\n")
   cat("               1000000 is our recommended number.\n")
   cat("  --res_path   Path of result files. Optional; \"./iRIGS_result/\" by default if not designated.\n")
   cat("  --res_pref   Prefix of the result file. Optional; SNP_file name by default if not designated.\n")
   cat("\nFor example Rscript Gibbs.R ./SNP_file/SCZ_108_loci 1000000 ./iRIGS_result_108_loci/ SCZ\n\n")
  } else
  {
   #source("extract_candidate_genes.R")
   #source("preparation.R")
   source("sampling_procedure_new.R")
   source("preparation_new.R")
   
   #source("extra_evi.R") 
  name <- paste(args[1], args[2], args[3], sep = "_")
   #args<-c("./SNP_file/SCZ_20_loci_test",1000000,"./iRIGS_result_108_loci/","SCZ")
  #sample command: BP Resnik BMA genes_in_ld_blocks.txt 100
  preparation(args)
  for (i in 1:args[5]) {
    output<-sampling_procedure(args[1], args[2], args[3], args[4], i)
  }
  
  calculate_ranking <- function(ont) {
    rank_table <- read.delim(paste("./", ont, "/", "rank_table_", ont, sep = ""))
    temp <- rank_table[, !names(rank_table) %in% c("gene", "peak_snp")]
    rank_table$average_rank <- rowMeans(temp)
    write.table(rank_table,paste("./", ont, "/", "rank_table_", ont, sep = ""),quote=F,row.names=F,sep="\t")
    
    return (rank_table)
  }
  
  rank_table <- calculate_ranking(args[1])
  
  select_genes <- function(table, ont) {
    final_result <- NULL
    split_result <- split(table, table$peak_snp)
    find_best <- function(one_ld) {
      if (is.null(final_result)) {
        final_result <<- one_ld[order(-one_ld$average_rank),][, c(1, 2)][1,]
      } else {
        final_result[nrow(final_result) + 1,] <<- one_ld[order(-one_ld$average_rank),][, c(1, 2)][1,]
      }
    }
    lapply(split_result, find_best)
    write.csv(final_result, paste(ont, "final_result.csv", sep = "_"))
  }
  
  select_genes(rank_table, name)

   cat("All analysis finished!\n")
  }
 }

 gibbs(args) 
 q("no")


