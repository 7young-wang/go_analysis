#data source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8183574/

tad_info <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)
split_tad <- split(tad_info, tad_info$Chr)

gene_info <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)
gene_info[1,]$chr

if_included <- function(df) {
  for (i in 1:nrow(gene_info)) {
    chr_num <- gene_info[i,]$chr
    if (df$Chr[1] == chr_num) {
      for (j in 1:nrow(df)) {
        #cat(paste("Comparing", gene_info[i,]$gene, "at", chr_num, gene_info[i,]$start, gene_info[i,]$end, "with", df$Chr[1], df[j,]$Start, df[j,]$End, "\n"))
        if (isTRUE(((df[j,]$Start <= gene_info[i,]$start) & (df[j,]$End >= gene_info[i,]$end)))) {
          cat(paste(gene_info[i,]$gene, "is found in TAD at", chr_num, df[j,]$Start, df[j,]$End, "\n"))
          return (1)
        }
      }
    }
  }
}

lapply(split_tad, if_included)


random_info <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)

#change the names a little bit to fit with the if_exist function
gene_info_temp <- gene_info
gene_info <- random_info
