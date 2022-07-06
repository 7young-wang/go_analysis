library(clusterProfiler)
bp_gene_list <- read.csv("bp_final.csv")$gene
temp_list <- read.csv("temp.csv")$gene
cc_gene_list <- read.csv("CC_final.csv")$gene
mf_gene_list <- read.csv("MF_final.csv")$gene
all_gene_list <- read.delim("./BP/rank_table_BP copy")
random_genes <- all_gene_list[sample(nrow(all_gene_list), 12), ]$gene
random_genes
enrichment_test <- function(gene, name) {
  ego <- enrichGO(gene = temp_list, OrgDb = "org.Dr.eg.db", ont = "ALL", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "SYMBOL")
  write.csv(ego, file = paste(name, ".csv", sep = ""))
}
temp_list
bp_gene_list
temp_list <- unname(temp_list)
enrichment_test(random_genes, "random_test")
enrichment_test(cc_gene_list, "cc_enrichment")
enrichment_test(mf_gene_list, "mf_enrichment")
enrichment_test(one_time, "one_time_bp_enrich")
enrichment_test(new_random, "random_enrich")
enrichment_test(temp_list, "temp_enrich")

ego <- enrichGO(gene = bp_gene_list, OrgDb = "org.Dr.eg.db", ont = "ALL", 
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "SYMBOL")
head(ego)
write.csv(ego, file = "no_pcutoff.csv")

write.csv(read.delim("SCZ_risk_genes_predicted_by_iRIGS"), file = "sample_test.csv")
sample_test <- read.csv(file = "sample_test.csv")
sample_test <- sample_test$gene
ego <- enrichGO(gene = small_sample, OrgDb = "org.Hs.eg.db", ont = "ALL", 
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "SYMBOL")
small_sample <- sample(sample_test, 12)

ego <- enrichGO(gene = temp_list, OrgDb = "org.Dr.eg.db", ont = "ALL", 
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 0.2, keyType = "SYMBOL")

one_time <- read.delim("./BP/test")
head(one_time)
write.csv(one_time, "one_time.csv")
one_time <- read.csv(file = "one_time.csv")
one_time <- one_time$gene

ggo <- groupGO(gene     = new_random,
               OrgDb    = "org.Dr.eg.db",
               ont      = "BP",
               level    = 3, keyType = "SYMBOL")
ggo_result <- ggo@result
write.csv(ggo_result, "random_just_group.csv")


split_all_gene_list <- split(all_gene_list, all_gene_list$peak_snp)
random_genes <- lapply(split_all_gene_list, function(x) x[sample(seq_len(nrow(x)), 1),])
random_genes <- unlist(random_genes)
names(random_genes)
new_random <- NULL
for (i in 1:length(random_genes)) {
  new_random <<- append(new_random, random_genes[[i]]$gene)
}
new_random
