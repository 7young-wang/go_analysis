
sampling_procedure<-function(ont, method, comb, origin, number_cycle)
{
  burnin_round<-3000          ### maximum sampling round allowed at burn in step
  after_burnin_round<-3000    ### maximum sampling round allowed at post burn in step
  exclude_samegene<-T
  
  ###======================= specify different network =============================###
  
  #transp <- paste("./supporting_files/", paste(ont, method, comb, "matrix", sep = "_"), ".RData", sep = "")
  transp <- paste("./supporting_files/", "BP_RWR", ".RData", sep = "")
  
  cat("Loading propagation probabilities...\n")
  load(transp); 
  target_matrix <- df
  #nodes<-colnames(target_matrix) #nodes are all the human genes in the matrix
  
  nodes<-colnames(target_matrix)
  
  #gene<-extract_candidate_genes(args) #get the target genes from SNP data; this is what I already have right? I don't actually need to use the extract_candidate_genes function
  #gene<-gene[(!is.na(gene$official_name)),] #get the official name of the genes, and exclude all the NA; in our case, genes
  #gene<-gene[is.element(gene$official_name,nodes),] #select the genes that are in our list and also in the matrix (we have their data)
  gene <- read.delim("genes_in_ld_blocks.txt") #convert the txt into a list object with columns genes, peak_snp, trait
  if (ont == "BP") {
    rank_table <- read.delim("./BP/rank_table_BP")
  } else if (ont == "MF") {
    rank_table <- read.delim("./MF/rank_table_MF")
    #target_matrix <- target_matrix_MF
  } else if (ont == "CC") {
    rank_table <- read.delim("./CC/rank_table_CC")
    #target_matrix <- target_matrix_CC
  } else {
    rank_table <- read.delim("./RWR_BP/0705RWR_table_RWR_BP")
  }
  gene <- gene[is.element(gene$genes,nodes),]
  region<-split(gene$genes,gene$peak_snp) #divide the genes by their LD block; we need to select only one gene from each block; group the genes by snp
  cat(paste(length(gene$genes),"genes from",length(region),"loci were found with propagation probability...\n"))
  
  #pro_p is matrix in this case; rename it so it's easier; these steps are just further modifying the matrix, we don't need these steps
  pro_p<-target_matrix[,(is.element(nodes,unique(gene$genes)))]
  #extra_weight<-extra_evi(gene) #??don't quite understand
  
  ####================= end of loading data ===============####
  
  ####================= burn in step ======================####
  
  #t0<-proc.time() 
  
  thres<-0.001; pickup<-0;
  num_region<-length(region); circle<-1; chosen<-NULL
  remaining<-unlist(lapply(region,function(x) sample(x,1)))
  num0<-rep(0,sum(unlist(lapply(region,length))))
  
  dif<-thres+1; dif_record<-NULL 
  #while(dif>thres && circle<100)
  while(dif>thres && circle<(burnin_round+1))
  {
    pickup<-pickup%%num_region+1  
    if(pickup==1)
      if(!is.null(chosen))
      {
        num1<-NULL
        for(j in 1:length(region))
          num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x)))) 
        num1<-num1+num0
        if(circle>1)
        {
          freq0<-num0/(num_region*(circle-1))
          freq1<-num1/(num_region*circle)
          dif<-(sum((freq0-freq1)^2))^0.5
          if( circle%%50==0 )
          {
            cat("Burnin sampling, sampling circle:",circle,"\n")
          }
          #dif_record<-c(dif_record,dif)
        }
        num0<-num1; chosen<-NULL; circle<-circle+1
      }
    
    pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
    pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]  
    
    if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
    {                                                                            ### conditional genes, exclude the same genes 
      pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
      if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
    }
    
    if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene  
    { pickup_p<-apply(pickup_p,1,sum);} #each row of pickup gene is a gene and we are calculating the sum here?  #pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p }
    
    if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0 
    
    remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
    chosen<-rbind(chosen,remaining)  
  }
  
  #proc.time()-t0 
  
  ###===================== end of burn in step ===============================###
  
  ###======================= post-burnin step ===================================###
  #t0<-proc.time()
  
  pickup<-0; num_region<-length(region); circle<-1; chosen<-NULL
  num0<-rep(0,sum(unlist(lapply(region,length))))
  
  joi_dis<-matrix(0,nrow=nrow(gene),ncol=nrow(gene))
  temp<-NULL;
  for(j in 1:length(region))
    temp<-c(temp,paste(names(region[j]),region[[j]],sep="_"))
  colnames(joi_dis)<-temp; rownames(joi_dis)<-temp
  
  thres<-0.001; dif<-thres+1
  while(dif>thres && circle<(after_burnin_round+1) )
  {
    pickup<-pickup%%num_region+1
    if(pickup==1)
      if(!is.null(chosen))  
      {
        ###================================= calculate frequency =========================###
        num1<-NULL
        for(j in 1:length(region))
          num1<-c(num1,unlist(lapply(region[[j]],function(x) sum(chosen[,j]==x))))
        num1<-num1+num0
        if(circle>1)
        {
          freq0<-num0/(num_region*(circle-1))
          freq1<-num1/(num_region*circle)
          dif<-(sum((freq0-freq1)^2))^0.5
          if( circle%%50==0 )
          {
            cat("Post-burnin sampling, sampling circle:",circle,"\n")
          }
          #dif_record<-c(dif_record,dif)
        }
        num0<-num1; circle<-circle+1; chosen<-NULL
        ###============================= end of calculating frequency =======================###
      }
    
    pickup_p<-pro_p[,is.element(colnames(pro_p),remaining[-pickup])]
    pickup_p<-pickup_p[is.element(rownames(pickup_p),unlist(region[pickup])),]
    
    if(exclude_samegene && !is.null(dim(pickup_p)))                              ### if there is overlap between candidate genes and
    {                                                                            ### conditional genes, exclude the same genes
      pickup_p<-pickup_p[,!is.element(colnames(pickup_p),rownames(pickup_p))]
      if( !is.null(dim(pickup_p)) && ncol(pickup_p)==0 )  stop("Error: no conditional genes!\n")
    }
    
    if(is.null(dim(pickup_p))) { pickup_p<-1; names(pickup_p)<-region[[pickup]] } else    ### when there is only one candiate gene
    { pickup_p<-apply(pickup_p,1,sum);} #pickup_p<-extra_weight[match(names(pickup_p),extra_weight$gene),]$extra_weight*pickup_p }
    
    if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i]<-1/length(pickup_p)     ### avoid probabilit=0
    
    remaining[pickup]<-sample(names(pickup_p),1,replace=T,prob=pickup_p)
    chosen<-rbind(chosen,remaining)
    
    ###=============================== calculating joint distribution  =====================###
    index_col<-match(paste(names(remaining[-pickup]),remaining[-pickup],sep="_"),colnames(joi_dis))  
    index_row<-match(paste(names(remaining[pickup]),remaining[pickup],sep="_"),colnames(joi_dis))  
    joi_dis[index_row,index_col]<-joi_dis[index_row,index_col]+1 
  }
  
  #proc.time()-t0
  
  ###=====================  end of post-burnin step  =======================================###
  
  ###===================== summarize and record the results ================================###
  freq<-cbind(unlist(region),freq1)
  
  region_indicator<-NULL
  gene_num<-as.numeric(lapply(region,length))
  for(i in 1:length(gene_num))
    region_indicator<-c(region_indicator,rep(names(region[i]),gene_num[i]))
  
  freq<-cbind(freq,region_indicator)
  colnames(freq)<-c("gene","post_prob","region")
  freq<-as.data.frame(freq,stringsAsFactors=F)
  freq[,2]<-as.numeric(freq[,2])
  
  output<-NULL                                      #### sort according to posterior probability
  ranking <- NULL
  for(i in unique(freq$region))
  {
    temp<-freq[freq$region==i,]
    output<-rbind(output,temp[order(temp$post_prob,decreasing=T),])
    ranking <- c(ranking, seq(1, nrow(temp)))
  }
  output <- cbind(output, ranking)
  colnames(output) <- c("gene","post_prob","region", "rank")
  rank_table <- merge(rank_table, output, by.x = c("gene", "peak_snp"), by.y = c("gene", "region"))
  rank_table <- rank_table[ , !names(rank_table) %in% c("region", "rank")]
  #rankings <- as.data.frame(output$gene)
  #colnames(rankings) <- c("gene")
  freq <- output
  
  
  cat("Recording the results!\n")
  
  res_path<-paste("./", ont, "/", sep = "")
  if(substr(res_path,nchar(res_path),nchar(res_path))!="/")  res_path<-paste(res_path,"/",sep="")
  if(!file.exists(res_path)) dir.create(res_path,showWarnings=T,recursive=T)
  
  if(is.na(number_cycle)) res_file<-strsplit(ont,"/")[[1]][length(strsplit(ont,"/")[[1]])]  else res_file<-number_cycle
  res_file<-paste(res_path,res_file,"",sep="")
  if (ont == "BP") {
    rank_file <- paste("./BP/","rank_table_BP","",sep="")
  } else if (ont == "MF") {
    rank_file <- paste("./MF/","rank_table_MF","",sep="")
  } else if (ont == "CC") {
    rank_file <- paste("./CC/","rank_table_CC","",sep="")
  } else {
    rank_file <- "./RWR_BP/0705RWR_table_RWR_BP"
  }
  write.table(freq,res_file,quote=F,row.names=F,sep="\t")
  write.table(rank_table,rank_file,quote=F,row.names=F,sep="\t")
  output
}








