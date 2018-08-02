phylo_genome <- function(sample_list){
  library(stringr)
  library(seqinr)
  sample_tab <- read.csv(sample_list,row.names = 1)
  spe <- sample_tab[,2]
  ids <- sample_tab[,1]
  sample_num <- length(ids)
  file1 <- as.character(sample_tab[,4])
  dir.create("genome_files_named_with_ids")
  dir.create("gene_prediction")
  dir.create("gene_prediction/protein_seqs")
  dir.create("gene_prediction/nucleotide_seqs")
  dir.create("gene_prediction/GFF")
  dir.create("gene_clusters")
  
  for(i in seq(length(file1))){
    if((file1[i] %in% dir())==FALSE){
      stop(str_c(i,"......",file1[i],"  does not exsit, please check sample list and sequence files!"))
    }
    file.copy(file1[i],"genome_files_named_with_ids/")
    file.rename(str_c("genome_files_named_with_ids/",file1[i]),str_c("genome_files_named_with_ids/",ids[i],".fasta"))
    seqs1 <- read.fasta(str_c("genome_files_named_with_ids/",ids[i],".fasta"))
    name1 <- as.list(paste(ids[i],seq(length(names(seqs1))),sep = "_"))
    write.fasta(sequences = seqs1,names = name1,file.out = str_c("genome_files_named_with_ids/",ids[i],".fasta"))
    print(str_c("Gene prediction[ ",file1[i], " ]...................",i,"/",length(file1)))
    system(str_c("prodigal -a gene_prediction/protein_seqs/aa_orf_",ids[i],".fasta -d gene_prediction/nucleotide_seqs/",ids[i],"_nucleotide.fasta -o gene_prediction/GFF/",ids[i],".gff -i genome_files_named_with_ids/", ids[i],".fasta"))
  }
  print("================================================================================================")
  print("Gene prediction=============================================================================DONE")
  print("================================================================================================")
  print("========================Begin clustering with CD-HIT algorithm==================================")
  print("This may take a little longer................................................................")
  print("Go and have a cup of coffee.............................................................................")
  system("cat gene_prediction/protein_seqs/*.fasta > gene_clusters/aa.fasta")
  system("cdhit -i gene_clusters/aa.fasta -o gene_clusters/cluster_results -c 0.5 -M 0 -T 8 -s 0.5 -n 2")
  res1 <- readLines("gene_clusters/cluster_results.clstr")
  seqall <- read.fasta("gene_clusters/aa.fasta")
  line1 <- grep("^>Cluster",res1)
  cluster1 <- NULL
  for(i in seq(from=1,to=length(line1)-1,by=1)){
    if((((line1[i+1]-1)-line1[i]))>(0.8*sample_num)){
      res2 <- res1[line1[i]:((line1[i+1]-1))]
      cluster_name <- str_replace(str_sub(res2[1],2,-1),pattern = " ","_")
      res3 <- res2[-1]
      res4 <- NULL
      for(j in seq(res3)){
        seq_name <- str_sub(str_split(res3[j],pattern = " ")[[1]][2],2,-4)
        ids1 <- str_split(seq_name,"_")[[1]][1]
        res4 <- rbind(res4,c(cluster_name,ids1,seq_name))
      }
      res5 <- res4[duplicated(res4[,2])==FALSE,]
    if(ncol(as.data.frame(res5))>1){
      if(length(res5[,2]) == sample_num){
        dir.create(str_c("gene_clusters/",cluster_name))
        write.fasta(sequences = seqall[res5[,3]],names = as.list(str_c(res5[,3],"_",res5[,1])),file.out = str_c("gene_clusters/",cluster_name,"/",cluster_name,".fasta"))
        system(str_c("muscle -in gene_clusters/",cluster_name,"/",cluster_name,".fasta -out gene_clusters/",cluster_name,"_msa.fasta"))
        system(str_c("qiime split_fasta_on_sample_ids.py -i gene_clusters/",cluster_name,"/",cluster_name,"_msa.fasta -o gene_clusters/",cluster_name,"/oneseq"))
        cluster1 <- c(cluster1,cluster_name)
      }else if((nrow(res5)>=(0.8*sample_num))&&(nrow(res5)!=sample_num)){
        dir.create(str_c("gene_clusters/",cluster_name))
        write.fasta(sequences = seqall[res5[,3]],names = as.list(str_c(res5[,3],"_",res5[,1])),file.out = str_c("gene_clusters/",cluster_name,"/",cluster_name,".fasta"))
        system(str_c("muscle -in gene_clusters/",cluster_name,"/",cluster_name,".fasta -out gene_clusters/",cluster_name,"/",cluster_name,"_msa.fasta"))
        system(str_c("qiime split_fasta_on_sample_ids.py -i gene_clusters/",cluster_name,"/",cluster_name,"_msa.fasta -o gene_clusters/",cluster_name,"/oneseq"))
        count_align <- nchar(read.alignment(str_c("gene_clusters/",cluster_name,"/",cluster_name,"_msa.fasta"),format = "fasta")$seq[[1]])
        no_seq <- ids[(ids %in% res5[,2])==FALSE]
        cluster1 <- c(cluster1,cluster_name)
        for(k in seq(no_seq)){
          seqs2 <- paste(rep("_",count_align),sep = "",collapse = "")
          write.fasta(sequences = seqs2,names = no_seq[k],file.out = str_c("gene_clusters/",cluster_name,"/oneseq/",no_seq[k],".fasta"))
        }
      }else{
        next()
      }
    }else{
       next()
     }
    }
  }
  dir.create("Results")
  dir.create("Results/sequences")
  for(i in seq(ids)){
    seqs4 <- NULL
    for(j in seq(cluster1)){
      file2 <- dir(str_c("gene_clusters/",cluster1[j],"/oneseq/"))
      file3 <- file2[grep(str_c("^",ids[i]),file2)]
      seqs3 <- readLines(str_c("gene_clusters/",cluster1[j],"/oneseq/",file3))[2]
      seqs4 <- paste(seqs4,seqs3,sep = "",collapse = "")
    }
    write.fasta(sequences = seqs4,names = ids[i],file.out = str_c("Results/sequences/",ids[i],".fasta"))
  }
  system(str_c("cat Results/sequences/*.fasta > Results/final.fasta"))
  print("Genome phyogenetic tree building..............................................")
  system(str_c("Gblocks Results/final.fasta -t=p -b3=10 -b4=5 -b5=H"))
  system(str_c("fasttree -wag Results/final.fasta-gb > Results/genome_tree.tre"))
  print("==========================All Done===================================")
}
