phylo_mlsa <- function(seq_dir,genes){
  library(stringr)
  library(seqinr)
  if(str_sub(seq_dir,-1,-1)!="/"){
    seq_dir <- str_c(seq_dir,"/")
  }
  # if("MLSA_Genes" %in% dir()){
  #   read
  # }else{
  #   print("###########No MLSA_Genes library###########")
  #   print("Program will search")
  #   dir.create("MLSA_Genes")
  #   str_c("wget -O  http://www.uniprot.org/uniprot/",P0A8V2,".fasta")
  # }
  seq_list <- dir(seq_dir)
  seq_name <- str_sub(seq_list,1,-7)
  sample_list <- str_sub(seq_list,8,-7)
  
  dir.create("phmmer_results_of_each_MLSA_genes",recursive = TRUE)
  dir.create("seqs_of_each_MLSA_genes",recursive = TRUE)
  dir.create("results",recursive = TRUE)
  dir.create("results/sequences",recursive = TRUE)
  dir.create("results/final_seq",recursive = TRUE)
  for(i in seq(genes)){
    for(j in seq(seq_list)){
      system(str_c("phmmer MLSA_Genes/",genes[i],".fasta ",seq_dir,seq_list[j]," > phmmer_results_of_each_MLSA_genes/",genes[i],"_",seq_name[j],"_results.txt"))
      print(str_c("analysing ",genes[i],"........................",i,"/",length(genes)))
    }
  }
  no_hit <- NULL
  for(i in seq(sample_list)){
    seqs <- read.fasta(str_c(seq_dir,seq_list[i]))
    for(j in seq(genes)){
      res <- readLines(str_c("phmmer_results_of_each_MLSA_genes/",genes[j],"_",seq_name[i],"_results.txt"))
      name <- str_split(res[grep(">>",res)][1]," ")[[1]][2]
      if(is.na(name)==TRUE){
        no_hit <- rbind(no_hit,c(sample_list[i],genes[j]))
      }else{
        res1 <- res[(grep(">>",res)[1]+3):(grep("==",res)[1]-3)]
        seqs2 <- NULL
        if(length(grep("!",res1))>0){
          for(k in seq(length(grep("!",res1)))){
            res2 <- res1[grep("!",res1)][k]
            line1 <- str_split(res2," ")[[1]][str_split(res2," ")[[1]]!=""]
            site <- as.numeric(line1[line1!=""][10:11])
            seqs1 <- paste(as.character(getFrag(seqs[[name]],site[1],site[2])),collapse = "")
            seqs2 <- paste(seqs2,seqs1)
            write.fasta(sequences = seqs2,names = str_c(sample_list[i],"_",genes[j]),file.out = str_c("seqs_of_each_MLSA_genes/",sample_list[i],"_",genes[j],".fasta"))
          }
        }else{
          for(k in seq(length(grep("?",res1)))){
            res2 <- res1[grep("?",res1)][k]
            line1 <- str_split(res2," ")[[1]][str_split(res2," ")[[1]]!=""]
            site <- as.numeric(line1[line1!=""][10:11])
            seqs1 <- paste(as.character(getFrag(seqs[[name]],site[1],site[2])),collapse = "")
            seqs2 <- paste(seqs2,seqs1,sep = "",collapse = "")
            write.fasta(sequences = seqs2,names = str_c(sample_list[i],"_",genes[j]),file.out = str_c("seqs_of_each_MLSA_genes/",sample_list[i],"_",genes[j],".fasta"))
          }
        }
      }
    }
  }
  
  for(j in seq(genes)){
    system(str_c("mkdir seqs_of_each_MLSA_genes/",genes[j]))
    system(str_c("mkdir seqs_of_each_MLSA_genes/",genes[j],"/msa"))
    system(str_c("cp seqs_of_each_MLSA_genes/*",genes[j],".fasta seqs_of_each_MLSA_genes/",genes[j]))
    system(str_c("cat seqs_of_each_MLSA_genes/",genes[j],"/*.fasta > seqs_of_each_MLSA_genes/",genes[j],"/",genes[j],"_seqs.fasta;"))
    system(str_c("muscle -in seqs_of_each_MLSA_genes/",genes[j],"/",genes[j],"_seqs.fasta -out seqs_of_each_MLSA_genes/",genes[j],"/",genes[j],"_msa.fasta"))
    system(str_c("qiime split_fasta_on_sample_ids.py -i seqs_of_each_MLSA_genes/",genes[j],"/",genes[j],"_msa.fasta -o seqs_of_each_MLSA_genes/",genes[j],"/msa/"))
    if(length(dir(str_c("seqs_of_each_MLSA_genes/",genes[j])))!=(length(sample_list)+3)){
      align_len <- length(read.fasta(str_c("seqs_of_each_MLSA_genes/",genes[j],"/",genes[j],"_msa.fasta"))[[1]])
      no_hit_seq <- str_c(sample_list,"_",genes[j],".fasta")[str_c(sample_list,"_",genes[j],".fasta") %in% dir(str_c("seqs_of_each_MLSA_genes/",genes[j],"/"))==FALSE]
      for(l in seq(no_hit_seq)){
        seqs3 <- paste(rep("-",align_len),sep="",collapse = "")
        write.fasta(sequences = seqs3,names = str_sub(no_hit_seq[l],1,-7),file.out = str_c("seqs_of_each_MLSA_genes/",genes[j],"/msa/",no_hit_seq[l]))
      }
    }
    system(str_c("cat seqs_of_each_MLSA_genes/",genes[j],"/msa/*.fasta > seqs_of_each_MLSA_genes/",genes[j],"/msa/",genes[j],"_msa.fasta"))
    system(str_c("Gblocks seqs_of_each_MLSA_genes/",genes[j],"/msa/",genes[j],"_msa.fasta -t=p -b3=10 -b4=5 -b5=H"))
    system(str_c("qiime split_fasta_on_sample_ids.py -i seqs_of_each_MLSA_genes/",genes[j],"/msa/",genes[j],"_msa.fasta-gb -o seqs_of_each_MLSA_genes/",genes[j],"/msa/oneseq_",genes[j],"/"))
  }
  
  for(i in seq(sample_list)){
    system(str_c("mkdir results/sequences/",sample_list[i]))
    seqs5 <- NULL
    for(j in seq(genes)){
      if(length(dir(str_c("seqs_of_each_MLSA_genes/",genes[j],"/msa/oneseq_",genes[j],"/")))==0){
        next
      }else{
        system(str_c("cp seqs_of_each_MLSA_genes/",genes[j],"/msa/oneseq_",genes[j],"/",sample_list[i],".fasta results/sequences/",sample_list[i],"/",sample_list[i],"_",genes[j],".fasta"))
      }
    }
    system(str_c("cat results/sequences/",sample_list[i],"/*.fasta > results/final_seq/",sample_list[i],".fasta"))
    len <- length(readLines(str_c("results/final_seq/",sample_list[i],".fasta")))
    seqs7 <- NULL
    for(k in seq(len)[seq(len) %%2 ==0 ]){
      seqs6 <- readLines(str_c("results/final_seq/",sample_list[i],".fasta"))[k]
      seqs7 <- paste0(seqs7,seqs6)
    }
    write.fasta(sequences = seqs7,names = sample_list[i],file.out = str_c("results/final_seq/oneseq_",sample_list[i],".fasta"))
  }
  system("cat results/final_seq/oneseq*.fasta > results/mlsa.fasta")
  #system(str_c("Gblocks results/mlsa.fasta -t=p -b3=10 -b4=5 -b5=H"))
  system("fasttree -wag results/mlsa.fasta > results/mlsa_tree.tre")
  print("MLSA gene phylogenetic analysis...........................................DONE")
}