phylo_ribo <- function(seq_dir){
  library(seqinr)
  library(stringr)
  if(str_sub(seq_dir,-1,-1)!="/"){
    seq_dir <- str_c(seq_dir,"/")
  }
  
  ##Ribosome protein family IDs
  ribo_list <- c("PF00687","PF00181","PF00297","PF00573","PF00281","PF17144","PF00347","PF00542","PF01281","PF00466",
                 "PF00298","PF00572","PF00238","PF00252","PF01196","PF00861","PF01245","PF00453","PF00829","PF00237",
                 "PF00276","PF17136","PF01386","PF01016","PF00830","PF00831","PF00327","PF01197","PF01783","PF00471",
                 "PF01632","PF00444","PF10501","PF00318","PF00189","PF00163","PF00333","PF01250","PF00177","PF00410",
                 "PF00380","PF00338","PF00411","PF00416","PF00253","PF00312","PF00886","PF00366","PF01084","PF00203",
                 "PF01649","PF01165","PF02482")
  seq_list <- dir(seq_dir)
  seq_name <- str_sub(seq_list,1,-7)
  sample_list <- str_sub(seq_list,8,-7)
  
  #creat main folders
  dir.create("hmmsearch_results_of_each_pfam",recursive = TRUE)
  dir.create("seqs_of_each_pfam",recursive = TRUE)
  dir.create("results",recursive = TRUE)
  dir.create("results/sequences",recursive = TRUE)
  dir.create("results/final_seq",recursive = TRUE)
  
  ##If no ribo hmm file, download
  if(("ribo_hmm_files" %in% dir())==FALSE){
    system("mkdir ribo_hmm_files")
    for(i in seq(ribo_list)){
      ##Download pfam hmm files from pfam websit
      system(str_c("wget -O ribo_hmm_files/",ribo_list[i],".hmm http://pfam.xfam.org/family/",ribo_list[i],"/hmm;"))
      print(str_c("Downloading Ribosome Protein Family Hmm File:",ribo_list[i],"..........",round(i/length(ribo_list),digits = 3)*100,"%"))
      ##sub folder of "hmmsearch_results_of_each_pfam"  named with pfam IDs were used for storing hmmsearch results of each pfam
    }
  }
  #search sequences based on hmm files
  for(i in seq(ribo_list)){
    ##sub folder of "hmmsearch_results_of_each_pfam"  named with pfam IDs were used for storing hmmsearch results of each pfam
    for(j in seq(seq_list)){
      system(str_c("hmmsearch ribo_hmm_files/",ribo_list[i],".hmm ",seq_dir,seq_list[j]," > hmmsearch_results_of_each_pfam/",ribo_list[i],"_",seq_name[j],"_results.txt"))
    }
    print(str_c("Searching and extraxting Ribosome Protein Sequences of ",ribo_list[i],"..........",round(i/length(ribo_list),digits = 3)*100,"%"))
  }
  #extract sequence from AA_orf sequences based on hmmsearch results of best hit
  no_hit <- NULL
  for(i in seq(sample_list)){
    seqs <- read.fasta(str_c(seq_dir,seq_list[i]))
    for(j in seq(ribo_list)){
      res <- readLines(str_c("hmmsearch_results_of_each_pfam/",ribo_list[j],"_",seq_name[i],"_results.txt"))
      name <- str_split(res[grep(">>",res)][1]," ")[[1]][2]
      if(is.na(name)==TRUE){
        no_hit <- rbind(no_hit,c(sample_list[i],ribo_list[j]))
      }else{
        res1 <- res[(grep(">>",res)[1]+3):(grep("==",res)[1]-3)]
        line2 <- NULL
        len4 <- NULL
        if(length(grep("!",res1))>0){
          for(k in seq(length(grep("!",res1)))){
            res2 <- res1[grep("!",res1)][k]
            line1 <- str_split(res2," ")[[1]][str_split(res2," ")[[1]]!=""]
            line2 <- rbind(line2,line1)
            len3 <- as.numeric(line1[11])-as.numeric(line1[10])
            len4 <- c(len4,len3)
          }
          res3 <- res1[grep("!",res1)][seq(nchar(len4))[len4==max(len4)]]
          line3 <- str_split(res3," ")[[1]][str_split(res3," ")[[1]]!=""]
          site <- as.numeric(line3[line3!=""][10:11])
          seqs1 <- paste(as.character(getFrag(seqs[[name]],site[1],site[2])),collapse = "")
          write.fasta(sequences = seqs1,names = str_c(sample_list[i],"_",ribo_list[j]),file.out = str_c("seqs_of_each_pfam/",sample_list[i],"_",ribo_list[j],".fasta"))
        }else{
          for(l in seq(length(grep("?",res1)))){
            res2 <- res1[grep("?",res1)][l]
            line1 <- str_split(res2," ")[[1]][str_split(res2," ")[[1]]!=""]
            line2 <- rbind(line2,line1)
            len3 <- as.numeric(line1[11])-as.numeric(line1[10])
            len4 <- c(len4,len3)
          }
          res3 <- res1[grep("?",res1)][seq(nchar(len4))[len4==max(len4)]]
          line3 <- str_split(res3," ")[[1]][str_split(res3," ")[[1]]!=""]
          site <- as.numeric(line3[line3!=""][10:11])
          seqs1 <- paste(as.character(getFrag(seqs[[name]],site[1],site[2])),collapse = "")
          write.fasta(sequences = seqs1,names = str_c(sample_list[i],"_",ribo_list[j]),file.out = str_c("seqs_of_each_pfam/",sample_list[i],"_",ribo_list[j],".fasta"))
        }
      }
    }
  }
  #move seqs to folder of each pfam
  #combine seqs and do MSA
  for(j in seq(ribo_list)){
    system(str_c("mkdir seqs_of_each_pfam/",ribo_list[j]))
    system(str_c("mkdir seqs_of_each_pfam/",ribo_list[j],"/msa"))
    system(str_c("cp seqs_of_each_pfam/*",ribo_list[j],".fasta seqs_of_each_pfam/",ribo_list[j]))
    system(str_c("cat seqs_of_each_pfam/",ribo_list[j],"/*.fasta > seqs_of_each_pfam/",ribo_list[j],"/",ribo_list[j],"_seqs.fasta;"))
    system(str_c("muscle -in seqs_of_each_pfam/",ribo_list[j],"/",ribo_list[j],"_seqs.fasta -out seqs_of_each_pfam/",ribo_list[j],"/",ribo_list[j],"_msa.fasta"))
    system(str_c("qiime split_fasta_on_sample_ids.py -i seqs_of_each_pfam/",ribo_list[j],"/",ribo_list[j],"_msa.fasta -o seqs_of_each_pfam/",ribo_list[j],"/msa/"))
    if(length(dir(str_c("seqs_of_each_pfam/",ribo_list[j])))!=(length(sample_list)+3)){
      align_len <- length(read.fasta(str_c("seqs_of_each_pfam/",ribo_list[j],"/",ribo_list[j],"_msa.fasta"))[[1]])
      no_hit_seq <- str_c(sample_list,"_",ribo_list[j],".fasta")[str_c(sample_list,"_",ribo_list[j],".fasta") %in% dir(str_c("seqs_of_each_pfam/",ribo_list[j],"/"))==FALSE]
      for(l in seq(no_hit_seq)){
        seqs3 <- paste(rep("-",align_len),sep="",collapse = "")
        write.fasta(sequences = seqs3,names = str_sub(no_hit_seq[l],1,-7),file.out = str_c("seqs_of_each_pfam/",ribo_list[j],"/msa/",no_hit_seq[l]))
      }
    }
    system(str_c("cat seqs_of_each_pfam/",ribo_list[j],"/msa/*.fasta > seqs_of_each_pfam/",ribo_list[j],"/msa/",ribo_list[j],"_msa.fasta"))
    system(str_c("qiime split_fasta_on_sample_ids.py -i seqs_of_each_pfam/",ribo_list[j],"/msa/",ribo_list[j],"_msa.fasta -o seqs_of_each_pfam/",ribo_list[j],"/msa/oneseq_",ribo_list[j],"/"))
  }
  
  for(i in seq(sample_list)){
    system(str_c("mkdir results/sequences/",sample_list[i]))
    seqs5 <- NULL
    for(j in seq(ribo_list)){
      if(length(dir(str_c("seqs_of_each_pfam/",ribo_list[j],"/msa/oneseq_",ribo_list[j],"/")))==0){
        next
      }else{
        system(str_c("cp seqs_of_each_pfam/",ribo_list[j],"/msa/oneseq_",ribo_list[j],"/",sample_list[i],".fasta results/sequences/",sample_list[i],"/",sample_list[i],"_",ribo_list[j],".fasta"))
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
  system("cat results/final_seq/oneseq*.fasta > results/ribo.fasta")
  system(str_c("Gblocks results/ribo.fasta -t=p -b3=10 -b4=5 -b5=H"))
  system("fasttree -wag results/ribo.fasta-gb > results/ribo_tree.tre")
  print("Ribosome gene phylogenetic analysis...........................................DONE")
}
