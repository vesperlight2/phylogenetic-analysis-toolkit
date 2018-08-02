pocp_cal <- function(sample_list){
  library(stringr)
  library(seqinr)
  sample_tab <- read.csv(sample_list,row.names = 1)
  list1 <- sample_tab[,1]
  dir.create("blast_res")
  for(i in seq(length(list1))){
    system(str_c("makeblastdb -in aa_orf_",list1[i],".fasta -dbtype prot"))
    for(j in seq(length(list1))){
      system(str_c("blastp -query aa_orf_", list1[j],".fasta -db aa_orf_",list1[i], ".fasta -outfmt 6 -evalue 1e-5 -num_threads 10 -out blast_res/", list1[i],".",list1[j],".txt"))
    }
    print(str_c(i,Sys.time()))
  }
  dis1 <- matrix(nrow=length(list1),ncol=length(list1))
  lenall <- NULL
  for(i in seq(list1)){
    seq1 <- read.fasta(str_c("aa_orf_",list1[i],".fasta"))
    len1 <- lapply(seq1,length)
    lenall <- c(lenall,length(seq1))
    for(j in seq(list1)){
      a <- read.table(str_c("blast_res/",list1[j],".",list1[i],".txt"))
      a1 <- as.data.frame(cbind(len1,names(len1)))
      a2 <- merge(x=a,y=a1,by.x = "V1",by.y="V2",all.y = T)
      a2 <- cbind(a2[,1:12],as.numeric(a2[,13]))
      a3 <- a2[a2[,3]>40 & a2[,4] >= (0.5*a2[,13]),]
      a4 <- a3[!duplicated(a3[,1]),]
      dis1[j,i] <- nrow(a4)
    }
    print(i)
  }
  
  lenall <- as.numeric(lenall)
  pocp <- matrix(nrow=length(list1),ncol=length(list1))
  for(i in seq(list1)){
    for(j in seq(list1)){
      pocp[i,j] <- (dis1[i,j]+dis1[j,i])/(lenall[i]+lenall[j])
    }
  }
  rownames(pocp) <- colnames(pocp) <- list1
  write.csv(pocp,"pocp.csv")
}