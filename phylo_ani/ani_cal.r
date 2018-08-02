ani_cal <- function(sample_list){
  library(stringr)
  sample_tab <- read.csv(sample_list,row.names = 1)
  f1 <- sample_tab[,1]
  a1 <- NULL
  dir.create("res")
  for(i in seq(length(f1))){
    for(j in seq(length(f1))){
      system(str_c("java -jar /home/sam/software/OAT_cmd.jar -fasta1 ",f1[i],".fasta -fasta2 ",f1[j],
      ".fasta -method ani -blastplus_dir /usr/bin | cat > res/",f1[i],"_",f1[j],".txt"))
      a <- readLines(str_c("res/",f1[i],"_",f1[j],".txt"))
      a <- a[length(a)-1]
      a1 <- rbind(a1,c(as.character(f1[i]),as.character(f1[j]),a))
    }
    print(str_c(i,Sys.time()))
  }  
  write.csv(a1,"ani_list.csv")
}
