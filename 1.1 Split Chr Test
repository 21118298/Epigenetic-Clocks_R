  for (i in 1:22){
  result <- split(x, x$chr)
  assign(paste0(file,"_chr",i,sep=""),clump_data(result[[i]], clump_p1 = 5e-08, clump_p2 = 5e-08)
  #write.table(result[[i]], file = paste0(file,"_chr",i,sep=""), quote = FALSE)
  #y <- read.table(file = paste0(file,"_chr",i,sep=""), header = T)
  #y <- clump_data(y, clump_p1 = 5e-08, clump_p2 = 5e-08)     
  }
