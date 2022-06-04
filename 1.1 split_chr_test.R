  for (i in 1:22){
  result <- split(x, x$chr)
  #assign(paste0(file,"_chr",i,sep=""),clump_data(result[[i]], clump_p1 = 5e-08, clump_p2 = 5e-08)
  y <- clump_data(result[[i]], clump_p1 = 5e-08, clump_p2 = 5e-08)
  write.table(y, file = paste0(file,"_chr",i,sep=""), quote = FALSE)
  }
  z <- rbind(...)
