read_data_overlap  <- function(working_dir){
  old_dir <- getwd()
  setwd(working_dir)
  df  <- data.frame()
  number  <- 0
  for (i in dir(path = working_dir, pattern = "*.overlap15")){
    if ((file.size(i, units = 'B'))[1] <= 1) { next
      cat(n)
    }
    temp_df  <- (read.table(i, stringsAsFactors = F))
    temp_df$V1 <- i
    df  <- rbind(df, temp_df)
    number <- number + 1
  }
  print(paste(as.character(number), "files were combined"))
  setwd(old_dir)
  return(df)
}
