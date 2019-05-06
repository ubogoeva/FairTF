make_deg <- function(my_data, colname='logFC_upd', lg=0) {
  colnames(my_data) <- c('GeneID', 'logFC', 'FDR')
  new_gens <- as.data.frame(my_data$GeneID)
  for (i in 1:nrow(my_data)) {
    if (is.na(my_data[i, 3]) | is.na(my_data[i, 2])) { next
      cat(n)
    }
    if (my_data[i, 3] < 0.05 & (my_data$logFC[i] > lg | my_data$logFC[i] < -lg)) {
      new_gens[i, 2] <- my_data$logFC[i]
    }
    else {
      new_gens[i, 2] <- 0
    }
  }
  colnames(new_gens) <- c('GeneID', colname)
  return(new_gens)
}