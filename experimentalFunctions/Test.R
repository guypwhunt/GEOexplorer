############################################
allColumns <- list ("GSM458594", "GSM458595", "GSM458596", "GSM458597", "GSM458598", "GSM458599", "GSM458600", "GSM458601") 
group1 <- list ("GSM458594", "GSM458595", "GSM458596", "GSM458597") 
group2 <- list ("GSM458598", "GSM458599", "GSM458600", "GSM458601") 
group3 <- list ()

calculateGsms <- function(columnNames,group1, group2, group3){
  lengthOfColumns <- sum(unlist(lapply(columnNames, length)))
  gsmsList <- vector(mode = "list", length = lengthOfColumns)
  i <- 1
  
  for (column in columnNames){
    if (column %in% group1) {
      gsmsList[[i]] <- 0
      i <- i+ 1
    } else if (column %in% group2) {
      gsmsList[[i]] <- 1
      i <- i+ 1
    } else if (column %in% group3) {
      gsmsList[[i]] <- 2
      i <- i+ 1
    }
  }
  gsms <- paste(gsmsList, collapse = '')
  return(gsms)
}

x <- calculateGsms(ex, group1, group2, group3)

print(x)


################################