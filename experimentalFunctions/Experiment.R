  # Input Values
  geoAccessionCode <- "GSE18380"
  platform <- "GPL4694"

    gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)

    platforms <- list()
    i <-1
    for(dataset in gset) {
      platforms[[i]] <- annotation(dataset)
      i <- i + 1
    }

    if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)

    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) {
      results <- "Yes"
    } else {
      results <- "Yes"
    }
    return(results)