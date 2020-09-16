interactions <- function(x, n = NULL) {
  
  meta:::chkclass(x, "netcomb")
  
  comps <- x$comps
  trts <- x$trts
  ##
  ints <- trts[!(trts %in% comps)]
  
  if (!is.null(n)) {
    if (!meta:::is.wholenumber(n))
      stop("Argument 'n' must be a whole number.")
    sel <- unlist(lapply(strsplit(ints, x$sep.comps, fixed = TRUE), length)) == n
    ints <- ints[sel]
  }
  
  ints
}


selcomp <- function(x, n = 1, sep = "+") {
  if (!meta:::is.wholenumber(n))
    stop("Argument 'n' must be a whole number.")
  ##
  res <- lapply(strsplit(x, sep, fixed = TRUE),
                gsub, pattern = " ", replacement = "")
  n.comps <- sapply(res, length)
  ##
  #if (any(n > n.comps))
  #  stop("Argument 'n' must be smaller equal than ", max(n.comps), ".")
  ##
  res <- sapply(res, "[", n)
  
  res
}


## Function for separating a treatment from the full NMA
## - x    a CNMA object
## - trt  a treatment name (character)
## - sep  a new separator for the treatment to be separated from all its combinations
##
## This function only works for combinations of two or three components
##
separate <- function(x, trt, sep = "&",
                     reference.group = x$reference.group) {
  
  meta:::chkclass(x, "netcomb")
  ##
  treat1 <- x$treat1
  treat2 <- x$treat2
  ##
  trts <- x$trts
  sep.comps <- x$sep.comps
  
  
  i2 <- interactions(x, 2)
  i3 <- interactions(x, 3)
  ##
  if (length(i2) > 0) {
    comp1 <- selcomp(i2, 1)
    comp2 <- selcomp(i2, 2)
    ##
    is.comp1 <- comp1 == trt
    is.comp2 <- comp2 == trt
    ##
    sel.trt <- is.comp1 | is.comp2
    ##
    sel1 <- treat1 %in% i2[sel.trt]
    sel2 <- treat2 %in% i2[sel.trt]
    ##
    if (any(sel1))
      treat1[sel1] <- gsub(sep.comps, sep, treat1[sel1], fixed = TRUE)
    if (any(sel2))
      treat2[sel2] <- gsub(sep.comps, sep, treat2[sel2], fixed = TRUE)
  }
  if (length(i3) > 0) {
    comp1 <- selcomp(i3, 1)
    comp2 <- selcomp(i3, 2)
    comp3 <- selcomp(i3, 3)
    ##
    is.comp1 <- comp1 == trt
    is.comp2 <- comp2 == trt
    is.comp3 <- comp3 == trt
    ##
    treat1.old <- treat1
    treat1.comp1 <- selcomp(treat1, 1)
    treat1.comp2 <- selcomp(treat1, 2)
    treat1.comp3 <- selcomp(treat1, 3)
    ##
    treat2.comp1 <- selcomp(treat2, 1)
    treat2.comp2 <- selcomp(treat2, 2)
    treat2.comp3 <- selcomp(treat2, 3)
    ##
    sel.trt <- is.comp1 | is.comp2 | is.comp3
    ##
    sel1 <- treat1 %in% i3[sel.trt]
    sel2 <- treat2 %in% i3[sel.trt]
    ##
    if (any(sel1)) {
      sel1.1 <- treat1 %in% i3[is.comp1]
      sel1.2 <- treat1 %in% i3[is.comp2]
      sel1.3 <- treat1 %in% i3[is.comp3]
      ##
      if (any(sel1.1))
        treat1[sel1.1] <- sub(sep.comps, sep, treat1[sel1.1], fixed = TRUE)
      ##
      if (any(sel1.3)) {
        treat1[sel1.3] <- gsub(sep.comps, sep, treat1[sel1.3], fixed = TRUE)
        treat1[sel1.3] <-  sub(sep, sep.comps, treat1[sel1.3], fixed = TRUE)
      }
      ##
      if (any(sel1.2))
        treat1[sel1.2] <- paste(treat1.comp1[sel1.2], sep.comps,
                                treat1.comp3[sel1.2], sep,
                                treat1.comp2[sel1.2])
    }
    ##
    if (any(sel2)) {
      sel2.1 <- treat2 %in% i3[is.comp1]
      sel2.2 <- treat2 %in% i3[is.comp2]
      sel2.3 <- treat2 %in% i3[is.comp3]
      ##
      if (any(sel2.1))
        treat2[sel2.1] <- sub(sep.comps, sep, treat2[sel2.1], fixed = TRUE)
      ##
      if (any(sel2.3)) {
        treat2[sel2.3] <- gsub(sep.comps, sep, treat2[sel2.3], fixed = TRUE)
        treat2[sel2.3] <-  sub(sep, sep.comps, treat2[sel2.3], fixed = TRUE)
      }
      ##
      if (any(sel2.2))
        treat2[sel2.2] <- paste(treat2.comp1[sel2.2], sep.comps,
                                treat2.comp3[sel2.2], sep,
                                treat2.comp2[sel2.2])
    }
  }
  ##
  res <- discomb(x$TE, x$seTE, treat1, treat2, x$studlab,
                 sep.comps = sep,
                 ref = reference.group, sm = x$sm)
  ##
  attr(res, "trt") <- trt
  attr(res, "rank.X") <- qr(res$X.matrix)$rank
  ##
  res
}
