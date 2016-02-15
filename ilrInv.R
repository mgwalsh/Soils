#' Nonorthonormal ilr inverses
#' J. Chen, January 2016

library(compositions)
library(raster)

# Inverse ilr transformation when the input is a dataframe of ilr transformed data. V is the ilr transformation matrix

ilrInv_nonorthonormal <- function(idata, V){
  n <- dim(idata)[2]
  LT <- (cbind(diag(n), rep(-1, n)))%*%V
  cdata_tmp <- as.matrix(idata)%*%solve(LT)
  cdata_tmp <- cbind(cdata_tmp, -rowSums(cdata_tmp))
  orig_data <- clrInv(cdata_tmp)
  return(orig_data)
}

# Inverse ilr transformation when the input data a raster stack of ilr transformed data. V is the ilr transformation matrix

ilrInv_nonorthonormal_stack <- function(idata_stack, V){
  idata_withcoords <- cbind(coordinates(idata_stack), getValues(idata_stack))
  idata_withcoords <- na.omit(idata_withcoords)
  idata <- idata_withcoords[,-c(1,2)]
  coords <- idata_withcoords[,1:2]
  n <- dim(idata)[2]
  LT <- (cbind(diag(n), rep(-1, n)))%*%V
  cdata_tmp <- idata%*%solve(LT)
  cdata_tmp <- cbind(cdata_tmp, -rowSums(cdata_tmp))
  orig_data <- clrInv(cdata_tmp)
  return(cbind(coords, orig_data))
}
