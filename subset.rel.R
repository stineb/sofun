subset.rel <- function( data, bycol, selection ){
  data_sub <- data[ is.element( as.character( data[[ bycol ]] ), as.character( selection ) ), ]
  return(data_sub)
}
