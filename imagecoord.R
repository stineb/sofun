imagecoord <- function( inpos, incoords ) {
  ## /////////////////////////////////////////////////////////////////////////
  ## Function returns x (y) coordinates as a function of pos0 given lon (lat)
  ## for use in image().
  ## Beni Stocker
  ## -------------------------------------------------------------------------
  
  dcoords <- incoords[2]-incoords[1]
  ncoords <- length(incoords)

  ## Shift lower margin to ursprung
  shiftcoords  <- incoords -incoords[1]+dcoords/2
  shiftpos     <- inpos    -incoords[1]+dcoords/2

  ## Normalise to [0,1]
  normcoords  <- shiftcoords /(shiftcoords[ncoords]+dcoords/2)
  normpos     <- shiftpos    /(shiftcoords[ncoords]+dcoords/2)

  return(normpos)

}



