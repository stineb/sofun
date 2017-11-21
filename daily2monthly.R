daily2monthly <- function( dval, method ){
  ##/////////////////////////////////////////////////////////////////////////
  ## Returns monthly values as a mean over daily values in each month.
  ## Arguments:
  ## - dval   : vector containing 365 (366 in case lapyear is TRUE) daily values
  ## - method : true of monthly values represent total of daily values in resp. month
  ##-------------------------------------------------------------------------

  ## parameters
  if (length(dval)==366) {
    ndaymonth <- c(31,29,31,30,31,30,31,31,30,31,30,31)    
  } else {
    ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  }
  ndayyear <- sum(ndaymonth)    
  nmonth <- 12

  istart <- cumsum(ndaymonth)-ndaymonth+1
  iend   <- cumsum(ndaymonth)

  ## loop over months and take sum/mean of days in that month
  mval <- rep( NA, nmonth )
  for (moy in 1:nmonth){
    
    if (method == "sum"){
      mval[moy] <- sum( dval[ istart[moy] : iend[moy] ])
    } else if (method == "mean") {
      mval[moy] <- sum( dval[ istart[moy] : iend[moy] ]) / ndaymonth[moy]
    } else {
      print( "DAILY2MONTHLY: select valid method (sum, mean)" )
    }
  }
  return( mval )
}
