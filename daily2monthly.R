daily2monthly <- function( dval, method, lapyear=FALSE ){
  ##/////////////////////////////////////////////////////////////////////////
  ## Returns monthly values as a mean over daily values in each month.
  ## Arguments:
  ## - dval   : vector containing 365 (366 in case lapyear is TRUE) daily values
  ## - method : true of monthly values represent total of daily values in resp. month
  ##-------------------------------------------------------------------------

  ## parameters
  if (lapyear){
    ndaymonth <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  } else {
    ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  }
  ndayyear <- sum(ndaymonth)    
  nmonths <- 12

  ## loop over months and take sum/mean of days in that month
  for (moy in 1:nmonth){
    if (method == "sum"){
      mval(moy) = sum( dval[ sum( ndaymonth[1:(moy-1)] )+1 : sum( ndaymonth[1:moy] ) ])
    } else if (method == "mean") {
      mval(moy) = sum( dval[ sum( ndaymonth[1:(moy-1)] ) +1 : sum( ndaymonth[1:moy] ) ] ) / ndaymonth[moy]
    } else {
      print( "DAILY2MONTHLY: select valid method (sum, mean)" )
    }
  }
  return( mval )
}
