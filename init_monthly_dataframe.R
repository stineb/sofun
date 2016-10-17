init_monthly_dataframe <- function( yrstart, yrend ){

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350)
  nyrs    <- length( yrstart:yrend )

  mdf <- data.frame( 
    year=rep(yrstart:yrend,each=12) , 
    moy=rep(1:12,nyrs), 
    doy=rep(middaymonth,nyrs)
    )
  mdf$date <- as.POSIXlt( as.Date( paste( as.character(mdf$year), "-01-01", sep="" ) ) + mdf$doy - 1 )
  mdf$year_dec <- mdf$year + ( mdf$doy - 1 ) / sum( ndaymonth )

  return( mdf)

}
