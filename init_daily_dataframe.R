init_daily_dataframe <- function( yrstart, yrend ){

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

  ndayyear <- sum(ndaymonth) 
  nmonth   <- length(ndaymonth)
  nyrs    <- length( yrstart:yrend )
  dm   <- rep( NA, sum(ndaymonth)*length(yrstart:yrend) )
  jdx <- 0
  for (yr in yrstart:yrend ){
    for (imoy in 1:nmonth){
      for (idm in 1:ndaymonth[imoy]){
        jdx <- jdx + 1 
        dm[jdx]   <- idm
      }
    }
  }
  ddf <- data.frame( 
    doy=rep( seq(ndayyear), nyrs ), 
    moy=rep( rep( seq(nmonth), times=ndaymonth ), times=nyrs ),
    dom=dm,
    year=rep( yrstart:yrend, each=ndayyear ) 
  )
  ddf$date <- as.POSIXlt( as.Date( paste( as.character(ddf$year), "-01-01", sep="" ) ) + ddf$doy - 1 )
  ddf$year_dec <- ddf$year + ( ddf$doy - 1 ) / sum( ndaymonth )

  return( ddf )
}
