get_pointdata_monthly_cru <- function( varnam, lon, lat, yr=NA ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  filn   <- paste( paste( myhome, "data/cru/ts_3.22/cru_ts3.22.1901.2013.", sep=""), varnam, ".dat.nc", sep="" )
  if (!is.na(yr)){
    istart <- max( 1     , (yr - 1901) * nmonth + 1 )
    iend   <- max( nmonth, (yr - 1901) * nmonth + nmonth )    
  }
  if ( file.exists( filn ) ){
    cmd <- paste( paste( myhome, "extract_pointdata_byfil.sh ", sep="" ), filn, varnam, "lon", "lat", lon, lat )
    print( paste( "executing command:", cmd ) )
    system( cmd )
    mdata_full <- read.table( paste( myhome, "sofun/getin/out.txt", sep="" ) )$V1
    if (!is.na(yr)){
      mdata <- mdata_full[ istart:iend ]
    } else {
      mdata <- mdata_full
    }
  } else {
    print( paste( "file", filn, "does not exist." ) )
    mdata <- NA
  }
  return( mdata )
}