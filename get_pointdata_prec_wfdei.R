get_pointdata_prec_wfdei <- function( lon, lat, mo, yr, ignore_leap=TRUE ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location). 
  ## Original data in kg/m2s, returns data in kg/m2/month.
  ##--------------------------------------------------------------------
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  
  ## rain
  filn <- paste( myhome, "data/watch_wfdei/Rainf_daily/Rainf_daily_WFDEI_CRU_", sprintf( "%4d", yr ), sprintf( "%02d", mo ), ".nc", sep="" )
  if ( file.exists( filn ) ){
    print( paste( "extracting from", filn ) )
    system( paste( paste( myhome, "sofun/getin/extract_pointdata_byfil.sh ", sep="" ), filn, "Rainf", "lon", "lat", sprintf("%.2f",lon), sprintf("%.2f",lat) ) )
    dprec <- read.table( paste( syshome, "/tmp/out.txt", sep="" ) )$V1
    dprec <- dprec*60*60*24 # kg/m2/s -> mm/day
  } else {
    # print( paste( "file", filn, "does not exist." ) )
    dprec <- rep( NA, ndaymonth[mo] )
  }
  # print( paste( "rain only: ", mprec))

  ## snow
  filn <- paste( myhome, "data/watch_wfdei/Snowf_daily/Snowf_daily_WFDEI_CRU_", sprintf( "%4d", yr ), sprintf( "%02d", mo ), ".nc", sep="" )
  if ( file.exists( filn ) ){
    print( paste( "extracting from", filn ) )
    system( paste( paste( myhome, "sofun/getin/extract_pointdata_byfil.sh ", sep="" ), filn, "Snowf", "lon", "lat", sprintf("%.2f",lon), sprintf("%.2f",lat) ) )
    dsnow <- read.table( paste( syshome, "/tmp/out.txt", sep="" ) )$V1
    dsnow <- dsnow*60*60*24 # kg/m2/s -> mm/day
    dprec <- dprec + dsnow
    # print( paste( "snow only: ", sum( dprec*60*60*24 )))
  }

  ## ignore leap years
  if (ignore_leap & length(dprec==29)){ dprec <- dprec[1:28] }

  return( dprec )
}
