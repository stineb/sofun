get_pointdata_prec_wfdei <- function( lon, lat, mo, yr, ignore_leap=TRUE ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location). 
  ## Original data in kg/m2s, returns data in kg/m2/month.
  ##--------------------------------------------------------------------
  
  ## rain
  filn <- paste( "../../data/watch_wfdei/Rainf_daily/Rainf_daily_WFDEI_CRU_", sprintf( "%4d", yr ), sprintf( "%02d", mo ), ".nc", sep="" )
  if ( file.exists( filn ) ){
    print( paste( "extracting from", filn ) )
    system( paste( "./extract_pointdata_byfil.sh", filn, "Rainf", "lon", "lat", lon, lat ) )
    dprec <- read.table( "out.txt" )$V1
    dprec <- dprec*60*60*24 # kg/m2/s -> mm/day
  } else {
    print( paste( "file", filn, "does not exist." ) )
  }
  # print( paste( "rain only: ", mprec))

  ## snow
  filn <- paste( "../../data/watch_wfdei/Snowf_daily/Snowf_daily_WFDEI_CRU_", sprintf( "%4d", yr ), sprintf( "%02d", mo ), ".nc", sep="" )
  if ( file.exists( filn ) ){
    print( paste( "extracting from", filn ) )
    system( paste( "./extract_pointdata_byfil.sh", filn, "Snowf", "lon", "lat", lon, lat ) )
    dsnow <- read.table( "out.txt" )$V1
    dsnow <- dsnow*60*60*24 # kg/m2/s -> mm/day
    dprec <- dprec + dsnow
    # print( paste( "snow only: ", sum( dprec*60*60*24 )))
  } else {
    # print( paste( "file", filn, "does not exist." ) )
    # print( paste( "snow only: ", 0.0 ) )
  }

  ## ignore leap years
  if (ignore_leap & length(dprec==29)){ dprec <- dprec[1:28] }

  return( dprec )
}
