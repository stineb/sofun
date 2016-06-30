get_pointdata_temp_wfdei <- function( lon, lat, mo, yr, ignore_leap=TRUE ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  filn <- paste( "../../data/watch_wfdei/Tair_daily/Tair_daily_WFDEI_", sprintf( "%4d", yr ), sprintf( "%02d", mo ), ".nc", sep="" )
  if ( file.exists( filn ) ){
    print( paste( "extracting from", filn ) )
    system( paste( "./extract_pointdata_byfil.sh", filn, "Tair", "lon", "lat", lon, lat ) )
    dtemp <- read.table( "out.txt" )$V1 - 273.15  # conversion from Kelving to Celsius
    if (ignore_leap & length(dtemp==29)){ dtemp <- dtemp[1:28] }
  }
  return( dtemp )
}
