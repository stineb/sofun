get_pointdata_fapar3g <- function( lon, lat ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  filn   <- "../../data/fAPAR/fAPAR3g/fAPAR3g_monthly_1982_2011.nc"

  if ( file.exists( filn ) ){
    
    cmd <- paste( "./extract_pointdata_byfil.sh ", filn, "fAPAR", "lon", "lat", lon, lat )
    print( paste( "executing command:", cmd ) )
    system( cmd )
    out <- read.table( "out.txt" )$V1

  } else {

    print( paste( "file", filn, "does not exist." ) )
    out <- NA

  }
  return( out )
}