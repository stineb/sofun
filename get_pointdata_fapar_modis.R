get_pointdata_fapar_modis <- function( lon, lat, yr ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  filn   <- paste( "../../data/fAPAR/monthly_0.5deg_MODIS-EVI-based/ISI-MIP_", as.character( yr ), "-fAPAR_0.5deg.nc", sep="" )

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