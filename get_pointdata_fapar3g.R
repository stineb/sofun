get_pointdata_fapar3g <- function( lon, lat ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  filn   <- paste( myhome, "data/fAPAR/fAPAR3g/fAPAR3g_monthly_1982_2011.nc", sep="" )

  if ( file.exists( filn ) ){
    
    cmd <- paste( paste( myhome, "sofun/getin/extract_pointdata_byfil.sh ", sep="" ), filn, "fAPAR", "lon", "lat", sprintf("%.2f",lon), sprintf("%.2f",lat) )
    print( paste( "executing command:", cmd ) )
    system( cmd )
    out <- read.table( paste( syshome, "/tmp/out.txt", sep="" ) )$V1

  } else {

    print( paste( "file", filn, "does not exist." ) )
    out <- NA

  }
  return( out )
}