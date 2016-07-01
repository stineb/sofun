get_pointdata_elv_watch <- function( lon, lat ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  filn   <- paste( myhome, "data/watch_wfdei/WFDEI-elevation.nc", sep="" )

  if ( file.exists( filn ) ){
    
    cmd <- paste( paste( myhome, "sofun/getin/extract_pointdata_byfil.sh "), filn, "elevation", "lon", "lat", lon, lat )
    system( cmd )
    out <- read.table( paste( syshome, "/tmp/out.txt", sep="" ) )$V1

  } else {

    print( paste( "file", filn, "does not exist." ) )
    out <- NA

  }
  return( out )
}