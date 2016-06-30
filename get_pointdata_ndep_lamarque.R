get_pointdata_ndep_lamarque <- function( lon, lat, yr=NA ){
  ##--------------------------------------------------------------------
  ## Extract monthly data from files for each year and attach to the 
  ## monthly dataframe (at the right location).
  ## Original data in K, returns data in K
  ##--------------------------------------------------------------------
  filn   <- "../../data/ndep_lamarque/Ndep_Lamarque11cc_historical_halfdeg_NEW.nc"
  if (!is.na(yr)){
    istart <- max( 1   , (yr - 1850) + 1 )
    iend   <- max( 2009, (yr - 1850) + 1 )    
  }
  if ( file.exists( filn ) ){
    
    ## TIME      
    cmd <- paste( "./extract_pointdata_byfil.sh ", filn, "TIME", "LON", "LAT", lon, lat )
    # print( paste( "executing command:", cmd ) )
    system( cmd )
    year <- read.table( "out.txt" )$V1

    ## NOy      
    cmd <- paste( "./extract_pointdata_byfil.sh ", filn, "NOy", "LON", "LAT", lon, lat )
    # print( paste( "executing command:", cmd ) )
    system( cmd )
    noy <- read.table( "out.txt" )$V1

    ## NHx      
    cmd <- paste( "./extract_pointdata_byfil.sh ", filn, "NHx", "LON", "LAT", lon, lat )
    # print( paste( "executing command:", cmd ) )
    system( cmd )
    nhx <- read.table( "out.txt" )$V1

    out_full <- data.frame( year=year, noy=noy, nhx=nhx )

    if (!is.na(yr)){
      out <- out_full[ istart:iend, ]
    } else {
      out <- out_full
    }

  } else {

    print( paste( "file", filn, "does not exist." ) )
    out <- NA

  }
  return( out )
}