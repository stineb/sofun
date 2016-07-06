library(plyr)
#library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/get_pointdata_ndep_lamarque.R", sep="" ) )
source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )

simsuite <- "fluxnet_cnmodel"

##--------------------------------------------------------------------
##--------------------------------------------------------------------
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

for (idx in seq(nsites)){
# for (idx in 1:1){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_ndep_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/ndep/", sitename, "/", sep="" )
  filnam_ndep_csv <- paste( dirnam_ndep_csv, "ndep_lamarque_", sitename, ".csv", sep="" )

  if (!file.exists(filnam_ndep_csv)){

    ##--------------------------------------------------------------------
    ## get annual N deposition data
    ##--------------------------------------------------------------------
    ndep_ann <- get_pointdata_ndep_lamarque( lon, lat )

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing ndep data frame into CSV file ", filnam_ndep_csv, "...") )
    system( paste( "mkdir -p", dirnam_ndep_csv ) )   
    write.csv( ndep_ann, file=filnam_ndep_csv, row.names=FALSE )

  } else {

    print( "daily climate data already available in CSV file ...")
    ndep_ann <- read.csv( filnam_ndep_csv )

  }

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )

  ## NOy in separate formatted file 
  filnam <- paste( dirnam_ndep_csv, "ndep_noy_lamarque_", sitename, ".txt", sep="" )
  write_sofunformatted( filnam, ndep_ann[,c(1,2)] )
  
  ## NHx in separate formatted file 
  filnam <- paste( dirnam_ndep_csv, "ndep_nhx_lamarque_", sitename, ".txt", sep="" )
  write_sofunformatted( filnam, ndep_ann[,c(1,3)] )

}

