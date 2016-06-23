library(plyr)
#library(dplyr)

source( "write_sofunformatted.R" )

overwrite <- FALSE

##--------------------------------------------------------------------
##--------------------------------------------------------------------
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

startyr <- 2000
endyr <- 2000

siteinfo <- read.csv( "../../input_splash_test_sofun/siteinfo_splash_test_sofun.csv" )
nsites <- dim(siteinfo)[1]

for (idx in seq(nsites)){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_clim_csv <- paste( "../../input_splash_test_sofun/sitedata/climate/", sitename, "/", sep="" )
  filnam_clim_csv <- paste( dirnam_clim_csv, "clim_daily_", sitename, ".csv", sep="" )

  system( paste( "mkdir -p", dirnam_clim_csv ) )   

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )
  for ( yr in startyr:endyr ){

    print( paste("... for year", yr))
    dirnam <- paste( "../../input_splash_test_sofun/sitedata/climate/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    dtemp <- read.table( "/alphadata01/bstocker/splash/data/daily_tair_2000_wfdei.txt" )
    filnam <- paste( dirnam, "dtemp_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, dtemp$V1 )
    
    dprec <- read.table( "/alphadata01/bstocker/splash/data/daily_pn_2000_wfdei.txt" )
    filnam <- paste( dirnam, "dprec_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, dprec$V1 )

    dfsun <- read.table( "/alphadata01/bstocker/splash/data/daily_sf_2000_cruts.txt" )
    filnam <- paste( dirnam, "dfsun_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, dfsun$V1 )

  }

}

