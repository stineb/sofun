library(plyr)
# library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/get_pointdata_monthly_cru.R", sep="" ) )
source( paste( myhome, "sofun/getin/get_pointdata_temp_wfdei.R", sep="" ) )
source( paste( myhome, "sofun/getin/get_pointdata_prec_wfdei.R", sep="" ) )
source( paste( myhome, "sofun/getin/find_nearest_cruland_by_lat.R", sep="" ) )
source( paste( myhome, "sofun/getin/get_daily_prec.R", sep="" ) )
source( paste( myhome, "sofun/getin/monthly2daily.R", sep="" ) )
source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/monthly2daily.R", sep="" ) )

overwrite <- FALSE

simsuite <- "swbm"

##--------------------------------------------------------------------
##--------------------------------------------------------------------
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
startyr_cru <- 1901
endyr_cru   <- 2013
nyrs_cru    <- length(startyr_cru:endyr_cru)
staryr_wfdei <- 1979
endyr_wfdei  <- 2012

## load meta data file for site simulation
siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ), as.is=TRUE )
nsites <- dim(siteinfo)[1]

for (idx in seq(nsites)){
# for (idx in 1:1){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_clim_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", sep="" )
  filnam_clim_csv <- paste( dirnam_clim_csv, "clim_daily_", sitename, ".csv", sep="" )

  if (overwrite || !file.exists(filnam_clim_csv)){

    ##--------------------------------------------------------------------
    ## get monthly data from CRU
    ##--------------------------------------------------------------------
    clim_cru_monthly <- data.frame( moy=rep( seq(nmonth), nyrs_cru ), year=rep( startyr_cru:endyr_cru, each=nmonth ) )

    ##--------------------------------------------------------------------
    ## air temperature
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$temp <- get_pointdata_monthly_cru( "tmp", lon, lat )
    #clim_cru_monthly$temp <- rep( 20*sin( seq(0, 2*pi, 2*pi/11)-0.5*pi), dim(clim_cru_monthly)[1]/12 )

    ## Check if data is provided, otherwise use nearest gridcell
    if (is.na(clim_cru_monthly$temp[1])) {
      lon_look <- find_nearest_cruland_by_lat( lon, lat )
      clim_cru_monthly$temp <- get_pointdata_monthly_cru( "tmp", lon_look, lat )
    } else {
      lon_look <- lon
    }

    ##--------------------------------------------------------------------
    ## precipitation
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$prec <- get_pointdata_monthly_cru( "pre", lon_look, lat )
    # clim_cru_monthly$prec <- rep( 50.0, dim(clim_cru_monthly)[1] )

    ##--------------------------------------------------------------------
    ## wet days
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$wetd <- get_pointdata_monthly_cru( "wet", lon_look, lat )
    # clim_cru_monthly$wetd <- rep( 5.0, dim(clim_cru_monthly)[1] )

    ##--------------------------------------------------------------------
    ## cloud cover
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$ccov <- get_pointdata_monthly_cru( "cld", lon_look, lat )
    # clim_cru_monthly$wetd <- rep( 5.0, dim(clim_cru_monthly)[1] )

    ##--------------------------------------------------------------------
    ## water vapour 
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$vapr <- get_pointdata_monthly_cru( "vap", lon_look, lat )
    # clim_cru_monthly$wetd <- rep( 5.0, dim(clim_cru_monthly)[1] )

    ##--------------------------------------------------------------------
    ## expanding to daily data (interpolating temp, and generating prec)
    ##--------------------------------------------------------------------
    print( paste( "expanding to daily data for station", sitename, "..." ) )

    dm   <- rep( NA, sum(ndaymonth)*length(startyr_cru:endyr_cru) )
    jdx <- 0
    for (yr in startyr_cru:endyr_cru ){
      for (imoy in 1:nmonth){
        for (idm in 1:ndaymonth[imoy]){
          jdx <- jdx + 1 
          dm[jdx]   <- idm
        }
      }
    }
    clim_daily <- data.frame( 
      doy=rep( seq(ndayyear), nyrs_cru ), 
      moy=rep( rep( seq(nmonth), times=ndaymonth ), times=nyrs_cru ),
      dom=dm,
      year=rep( startyr_cru:endyr_cru, each=ndayyear ) 
    )
    clim_daily$temp <- rep( NA, dim(clim_daily)[1] )
    clim_daily$prec <- rep( NA, dim(clim_daily)[1] )
    clim_daily$ccov <- rep( NA, dim(clim_daily)[1] )
    clim_daily$vapr <- rep( NA, dim(clim_daily)[1] )
    clim_daily$source <- rep( "cru", dim(clim_daily)[1] )

    imo <- 1
    for (iyr in seq( length(clim_cru_monthly$prec)/nmonth ) ){

      use_year <- clim_cru_monthly$year[imo]
      use_year_pvy <- max(startyr_cru,use_year-1)
      use_year_nxt <- min(endyr_cru,use_year+1)
      idxs <- which( clim_daily$year==use_year )
      
      ##--------------------------------------------------------------------
      ## air temperature: interpolate using polynomial
      ##--------------------------------------------------------------------
      #mtemp     <- select( filter( clim_cru_monthly, year==use_year ), temp )$temp
      mtemp <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$temp
      #mtemp_pvy <- select( filter( clim_cru_monthly, year==max(startyr_cru,use_year-1) ), temp )$temp
      mtemp_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$temp
      #mtemp_nxt <- select( filter( clim_cru_monthly, year==min(endyr_cru,use_year+1) ), temp )$temp
      mtemp_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$temp
      clim_daily$temp[ idxs ] <- monthly2daily( mtemp, "polynom", mtemp_pvy[nmonth], mtemp_nxt[1] )

      ##--------------------------------------------------------------------
      ## precipitation: interpolate weather generator
      ##--------------------------------------------------------------------
      #mprec <- select( filter( clim_cru_monthly, year==use_year ), prec )$prec
      mprec <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$prec
      #mwetd <- select( filter( clim_cru_monthly, year==use_year ), wetd )$wetd
      mwetd <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$wetd
      clim_daily$prec[ idxs ] <- get_daily_prec( mprec, mwetd )

      ##--------------------------------------------------------------------
      ## cloud cover: interpolate using polynomial
      ##--------------------------------------------------------------------
      mccov <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$ccov
      mccov_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$ccov
      mccov_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$ccov
      clim_daily$ccov[ idxs ] <- monthly2daily( mccov, "polynom", mccov_pvy[nmonth], mccov_nxt[1] )

      ##--------------------------------------------------------------------
      ## water vapour: interpolate using polynomial
      ##--------------------------------------------------------------------
      mvapr <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$vapr
      mvapr_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$vapr
      mvapr_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$vapr
      clim_daily$vapr[ idxs ] <- monthly2daily( mvapr, "polynom", mvapr_pvy[nmonth], mvapr_nxt[1] )

      imo <- imo + nmonth
    }

    ## Reduce CCOV to a maximum 100%
    clim_daily$ccov[ clim_daily$ccov>100.0 ] <- 100.0

    ##--------------------------------------------------------------------
    ## Get daily WFDEI data where available
    ##--------------------------------------------------------------------
    print( "getting daily data from WFDEI ...")
    for ( yr in staryr_wfdei:endyr_wfdei ){
      for ( moy in seq(nmonth) ){
        print( paste( "... found data for year and month:", yr, moy ) )
        tmp <- get_pointdata_temp_wfdei( lon, lat, moy, yr )
        if (!is.na(tmp[1])) { 
         useidx <- which( clim_daily$year==yr & clim_daily$moy==moy )
         clim_daily$temp[ useidx ]   <- tmp 
         clim_daily$source[ useidx ] <- "wfdei" 
        }
        tmp <- get_pointdata_prec_wfdei( lon, lat, moy, yr )
        if (!is.na(tmp[1])) { 
         clim_daily$prec[ which( clim_daily$year==yr & clim_daily$moy==moy ) ] <- tmp 
        }
      }
    }

    ##--------------------------------------------------------------------
    ## Get site-specific meto data for separate (site-specific) file
    ##--------------------------------------------------------------------
    if (!is.null(siteinfo$meteosource)){
      print("Using site-specific meteo data from separate file ...")
      filn <- paste( "../../", as.character(siteinfo$meteosource[idx] ), sep="" )
      print( paste( "file name", filn))
      meteo <- read.csv( filn )
      
      if ( is.null(meteo$moy) ) { meteo$moy <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mon + 1 }
      if ( is.null(meteo$dom) ) { meteo$dom <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mday }

      for (jdx in 1:dim(meteo)[1]){

        putjdx <- which( clim_daily$year==meteo$year[jdx] & clim_daily$moy==meteo$moy[jdx] & clim_daily$dom==meteo$dom[jdx] )

        if (length(putjdx)>1) { 
          print("PROBLEM: found multiple corresponding indeces for ...")
          print(paste("year =", meteo$year[jdx], "moy =", meteo$moy[jdx], "dom =", meteo$dom[jdx] ) ) 
        } else if (length(putjdx)==0){
          print("PROBLEM: found no corresponding indeces for ...")
          print(paste("year =", meteo$year[jdx], "moy =", meteo$moy[jdx], "dom =", meteo$dom[jdx] ) )           
        }
        if (!is.na(meteo$temp_mean[jdx])) { clim_daily$temp[ putjdx ] <- meteo$temp_mean[jdx] }
        if (!is.na(meteo$temp_mean[jdx])) { clim_daily$source[ putjdx ] <- "temp_sitedata" }

        if (!is.na(meteo$rainfall[jdx])) { clim_daily$prec[ putjdx ]   <- meteo$rainfall[jdx] }
        if (!is.na(meteo$rainfall[jdx])) { clim_daily$source[ putjdx ] <- "prec_sitedata" }

      }

    }

    ##--------------------------------------------------------------------
    ## Save (massive) daily and monthly climate data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing climate data frame into CSV file ", filnam_clim_csv, "..." ) )
    system( paste( "mkdir -p", dirnam_clim_csv ) )   
    write.csv( clim_daily, file=filnam_clim_csv, row.names=FALSE )

  } else {

    print( "daily climate data already available in CSV file ...")
    clim_daily <- read.csv( filnam_clim_csv )

  }

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )
  for ( yr in startyr_cru:endyr_cru ){

    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "dtemp_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, clim_daily$temp[ which( clim_daily$year==yr ) ] )
    
    filnam <- paste( dirnam, "dprec_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, clim_daily$prec[ which( clim_daily$year==yr ) ] )

    filnam <- paste( dirnam, "dfsun_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, ( 100.0 - clim_daily$ccov[ which( clim_daily$year==yr ) ] ) / 100.0 )

    filnam <- paste( dirnam, "dvapr_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, clim_daily$vapr[ which( clim_daily$year==yr ) ] )

  }

}

