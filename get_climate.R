##--------------------------------------------------------------------
## This script processes standard input data for SOFUN, reading from
## original data files sitting on Imperial's CX1.
## Reading from NetCDF maps uses NCO command 'ncks' (see extract_pointdata_byfil.sh)
## First, a monthly data frame is created holding all months given in 
## the CRU TS data (currently CRU TS 3.23), covering years 1901-2014
## Second, monthly data is interpolated using a weather generator for 
## precip and a polyonmial interpolation (conserving monthly averages)
## for other variables.
## Third actual daily data for temp. and precip. is read from WATCH-WFDEI 
## daily NetCDF data and overwrites interpolated data.
## Fourth, actual daily site-specific data is read if given and overwrites
## WATCH-WFDEI DATA.
## Fifth, the (massive) dataframe is written to a CSV file for each station.
## Finally, for FORTRAN SOFUN-read-in, text files for each station, year, and 
## variable are written.
##
## Standard (required) SOFUN input variables are:
## - air temerature (deg C)
## - precipitation (mm d-1)
## - cloud cover (%)
## - VPD (Pa)
##
## Optional SOFUN input variables (required if resp. simulation parameter is true):
## - PPFD (mol d-1 m-2)
## - net radiation (W m-2)
##--------------------------------------------------------------------
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
source( paste( myhome, "sofun/getin/calc_vpd.R", sep="" ) )

overwrite <- TRUE
ingest_meteodata <- TRUE

simsuite <- "fluxnet2015"

in_ppfd   <- FALSE
in_netrad <- FALSE

startyr_override <- 1980

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
startyr_cru <- 1901
endyr_cru   <- 2014
nyrs_cru    <- length(startyr_cru:endyr_cru)
staryr_wfdei <- 1979
endyr_wfdei  <- 2012

## load meta data file for site simulation
siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ), as.is=TRUE )
nsites <- dim(siteinfo)[1]
do.sites <- seq(nsites)
# do.sites <- 35:35

get_meteo_fluxnet2015 <- function( path ){
  ##--------------------------------------------------------------------
  ## Function returns a dataframe containing all the data of flux-derived
  ## GPP for the station implicitly given by path (argument).
  ## Specific for FLUXNET 2015 data
  ##--------------------------------------------------------------------
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  ndayyear <- sum(ndaymonth)
  nmonth <- length(ndaymonth)

  ## from flux to energy conversion, umol/J (Meek et al., 1984), same as used in SPLASH (see Eq.50 in spash_doc.pdf)
  kfFEC <- 2.04

  ## get daily meteo data
  meteo <- read.csv( path, na.strings="-9999" )  

  ## add three columns: year, month, day in month
  meteo$year <- as.numeric(substr(meteo$TIMESTAMP,start=1,stop=4))
  meteo$moy  <- as.numeric(substr(meteo$TIMESTAMP,start=5,stop=6))
  meteo$dom  <- as.numeric(substr(meteo$TIMESTAMP,start=7,stop=8))

  meteo$year_dec <- meteo$year + (meteo$moy-1)/nmonth + (meteo$dom-1)/sum(ndaymonth)

  ## rename variables (columns)
  meteo <- rename( meteo, c( "TA_F"="temp", "VPD_F"="vpd", "P_F"="prec", "NETRAD"="nrad" ) ) #, "LW_IN_F"="lwin", "SW_IN_F"="swin" ) )

  ## Take only net PPFD (=in - out)
  if (!is.null(meteo$PPFD_IN) && !is.null(meteo$PPFD_OUT)){
    meteo$ppfd <- meteo$PPFD_IN - meteo$PPFD_OUT    
  } else {
    meteo$ppfd <- rep( NA, dim(meteo)[1] )
  }

  ## convert units
  meteo$vpd  <- meteo$vpd  * 1e2  # given in hPa, required in Pa
  meteo$ppfd <- meteo$ppfd * 1.0e-6 * kfFEC * 60 * 60 * 24  # given in W m-2, required in mol m-2 d-1 
  meteo$nrad <- meteo$nrad * 60 * 60 * 24  # given in W m-2 (avg.), required in J m-2 (daily total)

  # meteo <- select( meteo, year, moy, dom, year_dec, temp, prec, vpd, ppfd, nrad ) 
  meteo <- subset( meteo, select=c( year, moy, dom, year_dec, temp, prec, vpd, ppfd, nrad ) )

  return( meteo )

}

for (idx in do.sites ){
# for (idx in 2:2){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  dirnam_clim_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", sep="" )
  filnam_clim_csv <- paste( dirnam_clim_csv, "clim_daily_", sitename, ".csv", sep="" )

  if (overwrite || !file.exists(filnam_clim_csv)){
  
    print( paste( "collecting monthly data for station", sitename, "..." ) )
  
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

    ##--------------------------------------------------------------------
    ## wet days
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$wetd <- get_pointdata_monthly_cru( "wet", lon_look, lat )

    ##--------------------------------------------------------------------
    ## cloud cover
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    clim_cru_monthly$ccov <- get_pointdata_monthly_cru( "cld", lon_look, lat )

    ##--------------------------------------------------------------------
    ## VPD 
    ## calculated as a function of vapour pressure and temperature, vapour
    ## pressure is given by CRU data.
    ##--------------------------------------------------------------------
    ## Extract data from monthly files and sum up each of these into one 
    ## monthly value which is then to be attached to the monthly dataframe 
    ## (at the right location).
    vap <- get_pointdata_monthly_cru( "vap", lon_look, lat )
    clim_cru_monthly$vpd <- rep( NA, dim(clim_cru_monthly)[1] )
    for (kdx in 1:dim(clim_cru_monthly)[1]){
      clim_cru_monthly$vpd[kdx] <- calc_vpd( vap[kdx], clim_cru_monthly$temp[kdx] )      
    }

    # ##--------------------------------------------------------------------
    # ## irradiance
    # ##--------------------------------------------------------------------
    # clim_cru_monthly$irad <- rep( NA, dim(clim_cru_monthly)[1] )

    ##--------------------------------------------------------------------
    ## expanding to daily data (interpolating temp, and generating prec)
    ##--------------------------------------------------------------------
    print( paste( "expanding to daily data for station", sitename, "..." ) )

    ## Reducing years
    startyr <- max( startyr_override, startyr_cru )
    nyrs    <- length( startyr:endyr_cru )
    clim_cru_monthly <- clim_cru_monthly[ clim_cru_monthly$year >= startyr, ]

    dm   <- rep( NA, sum(ndaymonth)*length(startyr:endyr_cru) )
    jdx <- 0
    for (yr in startyr:endyr_cru ){
      for (imoy in 1:nmonth){
        for (idm in 1:ndaymonth[imoy]){
          jdx <- jdx + 1 
          dm[jdx]   <- idm
        }
      }
    }
    clim_daily <- data.frame( 
      doy=rep( seq(ndayyear), nyrs ), 
      moy=rep( rep( seq(nmonth), times=ndaymonth ), times=nyrs ),
      dom=dm,
      year=rep( startyr:endyr_cru, each=ndayyear ) 
    )

    clim_daily$temp <- rep( NA, dim(clim_daily)[1] )
    clim_daily$prec <- rep( NA, dim(clim_daily)[1] )
    clim_daily$ccov <- rep( NA, dim(clim_daily)[1] )
    clim_daily$vpd  <- rep( NA, dim(clim_daily)[1] )
    # clim_daily$irad <- rep( NA, dim(clim_daily)[1] )

    clim_daily$source_temp <- rep( "NA", dim(clim_daily)[1] )
    clim_daily$source_prec <- rep( "NA", dim(clim_daily)[1] )
    clim_daily$source_ccov <- rep( "NA", dim(clim_daily)[1] )
    clim_daily$source_vpd  <- rep( "NA", dim(clim_daily)[1] )
    # clim_daily$source_irad <- rep( "NA", dim(clim_daily)[1] )

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
      clim_daily$source_temp[ idxs ] <- "CRU monthly interpolated (polynom)"

      ##--------------------------------------------------------------------
      ## precipitation: interpolate weather generator
      ##--------------------------------------------------------------------
      #mprec <- select( filter( clim_cru_monthly, year==use_year ), prec )$prec
      mprec <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$prec
      #mwetd <- select( filter( clim_cru_monthly, year==use_year ), wetd )$wetd
      mwetd <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$wetd
      clim_daily$prec[ idxs ] <- get_daily_prec( mprec, mwetd )
      clim_daily$source_prec[ idxs ] <- "CRU monthly (weather gen. from mon. precip., mon. wetd.)"

      ##--------------------------------------------------------------------
      ## cloud cover: interpolate using polynomial
      ##--------------------------------------------------------------------
      mccov     <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$ccov
      mccov_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$ccov
      mccov_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$ccov
      clim_daily$ccov[ idxs ] <- monthly2daily( mccov, "polynom", mccov_pvy[nmonth], mccov_nxt[1] )
      clim_daily$source_ccov[ idxs ] <- "CRU monthly interpolated (polynom)"

      ##--------------------------------------------------------------------
      ## VPD: interpolate using polynomial
      ##--------------------------------------------------------------------
      mvpd     <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$vpd
      mvpd_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$vpd
      mvpd_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$vpd
      clim_daily$vpd[ idxs ] <- monthly2daily( mvpd, "polynom", mvpd_pvy[nmonth], mvpd_nxt[1] )
      clim_daily$source_vpd[ idxs ] <- "function of CRU vap., monthly interpolated (polynom)"

      ##--------------------------------------------------------------------
      ## irradiance - nothing done yet. Is this identical to WATCH SW-DOWN (+LW-DOWN)? 
      ## (Long-wave downwards surface radiation W/m2 flux (average over previous 3 hours))
      ##--------------------------------------------------------------------

      imo <- imo + nmonth

    }

    ## Reduce CCOV to a maximum 100%
    clim_daily$ccov[ clim_daily$ccov>100.0 ] <- 100.0

    ##--------------------------------------------------------------------
    ## Get daily WFDEI data where available
    ##--------------------------------------------------------------------
    print( "getting daily data from WFDEI ...")
    startyr <- max( startyr_override, staryr_wfdei )
    for ( yr in startyr:endyr_wfdei ){
      for ( moy in seq(nmonth) ){
        print( paste( "... found data for year and month:", yr, moy ) )

        ## temperature
        tmp <- get_pointdata_temp_wfdei( lon, lat, moy, yr )
        if (!is.na(tmp[1])) { 
          useidx <- which( clim_daily$year==yr & clim_daily$moy==moy )
          clim_daily$temp[ useidx ]   <- tmp 
          clim_daily$source_temp[ useidx ] <- "WATCH-WFDEI" 
        }

        ## precipitation
        tmp <- get_pointdata_prec_wfdei( lon, lat, moy, yr )
        if (!is.na(tmp[1])) { 
          useidx <- which( clim_daily$year==yr & clim_daily$moy==moy )
          clim_daily$prec[ useidx ] <- tmp 
          clim_daily$source_prec[ useidx ] <- "WATCH-WFDEI" 
        }

        ## VPD to be calculated from Qair data available through WATCH-WFDEI

        ## PPFD to be calculated from SWdown data available through WATCH-WFDEI

      }
    }

    ##--------------------------------------------------------------------
    ## Save (massive) daily and monthly climate data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing climate data frame into CSV file ", filnam_clim_csv, "..." ) )
    system( paste( "mkdir -p", dirnam_clim_csv ) )   
    write.csv( clim_daily, file=filnam_clim_csv, row.names=FALSE )

  }

}

# for (idx in seq(nsites)){
for (idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_clim_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", sep="" )
  filnam_clim_csv <- paste( dirnam_clim_csv, "clim_daily_", sitename, ".csv", sep="" )

  if ( file.exists( filnam_clim_csv ) ) {
    ##--------------------------------------------------------------------
    ## Read daily data from CSV that may not contain site-specific meteo data yet
    ##--------------------------------------------------------------------
    print( "daily climate data already available in CSV file ...")
    clim_daily <- read.csv( filnam_clim_csv, as.is=TRUE )

    ## Reducing years
    startyr <- max( startyr_override, startyr_cru )
    nyrs    <- length( startyr:endyr_cru )
    clim_daily <- clim_daily[ clim_daily$year >= startyr, ]

    ## for some reason, DOM is not given in some files, add now
    if ( is.null( clim_daily$dom ) ){
      dm   <- rep( NA, dim(clim_daily)[1] )
      jdx <- 0
      for (yr in clim_daily$year[1]:clim_daily$year[dim(clim_daily)[1]] ){
        for (imoy in 1:nmonth){
          for (idm in 1:ndaymonth[imoy]){
            jdx <- jdx + 1 
            dm[jdx]   <- idm
          }
        }
      }
      clim_daily$dom <- dm
    }

    ## if data source column is not yet variable-specific, add info. (column 'source' should be avilable)
    print( "adding source information for each variable ...")
    if ( is.null( clim_daily$source_temp ) ) { clim_daily$source_temp <- clim_daily$source }
    if ( is.null( clim_daily$source_prec ) ) { clim_daily$source_prec <- clim_daily$source }
    if ( is.null( clim_daily$source_vpd  ) ) { clim_daily$source_vpd  <- clim_daily$source }

    ## previous version had 'vapr' and not 'vpd'. convert now.
    print( "calculating VPD ...")
    if ( is.null( clim_daily$vpd ) && !is.null( clim_daily$vapr )) { 
      clim_daily$vpd <- rep( NA, dim(clim_daily)[1] )
      for (kdx in 1:dim(clim_daily)[1]){
        clim_daily$vpd[kdx] <- calc_vpd( clim_daily$vapr[kdx], clim_daily$temp[kdx] )      
      }
    }      
  }

  if ( ingest_meteodata ){
    ##--------------------------------------------------------------------
    ## Get site-specific meto data for separate (site-specific) file
    ##--------------------------------------------------------------------
    print(paste("Ingesting meteo data"))
    if ( !is.null( siteinfo$meteosource ) || simsuite=="fluxnet2015" ){
      print("Using site-specific meteo data from separate file ...")
      if ( simsuite=="fluxnet2015" ){

        dirnam_obs <- paste( myhome, "data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/", sep="" )
        allfiles <- list.files( dirnam_obs )
        filnam_obs <- allfiles[ which( grepl( sitename, allfiles ) ) ]
        filn <- paste( dirnam_obs, filnam_obs, sep="" )
        if ( length(filnam_obs)>0 ){ 
          found <- TRUE 
          meteo <- try( get_meteo_fluxnet2015( filn ) ) 
        } else { 
          found <- FALSE 
        }

      } else {

        filn <- paste( myhome, as.character(siteinfo$meteosource[idx] ), sep="" )
        meteo <- read.csv( filn, as.is=TRUE )
        found <- TRUE

      }

      if ( found ){

        ## add MOY and DOM
        if ( is.null(meteo$moy) ) { meteo$moy <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mon + 1 }
        if ( is.null(meteo$dom) ) { meteo$dom <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mday }

        ## if site data goes beyond 2014, add rows
        if ( max(clim_daily$year) < max(meteo$year) ){

          print( paste( "adding years up to ", max(meteo$year), "..." ) )

          # colnames(clim_daily)[5:length(colnames(clim_daily))]

          ## Reducing years
          startyr <- max(clim_daily$year) + 1
          endyr   <- max(meteo$year)
          nyrs    <- length( startyr:endyr )

          dm  <- rep( NA, sum(ndaymonth)*length(startyr:endyr) )
          jdx <- 0
          for (yr in startyr:max(meteo$year) ){
            for (imoy in 1:nmonth){
              for (idm in 1:ndaymonth[imoy]){
                jdx <- jdx + 1 
                dm[jdx]   <- idm
              }
            }
          }
          addrows <- data.frame( 
            doy         =rep( seq(ndayyear), nyrs ), 
            moy         =rep( rep( seq(nmonth), times=ndaymonth ), times=nyrs ),
            dom         =dm,
            year        =rep( startyr:max(meteo$year), each=ndayyear ),
            temp        =rep( NA, nyrs*ndayyear ),
            prec        =rep( NA, nyrs*ndayyear ),
            ccov        =rep( NA, nyrs*ndayyear ),
            vapr        =rep( NA, nyrs*ndayyear ),
            source      =rep( NA, nyrs*ndayyear ),
            source_temp =rep( NA, nyrs*ndayyear ),
            source_prec =rep( NA, nyrs*ndayyear ),
            source_vpd  =rep( NA, nyrs*ndayyear ),
            vpd         =rep( NA, nyrs*ndayyear )
          )

          ## add rows to clim_daily
          addrows[,5:dim(clim_daily)[2]] <- NA

          ## add additional rows to clim_daily
          clim_daily <- rbind( clim_daily, addrows )

        }

        for (jdx in 1:dim(meteo)[1]){

          putjdx <- which( clim_daily$year==meteo$year[jdx] & clim_daily$moy==meteo$moy[jdx] & clim_daily$dom==meteo$dom[jdx] )

          if (length(putjdx)>1) { 
            print("PROBLEM: found multiple corresponding indeces for ...")
            print(paste("year =", meteo$year[jdx], "moy =", meteo$moy[jdx], "dom =", meteo$dom[jdx] ) ) 
          
          } else if (length(putjdx)==0){
            print("PROBLEM: found no corresponding indeces for ...")
            print(paste("year =", meteo$year[jdx], "moy =", meteo$moy[jdx], "dom =", meteo$dom[jdx] ) )           
          
          } else {

            ## temperature
            if (!is.na(meteo$temp[jdx])) { clim_daily$temp[ putjdx ] <- meteo$temp[jdx] }
            if (!is.na(meteo$temp[jdx])) { clim_daily$source_temp[ putjdx ] <- "temp. sitedata" }

            ## precipitation
            if (!is.na(meteo$prec[jdx])) { clim_daily$prec[ putjdx ] <- meteo$prec[jdx] }
            if (!is.na(meteo$prec[jdx])) { clim_daily$source_prec[ putjdx ] <- "prec. sitedata" }

            ## VPD
            if (!is.na(meteo$vpd[jdx])) { clim_daily$vpd[ putjdx ] <- meteo$vpd[jdx] }
            if (!is.na(meteo$vpd[jdx])) { clim_daily$source_vpd[ putjdx ] <- "VPD sitedata" }

          }

        }

      } else {

        print( paste( "WARNING: did not find file", filn ) )

      }

    }

  }

  ##--------------------------------------------------------------------
  ## Save (massive) daily and monthly climate data frames as CSV files
  ##--------------------------------------------------------------------
  print( paste( "writing climate data frame into CSV file ", filnam_clim_csv, "..." ) )
  system( paste( "mkdir -p", dirnam_clim_csv ) )   
  write.csv( clim_daily, file=filnam_clim_csv, row.names=FALSE )

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

    filnam <- paste( dirnam, "dvpd_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, clim_daily$vpd[ which( clim_daily$year==yr ) ] )

    # filnam <- paste( dirnam, "dirad_", sitename, "_", yr, ".txt", sep="" )
    # write_sofunformatted( filnam, clim_daily$irad[ which( clim_daily$year==yr ) ] )

  }

}

