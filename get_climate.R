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
## - net radiation (J m-2 d-1)
## - shortwave incoming radiation, swin: J m-2 d-1
##--------------------------------------------------------------------
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
source( paste( myhome, "sofun/getin/get_meteo_fluxnet2015.R", sep="" ) )
source( paste( myhome, "sofun/getin/get_meteo_swbm_meteoschweiz.R", sep="" ) )
source( paste( myhome, "sofun/getin/calc_netrad_orth.R", sep="" ) )
source( paste( myhome, "sofun/getin/init_daily_dataframe.R", sep="" ) )

overwrite <- TRUE
overwrite_byst <- TRUE

simsuite <- "fluxnet2015"

in_ppfd   <- TRUE
in_netrad <- FALSE

startyr_override <- 1990

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
startyr_cru <- 1901
endyr_cru   <- 2014
staryr_wfdei <- 1979
endyr_wfdei  <- 2012

## before knowing whether site data is available beyond CRU, assume...
endyr_act <- endyr_cru
startyr_act <- startyr_override

## load meta data file for site simulation
siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ), as.is=TRUE )
nsites <- dim(siteinfo)[1]
# do.sites <- seq(nsites)
# do.sites <- 21:nsites
do.sites <- 1:1

get_clim_cru_monthly <- function( lon, lat, startyr_cru, endyr_cru ){

  ##--------------------------------------------------------------------
  ## get monthly data from CRU
  ##--------------------------------------------------------------------
  nyrs_cru <- length(startyr_cru:endyr_cru)
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

  return( clim_cru_monthly )
}

expand_clim_cru_monthly <- function( clim_cru_monthly ){

  len <- range( clim_cru_monthly$year )
  dm   <- rep( NA, sum(ndaymonth)*length(range(clim_cru_monthly$year)[1]:range(clim_cru_monthly$year)[2]) )
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

  clim_daily$temp_cru_int <- rep( NA, dim(clim_daily)[1] )
  clim_daily$prec_cru_gen <- rep( NA, dim(clim_daily)[1] )
  clim_daily$ccov_cru_int <- rep( NA, dim(clim_daily)[1] )
  clim_daily$vpd_vap_cru_temp_cru_int  <- rep( NA, dim(clim_daily)[1] )
  # clim_daily$irad <- rep( NA, dim(clim_daily)[1] )

  imo <- 1
  for (iyr in seq(dim(clim_cru_monthly)[1]/nmonth) ){

    use_year <- clim_cru_monthly$year[imo]
    use_year_pvy <- max(startyr_cru,use_year-1)
    use_year_nxt <- min(endyr_cru,use_year+1)
    idxs <- which( clim_daily$year==use_year )
    
    ##--------------------------------------------------------------------
    ## air temperature: interpolate using polynomial
    ##--------------------------------------------------------------------
    mtemp <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$temp
    mtemp_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$temp
    mtemp_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$temp
    if (length(mtemp_pvy)==0){
      mtemp_pvy <- mtemp
    }
    if (length(mtemp_nxt)==0){
      mtemp_nxt <- mtemp
    }
    clim_daily$temp_cru_int[ idxs ] <- monthly2daily( mtemp, "polynom", mtemp_pvy[nmonth], mtemp_nxt[1] )

    ##--------------------------------------------------------------------
    ## precipitation: interpolate weather generator
    ##--------------------------------------------------------------------
    mprec <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$prec
    mwetd <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$wetd
    clim_daily$prec_cru_gen[ idxs ] <- get_daily_prec( mprec, mwetd )

    ##--------------------------------------------------------------------
    ## cloud cover: interpolate using polynomial
    ##--------------------------------------------------------------------
    mccov     <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$ccov
    mccov_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$ccov
    mccov_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$ccov
    if (length(mccov_pvy)==0){
      mccov_pvy <- mccov
    }
    if (length(mccov_nxt)==0){
      mccov_nxt <- mccov
    }
    clim_daily$ccov_cru_int[ idxs ] <- monthly2daily( mccov, "polynom", mccov_pvy[nmonth], mccov_nxt[1] )

    ## Reduce CCOV to a maximum 100%
    clim_daily$ccov_cru_int[ clim_daily$ccov_cru_int>100.0 ] <- 100.0

    ##--------------------------------------------------------------------
    ## VPD: interpolate using polynomial
    ##--------------------------------------------------------------------
    mvpd     <- clim_cru_monthly[ clim_cru_monthly$year==use_year, ]$vpd
    mvpd_pvy <- clim_cru_monthly[ clim_cru_monthly$year==use_year_pvy, ]$vpd
    mvpd_nxt <- clim_cru_monthly[ clim_cru_monthly$year==use_year_nxt, ]$vpd
    if (length(mvpd_pvy)==0){
      mvpd_pvy <- mvpd
    }
    if (length(mvpd_nxt)==0){
      mvpd_nxt <- mvpd
    }
    clim_daily$vpd_vap_cru_temp_cru_int[ idxs ] <- monthly2daily( mvpd, "polynom", mvpd_pvy[nmonth], mvpd_nxt[1] )

    ##--------------------------------------------------------------------
    ## irradiance - nothing done yet. Is this identical to WATCH SW-DOWN (+LW-DOWN)? 
    ## (Long-wave downwards surface radiation W/m2 flux (average over previous 3 hours))
    ##--------------------------------------------------------------------
    imo <- imo + nmonth

  }
  return( clim_daily )
}

add_watch <- function( clim_daily, lon, lat, startyr, endyr ){

  clim_daily$temp_watch <- rep( NA, dim(clim_daily)[1] )
  clim_daily$prec_watch <- rep( NA, dim(clim_daily)[1] )

  for ( yr in startyr:endyr ){
    for ( moy in seq(nmonth) ){
      print( paste( "... found data for year and month:", yr, moy ) )

      ## temperature
      tmp <- get_pointdata_temp_wfdei( lon, lat, moy, yr )
      if (!is.na(tmp[1])) { 
        useidx <- which( clim_daily$year==yr & clim_daily$moy==moy )
        clim_daily$temp_watch[ useidx ]   <- tmp 
      }

      ## precipitation
      tmp <- get_pointdata_prec_wfdei( lon, lat, moy, yr )
      if (!is.na(tmp[1])) { 
        useidx <- which( clim_daily$year==yr & clim_daily$moy==moy )
        clim_daily$prec_watch[ useidx ] <- tmp 
      }

      ## VPD to be calculated from Qair data available through WATCH-WFDEI

      ## PPFD to be calculated from SWdown data available through WATCH-WFDEI

    }
  }
  return( clim_daily )
}

addrows_to_clim_daily <- function( clim_daily, endyr_bysit ){

  print( paste( "adding years up to ", endyr_bysit, "..." ) )

  ## Reducing years
  startyr <- max(clim_daily$year) + 1
  endyr   <- endyr_bysit

  addrows             <- init_daily_dataframe <- function( startyr, endyr )

  data_in_clim_daily <- names( clim_daily )[5:ncol(clim_daily)]
  for (addcol in data_in_clim_daily){
    addrows[[ addcol ]] <- rep( NA, dim(addrows)[1] )
  }

  ## add additional rows to clim_daily
  clim_daily <- rbind( clim_daily, addrows )

  return( clim_daily )
}

ingest_meteodata <- function( clim_daily, meteo, in_netrad, in_ppfd ){

  ## add MOY and DOM to meteo data frame for association
  if ( is.null(meteo$moy) ) { meteo$moy <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mon + 1 }
  if ( is.null(meteo$dom) ) { meteo$dom <- as.POSIXlt( meteo$date, format="%d/%m/%Y" )$mday }

  ## add columns
  clim_daily$temp_meteo <- rep( NA, dim(clim_daily)[1]) 
  clim_daily$prec_meteo <- rep( NA, dim(clim_daily)[1]) 
  clim_daily$vpd_meteo  <- rep( NA, dim(clim_daily)[1]) 
  
  ## add column for netrad in dataframe clim_daily
  if ( in_netrad && !is.null( meteo$rad ) ){ 
    clim_daily$netrad_meteo <- rep( NA, dim(clim_daily)[1]) 
  }

  ## add column for ppfd in dataframe clim_daily
  if ( in_ppfd && !is.null( meteo$ppfd ) ){ 
    clim_daily$ppfd_sw_meteo <- rep( NA, dim(clim_daily)[1]) 
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
      if (!is.null(meteo$temp[jdx])) { clim_daily$temp_meteo[ putjdx ] <- meteo$temp[jdx] }

      ## precipitation
      if (!is.null(meteo$prec[jdx])) { clim_daily$prec_meteo[ putjdx ] <- meteo$prec[jdx] }

      ## VPD
      if (!is.null(meteo$vpd[jdx])) { clim_daily$vpd_meteo[ putjdx ] <- meteo$vpd[jdx] }

      ## PPFD
      if ( in_ppfd && !is.null( meteo$ppfd ) ){
        if (!is.null(meteo$ppfd[jdx])) { clim_daily$ppfd_sw_meteo[ putjdx ] <- meteo$ppfd[jdx] }
      }

      ## Net radiation
      if ( in_netrad && !is.null( meteo$rad ) ){
        if (!is.null(meteo$rad[jdx])) { clim_daily$netrad_meteo[ putjdx ] <- calc_netrad_orth(  meteo$rad[jdx] ) }
      }

    }

  }
  return( clim_daily )
}

fill_gaps_clim_daily <- function( clim_daily ){

  # print("filling gaps ...")
  
  cols <- names( clim_daily )[5:ncol( clim_daily )]

  for (icol in cols){

    if ( grepl( "prec", icol ) ){

      if (any(is.na(clim_daily[[ icol ]]))){
        clim_daily[[ icol ]][ is.na(clim_daily[[ icol ]]) ] <- 0.0
      }

    } else {

      if ( any(is.na(clim_daily[[ icol ]])) && any(!is.na(clim_daily[[ icol ]])) ){
        # print( paste( "found gaps for", icol ) )
        year_dec <- init_daily_dataframe( clim_daily$year[1], clim_daily$year[dim(clim_daily)[1]] )$year_dec
        if (dim(clim_daily)[1]==length(year_dec)){
          # print("using approx() function")
          clim_daily[[ icol ]] <- approx( year_dec, clim_daily[[ icol ]], xout=year_dec )$y
        }
      }

      if ( any(is.na(clim_daily[[ icol ]]))  && any(!is.na(clim_daily[[ icol ]])) ){
        # print( paste( "found gaps for", icol ) )
        ## if there are NAs at the end of the columns, fill with last non-NA value
        # print("XXXXX filling gaps at the tail")
        # print( icol )
        # print("before:")
        # print( head( clim_daily[[ icol ]] ) )
        for (idx in 1:nrow(clim_daily)){
          if ( any( is.na( tail( clim_daily[[ icol ]], n=idx ) ) ) && any( !is.na( tail( clim_daily[[ icol ]], n=(idx+1) ) ) ) ){
            clim_daily[[ icol ]][ (nrow(clim_daily)-idx+1):nrow(clim_daily) ] <- clim_daily[[ icol ]][ (nrow(clim_daily)-idx) ]
            break
          }
        }
        # print("after:")
        # print( head( clim_daily[[ icol ]] ) )
      }  

    }

  }
  return( clim_daily )
}

##--------------------------------------------------------------------
## First loop over stations: get monthly CRU data and daily WATCH-WFDEI
##--------------------------------------------------------------------
print("==============================================================")
print("GET_CLIMATE: getting data from NetCDF files ...")
print("==============================================================")
for ( idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste("site name", sitename))

  dirnam_clim_csv      <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", sep="" )
  filnam_clim_csv      <- paste( dirnam_clim_csv, "clim_daily_", sitename, ".csv", sep="" )  # NOT containing station-specific meteo data
  filnam_clim_byst_csv <- paste( dirnam_clim_csv, "clim_daily_byst_", sitename, ".csv", sep="" )

  if ( overwrite || ( !file.exists( filnam_clim_csv ) && !file.exists(filnam_clim_byst_csv) )){
    ##--------------------------------------------------------------------
    ## Get monthly data frame covering all CRU years
    ##--------------------------------------------------------------------
    print( paste( "collecting monthly data for station", sitename, "..." ) )
    clim_cru_monthly <- get_clim_cru_monthly( lon, lat, startyr_cru, endyr_cru )

    ##--------------------------------------------------------------------
    ## expanding to daily data (interpolating temp, and generating prec)
    ##--------------------------------------------------------------------
    print( paste( "expanding to daily data for station", sitename, "..." ) )

    ## Reducing years
    startyr <- max( startyr_override, startyr_cru )
    nyrs    <- length( startyr:endyr_cru )
    clim_cru_monthly <- clim_cru_monthly[ clim_cru_monthly$year >= startyr, ]

    clim_daily <- expand_clim_cru_monthly( clim_cru_monthly )

    ##--------------------------------------------------------------------
    ## Get daily WFDEI data where available
    ##--------------------------------------------------------------------
    print( "getting daily data from WFDEI ...")
    startyr <- max( startyr_override, staryr_wfdei )
    endyr   <- endyr_wfdei
    clim_daily <- add_watch( clim_daily, lon, lat, startyr, endyr )

    ##--------------------------------------------------------------------
    ## Save daily without station-specific data to CSV file
    ##--------------------------------------------------------------------
    print( paste( "writing climate data frame with meteo data into CSV file ", filnam_clim_csv, "..." ) )
    system( paste( "mkdir -p", dirnam_clim_csv ) )   
    write.csv( clim_daily, file=filnam_clim_csv, row.names=FALSE )

  } else {
    ##--------------------------------------------------------------------
    ## Read daily data from CSV 
    ##--------------------------------------------------------------------
    print("reading from available CSV file that has no site-specific data yet")
    if (file.exists( filnam_clim_csv )){
      clim_daily <- read.csv( filnam_clim_csv, as.is=TRUE )
    } else if (file.exists( filnam_clim_byst_csv )){
      clim_daily <- read.csv( filnam_clim_byst_csv, as.is=TRUE )      
    } else {
      print("I'm in the wrong place.")
    }

    # ## XXX das ist nur provisorisch>>>>>>>>>>
    # library( dplyr )
    # if (!is.null(clim_daily$temp) ){
    #   clim_daily$temp_watch <- clim_daily$temp
    # }
    # if ( !is.null(clim_daily$prec) ){
    #   clim_daily$prec_watch <- clim_daily$prec
    # }
    # if (!is.null(clim_daily$ccov) ){
    #   clim_daily$ccov_cru_int <- clim_daily$ccov
    # }
    # if ( !is.null(clim_daily$vpd) ){
    #   clim_daily$vpd_vap_cru_temp_cru_int <- clim_daily$vpd
    # }
    # clim_daily <- select( clim_daily, doy, moy, dom, year, temp_watch, prec_watch, ccov_cru_int, vpd_vap_cru_temp_cru_int )      
    # ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  if ( overwrite_byst || !file.exists( filnam_clim_byst_csv )){
    ##--------------------------------------------------------------------
    ## Get site-specific meto data for separate (site-specific) file
    ##--------------------------------------------------------------------
    if ( !is.null( siteinfo$meteosource ) || simsuite=="fluxnet2015" ){
  
      print("Using site-specific meteo data from separate file ...")
  
      ##--------------------------------------------------------------------
      ## Get separate data frame called 'meteo'
      ##--------------------------------------------------------------------
      if ( simsuite=="fluxnet2015" ){

        ## FLUXNET 2015 METEO DATA
        dirnam_obs <- paste( myhome, "data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/", sep="" )
        allfiles <- list.files( dirnam_obs )
        allfiles <- allfiles[ which( grepl( "FULLSET", allfiles ) ) ]
        filnam_obs <- allfiles[ which( grepl( sitename, allfiles ) ) ]
        filn <- paste( dirnam_obs, filnam_obs, sep="" )
        if ( length(filnam_obs)>0 ){ 
          found <- TRUE 
          meteo <- try( get_meteo_fluxnet2015( filn ) ) 
        } else { 
          found <- FALSE 
        }

      } else if ( simsuite=="swbm" ) {

        ## METEOSCHWEIZ DATA (for station Payerne, SWBM simsuite only)
        meteo <- get_meteo_swbm_meteoschweiz( paste( myhome, as.character(siteinfo$meteosource[idx] ), sep="" ) )
        found <- TRUE

      } else {

        ## OTHER METEO DATA
        filn <- paste( myhome, as.character(siteinfo$meteosource[idx] ), sep="" )
        meteo <- read.csv( filn, as.is=TRUE )
        found <- TRUE

      }

      ##--------------------------------------------------------------------
      ## Add data from 'meteo' to dataframe
      ##--------------------------------------------------------------------
      if ( found ){

        ## if site data goes beyond 2014 (maximum year of clim_daily), add rows
        endyr_bysit <- max( meteo$year )
        if ( max(clim_daily$year) < endyr_bysit ){
          clim_daily <- addrows_to_clim_daily( clim_daily, endyr_bysit )
        }

        ## Ingest meteodata        
        clim_daily <- ingest_meteodata( clim_daily, meteo, in_netrad, in_ppfd )

      } else {

        print( paste( "WARNING: did not find file", filn ) )

      }

      ##--------------------------------------------------------------------
      ## Save (massive) daily now containing station-specific data to CSV file
      ##--------------------------------------------------------------------
      print( paste( "writing climate data frame with meteo data into CSV file ", filnam_clim_byst_csv, "..." ) )
      system( paste( "mkdir -p", dirnam_clim_csv ) )   
      write.csv( clim_daily, file=filnam_clim_byst_csv, row.names=FALSE )

    } else {

      print("What meteo data should I use?")

    }

  } else {

    ##--------------------------------------------------------------------
    ## Read daily data from CSV that has station-specific meteo data
    ##--------------------------------------------------------------------
    clim_daily <- read.csv( filnam_clim_byst_csv, as.is=TRUE )

  }

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )
  startyr_act <- min( clim_daily$year )
  endyr_act   <- max( clim_daily$year )

  ## Get simulation years (equal to years in meteo data if available)
  if (found){
    simyears <- unique(meteo$year)
  } else {
    simyears <- startyr_act:endyr_act
  }

  ## select priority of sources for SOFUN-formatted files:
  out <- subset( clim_daily, select=c( doy, moy, dom, year ) )

  ## Delete directories that are not needed anymore
  for (yr in startyr_act:endyr_act){
    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "rm -rf", dirnam ) )
  }

  ## temperature: 1. meteo, 2. watch, 3. cru_int
  out$temp <- clim_daily$temp_meteo
  out$temp[ which( is.na( out$temp ) ) ] <- clim_daily$temp_watch[ which( is.na( out$temp ) ) ]
  if (!is.null(clim_daily$temp_cru_int)){
    out$temp[ which( is.na( out$temp ) ) ] <- clim_daily$temp_cru_int[ which( is.na( out$temp ) ) ]
  }

  ## precipitation: 1. meteo, 2. watch, 3. cru_gen
  out$prec <- clim_daily$prec_meteo
  out$prec[ which( is.na( out$prec ) ) ] <- clim_daily$prec_watch[ which( is.na( out$prec ) ) ]
  out$prec[ which( is.na( out$prec ) ) ] <- clim_daily$prec_cru_gen[ which( is.na( out$prec ) ) ]

  ## cloud cover fraction: cru_int
  out$ccov <- clim_daily$ccov_cru_int

  ## VPD
  out$vpd <- clim_daily$vpd_meteo
  out$vpd[ which( is.na( out$vpd ) ) ] <- clim_daily$vpd_vap_cru_temp_cru_int[ which( is.na( out$vpd ) ) ]

  ## PPFD
  if (in_ppfd){
    out$ppfd <- clim_daily$ppfd_sw_meteo 
    ## XXX use SWdown from watch data as second priority
  }

  ## Net radiation
  if (in_netrad){
    out$netrad <- clim_daily$netrad
  }

  ## reduce data frame
  out <- out[ is.element(out$year, simyears), ]

  ## write SOFUN-formatted input files
  for ( yr in simyears ){

    ## Subsetting and gap filling: interpolate linearly (temp and vpd) or set to zero (prec)
    tmp <- out[ which( out$year==yr ), ]
    tmp <- fill_gaps_clim_daily( tmp )

    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/climate/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "dtemp_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, tmp$temp )
    
    filnam <- paste( dirnam, "dprec_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, tmp$prec )

    filnam <- paste( dirnam, "dfsun_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, ( 100.0 - tmp$ccov ) / 100.0 )

    filnam <- paste( dirnam, "dvpd_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, tmp$vpd )

    if (in_netrad){
      filnam <- paste( dirnam, "dnetrad_", sitename, "_", yr, ".txt", sep="" )
      write_sofunformatted( filnam, tmp$netrad )
    }

    if (in_ppfd){
      filnam <- paste( dirnam, "dppfd_", sitename, "_", yr, ".txt", sep="" )
      write_sofunformatted( filnam, tmp$ppfd )
    }

  }

  print("... done")

}


