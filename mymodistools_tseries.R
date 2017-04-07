download_subset_modis_tseries <- function( lon, lat, bands, prod, sitename, TimeSeriesLength, end.date, savedir, overwrite, ignore_missing=FALSE ){

  # ######################
  # end.date <- Sys.Date()
  # savedir <- paste( "/alphadata01/bstocker/data/modis_gpp_fluxnet_cutouts_tseries/AR-SLu/", sep="" )
  # overwrite <- FALSE
  # # TimeSeriesLength <- 1

  # sitename <- "AR-SLu"
  # lon <- -66.4598
  # lat <- -33.4648
  # ######################

  require( MODISTools )

  ## Find file from which (crude) data is read
  filn <- list.files( path=savedir, pattern="*asc" )

  # test <- try( read.table(filn) )
  # if (class(test)=="try-error"){
  #   system( paste( "rm ", savedir, filn, sep="" ) )
  #   filn <- list.files( path=savedir, pattern="*asc" )    
  # }

  ## Account for different resolutions of the data
  if (prod=="MOD13Q1"){
    ## at 250 m, therefore use 1 km surrounding centre, that is 9x9=81 pixels
    size=c(1,1)
  } else if (prod=="MOD15A2"){
    ## at 1 km, therefore use 3 km surrounding centre, that is 9x9=81 pixels
    size=c(3,3)    
  }

  if ( (length(filn)==0||overwrite) && !ignore_missing ){

    print( paste( "deleting existing file", filn ) )
    if (length(filn)>0) { system( paste( "rm ", savedir, filn, sep="" ) ) }
    
    print( paste( "==========================="))
    print( paste( "DOWNLOADING MODIS DATA FOR:"))
    print( paste( "site :", sitename ) )
    print( paste( "lon  :", lon ) )
    print( paste( "lat  :", lat ) )
    print( paste( "---------------------------"))

    # system( paste( "mkdir ", savedir, sep="" ) )  

    modis.subset <- data.frame(  
      lat        = lat, 
      long       = lon, 
      end.date   = Sys.Date() #  end.date
      )

    ## Get subset for all pixels around centre (defined by input lon/lat), stretching 1km in lon and lat
    try( MODISSubsets( 
      LoadDat          = modis.subset,
      Products         = prod,
      Bands            = bands,
      Size             = size,    ## Get subset for all pixels around centre (defined by input lon/lat), stretching 1km in lon and lat (should be 9 in each direction = 81 in total)
      StartDate        = FALSE,
      SaveDir          = savedir,
      TimeSeriesLength = TimeSeriesLength       ## 26 years beyond first day of the year in 'end.date'. Will trigger reading longet possible automatially 
      # TimeSeriesLength = 1       ## 26 years beyond first day of the year in 'end.date'. Will trigger reading longet possible automatially 
      ) )

  } else {

    if (ignore_missing && length(filn)==0){

      print( paste( "==========================="))
      print( paste( "IGNORING MISSING DATA FOR:"))
      print( paste( "lon  :", lon) )
      print( paste( "lat  :", lat) )
      print( paste( "end  :", end.date ) )
      print( paste( "---------------------------"))


    } else {

      print( paste( "==========================="))
      print( paste( "FOUND DATA FOR:"))
      print( paste( "lon  :", lon) )
      print( paste( "lat  :", lat) )
      print( paste( "end  :", end.date ) )
      print( paste( "---------------------------"))

    }


  }

  return( filn[1] )

}

download_subset_modis <- function( lon, lat, bands, prod, start.date, end.date, savedir, overwrite, ignore_missing=FALSE ){
  #///////////////////////////////////////////////////////////////
  # Downloads a MODIS subset for the specified location (lon, lat)
  # and date (start.date, end.date)
  #---------------------------------------------------------------

  require( MODISTools )

  ## Find file from which (crude) data is read
  filn <- list.files( path=savedir, pattern="*asc" )

  # ## Deleting if file is erroneous
  # test <- try( read.table(filn) )
  # if (class(test)=="try-error"){
  #   system( paste( "rm ", savedir, filn, sep="" ) )
  #   filn <- list.files( path=savedir, pattern="*asc" )    
  # }

  ## Account for different resolutions of the data
  if (prod=="MOD13Q1"){
    ## at 250 m, therefore use 1 km surrounding centre, that is 9x9=81 pixels
    size=c(1,1)
  } else if (prod=="MOD15A2"){
    ## at 1 km, therefore use 3 km surrounding centre, that is 9x9=81 pixels
    size=c(3,3)    
  }

  if ( (length(filn)==0||overwrite) && !ignore_missing ){

    print( paste( "deleting existing file", filn ) )
    if (length(filn)>0) { system( paste( "rm ", savedir, filn, sep="" ) ) }

    print(paste("savedir", savedir))
    
    print( paste( "==========================="))
    print( paste( "DOWNLOADING MODIS DATA FOR:"))
    print( paste( "bands:", bands ) )
    print( paste( "prod :", prod  ) )
    print( paste( "site :", sitename ) )
    print( paste( "lon  :", lon ) )
    print( paste( "lat  :", lat ) )
    print( paste( "start:", start.date ) )
    print( paste( "end  :", end.date ) )
    print( paste( "---------------------------"))

    modis.subset <- data.frame(  
      lat        = lat, 
      long       = lon, 
      start.date = as.Date(start.date), 
      end.date   = as.Date(end.date)
      )

    MODISSubsets(
      LoadDat   = modis.subset, 
      Products  = prod,
      Bands     = bands,
      Size      = size,   # c(1,1) to get all pixels +/- 1 km around the centre = 81 pixels in total; c(0,0) to get only centre pixel
      StartDate = TRUE,
      SaveDir   = savedir
      )

    filn <- list.files( path=savedir, pattern="*asc" )
    print(paste("raw data file name:", filn))

  } else {

    if (ignore_missing && length(filn)==0){

      print( paste( "==========================="))
      print( paste( "IGNORING MISSING DATA FOR:"))
      print( paste( "lon  :", lon) )
      print( paste( "lat  :", lat) )
      print( paste( "start:", start.date ) )
      print( paste( "end  :", end.date ) )
      print( paste( "---------------------------"))


    } else {

      print( paste( "==========================="))
      print( paste( "FOUND DATA FOR:"))
      print( paste( "lon  :", lon) )
      print( paste( "lat  :", lat) )
      print( paste( "start:", start.date ) )
      print( paste( "end  :", end.date ) )
      print( paste( "---------------------------"))

    }


  }

  return( filn[1] )

}

read_crude_modis_tseries <- function( filn, savedir, band_var, band_qc=NA, expand_x, expand_y ){
  #///////////////////////////////////////////////////////////////
  # Reads MODIS data from downloaded ASCII file and returns
  # value, quality flag and associated date
  # arguments:
  # filn: file name of ASCII file holding MODIS "crude" data
  # savedir: directory, where to look for that file
  # expand_x : number of pixels to the right and left of centre
  # expand_y : number of pixels to the top and bottom of centre

  # Rank Key
  # Summary QA
  # Description
  # -1  Fill/No Data  Not Processed
  # 0 Good Data Use with confidence
  # 1 Marginal data Useful, but look at other QA information
  # 2 Snow/Ice  Target covered with snow/ice
  # 3 Cloudy  Target not visible, covered with cloud
  #---------------------------------------------------------------

  require( MODISTools )
  require( dplyr )

  # ######################
  # # for debugging:
  # savedir <- "/alphadata01/bstocker/data/modis_fapar_fluxnet_cutouts_tseries/FR-Pue/raw/"
  # filn    <- "Lat43.74140Lon3.59580Start1991-01-01End2017-02-21___MOD15A2.asc"
  # expand_x <- 1
  # expand_y <- 1
  # band_var <- "Fpar_1km"
  # band_qc  <- "FparLai_QC"
  # ######################

  if (band_var=="Fpar_1km"){
    ScaleFactor <- 1e-2
  } else if (band_var=="250m_16_days_EVI"){
    ScaleFactor <- 1e-4  # applied to output variable, in ascii file EVI value is multiplied by 1e4
  }
  ndayyear    <- 365

  ## Read dowloaded ASCII file
  path <- paste( savedir, filn, sep="" )
  print( paste( "reading file ", path ) )
 
  if (file.exists(path)){

    ## this is just read to get length of time series and dates
    tseries    <- as.data.frame( MODISTimeSeries( savedir, Band = band_var, Simplify=TRUE ) )   
    if (!is.na(band_qc)) { tseries_qc <- as.data.frame( MODISTimeSeries( savedir, Band = band_qc,  Simplify=TRUE ) )}

    ## determine centre pixel's column number in 'tseries'
    usecol <- ( ncol(tseries) - 1 ) / 2 + 1 
    npixels <- ncol( tseries )

    ## save centre pixel's info in separate column
    tseries$centre       <- tseries[,usecol] * ScaleFactor
    if (!is.na(band_qc)) { tseries_qc$centre_qc <- tseries_qc[,usecol] }

    ## Get mean of surrounding pixels and save in separate column. Take into account the quality flags of surrounding pixels.
    if (ncol(tseries)>2 && !is.na(band_qc)){

      ## Replace data points with quality flag = 2 (snow covered) by minimum of all pixel's values 
      tseries[ which( tseries_qc[,1:npixels] == 2 ), 1:npixels ] <- min( tseries[,1:npixels], na.rm=TRUE )
      tseries[ which( tseries_qc[,1:npixels] < 0.0 ) ] <- 0.0  ## remove negatives

      # ## Drop all data with quality flag != 0
      # tseries[ which( tseries_qc[,1:npixels] == 3 ),  1:npixels ] <- NA   # Target not visible, covered with cloud
      # tseries[ which( tseries_qc[,1:npixels] == -1 ), 1:npixels ] <- NA   # Not Processed

      tseries$centre_meansurr <- unname( apply( tseries[,1:npixels], 1, FUN=mean, na.rm=TRUE ) ) * ScaleFactor
    }

    ## add time information for tseries
    tmp              <- rownames( tseries )
    tseries$year     <- as.numeric( substr( tmp, start=2, stop=5 ))
    tseries$doy      <- as.numeric( substr( tmp, start=6, stop=8 ))
    tseries$date     <- as.POSIXct( as.Date( paste( as.character(tseries$year), "-01-01", sep="" ) ) + tseries$doy - 1 )
    tseries$year_dec <- tseries$year + ( tseries$doy - 1 ) / ndayyear

    ## add time information for tseries_qc
    if (!is.na(band_qc)){
      tmp                 <- rownames( tseries_qc )
      tseries_qc$year     <- as.numeric( substr( tmp, start=2, stop=5 ))
      tseries_qc$doy      <- as.numeric( substr( tmp, start=6, stop=8 ))
      tseries_qc$date     <- as.POSIXct( as.Date( paste( as.character(tseries_qc$year), "-01-01", sep="" ) ) + tseries_qc$doy - 1 )
      tseries_qc$year_dec <- tseries_qc$year + ( tseries_qc$doy - 1 ) / ndayyear

      ## Merge data frames
      tseries <- tseries %>% left_join( tseries_qc, by=c("date","year","doy","year_dec") )
      
    }

  } else {

    print("ERROR: raw file not found. Was looking for")
    print( path )

  }

  return( tseries )

}

read_crude_modis_tstep <- function( filn, savedir, band_var, band_qc, expand_x, expand_y ){
  #///////////////////////////////////////////////////////////////
  # Reads MODIS data from downloaded ASCII file and returns
  # value, quality flag and associated date
  # arguments:
  # filn: file name of ASCII file holding MODIS "crude" data
  # savedir: directory, where to look for that file
  # expand_x : number of pixels to the right and left of centre
  # expand_y : number of pixels to the top and bottom of centre
  # 
  # Quality flag codes:
  # -1  Fill/No Data  Not Processed
  # 0 Good Data Use with confidence
  # 1 Marginal data Useful, but look at other QA information
  # 2 Snow/Ice  Target covered with snow/ice
  # 3 Cloudy  Target not visible, covered with cloud
  #---------------------------------------------------------------

  # ######################
  # # for debugging:
  # end.date <- Sys.Date()
  # savedir <- paste( "/alphadata01/bstocker/data/modis_", varnam, "_fluxnet_cutouts/AT-Neu/data_AT-Neu_2000-02-18/", sep="" )
  # overwrite <- FALSE
  # sitename <- "AT-Neu"
  # lon <- 11.3175
  # lat <- 47.1167
  # expand_x <- 0
  # expand_y <- 0 
  # filn <- list.files( path=savedir, pattern="*asc" )
  # ######################

  require( MODISTools )
  require( dplyr )

  ScaleFactor <- 1e-4  # applied to output variable, in ascii file EVI value is multiplied by 1e4
  ndayyear    <- 365

  ## Read dowloaded ASCII file
  print( paste( "reading file ", paste( savedir, filn, sep="" ) ) )
  crude   <- read.csv( paste( savedir, filn, sep="" ), header = FALSE, as.is = TRUE )
  # crude   <- rename( crude, nrows=V1, ncols=V2, modislon_ll=V3, modislat_ll=V4, dxy_m=V5, id=V6, MODISprod=V7, yeardoy=V8, coord=V9, VMODISprocessdatetime=V10 )

  ## this is just read to get length of time series and dates
  tseries    <- MODISTimeSeries( savedir, Band = band_var )
  ntsteps    <- dim(tseries[[1]])[1]
  tmp        <- rownames( tseries[[1]] )
  time       <- data.frame( yr=as.numeric( substr( tmp, start=2, stop=5 )), doy=as.numeric( substr( tmp, start=6, stop=8 )) )
  time$dates <- as.POSIXct( as.Date( paste( as.character(time$yr), "-01-01", sep="" ) ) + time$doy - 1 )
  time$yr_dec<- time$yr + ( time$doy - 1 ) / ndayyear

  ## this may arise for some odd reason
  if (ntsteps!=1) {
    crude   <- crude[1,]
    ntsteps <-1
    time    <- time[1,]
  }

  ## get number of products for which data is in ascii file (not used)
  nprod <- dim(crude)[1] / ntsteps

  if ((dim(crude)[1]/nprod)!=ntsteps) { print("problem") }

  ## re-arrange data
  if ( dim(crude)[2]==11 && expand_x!=0 && expand_y!=0 ){

    # ---------------------------------------------------------------
    ## more data required than downloaded (additional pixels surrounding centre)
    #---------------------------------------------------------------
    print( "Only centre pixel downloaded, re-download data with surrounding pixels" )
    filn <- try( download_subset_modis( lon, lat, time$dates, time$dates, savedir, overwrite=TRUE ) )

    ## re-read crude date
    crude   <- try( read.csv( paste( savedir, filn, sep="" ), header = FALSE, as.is = TRUE ) )
    # crude   <- rename( crude, nrows=V1, ncols=V2, modislon_ll=V3, modislat_ll=V4, dxy_m=V5, id=V6, MODISprod=V7, yeardoy=V8, coord=V9, MODISprocessdatetime=V10 )

    ## this is just read to get length of time series and dates
    tseries    <- MODISTimeSeries( savedir, Band = band_var )
    ntsteps    <- dim(tseries[[1]])[1]
    tmp        <- rownames( tseries[[1]] )
    time       <- data.frame( yr=as.numeric( substr( tmp, start=2, stop=5 )), doy=as.numeric( substr( tmp, start=6, stop=8 )) )
    time$dates <- as.POSIXct( as.Date( paste( as.character(time$yr), "-01-01", sep="" ) ) + time$doy - 1 )
    time$yr_dec<- time$yr + ( time$doy - 1 ) / ndayyear

    ## get number of products for which data is in ascii file (not used)
    nprod <- dim(crude)[1] / ntsteps

    nice_all      <- ExtractTile( Data = crude[1:ntsteps,11:dim(crude)[2]] * ScaleFactor, Rows = c(crude$V1[1],expand_y), Cols = c(crude$V2[1],expand_x), Grid = TRUE )
    nice_qual_flg <- ExtractTile( Data = crude[(ntsteps+1):(2*ntsteps),11:dim(crude)[2]], Rows = c(crude$V1[1],expand_y), Cols = c(crude$V2[1],expand_x), Grid = TRUE )

  } else if ( dim(crude)[2]==11 && expand_x==0 && expand_y==0 ){

    # ---------------------------------------------------------------
    ## only one pixel downloaded
    # ---------------------------------------------------------------
    nice_all      <- as.matrix( crude$V11[1:ntsteps], dim(1,1,ntsteps) ) * ScaleFactor     ## EVI data
    nice_qual_flg <- as.matrix( crude$V11[(ntsteps+1):(2*ntsteps)], dim(1,1,ntsteps) )     ## pixel reliability data

  } else if ( dim(crude)[2]>11 ){
    
    # ---------------------------------------------------------------
    ## multiple pixels downloaded
    # ---------------------------------------------------------------
    # nice <- ExtractTile( Data = tseries, Rows = c(crude$V1,expand_y), Cols = c(crude$V2,expand_x), Grid = TRUE )    ## > is not working: applying ExtractTile to return of MODISTimeSeries
    nice_all      <- ExtractTile( Data = crude[1:ntsteps,11:dim(crude)[2]] * ScaleFactor, Rows = c(crude$V1[1],expand_y), Cols = c(crude$V2[1],expand_x), Grid = TRUE )
    nice_qual_flg <- ExtractTile( Data = crude[(ntsteps+1):(2*ntsteps),11:dim(crude)[2]], Rows = c(crude$V1[1],expand_y), Cols = c(crude$V2[1],expand_x), Grid = TRUE )

  } else {

    print( "Not sufficient data downloaded. Adjust expand_x and expand_y.")
    print( paste("dimension length of crude ", dim(crude) ) )

  }

  # ## Clean data for centre pixel: in case quality flag is not '0', use mean of all 8 surrounding pixels
  # if ( expand_x==1 && expand_y==1 ){
  #   # nice_centre <- nice_all[2,2,]
  #   # nice_centre[ which( nice_qual_flg[2,2,]!=0 ) ] <- apply( nice_all[,,which( nice_qual_flg[2,2,]!=0 )], c(3), FUN=mean )
  #   # for (idx in 1:ntsteps){
  #   #   if (nice_qual_flg[2,2,idx]!=0){
  #   #     nice_centre[idx] <- mean( nice_all[,,idx], na.rm=TRUE )
  #   #   }
  #   # }
  #   modis <- list( nice_all=nice_all, nice_centre=nice_centre, nice_qual_flg=nice_qual_flg, time=time )
  # } else {
  #   modis <- list( nice_all=nice_all, nice_qual_flg=nice_qual_flg, time=time )
  # }
  
  if (ntsteps==1 && length(dim(nice_all))==3 ) {
    nice_all <- nice_all[,,1]
    nice_qual_flg <- nice_qual_flg[,,1]
    time <- time[1,]
  }

  modis <- list( nice_all=nice_all, nice_qual_flg=nice_qual_flg, time=time )

  return( modis )

}

read_crude_modis <- function( varnam, sitename, lon, lat, band_var, band_qc, prod, expand_x, expand_y, overwrite, overwrite_dates=FALSE, ignore_missing=FALSE ){
  ##--------------------------------------
  ## Returns data frame containing EVI 
  ## (and year, moy, doy) for all available
  ## months. Interpolated to mid-months
  ## from original 16-daily data.
  ##--------------------------------------

  # ######################
  # ## for debugging:
  # varnam="evi"
  # overwrite=overwrite_raw
  # overwrite_dates=FALSE
  # # # overwrite <- FALSE
  # # # sitename <- "FR-Pue"
  # # # lon <- 3.5958
  # # # lat <- 43.7414
  # # varnam <- "evi"
  # # sitename <-  "AR-SLu"
  # # lon <- -66.4598
  # # lat <- -33.4648
  # # # sitename <- "AT-NEU"
  # # # lon <- 11.3175
  # # # lat <- 47.1167
  # # # sitename <- "AR-Vir"
  # # # lon <- -56.1886
  # # # lat <- -28.2395
  # # expand_x <- 1
  # # expand_y <- 1 
  # # overwrite <- FALSE
  # ######################

  require( MODISTools )
  require( dplyr )
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )
  source( paste( myhome, "sofun/getin/init_daily_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/init_monthly_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/monthly2daily.R", sep="" ) )

  dirnam_dates_csv <- paste( myhome, "data/modis_", varnam, "_fluxnet_cutouts/", sitename,"/", sep="" )
  system( paste( "mkdir -p ", dirnam_dates_csv, sep="" ) )  
  filnam_dates_csv <- paste( dirnam_dates_csv, "dates_MOD13Q1_", sitename, ".csv", sep="" )

  if ( !file.exists(filnam_dates_csv) || overwrite_dates ){
    ##--------------------------------------
    ## Get dates for which data is available
    ##--------------------------------------
    print( paste("lon", lon, "lat", lat))
    dates <- data.frame( dates=GetDates( Product = "MOD13Q1", Lat = lat, Long = lon ) )
    print("done")
    dates$yr  <- as.numeric( substr( dates$dates, start=2, stop=5 ))
    dates$doy <- as.numeric( substr( dates$dates, start=6, stop=8 ))
    dates$start <- as.POSIXct( as.Date( paste( as.character(dates$yr), "-01-01", sep="" ) ) + dates$doy - 1 )
    dates$end   <- as.POSIXct( as.Date( paste( as.character(dates$yr), "-01-01", sep="" ) ) + dates$doy + 10 )

    ## add absolute day (since 1. Jan. 2000)
    dates$ndayyear <- rep( 0, dim(dates)[1] )
    for (idx in 2:dim(dates)[1]){
      if (dates$yr[idx]>dates$yr[idx-1]){
        if ((dates$yr[idx-1]-2000)%%4==0){
          dates$ndayyear[idx] <- 366
        } else {
          dates$ndayyear[idx] <- 365
        }
      }
    }
    dates$absday <- cumsum(dates$ndayyear) + dates$doy

    ##--------------------------------------------------------------------
    ## Save dates as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing dates data frame into CSV file ", filnam_dates_csv, "...") )
    write.csv( dates, file=filnam_dates_csv, row.names=FALSE )

  } else {
    ##--------------------------------------------------------------------
    ## Read dates from CSV file
    ##--------------------------------------------------------------------
    print( "dates information already available in CSV file ...")
    dates <- read.csv( filnam_dates_csv, as.is=TRUE )

  }

  ##--------------------------------------
  ## Collect data for all available dates
  ##--------------------------------------
  modis <- subset( dates, select=c( yr, doy, start, end, absday ) )
  # print(dim(modis))
  modis$centre         <- rep( NA, dim(modis)[1] )
  modis$centre_meansurr<- rep( NA, dim(modis)[1] )
  modis$centre_qc      <- rep( NA, dim(modis)[1] )
  modis$yr_read        <- rep( NA, dim(modis)[1] )
  modis$doy_read       <- rep( NA, dim(modis)[1] )
  modis$date_read      <- rep( NA, dim(modis)[1] )
  modis$yr_dec_read    <- rep( NA, dim(modis)[1] )

  ## Loop over all dates and get data 
  for (idx in 1:dim(modis)[1]){
    
    savedir <- paste( myhome, "data/modis_", varnam, "_fluxnet_cutouts/", sitename,"/data_", sitename, "_", as.Date(modis$start[idx]), "/", sep="" )
    print(paste("Raw ascii data directory:", savedir))
    system( paste( "mkdir -p ", savedir, sep="" ) )  

    ##--------------------------------------------------------------------
    ## Download file with crude data if it's not there yet
    ##--------------------------------------------------------------------
    filn <- try( download_subset_modis( lon, lat, c( band_var, band_qc ), prod, modis$start[idx], modis$end[idx], savedir, overwrite=overwrite, ignore_missing=ignore_missing ) )
    ##--------------------------------------------------------------------

    if ( (ignore_missing && is.na(filn)) ){

      print( paste( "missing file, but ignoring it for site", sitename ) )

    } else {
      ##--------------------------------------------------------------------
      ## Read crude data file
      ##--------------------------------------------------------------------
      out <- read_crude_modis_tstep( filn, savedir, band_var, band_qc, expand_x=expand_x, expand_y=expand_y )
      ##--------------------------------------------------------------------

      if ( is.null( dim( out$nice_all ) ) && expand_x==0 && expand_y==0 ){

        modis$centre_qc[idx] <- out$nice_qual_flg
        modis$centre[idx]    <- out$nice_all

      } else if ( dim(out$nice_all)==c(3,3) ){

        ## get mean of surrounding pixels, first drop all data that is not quality flag 0 for getting mean across surrounding pixels
        surr <- out$nice_all
        # surr[ which( out$nice_qual_flg!=0 ) ] <- NA
        surr[ which( surr<0 ) ] <- NA
        modis$centre_meansurr[idx] <- mean( surr, na.rm=TRUE )
        if (is.nan(modis$centre_meansurr[idx])) { modis$centre_meansurr[idx] <- NA }

        ## save centre pixel
        modis$centre[idx] <- out$nice_all[2,2]

        ## save quality flag
        modis$centre_qc[idx] <- out$nice_qual_flg[2,2]
                
      } else if ( dim(out$nice_all)==c(1,1) && expand_x==0 && expand_y==0 ) { 

        modis$centre_qc[idx]    <- out$nice_qual_flg[1,1]
        modis$centre[idx]  <- out$nice_all[1,1]
     
      } else {

        print("weird...")
      
      }

      modis$yr_read[idx]      <- out$time$yr[1]
      modis$doy_read[idx]     <- out$time$doy[1]
      modis$date_read[idx]    <- out$time$dates[1]
      modis$yr_dec_read[idx]  <- out$time$yr_dec[1]
  
    }

  }

  ## Make nicer
  modis <- modis %>% 

    ## rename
    dplyr::rename( year_dec=yr_dec_read, year=yr, date=start ) %>% 

    ## select only relevant columns
    dplyr::select( year, doy, date, year_dec, absday, centre, centre_meansurr, centre_qc ) %>% 

    ## remove rows where dates were not read
    dplyr::filter( !is.na(year_dec) )


  return( modis )
    
}

interpolate_modis <- function( modis, sitename, lon, lat, prod, do_interpolate=TRUE ){
  ##--------------------------------------
  ## Returns data frame containing data 
  ## (and year, moy, doy) for all available
  ## months. Interpolated to mid-months
  ## from original 16-daily data.
  ##--------------------------------------

  # ######################
  # ## for debugging:
  # # overwrite <- FALSE
  # sitename        = "AR-SLu"
  # lon             = -66.4598
  # lat             = -33.4648
  # do_interpolate  = TRUE 
  #                # sitename = "AR-Vir"
  #                # lon = -56.1886
  #                # lat = -28.2395
  # ######################

  require( dplyr )
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )
  source( paste( myhome, "sofun/getin/init_daily_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/init_monthly_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/monthly2daily.R", sep="" ) )
  source( paste( myhome, "sofun/getin/get_consecutive.R", sep="" ) )
  source( paste( myhome, "sofun/getin/remove_outliers.R", sep="" ) )

  ##--------------------------------------
  ## CLEAN AND GAP-FILL
  ##--------------------------------------
  modis_gapfilled <- modis 
  error <- FALSE

  if (prod=="MOD13Q1"){
    ## MOD13Q1 contains EVI

    ##--------------------------------------
    ## data cleaning
    ##--------------------------------------
    ## Replace data points with quality flag = 2 (snow covered) by 0
    # modis_gapfilled$centre[ which(modis_gapfilled$centre_qc==2) ] <- max( min( modis_gapfilled$centre ), 0.0 )
    modis_gapfilled$raw <- modis_gapfilled$centre
    modis_gapfilled$centre[ which(modis_gapfilled$centre<0) ] <- NA

    if (!is.null(modis_gapfilled$centre_qc)){
      ## Drop all data with quality flag 3, 1 or -1
      modis_gapfilled$centre[ which(modis_gapfilled$centre_qc==3) ]  <- NA  # Target not visible, covered with cloud
      modis_gapfilled$centre[ which(modis_gapfilled$centre_qc==1) ]  <- NA  # Useful, but look at other QA information
      modis_gapfilled$centre[ which(modis_gapfilled$centre_qc==-1) ] <- NA  # Not Processed
    }

    ## open plot for illustrating gap-filling
    pdf( paste("fig/evi_fill_", sitename, ".pdf", sep="" ), width=10, height=6 )
    plot( modis_gapfilled$year_dec, modis_gapfilled$raw, pch=16, col='black', main=sitename, ylim=c(0,1), xlab="year", ylab="MODIS EVI 250 m", las=1 )
    left <- seq(2000, 2016, 2)
    right <- seq(2001, 2017, 2)
    rect( left, -99, right, 99, border=NA, col=rgb(0,0,0,0.2) )
    points( modis_gapfilled$year_dec, modis_gapfilled$centre, pch=16, col='red' )

    # ## Drop all data identified as outliers = lie outside 5*IQR
    # modis_gapfilled$centre <- remove_outliers( modis_gapfilled$centre, coef=5 ) ## maybe too dangerous - removes peaks

    ## add points to plot opened before
    points( modis_gapfilled$year_dec, modis_gapfilled$centre, pch=16, col='blue' )

    ##--------------------------------------
    ## get LOESS spline model for predicting daily values (below)
    ##--------------------------------------
    idxs    <- which(!is.na(modis_gapfilled$centre))
    myloess <- try( with( modis_gapfilled, loess( centre[idxs] ~ year_dec[idxs], span=0.01 ) ))
    i <- 0
    while (class(myloess)=="try-error" && i<50){
      i <- i + 1
      print(paste("i=",i))
      myloess <- try( with( modis_gapfilled, loess( centre[idxs] ~ year_dec[idxs], span=(0.01+0.002*(i-1)) ) ))
    }
    print("ok now...")

    ##--------------------------------------
    ## get spline model for predicting daily values (below)
    ##--------------------------------------
    spline <- try( with( modis_gapfilled, smooth.spline( year_dec[idxs], centre[idxs], spar=0.001 ) ) )

    ## aggregate by DOY
    agg <- aggregate( centre ~ doy, data=modis_gapfilled, FUN=mean, na.rm=TRUE )
    if (is.element("centre_meansurr", names(modis_gapfilled))){
      agg_meansurr <- aggregate( centre_meansurr ~ doy, data=modis_gapfilled, FUN=mean, na.rm=TRUE )
      agg <- agg %>% left_join( agg_meansurr ) %>% dplyr::rename( centre_meandoy=centre, centre_meansurr_meandoy=centre_meansurr )
    } else {
      agg <- agg %>% dplyr::rename( centre_meandoy=centre )      
    }
    modis_gapfilled <- modis_gapfilled %>% left_join( agg )

    ##--------------------------------------
    ## gap-fill with information from surrounding pixels 
    ##--------------------------------------
    idxs <- which( is.na(modis_gapfilled$centre) )
    if (is.element("centre_meansurr", names(modis_gapfilled))){
      ## get current anomaly of mean across surrounding pixels w.r.t. its mean annual cycle
      modis_gapfilled$anom_surr    <- modis_gapfilled$centre_meansurr / modis_gapfilled$centre_meansurr_meandoy
      modis_gapfilled$centre[idxs] <- modis_gapfilled$centre_meandoy[idxs] * modis_gapfilled$anom_surr[idxs]
    } else {
      modis_gapfilled$centre[idxs] <- modis_gapfilled$centre_meandoy[idxs]
    }
    with( modis_gapfilled, points( year_dec[idxs], centre[idxs], pch=16, col='green' ) )
    legend("topright", c("modis", "outliers", "after bad values dropped and outliers removed", "added from mean of surrounding" ), col=c("black", "red", "blue", "green" ), pch=16, bty="n" )
    # legend("topleft", c("R LOESS smoothing with span=0.01", "R smooth.spline"), col=c("red", "dodgerblue"), lty=1, bty="n" )


    # points( modis_gapfilled$yr_dec_read[idxs],  modis_gapfilled$centre[idxs],  pch=16 )
    # points( modis_gapfilled$yr_dec_read[-idxs], modis_gapfilled$centre[-idxs], pch=16, col='blue' )

    # ## Gap-fill remaining again by mean-by-DOY
    # idxs <- which( is.na(modis_gapfilled$centre) )
    # modis_gapfilled$centre[idxs] <- modis_gapfilled$centre_meandoy[idxs]
    # # points( modis_gapfilled$yr_dec_read[idxs], modis_gapfilled$centre[idxs], pch=16, col='red' )

    # ## Gap-fill still remaining by linear approximation
    # idxs <- which( is.na(modis_gapfilled$centre) )
    # if (length(idxs)>1){
    #   modis_gapfilled$centre <- approx( modis_gapfilled$year_dec[-idxs], modis_gapfilled$centre[-idxs], xout=modis_gapfilled$year_dec )$y
    # }

    # points( modis_gapfilled$yr_dec_read[idxs], modis_gapfilled$centre[idxs], pch=16, col='green' )
    # lines( modis_gapfilled$yr_dec_read, modis_gapfilled$centre )
    # dev.off()

  } else if (prod=="MOD15A2"){
    ## MOD15A2 contains fpar

    ##--------------------------------------
    ## data cleaning
    ##--------------------------------------
    # ## Drop all data with quality flag != 0
    # if (!is.null(modis_gapfilled$centre_qc)){
    #   modis_gapfilled$centre[ which(modis_gapfilled$centre_qc!=0) ] <- NA
    # }

    ## no quality info available for fpar from trevor

    ## remove values that are above 1
    modis_gapfilled$centre <- replace( modis_gapfilled$centre,  modis_gapfilled$centre>1.0, NA )

    ## open plot for illustrating gap-filling
    pdf( paste("fig/fpar_fill_", sitename, ".pdf", sep="" ), width=10, height=6 )
    plot( modis$year_dec, modis$centre, pch=16, col='black', main=sitename, ylim=c(0,1), xlab="year", ylab="MODIS FPAR 1 km", las=1 )
    left <- seq(2000, 2016, 2)
    right <- seq(2001, 2017, 2)
    rect( left, -99, right, 99, border=NA, col=rgb(0,0,0,0.2) )
    points( modis_gapfilled$year_dec, modis_gapfilled$centre, pch=16, col='red' )

    # ## Drop all data identified as outliers = lie outside 5*IQR
    # modis_gapfilled$centre <- remove_outliers( modis_gapfilled$centre, coef=5 ) ## maybe too dangerous - removes peaks

    ## add points to plot opened before
    points( modis_gapfilled$year_dec, modis_gapfilled$centre, pch=16, col='springgreen3' )

    ##--------------------------------------
    ## get LOESS spline model for predicting daily values (below)
    ##--------------------------------------
    idxs <- which(!is.na(modis_gapfilled$centre))
    if (length(idxs)>0){
      myloess <- try( with( modis_gapfilled, loess( centre[idxs] ~ year_dec[idxs], span=0.01 ) ))
      i <- 0
      while (class(myloess)=="try-error" && i<50){
        i <- i + 1
        print(paste("i=",i))
        myloess <- try( with( modis_gapfilled, loess( centre[idxs] ~ year_dec[idxs], span=(0.01+0.002*(i-1)) ) ))
      }
      print("ok now...")
    } else {
      missing <- TRUE
    }

    ##--------------------------------------
    ## get spline model for predicting daily values (below)
    ##--------------------------------------
    if (!missing){
      spline <- try( with( modis_gapfilled, smooth.spline( year_dec[idxs], centre[idxs], spar=0.001 ) ))
    }

    if (!missing){
      ## aggregate by DOY
      agg <- aggregate( centre ~ doy, data=modis_gapfilled, FUN=mean, na.rm=TRUE )
      if (is.element("centre_meansurr", names(modis_gapfilled))){
        agg_meansurr <- aggregate( centre_meansurr ~ doy, data=modis_gapfilled, FUN=mean, na.rm=TRUE )
        agg <- agg %>% left_join( agg_meansurr ) %>% dplyr::rename( centre_meandoy=centre, centre_meansurr_meandoy=centre_meansurr )
      } else {
        agg <- agg %>% dplyr::rename( centre_meandoy=centre )
      }
      modis_gapfilled <- modis_gapfilled %>% left_join( agg )

      ## get consecutive data gaps and fill only by mean seasonality if more than 3 consecutive dates are missing
      na_instances <- get_consecutive( is.na(modis_gapfilled$centre), leng_threshold=4, do_merge=FALSE )
      if (nrow(na_instances)>0){
        for (iinst in 1:nrow(na_instances)){
          idxs <- na_instances$idx_start[iinst]:(na_instances$idx_start[iinst]+na_instances$len[iinst]-1)
          if (is.element("centre_meansurr", names(modis_gapfilled))){
            ## get current anomaly of mean across surrounding pixels w.r.t. its mean annual cycle
            modis_gapfilled$anom_surr    <- modis_gapfilled$centre_meansurr / modis_gapfilled$centre_meansurr_meandoy
            modis_gapfilled$centre[idxs] <- modis_gapfilled$centre_meandoy[idxs] * modis_gapfilled$anom_surr[idxs]
          } else {
            modis_gapfilled$centre[idxs] <- modis_gapfilled$centre_meandoy[idxs]
          }
        }
      }

    }

    # ## Gap-fill still remaining by linear approximation
    # idxs <- which( is.na(modis_gapfilled$centre) )
    # if (length(idxs)>1){
    #   modis_gapfilled$centre <- approx( modis_gapfilled$year_dec[-idxs], modis_gapfilled$centre[-idxs], xout=modis_gapfilled$year_dec )$y
    # }

    # points( modis_gapfilled$year_dec[idxs], modis_gapfilled$centre[idxs], pch=16 )
    # points( modis_gapfilled$year_dec[-idxs], modis_gapfilled$centre[-idxs], pch=16, col='blue' )
    # lines(  modis_gapfilled$year_dec, modis_gapfilled$centre )

    # ## for pixels with low quality information, use mean of surroundings
    # if (is.element("centre_meansurr", names(modis_gapfilled))){
    #   for (idx in seq(nrow(modis_gapfilled))){
    #     if (modis_gapfilled_qc[idx,usecol]!=0) {
    #       modis_gapfilled$centre[idx] <- unname( apply( modis_gapfilled[idx,1:npixels], 1, FUN=mean, na.rm=TRUE ))
    #     }
    #   }
    # }

    # dev.off()

  }

  ##--------------------------------------
  ## MONTHLY DATAFRAME, Interpolate data to mid-months
  ##--------------------------------------
  yearstart  <- min( modis_gapfilled$year )
  yearend    <- max( modis_gapfilled$year )

  modis_monthly <- init_monthly_dataframe( yearstart, yearend )

  modis_monthly <- modis_monthly[ which( modis_monthly$year_dec>=modis_gapfilled$year_dec[1] ), ]
  modis_monthly <- modis_monthly[ which( modis_monthly$year_dec<=modis_gapfilled$year_dec[dim(modis)[1]] ), ]

  if (missing){

    modis_monthly$data <- rep( NA, nrow(modis_monthly) )

  } else {

    # modis_monthly$data <- approx( modis_gapfilled$year_dec, modis_gapfilled$centre, modis_monthly$year_dec )$y
    tmp <- try( with( modis_monthly, predict( myloess, year_dec ) ) )
    if (class(tmp)!="try-error"){
      modis_monthly$data <- tmp
    } else {
      tmp <- try( with( modis_monthly, predict( spline, year_dec ) )$y )
      if (class(tmp)!="try-error"){
        modis_monthly$data <- tmp
      } else {
        modis_monthly$data <- rep( NA, nrow(modis_monthly) )
      }
    }
    tmp <- try( with( modis_monthly, predict( spline, year_dec ) )$y)
    if (class(tmp)!="try-error"){
      modis_monthly$data_spline <- tmp
    } else {
      modis_monthly$data_spline <- rep( NA, nrow(modis_monthly) )
    }

    modis_monthly$data         <- replace( modis_monthly$data, modis_monthly$data<0, 0  )
    modis_monthly$data_spline  <- replace( modis_monthly$data_spline, modis_monthly$data_spline<0, 0  )
    
    modis_monthly$data         <- replace( modis_monthly$data, modis_monthly$data>1, 1  )
    modis_monthly$data_spline  <- replace( modis_monthly$data_spline, modis_monthly$data_spline>1, 1  )
  
  }

  ##--------------------------------------
  ## DAILY DATAFRAME
  ##--------------------------------------
  yearstart  <- min( modis_gapfilled$year )
  yearend    <- max( modis_gapfilled$year )

  modis_daily <- init_daily_dataframe( yearstart, yearend )
  
  modis_daily <- modis_daily[ which( modis_daily$year_dec>=modis_gapfilled$year_dec[1] ), ]
  modis_daily <- modis_daily[ which( modis_daily$year_dec<=modis_gapfilled$year_dec[dim(modis)[1]] ), ]

  if (missing){

    modis_daily$data <- rep( NA, nrow(modis_daily) )

  } else {

    if (do_interpolate){
      # modis_daily$data <- approx( modis_gapfilled$year_dec, modis_gapfilled$centre, modis_daily$year_dec )$y
      tmp <- try( with( modis_daily, predict( myloess, year_dec ) ) )
      if (class(tmp)!="try-error"){
        modis_daily$data <- tmp
      } else {
        tmp <- try( with( modis_daily, predict( spline, year_dec ) )$y )
        if (class(tmp)!="try-error"){
          modis_daily$data <- tmp
        } else {
          modis_daily$data <- rep( NA, nrow(modis_daily) )
        }      
      }
      tmp <- try( with( modis_daily, predict( spline, year_dec ) )$y)
      if (class(tmp)!="try-error"){
        modis_daily$data_spline <- tmp
      } else {
        modis_daily$data_spline <- rep( NA, nrow(modis_daily) )
      }

      ## savitzky golay filter
      modis_daily            <- modis_daily %>% left_join( select( modis_gapfilled, year, doy, centre ) )
      modis_daily$interpl    <- approx( modis_daily$year_dec, modis_daily$centre, xout=modis_daily$year_dec )$y 
      modis_daily$sgfiltered <- rep( NA, nrow(modis_daily) )
      idxs <- which(!is.na(modis_daily$interpl))
      modis_daily$sgfiltered[idxs] <- sgolayfilt( modis_daily$interpl[idxs], p=3, n=31 ) 
      
      # tmodis_daily <- sgolayfilt( p=3, n=5, modis$centre ) 

    } else {
      modis_daily$data <- rep( NA, nrow(modis_daily) )
      for (idx in seq(nrow(modis))){
        putidx <- which( modis_daily$year_dec==modis_gapfilled$year_dec[idx] )
        modis_daily$data[putidx] <- modis_gapfilled$data[idx]
      }
    }

    ## limit to within 0 and 1 (loess spline sometimes "explodes")
    modis_daily$data         <- replace( modis_daily$data, modis_daily$data<0, 0  )
    modis_daily$data_spline  <- replace( modis_daily$data_spline, modis_daily$data_spline<0, 0  )
    
    modis_daily$data         <- replace( modis_daily$data, modis_daily$data>1, 1  )
    modis_daily$data_spline  <- replace( modis_daily$data_spline, modis_daily$data_spline>1, 1  )

  }

  ## plot daily smoothed line and close plotting device
  if (!missing){
    # with( modis_daily, lines( year_dec, data_spline, col='darkgoldenrod', lwd=1 ) )
    # with( modis_daily, lines( year_dec, data, col='cyan', lwd=1 ) )
    # with( modis_monthly, lines( year_dec, data, col='blue', lwd=1 ) )
    with( modis_daily,   lines( year_dec, sgfiltered, col='red', lwd=1 ) )
  }
  dev.off()

  ##--------------------------------------
  ## select standard filtering method for daily data frames
  ##--------------------------------------
  if (prod=="MOD15A2" && !missing) { 
    modis_daily <- modis_daily %>% select( -data ) %>% dplyr::rename( data=sgfiltered ) 
  }

  print("done with interpolate_modis()")

  return( list( modis_gapfilled=modis_gapfilled, modis_daily=modis_daily, modis_monthly=modis_monthly ) )

}

