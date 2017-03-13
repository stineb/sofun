download_subset_modis <- function( lon, lat, bands, prod, sitename, TimeSeriesLength, end.date, savedir, overwrite, ignore_missing=FALSE ){

  # ######################
  # end.date <- Sys.Date()
  # savedir <- paste( "/alphadata01/bstocker/data/modis_gpp_fluxnet_cutouts_tseries/AR-SLu/", sep="" )
  # overwrite <- FALSE
  # # TimeSeriesLength <- 1

  # sitename <- "AR-SLu"
  # lon <- -66.4598
  # lat <- -33.4648
  # ######################

  library( MODISTools )

  ## Find file from which (crude) data is read
  filn <- list.files( path=savedir, pattern="*asc" )

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

read_crude_modis <- function( filn, savedir, band_var, band_qc, expand_x, expand_y ){
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

  source( paste( myhome, "sofun/getin/remove_outliers.R", sep="" ) )

  # ######################
  # # for debugging:
  # savedir <- "/alphadata01/bstocker/data/modis_fapar_fluxnet_cutouts_tseries/FR-Pue/raw/"
  # filn    <- "Lat43.74140Lon3.59580Start1991-01-01End2017-02-21___MOD15A2.asc"
  # expand_x <- 1
  # expand_y <- 1
  # band_var <- "Fpar_1km"
  # band_qc  <- "FparLai_QC"
  # ######################

  library( MODISTools )
  library( dplyr )

  ScaleFactor <- 1e-4  # applied to output variable, in ascii file EVI value is multiplied by 1e4
  ndayyear    <- 365

  ## Read dowloaded ASCII file
  path <- paste( savedir, filn, sep="" )
  print( paste( "reading file ", path ) )
 
  if (file.exists(path)){
    # crude   <- read.csv( paste( savedir, filn, sep="" ), header = FALSE, as.is = TRUE )
    # crude   <- rename( crude, nrows=V1, ncols=V2, modislon_ll=V3, modislat_ll=V4, dxy_m=V5, id=V6, MODISprod=V7, yeardoy=V8, coord=V9, VMODISprocessdatetime=V10 )

    ## this is just read to get length of time series and dates
    # tseries    <- MODISTimeSeries( savedir, Band = c("250m_16_days_EVI", "Gpp_1km") )
    tseries    <- as.data.frame( MODISTimeSeries( savedir, Band = band_var, Simplify=TRUE ) )   
    tseries_qc <- as.data.frame( MODISTimeSeries( savedir, Band = band_qc,  Simplify=TRUE ) )

    ## determine centre pixel's column number in 'tseries'
    usecol <- ( ncol(tseries) - 1 ) / 2 + 1 
    npixels <- ncol( tseries )
    tseries$data    <- tseries[,usecol]
    tseries_qc$data <- tseries_qc[,usecol]

    ntsteps          <- dim(tseries)[1]
    tmp              <- rownames( tseries )
    tseries$year     <- as.numeric( substr( tmp, start=2, stop=5 ))
    tseries$doy      <- as.numeric( substr( tmp, start=6, stop=8 ))
    tseries$date     <- as.POSIXct( as.Date( paste( as.character(tseries$year), "-01-01", sep="" ) ) + tseries$doy - 1 )
    tseries$year_dec <- tseries$year + ( tseries$doy - 1 ) / ndayyear

    tseries_qc <- cbind( tseries_qc, select( tseries, year, doy, date, year_dec ) )

    ## Drop all data with quality flag != 0
    tseries$data[ which(tseries_qc$data!=0) ] <- NA

    ## Drop all data identified as outliers = lie outside 3*IQR
    plot( tseries$year_dec, tseries$data*0.1, pch=16, col='red', xlim=c(2005,2008), main=savedir )
    tseries$data <- remove_outliers( tseries$data, coef=3.0 )

    ## Get mean of surrounding pixels and save in separate column
    if (npixels>1){
      tseries$data_meansurr <- unname( apply( tseries[,1:npixels], 1, FUN=mean, na.rm=TRUE ) )
    }

    ## aggregate by DOY
    agg <- aggregate( data ~ doy, data=tseries, FUN=mean, na.rm=TRUE )
    if (npixels>1){
      agg <- agg %>% rename( data_meandoy=data, data_meansurr_meandoy=data_meansurr )
    } else {
      agg <- agg %>% rename( data_meandoy=data )
    }
    tseries <- tseries %>% left_join( agg )

    idxs <- which( !is.na(tseries$data) )
    if (npixels>1){
      ## get current anomaly of mean across surrounding pixels w.r.t. its mean annual cycle
      tseries$anom_surr <- tseries$data_meansurr / tseries$data_meansurr_meandoy
      tseries$data[-idxs] <- tseries$data_meandoy[-idxs] * tseries$anom_surr[-idxs]
    } else {
      tseries$data[-idxs] <- tseries$data_meandoy[-idxs]
    }

    points( tseries$year_dec[idxs], tseries$data[idxs]*0.1, pch=16 )
    points( tseries$year_dec[-idxs], tseries$data[-idxs]*0.1, pch=16, col='blue' )
    lines( tseries$year_dec, tseries$data*0.1 )

    # ## for pixels with low quality information, use mean of surroundings
    # if (npixels>1){
    #   for (idx in seq(nrow(tseries))){
    #     if (tseries_qc[idx,usecol]!=0) {
    #       tseries$data[idx] <- unname( apply( tseries[idx,1:npixels], 1, FUN=mean, na.rm=TRUE ))
    #     }
    #   }
    # }

    tseries    <- select( tseries,    year_dec, year, doy, date, data )
    tseries_qc <- select( tseries_qc, year_dec, year, doy, date, data )


  } else {

    print("ERROR: raw file not found. Was looking for")
    print( path )

  }

  return( tseries )

}


interpolate_modis <- function( sitename, lon, lat, band_var, band_qc, prod, expand_x, expand_y, outdir, overwrite, ignore_missing=FALSE, do_interpolate=TRUE ){
  ##--------------------------------------
  ## Returns data frame containing data 
  ## (and year, moy, doy) for all available
  ## months. Interpolated to mid-months
  ## from original 16-daily data.
  ##--------------------------------------

  # ######################
  # ## for debugging:
  # # overwrite <- FALSE
  # sitename <- "AR-Vir"
  # lon <- -66.4598
  # lat <- -33.4648
  # # sitename <- "AR-SLu"
  # # sitename <- "AR-Vir"
  # # lon <- -56.1886
  # # lat <- -28.2395
  # expand_x <- 1
  # expand_y <- 1 
  # overwrite <- FALSE
  # ######################

  library( MODISTools )
  library( dplyr )
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )
  source( paste( myhome, "sofun/getin/init_daily_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/init_monthly_dataframe.R", sep="" ) )
  source( paste( myhome, "sofun/getin/monthly2daily.R", sep="" ) )

  ##--------------------------------------
  ## Collect full time series
  ##--------------------------------------
    end.date <- Sys.Date()
    # savedir <- paste( myhome, "data/modis_fluxnet_cutouts_tseries/", sitename, "_", end.date, "/", sep="" )
    savedir <- paste( outdir, "raw/", sep="" )
    system( paste( "mkdir -p ", savedir ) )
    print( savedir )

    ##--------------------------------------------------------------------
    # filn  <- download_subset_modis( lon, lat, , bands, prod, TimeSeriesLength=26, end.date=dates$end[dim(dates)[1]], savedir=savedir, overwrite=FALSE )
    filn  <- download_subset_modis( lon, lat, c( band_var, band_qc ), prod, sitename, TimeSeriesLength=26, end.date=Sys.Date(), savedir=savedir, overwrite=overwrite, ignore_missing=ignore_missing )
    ##--------------------------------------------------------------------
    modis <- read_crude_modis( filn, savedir, band_var, band_qc, expand_x=expand_x, expand_y=expand_y )
    ##--------------------------------------------------------------------

  ##--------------------------------------
  ## MONTHLY DATAFRAME, Interpolate data to mid-months
  ##--------------------------------------
    yearstart  <- min( modis$year )
    yearend    <- max( modis$year )

    modis_monthly <- init_monthly_dataframe( yearstart, yearend )

    modis_monthly <- modis_monthly[ which( modis_monthly$year_dec>=modis$year_dec[1] ), ]
    modis_monthly <- modis_monthly[ which( modis_monthly$year_dec<=modis$year_dec[dim(modis)[1]] ), ]

    modis_monthly$data <- approx( modis$year_dec, modis$nice_centre, modis_monthly$year_dec )$y


  ##--------------------------------------
  ## DAILY DATAFRAME
  ##--------------------------------------
    yearstart  <- min( modis$year )
    yearend    <- max( modis$year )

    modis_daily <- init_daily_dataframe( yearstart, yearend )
    
    modis_daily <- modis_daily[ which( modis_daily$year_dec>=modis$year_dec[1] ), ]
    modis_daily <- modis_daily[ which( modis_daily$year_dec<=modis$year_dec[dim(modis)[1]] ), ]

    if (do_interpolate){
      modis_daily$data <- approx( modis$year_dec, modis$nice_centre, modis_daily$year_dec )$y
    } else {
      modis_daily$data <- rep( NA, nrow(modis_daily) )
      for (idx in seq(nrow(modis))){
        putidx <- which( modis_daily$year_dec==modis$year_dec[idx] )
        modis_daily$data[putidx] <- modis$data[idx]
      }
    }

  return( list( modis=modis, modis_daily=modis_daily, modis_monthly=modis_monthly ) )

}

# savedir <- "/alphadata01/bstocker/data/modis_fluxnet_cutouts_tseries/AT-Neu_2016-10-14"
# overwrite <- FALSE
# end.date <- dates$end[dim(dates)[1]]
# TimeSeriesLength <- 1

# sitename <- "AT-Neu"
# lon <- 11.3175
# lat <- 47.1167

# out <- interpolate_modis( sitename, lon, lat, expand_x=1, expand_y=1  )

# plot(  out$modis$year_dec, out$modis$nice_all[2,2,], type="l" )
# lines( out$modis$year_dec, out$modis$nice_centre, col="red")

# lines( out$modis$year_dec, out$modis$nice_all[1,1,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[1,2,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[1,3,], col="grey70" )

# lines( out$modis$year_dec, out$modis$nice_all[2,1,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[2,2,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[2,3,], col="grey70" )

# lines( out$modis$year_dec, out$modis$nice_all[3,1,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[3,2,], col="grey70" )
# lines( out$modis$year_dec, out$modis$nice_all[3,3,], col="grey70" )

# lines( out$modis_daily$year_dec, out$modis_daily$data, col="blue", lwd=3 )
# lines( out$modis_monthly$year_dec, out$modis_monthly$data, col="green", lwd=3 )

