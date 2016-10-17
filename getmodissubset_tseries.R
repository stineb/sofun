download_subset_modis <- function( lon, lat, TimeSeriesLength, end.date, savedir, overwrite ){

  # ######################
  # end.date <- Sys.Date()
  # savedir <- paste( "/alphadata01/bstocker/data/modis_fluxnet_cutouts_tseries/AT-Neu_", end.date, "/", sep="" )
  # overwrite <- FALSE
  # TimeSeriesLength <- 1

  # sitename <- "AT-Neu"
  # lon <- 11.3175
  # lat <- 47.1167
  # ######################

  library( MODISTools )

  ## Find file from which (crude) data is read
  filn <- list.files( path=savedir, pattern="*asc" )

  if (length(filn)==0||overwrite){

    print( paste( "==========================="))
    print( paste( "DOWNLOADING MODIS DATA FOR:"))
    print( paste( "site :", sitename ) )
    print( paste( "lon  :", lon ) )
    print( paste( "lat  :", lat ) )
    print( paste( "---------------------------"))

    # system( paste( "mkdir ", savedir, sep="" ) )  

    bands <- c("250m_16_days_EVI", "250m_16_days_pixel_reliability")
    prod  <- "MOD13Q1"

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
      # Size             = c(1,1),  ## Get subset for all pixels around centre (defined by input lon/lat), stretching 1km in lon and lat (should be 9 in each direction = 81 in total)
      Size             = c(0,0),  ## Get subset for all pixels around centre (defined by input lon/lat), stretching 1km in lon and lat (should be 9 in each direction = 81 in total)
      StartDate        = FALSE,
      SaveDir          = savedir,
      TimeSeriesLength = TimeSeriesLength       ## 26 years beyond first day of the year in 'end.date'. Will trigger reading longet possible automatially 
      # TimeSeriesLength = 1       ## 26 years beyond first day of the year in 'end.date'. Will trigger reading longet possible automatially 
      ) )

  } else {

    print( paste( "==========================="))
    print( paste( "FOUND DATA FOR:" ) )
    print( paste( "file :", filn ) )
    print( paste( "site :", sitename ) )
    print( paste( "lon  :", lon ) )
    print( paste( "lat  :", lat ) )
    print( paste( "---------------------------"))

  }

  return( filn[1] )

}

read_crude_modis <- function( filn, savedir, expand_x, expand_y ){

  # arguments:
  # filn: file name of ASCII file holding MODIS "crude" data
  # savedir: directory, where to look for that file
  # expand_x : number of pixels to the right and left of centre
  # expand_y : number of pixels to the top and bottom of centre


  # ######################
  # end.date <- Sys.Date()
  # savedir <- paste( "/alphadata01/bstocker/data/modis_fluxnet_cutouts_tseries/AT-Neu_", end.date, "/", sep="" )
  # overwrite <- FALSE
  # TimeSeriesLength <- 1
  # sitename <- "AT-Neu"
  # lon <- 11.3175
  # lat <- 47.1167
  # expand_x <- 0
  # expand_y <- 0 
  # filn <- list.files( path=savedir, pattern="*asc" )
  # ######################

  # Rank Key
  # Summary QA
  # Description
  # -1  Fill/No Data  Not Processed
  # 0 Good Data Use with confidence
  # 1 Marginal data Useful, but look at other QA information
  # 2 Snow/Ice  Target covered with snow/ice
  # 3 Cloudy  Target not visible, covered with cloud

  library( MODISTools )

  ScaleFactor <- 0.0001  # applied to output variable
  ndayyear    <- 365

  ## Read dowloaded ASCII file
  print( paste( "reading file ", paste( savedir, filn, sep="" ) ) )
  crude   <- read.csv( paste( savedir, filn, sep="" ), header = FALSE, as.is = TRUE )
  crude   <- rename( crude, c( "V1"="nrows", "V2"="ncols", "V3"="modislon_ll", "V4"="modislat_ll", "V5"="dxy_m", "V6"="id", "V7"="MODISprod", "V8"="yeardoy", "V9"="coord", "V10"="MODISprocessdatetime" ) )

  ## this is just read to get length of time series and dates
  tseries    <- MODISTimeSeries( savedir, Band = "250m_16_days_EVI" )
  ntsteps    <- dim(tseries[[1]])[1]
  tmp        <- rownames( tseries[[1]] )
  time       <- data.frame( yr=as.numeric( substr( tmp, start=2, stop=5 )), doy=as.numeric( substr( tmp, start=6, stop=8 )) )
  time$dates <- as.POSIXlt( as.Date( paste( as.character(time$yr), "-01-01", sep="" ) ) + time$doy - 1 )
  time$yr_dec<- time$yr + ( time$doy - 1 ) / ndayyear

  ## get number of products for which data is in ascii file (not used)
  nprod <- dim(crude)[1] / ntsteps
  if ((dim(crude)[1]/nprod)!=ntsteps) { print("problem") }

  ## re-arrange data
  if ( dim(crude)[2]==11 && expand_x==0 && expand_y==0 ){
    ## only one pixel downloaded
    nice_all      <- as.matrix( crude$V11[1:ntsteps], dim(1,1,ntsteps) ) * ScaleFactor     ## EVI data
    nice_qual_flg <- as.matrix( crude$V11[(ntsteps+1):(2*ntsteps)], dim(1,1,ntsteps) )       ## pixel reliability data

  } else if ( dim(crude)[2]>11 ){
    ## multiple pixels downloaded
    # nice <- ExtractTile( Data = tseries, Rows = c(crude$nrows,expand_y), Cols = c(crude$ncols,expand_x), Grid = TRUE )    ## > is not working: applying ExtractTile to return of MODISTimeSeries
    nice_all      <- ExtractTile( Data = crude[1:ntsteps,11:dim(crude)[2]] * ScaleFactor, Rows = c(crude$nrows[1],expand_y), Cols = c(crude$ncols[1],expand_x), Grid = TRUE )
    nice_qual_flg <- ExtractTile( Data = crude[(ntsteps+1):(2*ntsteps),11:dim(crude)[2]], Rows = c(crude$nrows[1],expand_y), Cols = c(crude$ncols[1],expand_x), Grid = TRUE )

  } else {

    print( "Not sufficient data downloaded. Adjust expand_x and expand_y.")

  }

  ## Clean data for centre pixel: in case quality flag is not '0', use mean of all 8 surrounding pixels
  if ( expand_x==1 && expand_y==1 ){
    nice_centre <- nice_all[2,2,]
    nice_centre[ which( nice_qual_flg[2,2,]!=0 ) ] <- apply( nice_all[,,which( nice_qual_flg[2,2,]!=0 )], c(3), FUN=mean)
    # for (idx in 1:ntsteps){
    #   if (nice_qual_flg[2,2,idx]!=0){
    #     nice_centre[idx] <- mean( nice_all[,,idx], na.rm=TRUE )
    #   }
    # }
  }

  modis <- list( nice_all=nice_all, nice_centre=nice_centre, nice_qual_flg=nice_qual_flg, time=time )
  return( modis )

}


interpolate_modis <- function( sitename, lon, lat, expand_x, expand_y, overwrite ){
  ##--------------------------------------
  ## Returns data frame containing EVI 
  ## (and year, moy, doy) for all available
  ## months. Interpolated to mid-months
  ## from original 16-daily data.
  ##--------------------------------------
  library( MODISTools )
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
    savedir <- paste( myhome, "data/modis_fluxnet_cutouts_tseries/", sitename, "/", sep="" )
    system( paste( "mkdir -p ", savedir ) )
    print( savedir )

    ##--------------------------------------------------------------------
    # filn  <- download_subset_modis( lon, lat, TimeSeriesLength=26, end.date=dates$end[dim(dates)[1]], savedir=savedir, overwrite=FALSE )
    filn  <- download_subset_modis( lon, lat, TimeSeriesLength=26, end.date=Sys.Date(), savedir=savedir, overwrite=overwrite )
    # modis <- read_crude_modis( filn, savedir, expand_x=expand_x, expand_y=expand_y )
    ##--------------------------------------------------------------------

  ##--------------------------------------
  ## MONTHLY DATAFRAME, Interpolate data to mid-months
  ##--------------------------------------
    middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350)

    yrstart  <- min( modis$time$yr )
    yrend    <- max( modis$time$yr )

    modis_monthly <- init_monthly_dataframe( yrstart, yrend )

    modis_monthly <- modis_monthly[ which( modis_monthly$year_dec>=modis$time$yr_dec[1] ), ]
    modis_monthly <- modis_monthly[ which( modis_monthly$year_dec<=modis$time$yr_dec[dim(modis$time)[1]] ), ]

    modis_monthly$evi <- approx( modis$time$yr_dec, modis$nice_centre, modis_monthly$year_dec )$y


  ##--------------------------------------
  ## DAILY DATAFRAME
  ##--------------------------------------
    yrstart  <- min( modis$time$yr )
    yrend    <- max( modis$time$yr )

    modis_daily <- init_daily_dataframe( yrstart, yrend )
    
    modis_daily <- modis_daily[ which( modis_daily$year_dec>=modis$time$yr_dec[1] ), ]
    modis_daily <- modis_daily[ which( modis_daily$year_dec<=modis$time$yr_dec[dim(modis$time)[1]] ), ]

    modis_daily$evi <- approx( modis$time$yr_dec, modis$nice_centre, modis_daily$year_dec )$y


  # return( list( modis=modis, modis_daily=modis_daily, modis_monthly=modis_monthly ) )
  return( 9999 )

}

# savedir <- "/alphadata01/bstocker/data/modis_fluxnet_cutouts_tseries/AT-Neu_2016-10-14"
# overwrite <- FALSE
# end.date <- dates$end[dim(dates)[1]]
# TimeSeriesLength <- 1

# sitename <- "AT-Neu"
# lon <- 11.3175
# lat <- 47.1167

# out <- interpolate_modis( sitename, lon, lat, expand_x=1, expand_y=1  )

# plot(  out$modis$time$yr_dec, out$modis$nice_all[2,2,], type="l" )
# lines( out$modis$time$yr_dec, out$modis$nice_centre, col="red")

# lines( out$modis$time$yr_dec, out$modis$nice_all[1,1,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[1,2,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[1,3,], col="grey70" )

# lines( out$modis$time$yr_dec, out$modis$nice_all[2,1,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[2,2,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[2,3,], col="grey70" )

# lines( out$modis$time$yr_dec, out$modis$nice_all[3,1,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[3,2,], col="grey70" )
# lines( out$modis$time$yr_dec, out$modis$nice_all[3,3,], col="grey70" )

# lines( out$modis_daily$year_dec, out$modis_daily$evi, col="blue", lwd=3 )
# lines( out$modis_monthly$year_dec, out$modis_monthly$evi, col="green", lwd=3 )

