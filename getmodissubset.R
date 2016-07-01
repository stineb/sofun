getmodissubset_evi <- function( sitename, lon, lat, start.date, end.date, savedir ){

  library( MODISTools )

  if (!file.exists(savedir)){

    print( paste( "==========================="))
    print( paste( "DOWNLOADING MODIS DATA FOR:"))
    print( paste( "lon  :", lon) )
    print( paste( "lat  :", lat) )
    print( paste( "start:", start.date ) )
    print( paste( "end  :", end.date ) )
    print( paste( "---------------------------"))

    system( paste( "mkdir ", savedir, sep="" ) )  

    modis.subset <- data.frame(  
      lat=lat, 
      long=lon, 
      start.date=start.date, 
      end.date=end.date
      )

    MODISSubsets(
      LoadDat = modis.subset, 
      Products = "MOD13Q1",
      Bands = c("250m_16_days_EVI", "250m_16_days_pixel_reliability"),
      Size = c(0,0),
      StartDate = TRUE,
      SaveDir = savedir
      )

    MODISSummaries(
      LoadDat = modis.subset, 
      Dir = savedir,
      Product = "MOD13Q1", 
      Bands = "250m_16_days_EVI",
      ValidRange = c(-2000,10000), 
      NoDataFill = -3000, 
      StartDate = TRUE,
      Yield = TRUE,
      Interpolate = TRUE, 
      ScaleFactor = 0.0001,
      QualityScreen = TRUE, 
      QualityBand = "250m_16_days_pixel_reliability",
      QualityThreshold = 0
      )
    
  } else {

    print( paste( "==========================="))
    print( paste( "FOUND DATA FOR:"))
    print( paste( "lon  :", lon) )
    print( paste( "lat  :", lat) )
    print( paste( "start:", start.date ) )
    print( paste( "end  :", end.date ) )
    print( paste( "---------------------------"))

  }

  out <- read.csv( paste( savedir, list.files( path=savedir, pattern="MODIS_Summary*" ), sep="" ), as.is = TRUE )
  out$mystart <- start.date ## coincides with MODIS dates
  out$myend   <- end.date

  return(out)

}


get_evi_modis_250m <- function( sitename, lon, lat ){
  ##--------------------------------------
  ## Returns data frame containing EVI 
  ## (and year, moy, doy) for all available
  ## months. Interpolated to mid-months
  ## from original 16-daily data.
  ##--------------------------------------

  library( MODISTools )
  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  ##--------------------------------------
  ## Get dates for which data is available
  ##--------------------------------------
  dates <- data.frame( dates=GetDates( Product = "MOD13Q1", Lat = lat, Long = lon ) )
  dates$yr  <- as.numeric( substr( dates$dates, start=2, stop=5 ))
  dates$doy <- as.numeric( substr( dates$dates, start=6, stop=8 ))
  dates$start <- as.POSIXlt( as.Date( paste( as.character(dates$yr), "-01-01", sep="" ) ) + dates$doy - 1 )
  dates$end   <- as.POSIXlt( as.Date( paste( as.character(dates$yr), "-01-01", sep="" ) ) + dates$doy + 15 )

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

  ##--------------------------------------
  ## Collect data for all available dates
  ##--------------------------------------
  modis <- subset( dates, select=c(yr,doy,start,end,absday) )
  modis$evi <- rep( NA, dim(modis)[1] )

  for (idx in 1:dim(modis)[1]){
    
    tmp2 <- try(
                getmodissubset_evi( 
                                    sitename, lon, lat, modis$start[idx], modis$end[idx], 
                                    paste( 
                                            myhome,
                                            "data/modis_fluxnet_cutouts/data_",
                                            sitename, 
                                            "_", 
                                            as.Date(modis$start[idx]), 
                                            "/", 
                                            sep=""
                                            )
                                    )
                )

    if (class(tmp2)=="try-error") {
      modis$evi[idx] <- NA
    } else {
      modis$evi[idx] <- tmp2$mean.band
    } 
  
  }

  ##--------------------------------------
  ## Interpolate data to mid-months
  ##--------------------------------------
  middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350)

  yrstart <- dates$yr[1]
  yrend   <- dates$yr[dim(dates)[1]]

  dfm <- data.frame( 
    yr=rep(yrstart:yrend,each=12) , 
    moy=rep(1:12,length(yrstart:yrend)), 
    doy=rep(middaymonth,length(yrstart:yrend))
    )
  dfm$date <- as.POSIXlt( as.Date( paste( as.character(dfm$yr), "-01-01", sep="" ) ) + dfm$doy - 1 )

  dfm$ndayyear <- rep( 0, dim(dfm)[1] )
  for (idx in 2:dim(dfm)[1]){
    if (dfm$yr[idx]>dfm$yr[idx-1]){
      if ((dfm$yr[idx-1]-2000)%%4==0){
        dfm$ndayyear[idx] <- 366
      } else {
        dfm$ndayyear[idx] <- 365
      }
    }
  }
  dfm$absday <- cumsum(dfm$ndayyear) + dfm$doy

  dfm$evi <- approx( modis$absday, modis$evi, xout=dfm$absday )$y

  return(dfm)

}

