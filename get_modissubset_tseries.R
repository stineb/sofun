.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/mymodistools_tseries.R", sep="" ) )

##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
simsuite       = "fluxnet2015"
expand_x       = 1
expand_y       = 1
overwrite      = FALSE
ignore_missing = FALSE
do_interpolate = FALSE
bundle         = "fapar"
##--------------------------------------------------------------------


if (bundle=="fapar"){
  ##--------------------------------------------------------------------
  ## FAPAR
  ##--------------------------------------------------------------------
  band_var <- "Fpar_1km"
  band_qc  <- "FparLai_QC"
  prod     <- "MOD15A2"
  varnam   <- "fapar"

} else if (bundle=="lai"){
  ##--------------------------------------------------------------------
  ## LAI
  ##--------------------------------------------------------------------
  band_var <- "Lai_1km"
  band_qc  <- "FparLai_QC"
  prod     <- "MOD15A2"
  varnam   <- "lai"

} else if (bundle=="gpp"){
  ##--------------------------------------------------------------------
  ## GPP
  ##--------------------------------------------------------------------
  band_var <- "Gpp_1km"
  band_qc  <- "Psn_QC_1km"
  prod     <- "MOD17A2_51"
  varnam   <- "gpp"

} else if (bundle=="evi"){
  ##--------------------------------------------------------------------
  ## EVI
  ##--------------------------------------------------------------------
  band_var <- "250m_16_days_EVI"
  band_qc  <- "250m_16_days_pixel_reliability"
  prod   <- "MOD13Q1"
  varnam <- "evi"

}

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

year_start <- 1999
year_end   <- 2016
nyears <- year_end - year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]
do.sites <- seq(nsites)
do.sites <- c(88)

for (idx in do.sites){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_csv         <- paste( myhome, "data/modis_", varnam,"_fluxnet_cutouts_tseries/", sitename, "/", sep="" )
  filnam_modis_csv   <- paste( dirnam_csv, varnam,"_8d_modissubset_", sitename, ".csv", sep="" )
  filnam_monthly_csv <- paste( dirnam_csv, "m", varnam,"_modissubset_", sitename, ".csv", sep="" )
  filnam_daily_csv   <- paste( dirnam_csv, "d", varnam,"_modissubset_", sitename, ".csv", sep="" )

  print( paste( "outputs stored in", dirnam_csv ) )

  if (!file.exists(filnam_modis_csv)||overwrite){

    ##--------------------------------------------------------------------
    ## Download data from MODIS, read into a nice dataframe and interpolate
    ## to daily, monthly, and original time steps.
    ##--------------------------------------------------------------------
    out <- interpolate_modis( 
      sitename, 
      lon, 
      lat, 
      band_var,
      band_qc, 
      prod, 
      expand_x       = expand_x, 
      expand_y       = expand_y, 
      outdir         = dirnam_csv, 
      overwrite      = overwrite, 
      do_interpolate = do_interpolate, 
      ignore_missing = ignore_missing 
      )

    # ## original ASCII data is scaled and given in kgC m-2. Convert to gC m-2. 
    # if ()
    # tseries$data <- tseries$data * ScaleFactor * 1e3

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing data frame into CSV file ", filnam_modis_csv, "...") )
    system( paste( "mkdir -p", dirnam_csv ) )   
    write.csv( out$modis_monthly, file=filnam_monthly_csv, row.names=FALSE )
    write.csv( out$modis_daily,   file=filnam_daily_csv,   row.names=FALSE )
    write.csv( out$modis,         file=filnam_modis_csv,   row.names=FALSE )

  } else {

    ##--------------------------------------------------------------------
    ## Read CSV file directly for this site
    ##--------------------------------------------------------------------
    print( "monthly data already available in CSV file ...")
    out <- list()
    out$modis         <- read.csv( filnam_modis_csv )
    out$modis_daily   <- read.csv( filnam_daily_csv )
    out$modis_monthly <- read.csv( filnam_monthly_csv )

  }

  ##--------------------------------------------------------------------
  ## Check for some manually downloaded time series 
  ##--------------------------------------------------------------------
  # modis_manually <- read.csv( "/alphadata01/bstocker/data/modis_gpp_fluxnet_cutouts_tseries/FR-Pue/fromMODISwebsite/FR-Pue_MODIS.txt" )
  # modis_manually_gpp <- filter( modis_manually, Band=="Gpp_1km" )
  # modis_manually_qc <- filter( modis_manually, Band=="Psn_QC_1km" )

  # tmp <- modis_manually_gpp$Date
  # modis_manually_gpp$year <- as.numeric( substr( tmp, start=2, stop=5 ))
  # modis_manually_gpp$doy  <- as.numeric( substr( tmp, start=6, stop=8 ))
  # modis_manually_gpp$date <- as.POSIXlt( as.Date( paste( as.character(modis_manually_gpp$year), "-01-01", sep="" ) ) + modis_manually_gpp$doy - 1 )
  # modis_manually_gpp$year_dec <- modis_manually_gpp$year + ( modis_manually_gpp$doy - 1 ) / ndayyear

  # plot( out$modis$year_dec, out$modis$data, type='l', main=sitename )
  # lines( modis_manually_gpp$year_dec, modis_manually_gpp$X25*1e-1, col='red' )

  # ##--------------------------------------------------------------------
  # ## Write to Fortran-formatted output for each variable and year separately
  # ##--------------------------------------------------------------------
  # print( "writing formatted input files ..." )

  # ## in separate formatted file 
  # for (year in unique(df_monthly$year)){

  #   print( paste("... for year", year))
  #   dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(year), "/", sep="" )
  #   system( paste( "mkdir -p", dirnam ) )

  #   filnam <- paste( dirnam, "fapar_evi_modissubset_", sitename, "_", year, ".txt", sep="" )
  #   write_sofunformatted( filnam, df_monthly$evi[ which( df_monthly$year==year ) ] )
    
  # }

}

