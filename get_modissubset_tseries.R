.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/mymodistools_tseries.R", sep="" ) )

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

year_start <- 1999
year_end   <- 2016
nyears <- year_end - year_start + 1

##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
simsuite         = "fluxnet2015"
expand_x         = 1
expand_y         = 1
overwrite_raw    = FALSE
overwrite_rawcsv = TRUE
overwrite_csv    = TRUE
ignore_missing   = TRUE
do_interpolate   = TRUE
bundle           = "evi"
tseries_out      = TRUE
##--------------------------------------------------------------------

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]
# do.sites <- seq(nsites)
do.sites <- 108:nsites

if (bundle=="fapar"){
  ##--------------------------------------------------------------------
  ## FAPAR
  ##--------------------------------------------------------------------
  band_var <- "Fpar_1km"
  band_qc  <- "FparLai_QC"
  prod     <- "MOD15A2"
  varnam   <- "fapar"
  productnam <- "fpar"

} else if (bundle=="lai"){
  ##--------------------------------------------------------------------
  ## LAI
  ##--------------------------------------------------------------------
  band_var <- "Lai_1km"
  band_qc  <- "FparLai_QC"
  prod     <- "MOD15A2"
  varnam   <- "lai"
  productnam <- "lai"

} else if (bundle=="gpp"){
  ##--------------------------------------------------------------------
  ## GPP
  ##--------------------------------------------------------------------
  band_var <- "Gpp_1km"
  band_qc  <- "Psn_QC_1km"
  prod     <- "MOD17A2_51"
  varnam   <- "gpp"
  productnam <- "gpp"

} else if (bundle=="evi"){
  ##--------------------------------------------------------------------
  ## EVI
  ##--------------------------------------------------------------------
  band_var   <- "250m_16_days_EVI"
  band_qc    <- "250m_16_days_pixel_reliability"
  prod       <- "MOD13Q1"
  varnam     <- "fapar"
  productnam <- "evi"

}

for (idx in do.sites){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  ##--------------------------------------------------------------------
  ## Handle file names
  ##--------------------------------------------------------------------
  dirnam_csv_from       <- paste( myhome, "data/modis_", productnam, "_fluxnet_cutouts_tseries/", sitename, "/", sep="" )
  dirnam_csv_from_alt   <- paste( myhome, "data/modis_", productnam, "_fluxnet_cutouts/", sitename, "/", sep="" )
  dirnam_csv_to         <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", sep="" )

  ## This file is the link for _tseries and split scripts and should contain all the raw data, not gap-filled or interpolated!
  ## this file is written by 'interpolate_modis()'
  filnam_modis_raw_csv <- paste( dirnam_csv_from, varnam, "_", productnam, "_8d_RAW_modissubset_", sitename, ".csv", sep="" )

  ## These files are gap-filled and interpolated
  filnam_modis_csv   <- paste( dirnam_csv_to, varnam, "_", productnam, "_8d_modissubset_", sitename, ".csv", sep="" )
  filnam_monthly_csv <- paste( dirnam_csv_to, "m", varnam, "_", productnam, "_modissubset_", sitename, ".csv", sep="" )
  filnam_daily_csv   <- paste( dirnam_csv_to, "d", varnam, "_", productnam, "_modissubset_", sitename, ".csv", sep="" )

  print( paste( "outputs stored in", dirnam_csv_to ) )

  # ## Move some files
  # filnam_modis_csv_from   <- paste( dirnam_csv_from, varnam, "_", "fapar", "_8d_modissubset_", sitename, ".csv", sep="" )
  # filnam_monthly_csv_from <- paste( dirnam_csv_from, "m", varnam, "_", "fapar", "_modissubset_", sitename, ".csv", sep="" )
  # filnam_daily_csv_from   <- paste( dirnam_csv_from, "d", varnam, "_", "fapar", "_modissubset_", sitename, ".csv", sep="" )
  # if (file.exists(filnam_modis_csv_from)){
  #   system( paste( "mv", filnam_modis_csv_from, dirnam_csv_to ) )
  # }
  # if (file.exists(filnam_monthly_csv_from)){
  #   system( paste( "mv", filnam_monthly_csv_from, dirnam_csv_to ) )
  # }
  # if (file.exists(filnam_daily_csv_from)){
  #   system( paste( "mv", filnam_daily_csv_from, dirnam_csv_to ) )
  # }

  # ## if RAW-data CSV file is not in "data/modis_", varnam, "_", productnam, "_fluxnet_cutouts_tseries/" look for it in "data/modis_", varnam, "_", productnam, "_fluxnet_cutouts/"
  # if (!file.exists(filnam_modis_raw_csv)){
  #   filnam_modis_raw_csv <- paste( dirnam_csv_from_alt, varnam, "_", productnam, "_8d_RAW_modissubset_", sitename, ".csv", sep="" )
  # }

  if (!file.exists(filnam_modis_raw_csv)||overwrite_rawcsv){
    ##--------------------------------------------------------------------
    ## Check if full time series ascii files are available
    ##--------------------------------------------------------------------
    ## Find file from which (crude) data is read
    savedir <- paste( dirnam_csv_from, "raw/", sep="" )
    system( paste( "mkdir -p ", savedir ) )
    filn <- list.files( path=savedir, pattern="*asc" )

    ## If no raw text files are available, look in single time step directory
    if ( length(filn)==0 ){
      savedir <- paste( dirnam_csv_from_alt, "raw/", sep="" )
      filn <- list.files( path=savedir, pattern="*asc" )
      tseries_out <- FALSE
    }

    if (tseries_out){
      ##--------------------------------------------------------------------
      ## Get file name of raw data text file (and download it if they're not there yet)
      ##--------------------------------------------------------------------
      filn  <- download_subset_modis_tseries( 
                                              lon, 
                                              lat, 
                                              c( band_var, band_qc ), 
                                              prod, 
                                              sitename, 
                                              TimeSeriesLength = 26,  # number of years before end.date
                                              end.date         = Sys.Date(), 
                                              savedir          = savedir, 
                                              overwrite        = overwrite_raw, 
                                              ignore_missing   = ignore_missing 
                                              )

      ##--------------------------------------------------------------------
      ## Read raw data and create a data frame holding the complete time series
      ##--------------------------------------------------------------------
      modis <- read_crude_modis_tseries( 
                                        filn, 
                                        savedir, 
                                        band_var, 
                                        band_qc,
                                        prod,
                                        expand_x = expand_x,
                                        expand_y = expand_y 
                                        )

    } else {

      ##--------------------------------------------------------------------
      ## Get all at time steps individually. In this case 'read_crude_modis' is a wrapper
      ##--------------------------------------------------------------------
      modis <- read_crude_modis( 
                                productnam, 
                                sitename, 
                                lon, 
                                lat, 
                                band_var, 
                                band_qc, 
                                prod, 
                                expand_x        = expand_x,
                                expand_y        = expand_y, 
                                overwrite       = overwrite_raw, 
                                overwrite_dates = FALSE, 
                                ignore_missing  = ignore_missing 
                                )

    }

    ##--------------------------------------------------------------------
    ## Write full time series data to nice CSV file
    ##--------------------------------------------------------------------
    print("writing to raw CSV file...")
    write.csv( modis, file=filnam_modis_raw_csv, row.names=FALSE )


  } else {
    ##--------------------------------------------------------------------
    ## Read full data from nice CSV file
    ##--------------------------------------------------------------------
    print("reading from raw CSV file...")
    modis <- read.csv( filnam_modis_raw_csv )

  }

  if (!file.exists(filnam_modis_csv)||overwrite_csv){
    ##--------------------------------------------------------------------
    ## Clean (gapfill and interpolate) full time series data to 8-days, daily, and monthly
    ##--------------------------------------------------------------------
    print("interpolating...")
    out <- interpolate_modis(
                              modis,
                              sitename, 
                              lon, 
                              lat, 
                              prod,
                              do_interpolate = do_interpolate
                              )

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing data frame into CSV file ", filnam_modis_csv, "...") )
    system( paste( "mkdir -p", dirnam_csv_to ) )   
    write.csv( out$modis_monthly,   file=filnam_monthly_csv, row.names=FALSE )
    write.csv( out$modis_daily,     file=filnam_daily_csv,   row.names=FALSE )
    write.csv( out$modis_gapfilled, file=filnam_modis_csv,   row.names=FALSE )

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

  # #--------------------------------------------------------------------
  # # Check for some manually downloaded time series 
  # #--------------------------------------------------------------------
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

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )

  ## in separate formatted file 
  for (year in unique(out$modis_monthly$year)){

    print( paste("... for year", year))
    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(year), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "mfapar_evi_modissubset_", sitename, "_", year, ".txt", sep="" )
    write_sofunformatted( filnam, out$modis_monthly$data[ which( out$modis_monthly$year==year ) ] )

    filnam <- paste( dirnam, "dfapar_evi_modissubset_", sitename, "_", year, ".txt", sep="" )
    write_sofunformatted( filnam, out$modis_daily$data[ which( out$modis_daily$year==year ) ] )
    
  }

}

