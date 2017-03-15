library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/mymodistools.R", sep="" ) )
source( paste( myhome, "sofun/getin/init_monthly_dataframe.R", sep="" ) )
source( paste( myhome, "sofun/getin/init_daily_dataframe.R", sep="" ) )

##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
simsuite        = "fluxnet2015"
expand_x        = 1
expand_y        = 1
overwrite_csv  = TRUE
overwrite       = FALSE
overwrite_dates = FALSE
ignore_missing  = FALSE
bundle          = "fapar"
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

fapar_year_start <- 1999
fapar_year_end   <- 2016
nyears <- fapar_year_end - fapar_year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

## create data frame holding data for all sites
mdf_modis_allsites <- init_monthly_dataframe( fapar_year_start, fapar_year_end )
ddf_modis_allsites <- init_daily_dataframe( fapar_year_start, fapar_year_end )
path_mfapar_allsites_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/mfapar_", varnam, "_modissubset_allsites.csv", sep="")
path_dfapar_allsites_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/dfapar_", varnam, "_modissubset_allsites.csv", sep="")

overwrite_csv <- TRUE

do.sites <- seq(nsites)
do.sites <- 1:1
# do.sites <- 100

for (idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_csv_from       <- paste( myhome, "data/modis_", varnam,"_fluxnet_cutouts_tseries/", sitename, "/", sep="" )
  dirnam_csv_to         <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", sep="" )
  filnam_modisdates_csv <- paste( dirnam_csv, varnam,"_8d_modissubset_", sitename, ".csv", sep="" )
  filnam_monthly_csv    <- paste( dirnam_csv, "m", varnam,"_modissubset_", sitename, ".csv", sep="" )
  filnam_daily_csv      <- paste( dirnam_csv, "d", varnam,"_modissubset_", sitename, ".csv", sep="" )

  filnam_from <- paste( dirnam_csv_from, varnam,"_8d_modissubset_", sitename, ".csv", sep="" )
  filnam_to   <- paste( dirnam_csv_to,   varnam,"_8d_modissubset_", sitename, ".csv", sep="" )

  ## move files
  if ( file.exists( filnam_from ) ){
    system( paste( "mv",  filnam_from, dirnam_csv_to ) )
    system( paste( "ln -svf ", filnam_to, " ", dirnam_csv_from, ".", sep="" ) )
    print( paste( "mv",  filnam_from, dirnam_csv_to ) )
    print( paste( "ln -svf ", filnam_to, " ", dirnam_csv_from, ".", sep="" ) )
  }

  # if ( !file.exists(filnam_modisdates_csv)||overwrite_csv ){

  #   print(paste("file does not exist: ", filnam_modisdates_csv ))

  #   ##--------------------------------------------------------------------
  #   ## Download data from MODIS, read into a nice dataframe and interpolate
  #   ## to daily, monthly, and original time steps.
  #   ##--------------------------------------------------------------------
  #   out <- interpolate_modis(
  #       varnam,
  #       sitename,
  #       lon,
  #       lat,
  #       band_var,
  #       band_qc, 
  #       prod, 
  #       expand_x=expand_x,
  #       expand_y=expand_y,
  #       overwrite=overwrite,
  #       overwrite_dates=overwrite_dates,
  #       ignore_missing=ignore_missing  
  #       )
  #   ##--------------------------------------------------------------------
    
  #   if (!out$nodata){
      
  #     df_monthly <- out$modis_monthly
  #     df_daily   <- out$modis_daily

  #     ##--------------------------------------------------------------------
  #     ## Clean data frame that has MODIS-dates
  #     ##--------------------------------------------------------------------
  #     df_modis   <- out$modis
  #     df_modis$date_read <- unlist(df_modis$date_read)
  #     df_modis   <- as.data.frame( df_modis )
  #     df_modis$yr_read <- NULL
  #     df_modis$doy_read <- NULL
  #     df_modis$date_read <- NULL
  #     df_modis$yr_dec_read <- NULL
  #     df_modis$absday <- NULL

  #     ##--------------------------------------------------------------------
  #     ## add dummy year 1999 with median of each month in all subsequent years
  #     ##--------------------------------------------------------------------
  #     ## daily data frame
  #     dummy_1999 <- init_daily_dataframe( 1999, 1999 )
  #     dummy_1999$evi <- rep( NA, dim(dummy_1999)[1] )
  #     for (idx in 1:ndayyear){
  #       dummy_1999$evi[idx] <- median( df_daily$evi[ df_daily$doy==idx ], na.rm=TRUE )
  #     }
  #     df_daily <- rbind( dummy_1999, out$modis_daily )
      
  #     ## monthly data frame
  #     dummy_1999 <- init_monthly_dataframe( 1999, 1999 )
  #     dummy_1999$evi <- rep( NA, dim(dummy_1999)[1] )
  #     for (idx in 1:nmonth){
  #       dummy_1999$evi[idx] <- median( df_monthly$evi[ df_monthly$moy==idx ], na.rm=TRUE )
  #     }
  #     df_monthly <- rbind( dummy_1999, out$modis_monthly )
      
  #     ##--------------------------------------------------------------------
  #     ## Save data frames as CSV files
  #     ##--------------------------------------------------------------------
  #     print( paste( "writing fapar data frame into CSV file ", filnam_monthly_csv, "...") )
  #     system( paste( "mkdir -p", dirnam_fapar_csv ) )   
  #     # print("WARNING: NOT WRITING CSV FILES!!!")
  #     write.csv( df_monthly, file=filnam_monthly_csv,    row.names=FALSE )
  #     write.csv( df_daily,   file=filnam_daily_csv,      row.names=FALSE )
  #     write.csv( df_modis,   file=filnam_modisdates_csv, row.names=FALSE )
      
  #   }

  # } else {

  #   print(paste("file exists: ", filnam_modisdates_csv ))

  #   ##--------------------------------------------------------------------
  #   ## Read CSV file directly for this site
  #   ##--------------------------------------------------------------------
  #   print( "monthly fapar data already available in CSV file ...")
  #   df_monthly <- read.csv( filnam_monthly_csv )
  #   df_daily   <- read.csv( filnam_daily_csv )

  # }

  # # ##--------------------------------------------------------------------
  # # ## Attach to data frame for all sites
  # # ##--------------------------------------------------------------------
  # # mdf_modis_allsites[[ sitename ]] <- df_monthly$evi

  # ##--------------------------------------------------------------------
  # ## Write to Fortran-formatted output for each variable and year separately
  # ##--------------------------------------------------------------------
  # print( "writing formatted input files ..." )

  # ## in separate formatted file 
  # if (!out$nodata){
  #   for (yr in unique(df_monthly$year)){
      
  #     print( paste("... for year", yr))
  #     dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
  #     system( paste( "mkdir -p", dirnam ) )
      
  #     filnam <- paste( dirnam, "mfapar_", varnam, "_modissubset_", sitename, "_", yr, ".txt", sep="" )
  #     write_sofunformatted( filnam, df_monthly$evi[ which( df_monthly$year==yr ) ] )
      
  #     filnam <- paste( dirnam, "dfapar_", varnam, "_modissubset_", sitename, "_", yr, ".txt", sep="" )
  #     write_sofunformatted( filnam, df_daily$evi[ which( df_daily$year==yr ) ] )
      
  #   }
  # }
  
}

# ##--------------------------------------------------------------------
# ## Write CSV file for data frame for all sites
# ##--------------------------------------------------------------------
# write.csv( mdf_modis_allsites, file=path_mfapar_allsites_csv, row.names=FALSE )


