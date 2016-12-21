library(plyr)
#library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/mymodistools_tseries.R", sep="" ) )

##--------------------------------------------------------------------
## MANUAL SETTINGS
##--------------------------------------------------------------------
simsuite="fluxnet2015"
expand_x=1
expand_y=1
overwrite=FALSE
ignore_missing=FALSE
##--------------------------------------------------------------------

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

year_start <- 1999
year_end   <- 2016
nyears <- year_end - year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

# do.sites <- seq(nsites)
do.sites <- 3:3

for (idx in do.sites){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_csv <- paste( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/", sitename, "/", sep="" )
  filnam_modis_csv <- paste( dirnam_csv,   "gpp_8d_modissubset_", sitename, ".csv", sep="" )
  filnam_monthly_csv <- paste( dirnam_csv, "mgpp_modissubset_", sitename, ".csv", sep="" )
  filnam_daily_csv   <- paste( dirnam_csv, "dgpp_modissubset_", sitename, ".csv", sep="" )

  if (!file.exists(filnam_modis_csv)||overwrite){

    ##--------------------------------------------------------------------
    out <- interpolate_modis( sitename, lon, lat, expand_x=expand_x, expand_y=expand_y, outdir="data/modis_gpp_fluxnet_cutouts_tseries/", overwrite=overwrite, do_interpolate=FALSE )
    ##--------------------------------------------------------------------

    # ##--------------------------------------------------------------------
    # ## add dummy year 1999 with median of each month in all subsequent years
    # ##--------------------------------------------------------------------
    # ## daily data frame
    # dummy_1999 <- init_daily_dataframe( 1999, 1999 )
    # dummy_1999$evi <- rep( NA, dim(dummy_1999)[1] )
    # for (idx in 1:ndayyear){
    #   dummy_1999$evi[idx] <- median( out$modis_daily$evi[ out$modis_daily$doy==idx ], na.rm=TRUE )
    # }
    # out$modis_daily <- rbind( dummy_1999, out$modis_daily )

    # ## monthly data frame
    # dummy_1999 <- init_monthly_dataframe( 1999, 1999 )
    # dummy_1999$evi <- rep( NA, dim(dummy_1999)[1] )
    # for (idx in 1:nmonth){
    #   dummy_1999$evi[idx] <- median( out$modis_monthly$evi[ out$modis_monthly$moy==idx ], na.rm=TRUE )
    # }
    # out$modis_monthly <- rbind( dummy_1999, out$modis_monthly )

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
