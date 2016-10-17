library(plyr)
#library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/test.R", sep="" ) )

simsuite <- "fluxnet2015"

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

fapar_year_start <- 1999
fapar_year_end   <- 2016
nyears <- fapar_year_end - fapar_year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

## create data frame holding data for all sites
df_modis_allsites <- data.frame( year=rep( seq( fapar_year_start, fapar_year_end ), each=nmonth ), mo=rep(1:nmonth,nyears) )
path_fapar_allsites_csv <- paste( myhome, "sofun/input_", simsuite,"_sofun/sitedata/fapar/fapar_evi_modissubset_allsites.csv", sep="")

overwrite <- FALSE

do.sites <- seq(nsites)
# do.sites <- 3:10

for (idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_fapar_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", sep="" )
  filnam_monthly_fapar_csv <- paste( dirnam_fapar_csv, "mfapar_evi_modissubset_", sitename, ".csv", sep="" )
  filnam_daily_fapar_csv   <- paste( dirnam_fapar_csv, "dfapar_evi_modissubset_", sitename, ".csv", sep="" )

  if (!file.exists(filnam_daily_fapar_csv)||overwrite){

    ##--------------------------------------------------------------------
    out <- interpolate_modis( sitename, lon, lat, expand_x=1, expand_y=1, overwrite=overwrite  )
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

    # ##--------------------------------------------------------------------
    # ## Save data frames as CSV files
    # ##--------------------------------------------------------------------
    # print( paste( "writing fapar data frame into CSV file ", filnam_monthly_fapar_csv, "...") )
    # system( paste( "mkdir -p", dirnam_fapar_csv ) )   
    # write.csv( out$modis_monthly, file=filnam_monthly_fapar_csv, row.names=FALSE )
    # write.csv( out$modis_daily,   file=filnam_daily_fapar_csv,   row.names=FALSE )

  } else {

    # ##--------------------------------------------------------------------
    # ## Read CSV file directly for this site
    # ##--------------------------------------------------------------------
    # print( "monthly fapar data already available in CSV file ...")
    # df_monthly <- read.csv( filnam_fapar_csv )

  }

  # ##--------------------------------------------------------------------
  # ## Attach to data frame for all sites
  # ##--------------------------------------------------------------------
  # df_modis_allsites[[ sitename ]] <- df_monthly$evi

  # ##--------------------------------------------------------------------
  # ## Write to Fortran-formatted output for each variable and year separately
  # ##--------------------------------------------------------------------
  # print( "writing formatted input files ..." )

  # ## in separate formatted file 
  # for (yr in unique(df_monthly$yr)){

  #   print( paste("... for year", yr))
  #   dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
  #   system( paste( "mkdir -p", dirnam ) )

  #   filnam <- paste( dirnam, "fapar_evi_modissubset_", sitename, "_", yr, ".txt", sep="" )
  #   write_sofunformatted( filnam, df_monthly$evi[ which( df_monthly$yr==yr ) ] )
    
  # }
  
}

# ##--------------------------------------------------------------------
# ## Write CSV file for data frame for all sites
# ##--------------------------------------------------------------------
# write.csv( df_modis_allsites, file=path_fapar_allsites_csv, row.names=FALSE )

