library(plyr)
#library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/getmodissubset.R", sep="" ) )

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

for (idx in seq(nsites)){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]

  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_fapar_csv <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", sep="" )
  filnam_fapar_csv <- paste( dirnam_fapar_csv, "mfapar_evi_modissubset_", sitename, ".csv", sep="" )

  if (!file.exists(filnam_fapar_csv)||overwrite){

    ##--------------------------------------------------------------------
    out_get_evi_modis_250m <- get_evi_modis_250m( sitename, lon, lat )
    ##--------------------------------------------------------------------

    ##--------------------------------------------------------------------
    # get dataframe with data interpolated to monthly values
    ##--------------------------------------------------------------------
    df_monthly   <- out_get_evi_modis_250m$dfm
    df_origdates <- out_get_evi_modis_250m$df_origdates

    ##--------------------------------------------------------------------
    ## gap-fill data with median of corresponding months
    ##--------------------------------------------------------------------
    df_monthly$flag <- rep(0,dim(df_monthly)[1])
    for (idx in 1:dim(df_monthly)[1]){
      if (is.na(df_monthly$evi[idx])){
        df_monthly$evi[idx] <- median( df_monthly$evi[ df_monthly$moy==df_monthly$moy[idx] ], na.rm=TRUE )
        df_monthly$flag[idx] <- 1
      }
    }

    ##--------------------------------------------------------------------
    ## add dummy year 1999 with median of each month in all subsequent years
    ##--------------------------------------------------------------------
    dummy_1999 <- df_monthly[ 1:12, ]
    for (imo in 1:nmonth){
      dummy_1999$evi[imo] <- median( df_monthly$evi[ df_monthly$moy==imo ] )
      dummy_1999$yr[imo] <- 1999 
      dummy_1999$flag[imo] <- 1
    }
    df_monthly <- rbind( dummy_1999, df_monthly )

    ##--------------------------------------------------------------------
    ## Reduce data
    ##--------------------------------------------------------------------
    df_monthly <- subset( df_monthly, select=c( yr, moy, evi ) )

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing fapar data frame into CSV file ", filnam_fapar_csv, "...") )
    system( paste( "mkdir -p", dirnam_fapar_csv ) )   
    write.csv( df_monthly, file=filnam_fapar_csv, row.names=FALSE )

  } else {

    ##--------------------------------------------------------------------
    ## Read CSV file directly for this site
    ##--------------------------------------------------------------------
    print( "monthly fapar data already available in CSV file ...")
    df_monthly <- read.csv( filnam_fapar_csv )

  }

  ##--------------------------------------------------------------------
  ## Attach to data frame for all sites
  ##--------------------------------------------------------------------
  df_modis_allsites[[ sitename ]] <- df_monthly$evi

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )

  ## in separate formatted file 
  for (yr in unique(df_monthly$yr)){

    print( paste("... for year", yr))
    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "fapar_evi_modissubset_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, df_monthly$evi[ which( df_monthly$yr==yr ) ] )
    
  }
  
}

##--------------------------------------------------------------------
## Write CSV file for data frame for all sites
##--------------------------------------------------------------------
write.csv( df_modis_allsites, file=path_fapar_allsites_csv, row.names=FALSE )

