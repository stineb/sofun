library(plyr)
#library(dplyr)

# install.packages("MODISTools")

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )
source( paste( myhome, "sofun/getin/getmodissubset.R", sep="" ) )

simsuite <- "fluxnet"

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
  filnam_fapar_csv <- paste( dirnam_fapar_csv, "fapar_evi_modissubset_", sitename, ".csv", sep="" )

  if (!file.exists(filnam_fapar_csv)||overwrite){

    ##--------------------------------------------------------------------
    df_modis <- get_evi_modis_250m( sitename, lon, lat )
    ##--------------------------------------------------------------------

    ##--------------------------------------------------------------------
    ## gap-fill data with median of corresponding months
    ##--------------------------------------------------------------------
    df_modis$flag <- rep(0,dim(df_modis)[1])
    for (idx in 1:dim(df_modis)[1]){
      if (is.na(df_modis$evi[idx])){
        df_modis$evi[idx] <- median( df_modis$evi[ df_modis$moy==df_modis$moy[idx] ], na.rm=TRUE )
        df_modis$flag[idx] <- 1
      }
    }

    ##--------------------------------------------------------------------
    ## add dummy year 1999 with median of each month in all subsequent years
    ##--------------------------------------------------------------------
    dummy_1999 <- df_modis[ 1:12, ]
    for (imo in 1:nmonth){
      dummy_1999$evi[imo] <- median( df_modis$evi[ df_modis$moy==imo ] )
      dummy_1999$yr[imo] <- 1999 
      dummy_1999$flag[imo] <- 1
    }
    df_modis <- rbind( dummy_1999, df_modis )

    ##--------------------------------------------------------------------
    ## Reduce data
    ##--------------------------------------------------------------------
    df_modis <- subset( df_modis, select=c(yr,moy,evi))

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing fapar data frame into CSV file ", filnam_fapar_csv, "...") )
    system( paste( "mkdir -p", dirnam_fapar_csv ) )   
    write.csv( df_modis, file=filnam_fapar_csv, row.names=FALSE )

  } else {

    ##--------------------------------------------------------------------
    ## Read CSV file directly for this site
    ##--------------------------------------------------------------------
    print( "monthly fapar data already available in CSV file ...")
    df_modis <- read.csv( filnam_fapar_csv )

  }

  ##--------------------------------------------------------------------
  ## Attach to data frame for all sites
  ##--------------------------------------------------------------------
  df_modis_allsites[[ sitename ]] <- df_modis$evi

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )

  ## in separate formatted file 
  for (yr in unique(df_modis$yr)){

    print( paste("... for year", yr))
    dirnam <- paste( myhome, "sofun/input_", simsuite, "_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "fapar_evi_modissubset_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, df_modis$evi[ which( df_modis$yr==yr ) ] )
    
  }
  
}

##--------------------------------------------------------------------
## Write CSV file for data frame for all sites
##--------------------------------------------------------------------
write.csv( df_modis_allsites, file=path_fapar_allsites_csv, row.names=FALSE )

