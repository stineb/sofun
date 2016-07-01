library(plyr)
#library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( "get_pointdata_fapar_modis.R" )
source( "write_sofunformatted.R" )

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

fapar_year_start <- 2001
fapar_year_end   <- 2013
nyears <- fapar_year_end - fapar_year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet_sofun/siteinfo_fluxnet_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

## create data frame holding data for all sites
df_fapar_allsites <- data.frame( year=rep( seq( fapar_year_start, fapar_year_end ), each=nmonth ), mo=rep(1:nmonth,nyears) )
path_fapar_allsites_csv <- paste( myhome, "sofun/input_fluxnet_sofun/sitedata/fapar/fapar_modis_allsites.csv", sep="" )

overwrite <- FALSE

for (idx in seq(nsites)){
# for (idx in 1:2){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_fapar_csv <- paste( myhome, "sofun/input_fluxnet_sofun/sitedata/fapar/", sitename, "/", sep="" )
  filnam_fapar_csv <- paste( dirnam_fapar_csv, "fapar_modis_", sitename, ".csv", sep="" )

  df_fapar <- data.frame( year=rep( seq( fapar_year_start, fapar_year_end ), each=nmonth ), mo=rep(1:nmonth,nyears), fapar=rep(NA,nmonth*nyears) )

  if (!file.exists(filnam_fapar_csv)||overwrite){

    ##--------------------------------------------------------------------
    ## get annual MODIS fAPAR data
    ##--------------------------------------------------------------------
    idx <- 0
    for (yr in fapar_year_start:fapar_year_end){
      df_fapar$fapar[(idx+1):(idx+nmonth)] <- get_pointdata_fapar_modis( lon, lat,  yr )
      idx <- idx + nmonth
    }

    ##--------------------------------------------------------------------
    ## For MODIS data, add dummy years 2000 and 2014, use median of all 
    ## available years by month.
    ##--------------------------------------------------------------------
    tmp <- rep(NA,nmonth)
    for (imo in 1:nmonth){
      tmp[imo] <- median( df_fapar$fapar[ df_fapar$mo==imo ] )
    }
    dummy_2000 <- data.frame( year=2000, mo=1:nmonth, fapar=tmp )
    dummy_2014 <- data.frame( year=2014, mo=1:nmonth, fapar=tmp )
    df_fapar <- rbind( dummy_2000, df_fapar, dummy_2014 )

    ##--------------------------------------------------------------------
    ## Save data frames as CSV files
    ##--------------------------------------------------------------------
    print( paste( "writing fapar data frame into CSV file ", filnam_fapar_csv, "...") )
    system( paste( "mkdir -p", dirnam_fapar_csv ) )   
    write.csv( df_fapar, file=filnam_fapar_csv, row.names=FALSE )

  } else {

    print( "monthly fapar data already available in CSV file ...")
    df_fapar <- read.csv( filnam_fapar_csv )

  }

  ##--------------------------------------------------------------------
  ## Attach to data frame for all sites
  ##--------------------------------------------------------------------
  df_fapar_allsites[[ sitename ]] <- df_fapar$fapar

  ##--------------------------------------------------------------------
  ## Write to Fortran-formatted output for each variable and year separately
  ##--------------------------------------------------------------------
  print( "writing formatted input files ..." )

  ## in separate formatted file 
  for (yr in unique(df_fapar$year)){

    print( paste("... for year", yr))
    dirnam <- paste( myhome, "sofun/input_fluxnet_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "fapar_modis_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, df_fapar$fapar[ which( df_fapar$year==yr ) ] )
    
  }
  
}

##--------------------------------------------------------------------
## Write CSV file for data frame for all sites
##--------------------------------------------------------------------
write.csv( df_fapar_allsites, file=path_fapar_allsites_csv, row.names=FALSE )

