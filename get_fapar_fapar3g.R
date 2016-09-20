library(plyr)
#library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/getin/get_pointdata_fapar3g.R", sep="" ) )
source( paste( myhome, "sofun/getin/write_sofunformatted.R", sep="" ) )

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)

simsuite <- "atkinfull"

fapar_year_start <- 1982
fapar_year_end   <- 2011
nyears <- fapar_year_end - fapar_year_start + 1

siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite,"_sofun/siteinfo_", simsuite,"_sofun.csv", sep="" ) )
nsites <- dim(siteinfo)[1]

## create data frame holding data for all sites
df_fapar_allsites <- data.frame( year=rep( seq( fapar_year_start, fapar_year_end ), each=nmonth ), mo=rep(1:nmonth,nyears) )
path_fapar_allsites_csv <- paste( myhome, "sofun/input_", simsuite,"_sofun/sitedata/fapar/fapar_fapar3g_allsites.csv", sep="" )

overwrite <- TRUE

for (idx in seq(nsites)){
# for (idx in 1:1){

  sitename <- as.character(siteinfo$mysitename[idx])
  lon      <- siteinfo$lon[idx]
  lat      <- siteinfo$lat[idx]
  print( paste( "collecting monthly data for station", sitename, "..." ) )

  dirnam_fapar_csv <- paste( myhome, "sofun/input_", simsuite,"_sofun/sitedata/fapar/", sitename, "/", sep="" )
  filnam_fapar_csv <- paste( dirnam_fapar_csv, "fapar_fapar3g_", sitename, ".csv", sep="" )

  df_fapar <- data.frame( year=rep( seq( fapar_year_start, fapar_year_end ), each=nmonth ), mo=rep(1:nmonth,nyears) )

  if (!file.exists(filnam_fapar_csv)||overwrite){

    ##--------------------------------------------------------------------
    ## get data for all years and months (360 months covering 1982-2011)
    ##--------------------------------------------------------------------
    df_fapar$fapar <- get_pointdata_fapar3g( lon, lat )

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
    dirnam <- paste( myhome, "sofun/input_", simsuite,"_sofun/sitedata/fapar/", sitename, "/", as.character(yr), "/", sep="" )
    system( paste( "mkdir -p", dirnam ) )

    filnam <- paste( dirnam, "fapar_fapar3g_", sitename, "_", yr, ".txt", sep="" )
    write_sofunformatted( filnam, df_fapar$fapar[ which( df_fapar$year==yr ) ] )
    
  }
  
}

##--------------------------------------------------------------------
## Write CSV file for data frame for all sites
##--------------------------------------------------------------------
write.csv( df_fapar_allsites, file=path_fapar_allsites_csv, row.names=FALSE )

