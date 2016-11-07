## This is to check whether input data for PPFD corresponds to original data in the files and how it compares across different sources for original data

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

simsuite <- "fluxnet2015"

## load meta data file for site simulation
siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ), as.is=TRUE )
nsites <- dim(siteinfo)[1]
# do.sites <- seq(nsites)
do.sites <- 49:49

for ( idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])

  ## prepared input
  meteo <- read.csv( paste( "../input_fluxnet2015_sofun/sitedata/climate/", sitename, "/clim_daily_byst_", sitename, ".csv", sep="" ), header=TRUE )
  meteo$year_dec <- meteo$year + (meteo$doy - 1)/365

  ## read output from SPLASH
  splash <- read.table( paste( "../trunk/output/", sitename, ".d.ppfd.out", sep="" ), col.names = c("year_dec","ppfd"))

  ## Plot PPFD from FLUXNET data, derived from SWin
  plot( meteo$year_dec[ !is.na(meteo$ppfd) ], meteo$ppfd[ !is.na(meteo$ppfd) ], type = "l" )
  lines( splash$year_dec, splash$ppfd, col="red", lwd=2 )
  title( sitename )

}