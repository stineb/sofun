## This is to check whether input data for PPFD corresponds to original data in the files and how it compares across different sources for original data
library( dplyr )
syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )
source( "init_daily_dataframe.R" )
source( "add_year_dec.R" )

simsuite <- "fluxnet2015"

## load meta data file for site simulation
siteinfo <- read.csv( paste( myhome, "sofun/input_", simsuite, "_sofun/siteinfo_", simsuite, "_sofun.csv", sep="" ), as.is=TRUE )
nsites <- dim(siteinfo)[1]
do.sites <- seq(nsites)
# do.sites <- 1:10


pdf( "compare_evi.pdf", width=15, height=6 )
for ( idx in do.sites ){

  sitename <- as.character(siteinfo$mysitename[idx])

  devi <- init_daily_dataframe( 1999, 2016 )

  ## old monthly EVI input (as before)
  mevi_old <- read.csv( paste( "../input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/fapar_evi_modissubset_", sitename, ".csv", sep="" ), header=TRUE )
  devi$mold <- rep( NA, dim(devi)[1] )
  for (jdx in 1:dim(mevi_old)[1]){
    devi$mold[ which( devi$year == mevi_old$yr[jdx] & devi$moy == mevi_old$moy[jdx] ) ] <- mevi_old$evi[jdx]
  }

  ## new monthly EVI
  mevi_new <- read.csv( paste( "../input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/mfapar_evi_modissubset_", sitename, ".csv", sep="" ), header=TRUE )
  devi$mnew <- rep( NA, dim(devi)[1] )
  for (jdx in 1:dim(mevi_new)[1]){
    devi$mnew[ which( devi$year == mevi_new$year[jdx] & devi$moy == mevi_new$moy[jdx] ) ] <- mevi_new$evi[jdx]
  }

  ## new daily EVI
  devi_new <- read.csv( paste( "../input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/dfapar_evi_modissubset_", sitename, ".csv", sep="" ), header=TRUE )
  devi$dnew <- rep( NA, dim(devi)[1] )
  for (jdx in 1:dim(devi_new)[1]){
    devi$dnew[ which( devi$year == devi_new$year[jdx] & devi$moy == devi_new$moy[jdx] & devi$doy == devi_new$doy[jdx] ) ] <- devi_new$evi[jdx]
  }

  ## monthly fAPAR3g
  mevi_fapar3g <- try( read.csv( paste( "../input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/fapar_fapar3g_", sitename, ".csv", sep="" ), header=TRUE ) )
  if (class(mevi_fapar3g)!="try-error"){
    devi$mevi_fapar3g <- rep( NA, dim(devi)[1] )
    for (jdx in 1:dim(mevi_fapar3g)[1]){
      devi$mevi_fapar3g[ which( devi$year == mevi_fapar3g$year[jdx] & devi$moy == mevi_fapar3g$mo[jdx] ) ] <- mevi_fapar3g$fapar[jdx]
    }
  }

  ## monthly MODIS-fAPAR
  mevi_faparmodis <- try( read.csv( paste( "../input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/fapar_modis_", sitename, ".csv", sep="" ), header=TRUE ))
  if (class(mevi_faparmodis)!="try-error"){
    devi$mevi_faparmodis <- rep( NA, dim(devi)[1] )
    for (jdx in 1:dim(mevi_faparmodis)[1]){
      devi$mevi_faparmodis[ which( devi$year == mevi_faparmodis$year[jdx] & devi$moy == mevi_faparmodis$mo[jdx] ) ] <- mevi_faparmodis$fapar[jdx]
    }
  }

  ## Plot PPFD from FLUXNET data, derived from SWin
  par( las=1 )
  plot( devi$year_dec, devi$mold, type = "l", ylim=c(0,max(devi$dnew, na.rm=TRUE)), xlab="year", ylab="EVI" )
  if (class(mevi_fapar3g)!="try-error"){
    lines( devi$year_dec, devi$mevi_fapar3g, col="blue", lwd=1 )
    legend( "bottom", c("FAPAR3G"), col=c("blue"), lty=1, bty="n" )
  }
  if (class(mevi_faparmodis)!="try-error"){
    lines( devi$year_dec, devi$mevi_faparmodis, col="tomato", lwd=1 )
    legend( "bottomright", c("MODIS 1deg"), col=c("tomato"), lty=1, bty="n" )
  }

  lines( devi$year_dec, devi$mnew, col="red", lwd=1 )
  lines( devi$year_dec, devi$dnew, col="green", lwd=2 )
  title( sitename )
  legend( "bottomleft", c("monthly EVI, 250m, old", "monthly EVI, 250m, new", "daily EVI, 250m"), col=c("black", "red", "green"), lty=1, bty="n", lwd=c(1,1,2) )

}
dev.off()
