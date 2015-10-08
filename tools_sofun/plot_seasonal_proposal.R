runname <- "test"
outdir <- "/alphadata01/bstocker/sofun/trunk/output/"
dvars <- c(
  "gpp","npp","cex","nup","nup_pas","nup_act","nup_fix","nup_ret","nfixfree",
  "cleaf","netmin","ninorg","ccost","cleaf","croot","clabl","nlabl","clitt","nlitt","csoil","nsoil","netmin_litt",
  "netmin_soil","nloss","nvol","denitr","nitr","nleach","soiltemp","temp"
  )
avars <- c("calloc","nalloc","clit2soil","nlit2soil","nreq","cveg2lit","nveg2lit")

plotyear <- 2000

nmonths <- 12
ndayyear <- 365
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndvars <- length(dvars)
navars <- length(avars)

##--------------------------------------
## read data from "imbalance analysis"
##--------------------------------------
imbal <- read.csv( "imbalance_analysis.csv")
## try:
imbal$labile.C[7] <- 23.0

##--------------------------------------
## construct data frame for daily values
##--------------------------------------
  ## get time data
  filn <- paste(outdir,runname,".d.",dvars[1],".out",sep="")
  tmp <- read.table( filn )
  time <- tmp$V1
  nyears <- length(time)/ndayyear

  daily <- data.frame( time = time )
  print("reading daily output files ...")
  for (ivar in seq(ndvars)) {
    filn <- paste(outdir,runname,".d.",dvars[ivar],".out",sep="")
    print(filn)
    daily[[ dvars[ivar] ]] <- read.table( filn )[,2]
  }
  daily[[ "year" ]] <- as.integer(time)
  daily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
  daily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

  ## remove dummy value -9999 with NA
  daily[daily==-9999] <- NA

  ## take subset of one year and normalise all variables to [0,1]
  daily_sub <- daily[daily$year==plotyear,]
  for (ivar in seq(ndvars)) {
    daily_sub[[ paste(dvars[ivar],".norm", sep="") ]] <- daily_sub[[ dvars[ivar] ]] / max( daily_sub[[ dvars[ivar] ]], na.rm=TRUE )
  }

##--------------------------------------
## construct data frame for annual values
##--------------------------------------
  ## get time data
  tmp <- read.table( paste(outdir,runname,".a.",avars[1],".out",sep=""))
  time <- tmp$V1
  nyears <- length(time)

  annual <- data.frame( time = time )
  print("reading annual output files ...")
  for (ivar in seq(navars)) {
    filn <- paste(outdir,runname,".a.",avars[ivar],".out",sep="")
    print(filn)
    annual[[ avars[ivar] ]] <- read.table( filn )[,2]
  }
  annual[[ "year" ]] <- as.integer(time)

  ## derived efficiency of microbial decomposition
  annual$eff <- annual$clit2soil / annual$cveg2lit

  annual_sub <- annual[which(annual$year==plotyear),]

##--------------------------------------
## read GPP data from observation (CH-Oe1_2002)
##--------------------------------------
  filn <- paste( outdir, 'save/CH-Oe1_2002.d.gpp.out', sep="" )
  gpp_obs_daily <- read.table( filn, col.names=c('time','gpp') )
  nyears <- length(gpp_obs_daily$time)/ndayyear

  obsdaily <- gpp_obs_daily
  print("reading daily observation file ...")
  obsdaily[[ "year" ]] <- as.integer(gpp_obs_daily$time)
  obsdaily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
  obsdaily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

  ## take subset of one year and normalise all variables to [0,1]
  obsdaily_sub <- obsdaily[obsdaily$year==plotyear,]


##--------------------------------------
## Plot parameters
##--------------------------------------
aspect <- 0.9
magn <- 3
ncols <- 3
nrows <- 1
widths <- rep(magn,ncols)
widths[ncols] <- widths[1]*1.2
heights <- rep(aspect*magn,nrows)
# heights[nrows] <- 0.3*magn


  # filn <- paste( "SOFUN_overview_for_proposal_small", runname, ".pdf", sep="" )
  # pdf( filn, width=4, height=2.8 ) 

  filn <- paste( "SOFUN_overview_for_proposal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)

  ##--------------------------------------
  ## C FLUXES
  ##--------------------------------------
  par(mar=c(3,4,2,0), las=1)
  plot( obsdaily_sub$doy, obsdaily_sub$gpp, type="l", xlab="", ylab="", ylim=c(0,1.0), col="grey50" )
  mtext( "day of year", side=1, line=2, las=1, cex=0.7 )
  mtext( expression( "gC" ~ m^{-2} ~ d^{-1} ), side=2, line=2.5, las=0, cex=0.7 )
  lines( daily_sub$doy, daily_sub$gpp, col="red" )
  lines( daily_sub$doy, daily_sub$npp, col="magenta" )
  lines( daily_sub$doy, daily_sub$npp-daily_sub$cex, col="green" )
  lines( daily_sub$doy, daily_sub$cex, col="steelblue3" )
  legend( "topleft", c("GPP observations","GPP","NPP","BP",expression(paste("C"[ex]))), col=c("grey50","red","magenta","green","steelblue3"), bty="n", lty=1 )
  title("Carbon fluxes")

  ##--------------------------------------
  ## N FLUXES
  ##--------------------------------------
  par(mar=c(3,4,2,0), las=1 )
  plot( daily_sub$doy, daily_sub$nup*1000, type="l", xlab="", ylab="", ylim=range(c(daily_sub$nup, daily_sub$netmin_litt, daily_sub$netmin_soil, daily_sub$nfixfree) )*1000, col="green" )
  mtext( "day of year", side=1, line=2, las=1, cex=0.7 )
  mtext( expression( "mgN" ~ m^{-2} ~ d^{-1} ), side=2, line=2, las=0, cex=0.7 )
  lines( daily_sub$doy, daily_sub$netmin*1000, col="steelblue3" )
  lines( daily_sub$doy, daily_sub$nfixfree*1000, col="red" )
  lines( daily_sub$doy, daily_sub$nloss*1000, col="magenta" )
  legend( "topleft", c("N uptake","net N mineralisation","N deposition","N loss"), col=c("green","steelblue3","red","magenta"),lty=1, bty="n" )
  title("Nitrogen fluxes")

  ##--------------------------------------
  ## IMBALANCE ANALYSIS
  ##--------------------------------------
  par(mar=c(3,4,2,4), las=1 )

  plot( imbal$f_shoot, imbal$labile.N, type="l", col="steelblue3", axes=FALSE, xlab="", ylab="" )
  abline( v=0.655, col="red")
  text( 0.71, 0.6, "optimum", col="red")
  mtext( "allocation fraction to leaves", side=1, line=2, las=1, cex=0.7 )
  box()
  axis(1)
  axis(2, col="steelblue3", xlab="", col.axis="steelblue3", las=1, ylab="labile C" )
  # mtext( "labile C", side=2, col="green", line=4 ) 
  mtext( expression("labile N (gN" ~ m^{-2} ~ ")"), side=2, line=2.5, las=0, cex=0.7, col="steelblue3" )
  par( new=TRUE )

  plot( imbal$f_shoot, imbal$labile.C, type="l", col="green", axes=FALSE, xlab="", ylab="", xlim=c(0.4,0.8))
  axis(4, col="green", ylab="labile N", col.axis="green", col.lab="green",las=1 )
  mtext( expression("labile C (gC" ~ m^{-2} ~ ")"), side=4, line=2.5, las=0, cex=0.7, col="green" )
  title("Labile pool size at end of growing season", cex=0.7)

  # ##--------------------------------------
  # ## IMBALANCE ANALYSIS
  # ##--------------------------------------
  # par(mar=c(3,4,1,4), las=1 )

  # plot( imbal$f_shoot, imbal$labile.N, type="l", col="steelblue3", axes=FALSE, xlab="", ylab="", lwd=2 )
  # abline( v=0.655, col="red", lwd=2)
  # text( 0.73, 0.64, "optimum", col="red")
  # mtext( "allocation fraction to leaves", side=1, line=2, las=1, cex=1.0 )
  # box()
  # axis(1)
  # axis(2, col="steelblue3", xlab="", col.axis="steelblue3", las=1, ylab="labile C" )
  # # mtext( "labile C", side=2, col="springgreen3", line=4 ) 
  # mtext( expression("labile N (gN" ~ m^{-2} ~ ")"), side=2, line=2.5, las=0, cex=1.0, col="steelblue3" )
  # par( new=TRUE )

  # plot( imbal$f_shoot, imbal$labile.C, type="l", col="springgreen3", axes=FALSE, xlab="", ylab="", xlim=c(0.4,0.8), lwd=2)
  # axis(4, col="springgreen3", ylab="labile N", col.axis="springgreen3", col.lab="springgreen3",las=1 )
  # mtext( expression("labile C (gC" ~ m^{-2} ~ ")"), side=4, line=2.5, las=0, cex=1.0, col="springgreen3" )
  # # title("Labile pool size at end of growing season", font.main=1)

  dev.off()
