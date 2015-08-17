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

# ##--------------------------------------
# ## construct data frame for daily values
# ##--------------------------------------
#   ## get time data
#   filn <- paste(outdir,runname,".d.",dvars[1],".out",sep="")
#   tmp <- read.table( filn )
#   time <- tmp$V1
#   nyears <- length(time)/ndayyear

#   daily <- data.frame( time = time )
#   print("reading daily output files ...")
#   for (ivar in seq(ndvars)) {
#     filn <- paste(outdir,runname,".d.",dvars[ivar],".out",sep="")
#     print(filn)
#     daily[[ dvars[ivar] ]] <- read.table( filn )[,2]
#   }
#   daily[[ "year" ]] <- as.integer(time)
#   daily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
#   daily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

#   ## remove dummy value -9999 with NA
#   daily[daily==-9999] <- NA

#   ## take subset of one year and normalise all variables to [0,1]
#   daily_sub <- daily[daily$year==plotyear,]
#   for (ivar in seq(ndvars)) {
#     daily_sub[[ paste(dvars[ivar],".norm", sep="") ]] <- daily_sub[[ dvars[ivar] ]] / max( daily_sub[[ dvars[ivar] ]], na.rm=TRUE )
#   }

# ##--------------------------------------
# ## construct data frame for annual values
# ##--------------------------------------
#   ## get time data
#   tmp <- read.table( paste(outdir,runname,".a.",avars[1],".out",sep=""))
#   time <- tmp$V1
#   nyears <- length(time)

#   annual <- data.frame( time = time )
#   print("reading annual output files ...")
#   for (ivar in seq(navars)) {
#     filn <- paste(outdir,runname,".a.",avars[ivar],".out",sep="")
#     print(filn)
#     annual[[ avars[ivar] ]] <- read.table( filn )[,2]
#   }
#   annual[[ "year" ]] <- as.integer(time)

#   ## derived efficiency of microbial decomposition
#   annual$eff <- annual$clit2soil / annual$cveg2lit

#   annual_sub <- annual[which(annual$year==plotyear),]

# ##--------------------------------------
# ## read GPP data from observation (CH-Oe1_2002)
# ##--------------------------------------
#   filn <- paste( outdir, 'save/CH-Oe1_2002.d.gpp.out', sep="" )
#   gpp_obs_daily <- read.table( filn, col.names=c('time','gpp') )
#   nyears <- length(gpp_obs_daily$time)/ndayyear

#   obsdaily <- gpp_obs_daily
#   print("reading daily observation file ...")
#   obsdaily[[ "year" ]] <- as.integer(gpp_obs_daily$time)
#   obsdaily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
#   obsdaily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

#   ## take subset of one year and normalise all variables to [0,1]
#   obsdaily_sub <- obsdaily[obsdaily$year==plotyear,]


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
  lines( daily_sub$doy, daily_sub$cex, col="blue" )
  legend( "topleft", c("GPP observations","GPP","NPP","BP",expression(paste("C"[ex]))), col=c("grey50","red","magenta","green","blue"), bty="n", lty=1 )
  title("Carbon fluxes")

  ##--------------------------------------
  ## N FLUXES
  ##--------------------------------------
  par(mar=c(3,4,2,0), las=1 )
  plot( daily_sub$doy, daily_sub$nup*1000, type="l", xlab="", ylab="", ylim=range(c(daily_sub$nup, daily_sub$netmin_litt, daily_sub$netmin_soil, daily_sub$nfixfree) )*1000, col="green" )
  mtext( "day of year", side=1, line=2, las=1, cex=0.7 )
  mtext( expression( "mgN" ~ m^{-2} ~ d^{-1} ), side=2, line=2, las=0, cex=0.7 )
  lines( daily_sub$doy, daily_sub$netmin*1000, col="blue" )
  lines( daily_sub$doy, daily_sub$nfixfree*1000, col="red" )
  lines( daily_sub$doy, daily_sub$nloss*1000, col="magenta" )
  legend( "topleft", c("N uptake","net N mineralisation","N deposition","N loss"), col=c("green","blue","red","magenta"),lty=1, bty="n" )
  title("Nitrogen fluxes")

  ##--------------------------------------
  ## IMBALANCE ANALYSIS
  ##--------------------------------------
  par(mar=c(3,5,2,4), las=1 )

  plot( imbal$f_shoot, imbal$labile.N, type="l", col="blue", axes=FALSE, xlab="", ylab="" )
  mtext( "allocation fraction to leaves", side=1, line=2, las=1, cex=0.7 )
  box()
  axis(1)
  axis(2, col="blue", xlab="", col.axis="blue", las=1, ylab="labile C" )
  # mtext( "labile C", side=2, col="green", line=4 ) 
  mtext( expression("labile N (gN" ~ m^{-2} ~ ")"), side=2, line=2.5, las=0, cex=0.7, col="blue" )
  par( new=TRUE )

  plot( imbal$f_shoot, imbal$labile.C, type="l", col="green", axes=FALSE, xlab="", ylab="", xlim=c(0.4,0.8))
  axis(4, col="green", ylab="labile N", col.axis="green", col.lab="green",las=1 )
  mtext( expression("labile C (gC" ~ m^{-2} ~ ")"), side=4, line=2.5, las=0, cex=0.7, col="green" )

  title("Labile pool size at end of growing season")

  dev.off()

# ##--------------------------------------
# ## C POOLS
# ##--------------------------------------
#   filn <- paste( "C_pool_overview_sofun_seasonal_", runname, ".pdf", sep="" )
#   pdf( filn, width=sum(widths), height=sum(heights) )
#   panel <- layout(
#                   matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
#                   widths=widths,
#                   heights=heights,
#                   TRUE
#                   )
#   # layout.show(panel)

#   par(mar=c(4,4,2,1), las=1)
#   plot( daily_sub$doy, daily_sub$cleaf, type="l", xlab="", ylab="gC/m2/d", ylim=range(c(daily_sub$cleaf,daily_sub$croot)) )
#   lines( daily_sub$doy, daily_sub$croot, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
#   lines( daily_sub$doy, daily_sub$clabl, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
#   legend( "bottomleft", c("leaf C","root C","labile C"), col=c("black","green","blue"), bty="n", lty=1 )

#   # par(mar=c(0,0,2,0))
#   # plot( c(0,1), c(0,1), type="n", axes=FALSE )
#   # lab <- expression(paste("(gC m"^-2," yr"^-1,")", sep=""))
#   # text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
#   # text( 0, 0.93, "GPP", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
#   # text( 0, 0.88, "NPP", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(sum(daily_sub$npp),digits=1,format="f")) , adj=c(1,0))
#   # text( 0, 0.83, "CEX ", adj=c(0,0));         text( 0.5, 0.83, as.character(formatC(sum(daily_sub$cex),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.83, paste( as.character(formatC(sum(daily_sub$cex)/sum(daily_sub$npp)*100,digits=1,format="f")),"% of NPP",sep=""), adj=c(0,0))
#   # text( 0, 0.78, "C -> veg ", adj=c(0,0));    text( 0.5, 0.78, as.character(formatC(annual_sub$calloc,digits=1,format="f")), adj=c(1,0))
#   # text( 0, 0.73, "C veg -> lit", adj=c(0,0)); text( 0.5, 0.73, as.character(formatC(annual_sub$cveg2lit,digits=1,format="f")), adj=c(1,0))
#   # text( 0, 0.68, "C lit -> soil", adj=c(0,0));text( 0.5, 0.68, as.character(formatC(annual_sub$clit2soil,digits=1,format="f")), adj=c(1,0)); text( 0.55, 0.68, paste(as.character(formatC(annual_sub$eff*100,digits=1,format="f")),"% efficiency",sep=""), adj=c(0,0))

#   dev.off()
