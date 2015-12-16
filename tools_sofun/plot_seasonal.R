# runname <- "CH-Oe1_2002"
namstat <- "CH-Oe1"
runname <- "test"
outdir <- "/alphadata01/bstocker/sofun/trunk/output/constalloc"
dvars <- c(
  "gpp","npp","cex","nup","nup_pas","nup_act","nup_fix","nup_ret","nfixfree",
  "cleaf","netmin","ninorg","ccost","cleaf","croot","clabl","clitt","nlitt","csoil","nsoil","netmin_litt",
  "netmin_soil","nloss","nvol","denitr","nitr","nleach","soiltemp","temp","lai","transp","ea_n"
  )
avars <- c("calclm","nalclm","clit2soil","nlit2soil","nreq","cveg2lit","nveg2lit")

plotyear <- 2000

nmonths <- 12
ndayyear <- 365
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndvars <- length(dvars)
navars <- length(avars)

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

# ##--------------------------------------
# ## read GPP data from observation (CH-Oe1_2002)
# ##--------------------------------------
#   # filn <- paste( outdir, 'save/CH-Oe1_2002.d.gpp.out', sep="" )
#   # gpp_obs_daily <- read.table( filn, col.names=c('time','gpp') )
#   gpp_obs_daily <- read.csv( paste( "../data/gepisat/gpp_daily/",namstat,"_daily_gpp.txt", sep="" ) ) 
#   nyears <- length(gpp_obs_daily$time)/ndayyear

#   obsdaily <- gpp_obs_daily
#   print("reading daily observation file ...")
#   obsdaily[[ "year" ]] <- as.integer(gpp_obs_daily$time)
#   obsdaily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
#   obsdaily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

#   ## take subset of one year and normalise all variables to [0,1]
#   obsdaily_sub <- obsdaily[obsdaily$year==plotyear,]


##--------------------------------------
## Plot absolute C variables
##--------------------------------------
aspect <- 0.5
magn <- 8
ncols <- 2
nrows <- 1
widths <- rep(magn,ncols)
widths[2] <- 0.4*magn
heights <- rep(aspect*magn,nrows)
# heights[nrows] <- 0.3*magn

##--------------------------------------
## C FLUXES
##--------------------------------------
  filn <- paste( "C_flux_overview_sofun_seasonal_", runname, ".pdf", sep="" )
  ylim <- c(0,max(daily_sub$gpp))
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)

  par(mar=c(4,4,2,1), las=1)
  # plot( obsdaily_sub$doy, obsdaily_sub$gpp, type="n", xlab="", ylab="gC/m2/d", ylim=c(0,max(daily_sub$gpp)) )
  plot( c(1,365), ylim, type="n", xlab="", ylab="gC/m2/d", ylim=ylim )
  lines( daily_sub$doy, daily_sub$gpp, col="red" )
  lines( daily_sub$doy, daily_sub$npp, col="green" )
  lines( daily_sub$doy, daily_sub$cex, col="blue" )
  legend( "bottomleft", c("GPP obs.","GPP","NPP","EXU"), col=c("black","red","green","blue"), bty="n", lty=1 )

  par(mar=c(0,0,2,0))
  plot( c(0,1), c(0,1), type="n", axes=FALSE )
  lab <- expression(paste("(gC m"^-2," yr"^-1,")", sep=""))
  text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
  # text( 0, 0.93, "GPP obs.", adj=c(0,0));     text( 0.5, 0.93, as.character(formatC(sum(obsdaily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
  text( 0, 0.88, "GPP", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(sum(daily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
  text( 0, 0.83, "NPP", adj=c(0,0));          text( 0.5, 0.83, as.character(formatC(sum(daily_sub$npp),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.83, paste( as.character(formatC(sum(daily_sub$npp)/sum(daily_sub$gpp)*100,digits=1,format="f")),"% of GPP",sep=""), adj=c(0,0))
  text( 0, 0.78, "CEX ", adj=c(0,0));         text( 0.5, 0.78, as.character(formatC(sum(daily_sub$cex),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.78, paste( as.character(formatC(sum(daily_sub$cex)/sum(daily_sub$npp)*100,digits=1,format="f")),"% of NPP",sep=""), adj=c(0,0))
  text( 0, 0.73, "C -> veg ", adj=c(0,0));    text( 0.5, 0.73, as.character(formatC(annual_sub$calloc,digits=1,format="f")), adj=c(1,0))
  text( 0, 0.68, "C veg -> lit", adj=c(0,0)); text( 0.5, 0.68, as.character(formatC(annual_sub$cveg2lit,digits=1,format="f")), adj=c(1,0))
  text( 0, 0.63, "C lit -> soil", adj=c(0,0));text( 0.5, 0.63, as.character(formatC(annual_sub$clit2soil,digits=1,format="f")), adj=c(1,0)); text( 0.55, 0.63, paste(as.character(formatC(annual_sub$eff*100,digits=1,format="f")),"% efficiency",sep=""), adj=c(0,0))

  dev.off()


##--------------------------------------
## N FLUXES
##--------------------------------------
  filn <- paste( "N_flux_overview_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  par(mar=c(4,4,2,1), las=1 )
  plot( daily_sub$doy, daily_sub$nup, type="l", xlab="day of year", ylab="gC/m2/d", ylim=range(c(daily_sub$nup, daily_sub$netmin_litt, daily_sub$netmin_soil, daily_sub$nfixfree) ), col="green" )
  lines( daily_sub$doy, daily_sub$netmin, col="blue" )
  lines( daily_sub$doy, daily_sub$netmin_litt, col="blue", lty=2 )
  lines( daily_sub$doy, daily_sub$netmin_soil, col="blue", lty=3 )
  lines( daily_sub$doy, daily_sub$nup_fix, col="cyan" )
  lines( daily_sub$doy, daily_sub$nfixfree, col="red" )
  lines( daily_sub$doy, daily_sub$nloss, col="magenta" )
  lines( daily_sub$doy, daily_sub$nvol, col="magenta", lty=2 )
  lines( daily_sub$doy, daily_sub$denitr, col="magenta", lty=3 )
  lines( daily_sub$doy, daily_sub$nitr, col="magenta", lty=4 )
  lines( daily_sub$doy, daily_sub$nleach, col="magenta", lty=5 )
  abline( h=0.0, col=rgb(0,0,0,0.3))
  legend( "topleft", c("Nup","net N min.","net N min., litter","net N min., soil","N fix., symbiotic","N fix., free-living","N loss tot.","N vol.","denitr.","nitr.","N leach"), col=c("green","blue","blue","blue","cyan","red","magenta","magenta","magenta","magenta","magenta"),lty=c(1,1,2,3,1,1,2,3,4,5), bty="n" )

  par(mar=c(0,0,2,0))
  plot( c(0,1), c(0,1), type="n", axes=FALSE )
  lab <- expression(paste("(gN m"^-2," yr"^-1,")", sep=""))
  text( 0.10, 0.98, "ANNUAL TOTALS" , adj=c(0,0),font=2); text( 0.60, 0.98, lab , adj=c(0,0),font=2, cex=0.7)
  text( 0.10, 0.93, "N uptake", adj=c(0,0));              text( 0.60, 0.93, as.character(formatC(sum(daily_sub$nup),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.88, "N up. pas.", adj=c(0,0), cex=1.0);            text( 0.60, 0.88, as.character(formatC(sum(daily_sub$nup_pas),digits=2,format="f")) , adj=c(1,0), cex=1.0)
  text( 0.10, 0.83, "N up. act.", adj=c(0,0), cex=1.0);            text( 0.60, 0.83, as.character(formatC(sum(daily_sub$nup_act),digits=2,format="f")) , adj=c(1,0), cex=1.0)
  text( 0.10, 0.78, "N up. fix.", adj=c(0,0), cex=1.0);            text( 0.60, 0.78, as.character(formatC(sum(daily_sub$nup_fix),digits=2,format="f")) , adj=c(1,0), cex=1.0)
  text( 0.10, 0.73, "N up. ret.", adj=c(0,0), cex=1.0);            text( 0.60, 0.73, as.character(formatC(sum(daily_sub$nup_ret),digits=2,format="f")) , adj=c(1,0), cex=1.0)


  text( 0.10, 0.68, "N -> veg", adj=c(0,0));              text( 0.60, 0.68, as.character(formatC(annual_sub$nalloc,digits=2,format="f")) , adj=c(1,0));     text( 0.65, 0.68, paste(as.character(formatC(annual_sub$calloc/annual_sub$nalloc,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0.10, 0.63, "N veg -> lit", adj=c(0,0));          text( 0.60, 0.63, as.character(formatC(annual_sub$nveg2lit,digits=2,format="f")) , adj=c(1,0));   text( 0.65, 0.63, paste(as.character(formatC(annual_sub$cveg2lit/annual_sub$nveg2lit,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0.10, 0.58, "N lit -> soil", adj=c(0,0));         text( 0.60, 0.58, as.character(formatC(annual_sub$nlit2soil,digits=2,format="f")) , adj=c(1,0));  text( 0.65, 0.58, paste(as.character(formatC(annual_sub$clit2soil/annual_sub$nlit2soil,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0.10, 0.53, "N -> soil req.", adj=c(0,0));        text( 0.60, 0.53, as.character(formatC(annual_sub$nreq,digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.48, "N fix.", adj=c(0,0));                text( 0.60, 0.48, as.character(formatC(sum(daily_sub$nfixfree),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.43, "N net min.", adj=c(0,0));            text( 0.60, 0.43, as.character(formatC(sum(daily_sub$netmin),digits=2,format="f")) , adj=c(1,0))

  text( 0.10, 0.38, "N loss tot.", adj=c(0,0));          text( 0.60, 0.38, as.character(formatC(sum(daily_sub$nloss),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.33, "N vol", adj=c(0,0));                text( 0.60, 0.33, as.character(formatC(sum(daily_sub$nvol),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.28, "N denitr.", adj=c(0,0));            text( 0.60, 0.28, as.character(formatC(sum(daily_sub$denitr),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.23, "N nitr.", adj=c(0,0));              text( 0.60, 0.23, as.character(formatC(sum(daily_sub$nitr),digits=2,format="f")) , adj=c(1,0))
  text( 0.10, 0.18, "N leach.", adj=c(0,0));             text( 0.60, 0.18, as.character(formatC(sum(daily_sub$nleach),digits=2,format="f")) , adj=c(1,0))

  dev.off()

##--------------------------------------
## C FLUXES
##--------------------------------------
  filn <- paste( "water_overview_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)
  xlim <- range( daily_sub$doy )
  ylim <- c( 0.0, max( daily_sub$transp*1e-3, daily_sub$ea_n ) )
  par(mar=c(4,4,2,1), las=1)
  plot( xlim, ylim, type="n", xlab="DOY", ylab="mm/d", ylim=ylim )

  lines( daily_sub$doy, daily_sub$transp * 1e-3, col="red" )
  lines( daily_sub$doy, daily_sub$ea_n, col="blue" ) # ea_n is given in mm, conversion to gH2O/m2

  legend( "topleft", c("T (P-model)","AET (SPLASH)"), col=c("red","blue"), bty="n", lty=1 )

  par(mar=c(0,0,2,0))
  plot( c(0,1), c(0,1), type="n", axes=FALSE )
  lab <- expression(paste("(mm"," yr"^-1,")", sep=""))
  text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
  text( 0, 0.93, "T (P-model) ", adj=c(0,0) );        text( 0.7, 0.93, as.character(formatC(sum(daily_sub$transp*1e-3),digits=1,format="f")) , adj=c(1,0))
  text( 0, 0.88, "AET (SPLASH)", adj=c(0,0) );        text( 0.7, 0.88, as.character(formatC(sum(daily_sub$ea_n),digits=1,format="f")) , adj=c(1,0))

  dev.off()


##--------------------------------------
## N COST ANALYSIS
##--------------------------------------
  filn <- paste( "CostOverview_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  par(mar=c(4,4,2,1), las=1 )
  plot( daily_sub$doy, daily_sub$ccost, xlab="", ylab="", type="l", col="magenta" )
  par(new=TRUE)
  plot(  daily_sub$doy, daily_sub$nup/daily_sub$ninorg, axes=FALSE, type="l", col="green" )
  lines( daily_sub$doy, daily_sub$nloss/daily_sub$ninorg, col="red" )
  axis( 4, col="green" )
  legend( "bottomleft", c("C price of N", "fraction N taken up", "fraction N lost"), col=c("magenta","green","red"), bty="n", lty=1 )

  par(mar=c(0,0,2,0))
  dev.off()


##--------------------------------------
## CUMULATIVE N FLUXES
##--------------------------------------
  filn <- paste( "N_flux_overview_cumulative_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  par(mar=c(4,4,2,1), las=1 )
  plot( daily_sub$doy,  daily_sub$ninorg, type="l", xlab="day of year", ylab="gN/m2/d", ylim=c(0,max(cumsum(daily_sub$netmin))) )
  lines( daily_sub$doy, cumsum(daily_sub$nup), type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
  lines( daily_sub$doy, cumsum(daily_sub$netmin), type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
  lines( daily_sub$doy, cumsum(daily_sub$nfixfree), type="l", xlab="day of year", ylab="gC/m2/d", col="red" )
  legend( "bottomleft", c("Ninorg","Nup","net N min.","free-living BNF"), col=c("black","green","blue","red"), bty="n", lty=1 )

  par(mar=c(0,0,2,0))
  plot( c(0,1), c(0,1), type="n", axes=FALSE )
  lab <- expression(paste("(gN m"^-2," yr"^-1,")", sep=""))
  text( 0, 0.98, "ANNUAL TOTALS" , adj=c(0,0),font=2); text( 0.5, 0.98, lab , adj=c(0,0),font=2, cex=0.7)
  text( 0, 0.93, "N uptake", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$nup),digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.88, "N up. pas.", adj=c(0,0));        text( 0.5, 0.88, as.character(formatC(sum(daily_sub$nup_pas),digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.83, "N up. act.", adj=c(0,0));        text( 0.5, 0.83, as.character(formatC(sum(daily_sub$nup_act),digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.78, "N up. fix.", adj=c(0,0));        text( 0.5, 0.78, as.character(formatC(sum(daily_sub$nup_fix),digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.73, "N up. ret.", adj=c(0,0));        text( 0.5, 0.73, as.character(formatC(sum(daily_sub$nup_ret),digits=2,format="f")) , adj=c(1,0))


  text( 0, 0.68, "N -> veg", adj=c(0,0));          text( 0.5, 0.68, as.character(formatC(annual_sub$nalloc,digits=2,format="f")) , adj=c(1,0));     text( 0.55, 0.68, paste(as.character(formatC(annual_sub$calloc/annual_sub$nalloc,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0, 0.63, "N veg -> lit", adj=c(0,0));      text( 0.5, 0.63, as.character(formatC(annual_sub$nveg2lit,digits=2,format="f")) , adj=c(1,0));   text( 0.55, 0.63, paste(as.character(formatC(annual_sub$cveg2lit/annual_sub$nveg2lit,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0, 0.58, "N lit -> soil", adj=c(0,0));     text( 0.5, 0.58, as.character(formatC(annual_sub$nlit2soil,digits=2,format="f")) , adj=c(1,0));  text( 0.55, 0.58, paste(as.character(formatC(annual_sub$clit2soil/annual_sub$nlit2soil,digits=2,format="f")),"C:N") , adj=c(0,0))
  text( 0, 0.53, "N -> soil req.", adj=c(0,0));    text( 0.5, 0.53, as.character(formatC(annual_sub$nreq,digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.48, "N fBNF", adj=c(0,0));            text( 0.5, 0.48, as.character(formatC(sum(daily_sub$nfixfree),digits=2,format="f")) , adj=c(1,0))
  text( 0, 0.43, "N net min.", adj=c(0,0));        text( 0.5, 0.43, as.character(formatC(sum(daily_sub$netmin),digits=2,format="f")) , adj=c(1,0))
  dev.off()


##--------------------------------------
## C POOLS
##--------------------------------------
  filn <- paste( "C_pool_overview_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)

  par(mar=c(4,4,2,1), las=1)
  plot( daily_sub$doy, daily_sub$cleaf, type="l", xlab="", ylab="gC/m2", ylim=range(c(daily_sub$cleaf,daily_sub$croot)) )
  lines( daily_sub$doy, daily_sub$croot, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
  lines( daily_sub$doy, daily_sub$clabl, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
  legend( "bottomleft", c("leaf C","root C","labile C"), col=c("black","green","blue"), bty="n", lty=1 )

  # par(mar=c(0,0,2,0))
  # plot( c(0,1), c(0,1), type="n", axes=FALSE )
  # lab <- expression(paste("(gC m"^-2," yr"^-1,")", sep=""))
  # text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
  # text( 0, 0.93, "GPP", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
  # text( 0, 0.88, "NPP", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(sum(daily_sub$npp),digits=1,format="f")) , adj=c(1,0))
  # text( 0, 0.83, "CEX ", adj=c(0,0));         text( 0.5, 0.83, as.character(formatC(sum(daily_sub$cex),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.83, paste( as.character(formatC(sum(daily_sub$cex)/sum(daily_sub$npp)*100,digits=1,format="f")),"% of NPP",sep=""), adj=c(0,0))
  # text( 0, 0.78, "C -> veg ", adj=c(0,0));    text( 0.5, 0.78, as.character(formatC(annual_sub$calloc,digits=1,format="f")), adj=c(1,0))
  # text( 0, 0.73, "C veg -> lit", adj=c(0,0)); text( 0.5, 0.73, as.character(formatC(annual_sub$cveg2lit,digits=1,format="f")), adj=c(1,0))
  # text( 0, 0.68, "C lit -> soil", adj=c(0,0));text( 0.5, 0.68, as.character(formatC(annual_sub$clit2soil,digits=1,format="f")), adj=c(1,0)); text( 0.55, 0.68, paste(as.character(formatC(annual_sub$eff*100,digits=1,format="f")),"% efficiency",sep=""), adj=c(0,0))

  dev.off()

##--------------------------------------
## CLIMATE
##--------------------------------------
  filn <- paste( "climate_overview_sofun_seasonal_", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)

  par(mar=c(4,4,2,1), las=1)
  plot( daily_sub$doy, daily_sub$temp, type="l", xlab="", ylab="gC/m2/d", ylim=range( daily_sub$temp ) )
  lines( daily_sub$doy, daily_sub$soiltemp, col="red" )
  legend( "bottomleft", c("air temp.","soiltemp"), col=c("black","red"), bty="n", lty=1 )

  par(mar=c(0,0,2,0))
  plot( c(0,1), c(0,1), type="n", axes=FALSE )
  # lab <- expression(paste("(gC m"^-2," yr"^-1,")", sep=""))
  # text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
  # text( 0, 0.93, "GPP obs.", adj=c(0,0));     text( 0.5, 0.93, as.character(formatC(sum(obsdaily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
  
  dev.off()


##--------------------------------------
## LITTER/SOIL C and N POOLS
##--------------------------------------
  filn <- paste( "cn_littersom_overview_sofun_seasonal", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  par(mar=c(4,4,2,1), las=1 )
  plot(  daily_sub$doy,  daily_sub$ninorg, type="l", xlab="day of year", ylab="gN/m2", ylim=c(0,max(cumsum(daily_sub$netmin))), col="blue" )
  lines( daily_sub$doy,  daily_sub$nlitt - min(daily_sub$nlitt), col="brown" )
  lines( daily_sub$doy,  daily_sub$nsoil - min(daily_sub$nsoil), col="red" )

  legend( "topleft", c("N inorg","N litt. change","N soil change"), col=c("blue","brown","red"), bty="n", lty=1 )

  dev.off()


##--------------------------------------
## LITTER/SOIL C and N POOLS
##--------------------------------------
  filn <- paste( "lai_overview_sofun_seasonal", runname, ".pdf", sep="" )
  pdf( filn, width=sum(widths), height=sum(heights) )
  panel <- layout(
                  matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  par(mar=c(4,4,2,1), las=1 )
  plot(  daily_sub$doy,  daily_sub$lai, type="l", xlab="day of year", ylab="LAI", ylim=c(0,max(daily_sub$lai)), col="red" )

  legend( "topleft", c("LAI"), col=c("red"), bty="n", lty=1 )

  dev.off()


