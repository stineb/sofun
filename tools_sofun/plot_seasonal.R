runname <- "RUNNAME"
outdir <- "/alphadata01/bstocker/sofun/trunk/output/"
dvars <- c("cex","cleaf","gpp","netmin","ninorg","npp","nup","nfixfree","ccost")
avars <- c("calloc","nalloc","clit2soil","nlit2soil","nreq","cveg2lit","nveg2lit")

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
for (ivar in seq(ndvars)) {
  filn <- paste(outdir,runname,".d.",dvars[ivar],".out",sep="")
  print(filn)
  daily[[ dvars[ivar] ]] <- read.table( filn )[,2]
}
daily[[ "year" ]] <- as.integer(time)
daily[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
daily[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

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
for (ivar in seq(navars)) {
  annual[[ avars[ivar] ]] <- read.table( paste(outdir,runname,".a.",avars[ivar],".out",sep=""))[,2]
}
annual[[ "year" ]] <- as.integer(time)

## derived efficiency of microbial decomposition
annual$eff <- annual$clit2soil / annual$cveg2lit

annual_sub <- annual[which(annual$year==plotyear),]

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

pdf(paste("Coverview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
# layout.show(panel)

par(mar=c(4,4,2,1), las=1)
plot( daily_sub$doy, daily_sub$gpp, type="l", xlab="", ylab="gC/m2/d", ylim=c(0,max(daily_sub$gpp)) )
lines( daily_sub$doy, daily_sub$npp, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
lines( daily_sub$doy, daily_sub$cex, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
legend( "bottomleft", c("GPP","NPP","EXU"), col=c("black","green","blue"), bty="n", lty=1 )

par(mar=c(0,0,2,0))
plot( c(0,1), c(0,1), type="n", axes=FALSE )
text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2)
text( 0, 0.93, "GPP", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
text( 0, 0.88, "NPP", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(sum(daily_sub$npp),digits=1,format="f")) , adj=c(1,0))
text( 0, 0.83, "CEX ", adj=c(0,0));         text( 0.5, 0.83, as.character(formatC(sum(daily_sub$cex),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.83, paste( as.character(formatC(sum(daily_sub$cex)/sum(daily_sub$npp)*100,digits=1,format="f")),"% of NPP",sep=""), adj=c(0,0))
text( 0, 0.78, "C -> veg ", adj=c(0,0));    text( 0.5, 0.78, as.character(formatC(annual_sub$calloc,digits=1,format="f")), adj=c(1,0))
text( 0, 0.73, "C veg -> lit", adj=c(0,0)); text( 0.5, 0.73, as.character(formatC(annual_sub$cveg2lit,digits=1,format="f")), adj=c(1,0))
text( 0, 0.68, "C lit -> soil", adj=c(0,0));text( 0.5, 0.68, as.character(formatC(annual_sub$clit2soil,digits=1,format="f")), adj=c(1,0)); text( 0.55, 0.68, paste(as.character(formatC(annual_sub$eff*100,digits=1,format="f")),"%",sep=""), adj=c(0,0))

dev.off()


pdf(paste("Noverview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
par(mar=c(4,4,2,1), las=1 )
plot( daily_sub$doy, daily_sub$ninorg, type="l", xlab="day of year", ylab="gN/m2/d", ylim=c(0,max(daily_sub$ninorg)) )
lines( daily_sub$doy, daily_sub$nup, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
lines( daily_sub$doy, daily_sub$netmin, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
lines( daily_sub$doy, daily_sub$nfixfree, type="l", xlab="day of year", ylab="gC/m2/d", col="red" )
legend( "bottomleft", c("Ninorg","Nup","net N min.","free-living BNF"), col=c("black","green","blue","red"), bty="n", lty=1 )

par(mar=c(0,0,2,0))
plot( c(0,1), c(0,1), type="n", axes=FALSE )
text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2)
text( 0, 0.93, "N uptake", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$nup),digits=2,format="f")) , adj=c(1,0))
text( 0, 0.88, "N -> veg", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(annual_sub$nalloc,digits=2,format="f")) , adj=c(1,0))
text( 0, 0.83, "N veg -> lit", adj=c(0,0));      text( 0.5, 0.83, as.character(formatC(annual_sub$nveg2lit,digits=2,format="f")) , adj=c(1,0))
dev.off()

