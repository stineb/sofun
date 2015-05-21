# runname <- "CH-Oe1_2002"
runname <- "test_lastworking"
outdir <- "/alphadata01/bstocker/sofun/trunk/output/"
dvars <- c("gpp","npp","cex","nup","nup_pas","nup_act","nup_fix","nup_ret","nfixfree","cleaf","netmin","ninorg","ccost","cleaf","croot","clabl","clitt")
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
pdf(paste("C_flux_overview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
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
lab <- expression(paste("(gC m"^-2," yr"^-1,")", sep=""))
text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2); text( 0.5, 0.98, lab, adj=c(0,0),font=2, cex=0.7 )
text( 0, 0.93, "GPP", adj=c(0,0));          text( 0.5, 0.93, as.character(formatC(sum(daily_sub$gpp),digits=1,format="f")) , adj=c(1,0))
text( 0, 0.88, "NPP", adj=c(0,0));          text( 0.5, 0.88, as.character(formatC(sum(daily_sub$npp),digits=1,format="f")) , adj=c(1,0))
text( 0, 0.83, "CEX ", adj=c(0,0));         text( 0.5, 0.83, as.character(formatC(sum(daily_sub$cex),digits=1,format="f")) , adj=c(1,0)); text( 0.55, 0.83, paste( as.character(formatC(sum(daily_sub$cex)/sum(daily_sub$npp)*100,digits=1,format="f")),"% of NPP",sep=""), adj=c(0,0))
text( 0, 0.78, "C -> veg ", adj=c(0,0));    text( 0.5, 0.78, as.character(formatC(annual_sub$calloc,digits=1,format="f")), adj=c(1,0))
text( 0, 0.73, "C veg -> lit", adj=c(0,0)); text( 0.5, 0.73, as.character(formatC(annual_sub$cveg2lit,digits=1,format="f")), adj=c(1,0))
text( 0, 0.68, "C lit -> soil", adj=c(0,0));text( 0.5, 0.68, as.character(formatC(annual_sub$clit2soil,digits=1,format="f")), adj=c(1,0)); text( 0.55, 0.68, paste(as.character(formatC(annual_sub$eff*100,digits=1,format="f")),"% efficiency",sep=""), adj=c(0,0))

dev.off()


##--------------------------------------
## N FLUXES
##--------------------------------------
pdf(paste("N_flux_overview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
par(mar=c(4,4,2,1), las=1 )
plot( daily_sub$doy, daily_sub$nup,    type="l", xlab="day of year", ylab="gC/m2/d", ylim=c(0,max(daily_sub$nup)), col="green" )
lines( daily_sub$doy, daily_sub$netmin, col="blue" )
lines( daily_sub$doy, daily_sub$nfixfree, col="red" )
legend( "bottomleft", c("Ninorg","Nup","net N min.","free-living BNF"), col=c("black","green","blue","red"), bty="n", lty=1 )

par(mar=c(0,0,2,0))
plot( c(0,1), c(0,1), type="n", axes=FALSE )
lab <- expression(paste("(gN m"^-2," yr"^-1,")", sep=""))
text( 0.10, 0.98, "ANNUAL TOTALS" , adj=c(0,0),font=2); text( 0.60, 0.98, lab , adj=c(0,0),font=2, cex=0.7)
text( 0.10, 0.93, "N uptake", adj=c(0,0));              text( 0.60, 0.93, as.character(formatC(sum(daily_sub$nup),digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.88, "N up. pas.", adj=c(0,0));            text( 0.60, 0.88, as.character(formatC(sum(daily_sub$nup_pas),digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.83, "N up. act.", adj=c(0,0));            text( 0.60, 0.83, as.character(formatC(sum(daily_sub$nup_act),digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.78, "N up. fix.", adj=c(0,0));            text( 0.60, 0.78, as.character(formatC(sum(daily_sub$nup_fix),digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.73, "N up. ret.", adj=c(0,0));            text( 0.60, 0.73, as.character(formatC(sum(daily_sub$nup_ret),digits=2,format="f")) , adj=c(1,0))


text( 0.10, 0.68, "N -> veg", adj=c(0,0));              text( 0.60, 0.68, as.character(formatC(annual_sub$nalloc,digits=2,format="f")) , adj=c(1,0));     text( 0.65, 0.68, paste(as.character(formatC(annual_sub$calloc/annual_sub$nalloc,digits=2,format="f")),"C:N") , adj=c(0,0))
text( 0.10, 0.63, "N veg -> lit", adj=c(0,0));          text( 0.60, 0.63, as.character(formatC(annual_sub$nveg2lit,digits=2,format="f")) , adj=c(1,0));   text( 0.65, 0.63, paste(as.character(formatC(annual_sub$cveg2lit/annual_sub$nveg2lit,digits=2,format="f")),"C:N") , adj=c(0,0))
text( 0.10, 0.58, "N lit -> soil", adj=c(0,0));         text( 0.60, 0.58, as.character(formatC(annual_sub$nlit2soil,digits=2,format="f")) , adj=c(1,0));  text( 0.65, 0.58, paste(as.character(formatC(annual_sub$clit2soil/annual_sub$nlit2soil,digits=2,format="f")),"C:N") , adj=c(0,0))
text( 0.10, 0.53, "N -> soil req.", adj=c(0,0));        text( 0.60, 0.53, as.character(formatC(annual_sub$nreq,digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.48, "N fBNF", adj=c(0,0));                text( 0.60, 0.48, as.character(formatC(sum(daily_sub$nfixfree),digits=2,format="f")) , adj=c(1,0))
text( 0.10, 0.43, "N net min.", adj=c(0,0));            text( 0.60, 0.43, as.character(formatC(sum(daily_sub$netmin),digits=2,format="f")) , adj=c(1,0))
dev.off()


##--------------------------------------
## N COST ANALYSIS
##--------------------------------------
pdf(paste("CostOverview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
par(mar=c(4,4,2,1), las=1 )
plot( daily_sub$doy, daily_sub$ccost, xlab="", ylab="", type="l", col="magenta" )
par(new=TRUE)
plot( daily_sub$doy, daily_sub$nup/daily_sub$ninorg, axes=FALSE, type="l", col="green" )
axis( 4, col="green" )
legend( "bottomleft", c("C price of N", "fraction N taken up"), col=c("magenta","green"), bty="n", lty=1 )

par(mar=c(0,0,2,0))
dev.off()


##--------------------------------------
## CUMULATIVE N FLUXES
##--------------------------------------
pdf(paste("N_flux_overview_cumulative_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
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

# ##--------------------------------------
# ## GPP model - obs.
# ##--------------------------------------
# pdf(paste("gpp_model_obs_seasonal.pdf"),width=sum(widths),height=sum(heights))
# panel <- layout(
#                 matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
#                 widths=widths,
#                 heights=heights,
#                 TRUE
#                 )
# # layout.show(panel)

# par(mar=c(4,4,2,1), las=1)
# plot( daily_sub$doy, daily_sub$gpp, type="l", xlab="", ylab="gC/m2/d", ylim=c(0,max(daily_sub$gpp)) )
# legend( "bottomleft", c("simulated GPP"), col=c("black"), bty="n", lty=1 )

# dev.off()

##--------------------------------------
## C POOLS
##--------------------------------------
pdf(paste("C_pool_overview_sofun_seasonal.pdf"),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
# layout.show(panel)

par(mar=c(4,4,2,1), las=1)
plot( daily_sub$doy, daily_sub$cleaf, type="l", xlab="", ylab="gC/m2/d", ylim=c(0,max(daily_sub$cleaf)) )
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


