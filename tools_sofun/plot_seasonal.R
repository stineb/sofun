runname <- "RUNNAME"
outdir <- "/alphadata01/bstocker/sofun/trunk/output/"
vars <- c("cex","cleaf","gpp","netmin","ninorg","npp","nup","nfixfree")

plotyear <- 2000

nmonths <- 12
ndayyear <- 365
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
nvars <- length(vars)

## get time data
tmp <- read.table( paste(outdir,runname,".d.",vars[1],".dat",sep=""))
time <- tmp$V1
nyears <- length(time)/ndayyear

## construct data frame
df <- data.frame( time = time )
for (ivar in seq(nvars)) {
  df[[ vars[ivar] ]] <- read.table( paste(outdir,runname,".d.",vars[ivar],".dat",sep=""))[,2]
}
df[[ "year" ]] <- as.integer(time)
df[[ "moy" ]]  <- rep(rep(seq(12),times=ndaymonth),nyears)
df[[ "doy" ]]  <- rep(seq(ndayyear),nyears)

## take subset of one year and normalise all variables to [0,1]
df_sub <- df[df$year==plotyear,]
for (ivar in seq(nvars)) {
  df_sub[[ paste(vars[ivar],".norm", sep="") ]] <- df_sub[[ vars[ivar] ]] / max( df_sub[[ vars[ivar] ]], na.rm=TRUE )
}

## Plot absolute C variables
aspect <- 0.5
magn <- 4
ncols <- 2
nrows <- 2
widths <- rep(magn,ncols)
widths[2] <- 0.4*magn
heights <- rep(aspect*magn,nrows)
# heights[nrows] <- 0.3*magn

# pdf(paste("fig/plf_maps_global_",res,".pdf",sep=""),width=sum(widths),height=sum(heights))
panel <- layout(
                matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                widths=widths,
                heights=heights,
                TRUE
                )
# layout.show(panel)

par(mar=c(4,4,0,1))
plot( df_sub$doy, df_sub$gpp, type="l", xlab="", ylab="gC/m2/d", ylim=c(0,max(df_sub$gpp)) )
lines( df_sub$doy, df_sub$npp, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
lines( df_sub$doy, df_sub$cex, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
legend( "bottomleft", c("GPP","NPP","EXU"), col=c("black","green","blue"), bty="n", lty=1 )

par(mar=c(0,0,0,0))
plot( c(0,1), c(0,1), type="n", axes=FALSE )
text( 0, 0.98, "ANNUAL TOTALS", adj=c(0,0),font=2)
text( 0, 0.93, paste("NPP",as.character(format(sum(df_sub$npp),digits=2))) , adj=c(0,0))
text( 0, 0.88, paste("CEX ",as.character(format(sum(df_sub$cex),digits=2))," (",as.character(format(sum(df_sub$cex)/sum(df_sub$npp)*100,digits=2)),"% of NPP)",sep=""), adj=c(0,0))


par(mar=c(4,4,0,1))
plot( df_sub$doy, df_sub$ninorg, type="l", xlab="day of year", ylab="gN/m2/d", ylim=c(0,max(df_sub$ninorg)) )
lines( df_sub$doy, df_sub$nup, type="l", xlab="day of year", ylab="gC/m2/d", col="green" )
lines( df_sub$doy, df_sub$netmin, type="l", xlab="day of year", ylab="gC/m2/d", col="blue" )
lines( df_sub$doy, df_sub$nfixfree, type="l", xlab="day of year", ylab="gC/m2/d", col="red" )
legend( "bottomleft", c("Ninorg","Nup","net N min.","free-living BNF"), col=c("black","green","blue","red"), bty="n", lty=1 )
