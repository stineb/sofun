daily2monthly <- function( dvals ){
  ## /////////////////////////////////////
  ## get monthly values
  ## dvals is vector of daily values
  ## -------------------------------------

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  ndayyear <- sum(ndaymonth)
  nmonth <- length(ndaymonth)

  ## can only treat even years, i.e. 365 days in a year
  if (length(dvals)%%ndayyear!=0) {stop}

  ## get number of years covered by dvals
  nyears <- length(dvals)/ndayyear
  moy <- rep(rep(seq(12),times=ndaymonth),nyears)
  yearno <- rep(seq(nyears),each=ndayyear)

  ## average over daily values of each month
  mvals <- rep(NA,nmonth*nyears)
  for (yr in 1:nyears){
    for (m in 1:nmonth){
      istart <- (yr-1)*ndayyear+1
      iend   <- yr*ndayyear
      # print(paste(istart,iend))
      tmp <- dvals[istart:iend]
      # print(paste('writing to index',(yr-1)*nmonth+m))
      # print(tmp[which(moy[istart:iend]==m)])
      mvals[(yr-1)*nmonth+m] <- mean(tmp[which(moy[istart:iend]==m)])
    }
  }
  return(mvals)
}

monthly2daily <- function( mvals ){

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  ndayyear <- sum(ndaymonth)
  nmonth <- length(ndaymonth)

  ## can only treat even years, i.e. 365 days in a year
  if (length(mvals)%%nmonth!=0) {stop}

  

}


## get soil moisture and temperature
sofun.out <- read.table( '/alphadata01/bstocker/sofun/trunk/output/RUNNAME.d.wn.out', col.names=c("year","dwtot") )
tmp <- read.table( '/alphadata01/bstocker/sofun/trunk/output/RUNNAME.d.ea_n.out')
sofun.out$daet <- tmp[,2]
tmp <- read.table( '/alphadata01/bstocker/sofun/trunk/output/RUNNAME.d.eq_n.out')
sofun.out$deet <- tmp[,2]
sofun.out$dcpa <- sofun.out$daet/sofun.out$deet # cramer-prentice alpha
df.temp <- read.table( '/alphadata01/bstocker/sofun/trunk/input/dtemp_CH-Oe1_2000_365.txt', col.names="dtemp" )

nyears <- dim(sofun.out)[1]/ndayyear
mtemp <- daily2monthly( rep(df.temp$dtemp,nyears) )
mwtot <- daily2monthly( sofun.out$dcpa )

mtemp_soil <- rep(NA,nmonth*nyears)
dtemp_soil <- rep(NA,ndayyear*nyears)

# mtemp <- rep(18*sin(pi/6*(1:12+9))+1, nyears)

## for soil parameters, use LPJ-soil code 2: medium, silty clay loam 
diffus_wp    <- 0.2        # soilpar(5,jpngr) = store(3,soilcode(jpngr))
diffus_15whc <- 0.65       # soilpar(6,jpngr) = store(4,soilcode(jpngr))
diffus_fc    <- 0.4        # soilpar(7,jpngr) = store(5,soilcode(jpngr))

## Average over preceeding 365 days
avetemp <- mean(df.temp$dtemp)
meanw1  <- mean(sofun.out$dcpa)
    
## Interpolate thermal diffusivity function against soil water content
if (meanw1<0.15) {
  diffus <- (diffus_15whc-diffus_wp)/0.150*meanw1+diffus_wp
} else {
  diffus <- (diffus_fc-diffus_15whc)/0.850*(meanw1-0.150)+diffus_15whc
}

## Convert diffusivity from mm2/s to m2/month
## multiplication by 1e-6 (-> m2/s) * 2.628e6 (s/month)  <-  2.628
mdiffus <- diffus*2.6280
ddiffus <- diffus*0.0864
    
## Calculate mamplitude fraction and lag at soil depth 0.25 m
malag <- 0.250/sqrt(12.00*mdiffus/pi)
mamp <- exp(-malag)
mlag <- malag*(6.00/pi)                                 ##convert lag from angular units to months

# dalag <- 0.250/sqrt(12.00*ddiffus/pi)
# damp <- exp(-dalag)
# dlag <- dalag*(6.00/pi)                                 ##convert lag from angular units to months

dalag <- 0.250/sqrt(365*ddiffus/pi)
damp <- exp(-dalag)
dlag <- dalag*(365/(2*pi))                                 ##convert lag from angular units to days
    
## Calculate monthly soil temperatures for this year.  For each month,
## calculate average air temp for preceding 12 months (including this one)
    
## Estimate air temperature "lag" months ago by linear interpolation
## between air temperatures for this and last month
for (yr in 1: nyears) {
  for (m in 1:nmonth) {
    if (m==1) {
      pm <- nmonth
    } else {
      pm <- m-1    
    }
    tempthismonth <- mtemp[m]
    templastmonth <- mtemp[pm]
    lagtemp <- (tempthismonth-templastmonth)*(1.0-mlag)+templastmonth
        
    ## Adjust mamplitude of lagged air temp to give estimated soil temp
    mtemp_soil[(yr-1)*nmonth+m] <- avetemp+mamp*(lagtemp-avetemp)
  }
} 

for (yr in 1:nyears){
  day<-0
  for (m in 1:nmonth){
    for (dm in 1:ndaymonth[m]){
      day<-day+1
      dtemp_soil[(yr-1)*ndayyear+day] <- mtemp_soil[(yr-1)*nmonth+m]
    }
  }
}

# for (yr in 1:nyears){
#   for (day in 1:ndayyear){
#     if (day==1){
#       pd <- ndayyear
#     } else {
#       pd <- day - 1
#     }
#     tempthisday <- df.temp$dtemp[day]
#     templastday <- df.temp$dtemp[pd]
#     lagtemp <- (tempthisday-templastday)*(1.0-dlag)+templastday
        
#     ## Adjust mamplitude of lagged air temp to give estimated soil temp
#     dtemp_soil[(yr-1)*ndayyear+day] <- avetemp+damp*(lagtemp-avetemp)
#   }
# }


## plotting
# plot( 1:(3*nmonth), mtemp[1:(3*nmonth)], type="l" )
# lines( 1:(3*nmonth), mtemp_soil[1:(3*nmonth)], col="red" )
# lines( 1:(3*nmonth), daily2monthly(dtemp_soil[1:(3*ndayyear)]), col="blue" )

plot( 1:(ndayyear), df.temp$dtemp, type="l")
lines( 1:(ndayyear), dtemp_soil[1:ndayyear], type="l", col="blue")
