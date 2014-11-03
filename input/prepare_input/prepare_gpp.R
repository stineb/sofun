## ////////////////////////////////////////////////////////////////////
## Prepare GPP input from Tyler's fluxnet data for each day available 
## available in the original free-and-fair-use FLUXNET dataset.
## File created here is formatted like this:
## '2002  1  1 0.496632 0.054053', which represents 
## 'YYYY MM DM      GPP GPP err.'.
## --------------------------------------------------------------------
data.in <- read.csv('/alphadata01/bstocker/data/gepisat/gpp_daily/AU-How_daily_gpp.txt',header=TRUE)
data.out <- array(NA,dim=c(dim(data.in)[1],5)) # five columns: year, month, day in month, gpp, gpp-error

len <- dim(data.in)[1]

## read in array from original file
for (idx in seq(len)){
  data.out[idx,1] <- as.numeric(substr(data.in$Timestamp[idx],start=1,stop=4))
  data.out[idx,2] <- as.numeric(substr(data.in$Timestamp[idx],start=6,stop=7))
  data.out[idx,3] <- as.numeric(substr(data.in$Timestamp[idx],start=9,stop=10))
  data.out[idx,4] <- data.in[idx,2]
  data.out[idx,5] <- data.in[idx,3]
}

## write output with standard formatting: (year, mo, dm, gpp, gpp-err)
formatted.out <- vector("character",len)
for (i in 1:len){
  formatted.out[i] <- sprintf("%4i %2i %2i %6f %6f",data.out[i,1],data.out[i,2],data.out[i,3],data.out[i,4],data.out[i,5])
  print(formatted.out[i])
}
writeLines(formatted.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_daily_gpp_STANDARD.txt")

## plot GPP for each day in the orignial time series
plot( 1:len, data.in$GPP_mol_m2, type="l", ylim=c(0,1) )


## ////////////////////////////////////////////////////////////////////
## Prepare GPP input from Tyler's fluxnet data containing the median
## seasonality GPP over all available years in the multi-year the 
## original free-and-fair-use FLUXNET dataset.
## File created here is formatted like this:
## ' 1  1 0.496632 0.054053', which represents 
## 'MM DM      GPP GPP err.'.
## Below code uses 'data.out' defined above!
## --------------------------------------------------------------------

## average over each day of the year for all years available in the 
## original data set.
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
nmonth <- 12
data.med <- array( NA, dim=c(sum(ndaymonth),4 ))
day <- 0
for (mo in 1:nmonth){
  for (dm in 1:ndaymonth[mo]){
    day <- day+1
    ## find values in entire time series that corresponds to this day
    tmp <- NA
    tmp.err <- NA
    for (i in 1:len){
      if (data.out[i,2]==mo && data.out[i,3]==dm){
        tmp <- c(tmp,data.out[i,4])      
        tmp.err <- c(tmp.err,data.out[i,5])      
      }
    }
    data.med[day,1] <- mo
    data.med[day,2] <- dm
    data.med[day,3] <- median(tmp,na.rm=TRUE)
    data.med[day,4] <- median(tmp.err,na.rm=TRUE)
  }
}

## write output with standard formatting: (mo, dm, gpp, gpp-err)
formatted.out <- vector("character",sum(ndaymonth))
for (i in 1:sum(ndaymonth)){
  formatted.out[i] <- sprintf("%2i %2i %6f %6f",data.med[i,1],data.med[i,2],data.med[i,3],data.med[i,4])
  print(formatted.out[i])
}
writeLines(formatted.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_daily_gpp_med_STANDARD.txt")

## plot median-seasonality GPP
plot( 1:sum(ndaymonth), data.med[,3], type="l" )

