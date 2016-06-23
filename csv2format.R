## ////////////////////////////////////////////////////////////////////
## Prepare climate input from Tyler's fluxnet data containing 
## fsun, tair, pre
## File created here is formatted like this:
## ' 1  1 0.496632 0.054053', which represents 
## 'MM DM      GPP GPP err.'.
## Below code uses 'dtemp.out' defined above!

##////////////////////////////////////////////////////////////////
##  SR reads one (annual) value corresponding to the given year 
##  from a time series ascii file. File has to be located in 
##  ./input/ and has to contain only rows formatted like
##  '2002  1  1 0.496632 0.054053', which represents 
##  'YYYY MM DM      GPP GPP err.'. DM is the day within the month.
##  If 'realyear' is dummy (-9999), then it's interpreted as to 
##  represent a mean climatology for the course of a year.
##----------------------------------------------------------------


## --------------------------------------------------------------------
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- 12

data.in <- read.csv('/alphadata01/bstocker/sofun/trunk/input/data_san-fran_2000-day.csv',header=TRUE)

# dtemp.out <- array( NA, dim=c(sum(ndaymonth), 5 ))
# dprec.out <- array( NA, dim=c(sum(ndaymonth), 5 ))
# dfsun.out <- array( NA, dim=c(sum(ndaymonth), 5 ))

# ## read in array from original file
# dm <- 0
# moy <- 1
# for (idx in seq(ndayyear)){
  
#   ## define 'dm' (day of month) and 'moy' (month of year)
#   dm <- dm+1
#   if (dm>ndaymonth[moy]) {
#     dm <- 1
#     moy <- moy+1 
#   }

#   dtemp.out[idx,1] <- -9999
#   dtemp.out[idx,2] <- moy
#   dtemp.out[idx,3] <- dm
#   dtemp.out[idx,4] <- data.in[idx,2]
#   dtemp.out[idx,5] <- -9999

#   dprec.out[idx,1] <- -9999
#   dprec.out[idx,2] <- moy
#   dprec.out[idx,3] <- dm
#   dprec.out[idx,4] <- data.in[idx,3]
#   dprec.out[idx,5] <- -9999

#   dfsun.out[idx,1] <- -9999
#   dfsun.out[idx,2] <- moy
#   dfsun.out[idx,3] <- dm
#   dfsun.out[idx,4] <- data.in[idx,1]
#   dfsun.out[idx,5] <- -9999

# }

## write output with standard formatting: (year, mo, dm, gpp, gpp-err)
formatted.dtemp.out <- vector("character",ndayyear)
formatted.dprec.out <- vector("character",ndayyear)
formatted.dfsun.out <- vector("character",ndayyear)

for (i in 1:ndayyear){
  # formatted.dtemp.out[i] <- sprintf("%4i %2i %2i %6f %6f",dtemp.out[i,1],dtemp.out[i,2],dtemp.out[i,3],dtemp.out[i,4],dtemp.out[i,5])
  # formatted.dprec.out[i] <- sprintf("%4i %2i %2i %6f %6f",dprec.out[i,1],dprec.out[i,2],dprec.out[i,3],dprec.out[i,4],dprec.out[i,5])
  # formatted.dfsun.out[i] <- sprintf("%4i %2i %2i %6f %6f",dfsun.out[i,1],dfsun.out[i,2],dfsun.out[i,3],dfsun.out[i,4],dfsun.out[i,5])
  # print(formatted.dtemp.out[i])
  formatted.dtemp.out[i] <- sprintf("%6f",data.in$tair[i])
  formatted.dprec.out[i] <- sprintf("%6f",data.in$pre[i])
  formatted.dfsun.out[i] <- sprintf("%6f",data.in$fsun[i])
  print(formatted.dtemp.out[i])
}

writeLines(formatted.dtemp.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_dtemp_2000_STANDARD.txt")
writeLines(formatted.dprec.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_dprec_2000_STANDARD.txt")
writeLines(formatted.dfsun.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_dfsun_2000_STANDARD.txt")
