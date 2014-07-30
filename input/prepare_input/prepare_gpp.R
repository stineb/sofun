## prepare GPP input from Tyler's fluxnet data
data.in <- read.csv('/alphadata01/bstocker/data/gepisat/gpp_daily/AU-How_daily_gpp.txt',header=TRUE)
data.out <- array(NA,dim=c(dim(data.in)[1],5)) # five columns: year, month, day in month, gpp, gpp-error

len <- dim(data.in)[1]

for (idx in seq(len)){
  data.out[idx,1] <- as.numeric(substr(data.in$Timestamp[idx],start=1,stop=4))
  data.out[idx,2] <- as.numeric(substr(data.in$Timestamp[idx],start=6,stop=7))
  data.out[idx,3] <- as.numeric(substr(data.in$Timestamp[idx],start=9,stop=10))
  data.out[idx,4] <- data.in[idx,2]
  data.out[idx,5] <- data.in[idx,3]
}

formatted.out <- vector("character",len)
for (i in 1:len){
#  formatted.out[i] <- sprintf("%4i %2i %2i %6f %6f",data.out[i,1],data.out[i,2],data.out[i,3],data.out[i,4],data.out[i,5])
  formatted.out[i] <- sprintf("%4i %2i %2i %6f %6f",data.out[i,1],data.out[i,2],data.out[i,3],data.out[i,4],data.out[i,5])
  print(formatted.out[i])
}

writeLines(formatted.out,"/alphadata01/bstocker/sofun/trunk/input/AU-How_daily_gpp_STANDARD.txt")

# write.table(
#   data.out,
#   '/alphadata01/bstocker/data/gepisat/gpp_daily/AU-How_daily_gpp_STANDARD.txt',
#   col.names=c('yr','mo','dm','gpp_mol_m2','gpp_err'),
#   row.names=FALSE
#   )