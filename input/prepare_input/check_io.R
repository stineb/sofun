do.gpp <- FALSE
do.stash <- TRUE

if (do.gpp){
  indata  <- read.table("/alphadata01/bstocker/sofun/trunk/input/AU-How_daily_gpp_med_STANDARD.txt", col.names=c("mo","dm","gpp","gpp-err"))
  outdata <- read.table("/alphadata01/bstocker/sofun/trunk/output/RUNNAME.d.gpp.dat", col.names=c("itime","gpp"))

  ## input data
  plot( 1:dim(indata)[1], indata$gpp, type="l" )

  ## output data
  lines( 1:dim(indata)[1], outdata$gpp[1:dim(indata)[1]], col="red" )  
}

if (do.stash){
  varnam <- "wn"

  sofun_out <- read.table(paste("/alphadata01/bstocker/sofun/trunk/output/RUNNAME.d.",varnam,".out", sep=""), col.names=c("itime",varnam))
  sofun_out <- sofun_out[as.integer(sofun_out$itime)==2000,]
  stash_out <- read.table(paste("/alphadata01/bstocker/stash/f90_version/output/",varnam,".d.out", sep=""), col.names=c(varnam))

  ylim <- c(min(stash_out[[ varnam ]], sofun_out[[ varnam ]]),max( stash_out[[ varnam ]], sofun_out[[ varnam ]] )) 

  ## stash data
  plot( 1:dim(stash_out)[1], stash_out[[ varnam ]], type="l", ylim=ylim )

  ## sofun data
  lines( 1:dim(sofun_out)[1], sofun_out[[ varnam ]], col="red" )

}