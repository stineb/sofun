library(zoo)

df <- read.table( "/alphadata01/bstocker/data/co2/cCO2_rcp85_const850-1765.dat", col.names=c("year","co2") )

len <- length(df[,1])

## constant at year 1993 value (although simulation starts in 1982)
# idx <- which( df$year == 1993 )
# amb <- df$co2[idx]
amb <- 350.0
elv <- 600.0

print(paste("Constant at", amb))
out.ctrl <- df
out.ctrl$co2[1:len] <- amb

print(paste("Step change to", elv))
out.step <- df
out.step$co2[1:idx]   <- amb
out.step$co2[(idx+1):len] <- elv

## gradual increase from 1993 to 2073
idx <- which( df$year == 1993 )
jdx <- which( df$year == 2073 )
lev0 <- amb
lev1 <- elv
print( paste( "Linear increase from", lev0, "to", lev1 ) )
out.grad <- out.step
out.grad$co2[]        <- NA
out.grad$co2[1:idx]   <- lev0
out.grad$co2[jdx:len] <- lev1
out.grad$co2          <- na.approx( out.grad$co2 )

## test plot
plot( out.grad$year,  out.grad$co2,type="l", xlim=c(1982,2100) )
lines( out.grad$year, out.step$co2,col="red" )

## write output with standard formatting: (year, co2)
formatted.ctrl.out <- vector("character",len)
formatted.step.out <- vector("character",len)
formatted.grad.out <- vector("character",len)

for (i in 1:len){
  formatted.ctrl.out[i] <- sprintf("%16.6f    %f", out.ctrl$year[i], out.ctrl$co2[i])
  formatted.step.out[i] <- sprintf("%16.6f    %f", out.step$year[i], out.step$co2[i])
  formatted.grad.out[i] <- sprintf("%16.6f    %f", out.grad$year[i], out.grad$co2[i])
  print(formatted.step.out[i])
}

writeLines(formatted.ctrl.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/co2/SwissFACE/cCO2_ctrl.dat")
writeLines(formatted.step.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/co2/SwissFACE/cCO2_step.dat")
writeLines(formatted.grad.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/co2/SwissFACE/cCO2_grad.dat")
