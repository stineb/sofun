library(zoo)

years <- 850:2500
df <- data.frame( year=years, ndep=rep(1.0,length(years)))

len <- length(df[,1])

## constant at year 2002 value (although simulation starts in 1982)
# idx <- which( df$year == 2002 )
# amb <- df$ndep[idx]

## levels of NHx and NOy from SwissFACE experiment
amb <- 7.0
elv <- 28.0

print(paste("Constant at", amb))
## nhx
out_nhx.ctrl <- df
out_nhx.ctrl$ndep[1:len] <- amb
## noy
out_noy.ctrl <- df
out_noy.ctrl$ndep[1:len] <- amb

## step change after 20 years
idx <- which( df$year == 1993 )
print(paste("Step change to", elv))
## nhx
out_nhx.step <- df
out_nhx.step$ndep[1:idx]   <- amb
out_nhx.step$ndep[(idx+1):len] <- elv
## noy
out_noy.step <- df
out_noy.step$ndep[1:idx]   <- amb
out_noy.step$ndep[(idx+1):len] <- elv

## gradual increase from 1993 to 2073
idx <- which( df$year == 1993 )
jdx <- which( df$year == 2073 )
lev0 <- amb
lev1 <- elv
print( paste( "Linear increase from", lev0, "to", lev1 ) )
## nhx
out_nhx.grad <- out_nhx.step
out_nhx.grad$ndep[]        <- NA
out_nhx.grad$ndep[1:idx]   <- lev0
out_nhx.grad$ndep[jdx:len] <- lev1
out_nhx.grad$ndep          <- na.approx( out_nhx.grad$ndep )
## noy
out_noy.grad <- out_noy.step
out_noy.grad$ndep[]        <- NA
out_noy.grad$ndep[1:idx]   <- lev0
out_noy.grad$ndep[jdx:len] <- lev1
out_noy.grad$ndep          <- na.approx( out_noy.grad$ndep )


## test plot
plot( out_nhx.grad$year,  out_nhx.grad$ndep,type="l", xlim=c(1982,2100) )
lines( out_nhx.grad$year, out_nhx.step$ndep,col="red" )

plot( out_noy.grad$year,  out_noy.grad$ndep,type="l", xlim=c(1982,2100) )
lines( out_noy.grad$year, out_noy.step$ndep,col="red" )

## write output with standard formatting: (year, co2)
formatted_noy.ctrl.out <- vector("character",len)
formatted_noy.step.out <- vector("character",len)
formatted_noy.grad.out <- vector("character",len)

formatted_nhx.ctrl.out <- vector("character",len)
formatted_nhx.step.out <- vector("character",len)
formatted_nhx.grad.out <- vector("character",len)

for (i in 1:len){
  formatted_noy.ctrl.out[i] <- sprintf("%16.6f    %f", out_noy.ctrl$year[i], out_noy.ctrl$ndep[i])
  formatted_noy.step.out[i] <- sprintf("%16.6f    %f", out_noy.step$year[i], out_noy.step$ndep[i])
  formatted_noy.grad.out[i] <- sprintf("%16.6f    %f", out_noy.grad$year[i], out_noy.grad$ndep[i])

  formatted_nhx.ctrl.out[i] <- sprintf("%16.6f    %f", out_nhx.ctrl$year[i], out_nhx.ctrl$ndep[i])
  formatted_nhx.step.out[i] <- sprintf("%16.6f    %f", out_nhx.step$year[i], out_nhx.step$ndep[i])
  formatted_nhx.grad.out[i] <- sprintf("%16.6f    %f", out_nhx.grad$year[i], out_nhx.grad$ndep[i])
}

writeLines(formatted_nhx.ctrl.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_nhx_ctrl.dat")
writeLines(formatted_nhx.step.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_nhx_step.dat")
writeLines(formatted_nhx.grad.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_nhx_grad.dat")

writeLines(formatted_noy.ctrl.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_noy_ctrl.dat")
writeLines(formatted_noy.step.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_noy_step.dat")
writeLines(formatted_noy.grad.out,"/alphadata01/bstocker/sofun/trunk/input/sitedata/ndep/SwissFACE/ndep_noy_grad.dat")
