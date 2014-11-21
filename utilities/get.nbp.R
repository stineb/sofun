get.nbp <- function( name, dir,
                    lu=FALSE,
                    netcdf=FALSE,
                    monthly=FALSE,
                    tstart=NA, tend=NA,
                    mask.id=NA
                    ) {
  
  if (netcdf){
    
    ## ////////////////////////////////////////////
    ## use NetCDF outputs
    ## --------------------------------------------
    library(RNetCDF)
    
    ## read NetCDF files
    print("reading LPX NetCDF output ...")
    if (monthly) {
      filn <- paste(dir,name,"_m.cdf",sep="" )
    } else {
      filn <- paste(dir,name,".cdf",sep="" )
    }
    print(paste("opening file:",filn))
    nc <- open.nc( filn )
    print("reading lon, lat, and time ...")
    lon <- var.get.nc(nc,"LONGITUDE")
    lat <- var.get.nc(nc,"LATITUDE")
    time <- var.get.nc(nc,"TIME")
    print("reading example variable ...")
    nep <- var.get.nc(nc,"nep")

    ## ## determine whether no-landuse simulation has k dimension
    ## noludim <- FALSE
    ## if (length(dim(nep))==3){
    ##   print("no LU dimension in file without LUC")
    ##   noludim <- TRUE
    ##   count[3] <- 1
    ## } else {
    ##   print("full LU dimension in file without LUC")
    ##   llu <- dim(nep.nolu)[3]
    ## }
    
    ## Determine domain to be read (time)
    if (!is.na(tstart)&&!is.na(tend)){
      print(paste("cutting to years ",tstart,tend))
      start <- dim(nep)
      start[] <- 1
      start[length(start)] <- which.min(abs(time-tstart))
      count <- dim(nep)
      count[length(count)] <- which.min(abs(time-tend))-which.min(abs(time-tstart))+1
    } else {
      print(paste("reading all time steps "))
      start <- dim(nep)
      start[] <- 1
      count <- dim(nep)
    }
    rm(nep)
    gc() # garbage collection

    ## Cut time vector to desired length
    if (!is.na(tstart)&&!is.na(tend)){
      time <- time[which.min(abs(time-tstart)):which.min(abs(time-tend))]
    }
    
    ## test
    if (count[4]!=length(time)) {
      print("length problem")
      stop
    }
    
    ## Read NetCDF
    print("reading NEP ...")
    nep <- var.get.nc(nc,"nep"
                         ,start=start,count=count
                         )
    print("reading time ...")
    time <- var.get.nc(nc,"TIME")
    print("reading LU_AREA ...")
    luarea <- var.get.nc(nc,"lu_area"
                         ,start=start,count=count
                         )
    if (lu) {
      ## Read only if output is from a simulation with landuse (product pool decay)
      print(paste("simulation with landuse:",name))
      print("reading PRODUCT C FLUX ...")
      cflux.prod <- var.get.nc(nc,"acflux_products"
                               ,start=start,count=count
                               )
    }
    close.nc(nc)
    gc() # garbage collection

    ## get gridcell area
    nc <- open.nc( paste(dir,name,".cdf",sep="") )
    area <- var.get.nc(nc,"area")
    close.nc(nc)
    gc() # garbage collection
    
    ## ## Read no-LU file
    ## print(paste("no-LU sim:",name.nolu))
    ## if (monthly) {
    ##   nc <- open.nc( paste(dir,name.nolu,"_m.cdf",sep="") )
    ## } else {
    ##   nc <- open.nc( paste(dir,name.nolu,".cdf",sep="") )
    ## }
    ## print("reading without-LU-variable (example) ...")
    ## nep.nolu <- var.get.nc(nc,"nep")

    ## Determine domain to be read (time)
    if (!is.na(tstart)&&!is.na(tend)){
      print(paste("cutting to years ",tstart,tend))
    } else {
      print(paste("reading all time steps "))
    }

    ## Cut time vector to desired length
    if (!is.na(tstart)&&!is.na(tend)){
      time <- time[which.min(abs(time-tstart)):which.min(abs(time-tend))]
    }
      
    ## Apply spatial mask (continents)
    if (!is.na(mask.id)) {
      print("reading continents file ...")
      nconts <- 8
      nc <- open.nc( "/card/forcings/lpx/regionmasks/regmask_1x1deg.nc" )
      ## MASK[k=1]:MASK[k=8] are in this order for:
      ## latin america, africa, south/southeast asia, china, europe,
      ## australia and oceania, north america, russia
      mask <- var.get.nc(nc,"mask")
    }

    ## determine dimension lenghts
    llon  <- dim(luarea)[1]
    llat  <- dim(luarea)[2]
    llu   <- dim(luarea)[3]
    ltime <- dim(luarea)[4]
            
    ## element-wise multiplication to weigh flux densities with
    ## land use area fractions [-> gC/m2/yr per gridcell and land use class]
    print("multiplying by luarea ...")
    lunep <- nep*luarea
    if (lu) {
      lucflux.prod <- cflux.prod
    }

    ## remove flux from products from NEP
    if (lu) {
      lunep <- lunep - lucflux.prod
    }
    rm(cflux.prod,lucflux.prod)
    gc()
    
    
    ## Reduce dimension: sum over land use class
    ## 'apply' function did not work for some reason [-> gC/m2/yr per gridcell]
    print("sum over luarea dimension ...")

    lunep.lusum        <- array(NA,dim=c(dim(lunep)[1:2],dim(lunep)[4]))

    for (ilon in seq(llon)){
      for (ilat in seq(llat)){
        if (!is.na(luarea[ilon,ilat,1,1])){
          for (itime in seq(ltime)){
            lunep.lusum[ilon,ilat,itime] <- sum(lunep[ilon,ilat,,itime], na.rm=TRUE)
          }
        }
      }
    }
    rm(lunep)
    gc()
    
    ## multiply with grid cell area to get absolute fluxes
    ## [-> gC/yr per gridcell]
    print("multiply field with grid cell area to get absolute values ...")

    lunep.abs <- lunep.lusum*NA

    for (itime in seq(ltime)){
      lunep.abs[,,itime] <- lunep.lusum[,,itime]*area
    }
    
    ## sum over all gridcells for each time step individually
    ## [-> PgC/yr global total]
    print("sum over all gridcells to get global total ...")

    nbp <- rep(NA,ltime)
    
    for (itime in seq(ltime)){
      nbp[itime] <- sum( lunep.abs[,,itime], na.rm=TRUE )/1e15
    }
    
    ## return land use flux time series 'f.luc.tseries', and spatially resolved field
    ## 'f.luc.field.out' attached to the list 'out.f.luc'.
    ## this tseries gives identical results as tseries when computed with ascii output
    ## (within numerical imprecision of ~1e-8 GtC/yr)
    print("attach output to list 'out.nbp' ...")
    out.nbp <- list()
    nbp.tseries <- data.frame( time=time, nbp=nbp, cum.nbp=cumsum(nbp) )
    out.nbp$tseries <- nbp.tseries
    out.nbp$field <- lunep.lusum
    out.nbp$lon <- lon
    out.nbp$lat <- lat
    
    return(out.nbp)
        
  } else {
    ## ////////////////////////////////////////////
    ## use ASCII outputs
    ## --------------------------------------------
    print("reading LPX ASCII output ...")
    
    ## --------------------------------------------
    ## compute by flux (default)
    ## --------------------------------------------
    
    if (lu) {
      tmp      <- read.table( paste( dir, "trans_", name,   ".nep.out", sep="") )
      nep   <- tmp[,2]
      time  <- tmp[,1]
      tmp      <- read.table( paste( dir, "trans_", name, ".cflux_prod.out", sep="") )
      cflux.prod <- tmp[,2]
    } else {
      tmp      <- read.table( paste( dir, "trans_", name, ".nep.out", sep="") )
      print(paste( dir, "trans_", name, ".nep.out", sep="" ))
      time<- tmp[,1]
      nep <- tmp[,2]
    }
    
    if ( is.na(tstart) && is.na(tend) ){
      if (length(time)!=length(time)) {
        print("selecting subset of time series available in both datasets")
        tstart <- head(time,n=1)
        tend <- tail(time,n=1)
        tstart <- head(time,n=1)
        tend <- tail(time,n=1)
        if (tstart < tstart) {
          tstart <- tstart
        } else {
          tstart <- tstart
        }
        if (tend < tend){
          tend <- tend
        } else {
          tend <- tend
        }
        nep  <- nep[which(time==tstart):which(time==tend)]
        nep <- nep[which(time==tstart):which(time==tend)]
        cflux.prod <- cflux.prod[which(time==tstart):which(time==tend)]
        time <- time[which(time==tstart):which(time==tend)]
        time <- time[which(time==tstart):which(time==tend)]
        time <- time[which(time==tstart):which(time==tend)]
      }
    } else {
      nep  <- nep[which(time==tstart):which(time==tend)]
      nep    <- nep[which(time==tstart):which(time==tend)]
      cflux.prod <- cflux.prod[which(time==tstart):which(time==tend)]
      time <- time[which(time==tstart):which(time==tend)]
      time   <- time[which(time==tstart):which(time==tend)]
    }

    out.nbp <- list()
    if (lu) {
      nbp.tseries <- data.frame( time=time, nbp=nep-cflux.prod, cum.nbp=cumsum(nep-cflux.prod) )
    } else {
      nbp.tseries <- data.frame( time=time, nbp=nep, cum.nbp=cumsum(nep) )
    }
    out.nbp$tseries <- nbp.tseries
    return(out.nbp)
    
  }
  
}
