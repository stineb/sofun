get.f.luc <- function( name.nolu, name.lu, dir,
                      netcdf=FALSE,
                      monthly=FALSE,
                      tstart=NA, tend=NA,
                      grossfluxes=FALSE,
                      byflux=TRUE,
                      mask.id=NA,
                      offset=FALSE
                      ) {
  
  if (netcdf){
    
    ## ////////////////////////////////////////////
    ## use NetCDF outputs
    ## --------------------------------------------
    library(RNetCDF)
        
    ## read NetCDF files
    print("reading LPX NetCDF output ...")
    print(paste("LU sim:",name.lu))
    if (monthly) {
      nc <- open.nc( paste(dir,name.lu,"_m.cdf",sep="") )
    } else {
      nc <- open.nc( paste(dir,name.lu,".cdf",sep="") )
    }
    print("reading lon, lat, and time ...")
    lon <- var.get.nc(nc,"LONGITUDE")
    lat <- var.get.nc(nc,"LATITUDE")
    time <- var.get.nc(nc,"TIME")
    print("reading with-LU-variable (example) ...")
    nep.lu <- var.get.nc(nc,"nep")

    ## Determine domain to be read (time)
    if (!is.na(tstart)&&!is.na(tend)){
      print(paste("cutting to years ",tstart,tend))
      start.lu <- dim(nep.lu)
      start.lu[] <- 1
      start.lu[length(start.lu)] <- which.min(abs(time-tstart))
      count.lu <- dim(nep.lu)
      count.lu[length(count.lu)] <- which.min(abs(time-tend))-which.min(abs(time-tstart))+1
    } else {
      print(paste("reading all time steps "))
      start.lu <- dim(nep.lu)
      start.lu[] <- 1
      count.lu <- dim(nep.lu)
    }
    rm(nep.lu)
    gc() # garbage collection

    ## Cut time vector to desired length
    if (!is.na(tstart)&&!is.na(tend)){
      time <- time[which.min(abs(time-tstart)):which.min(abs(time-tend))]
    }
    
    ## test
    if (count.lu[4]!=length(time)) {
      print("length problem")
      stop
    }
    
    ## Read with-LU file
    print(paste("LU sim:",name.lu))
    print("reading NEP ...")
    nep.lu <- var.get.nc(nc,"nep"
                         ,start=start.lu,count=count.lu
                         )
    print("reading PRODUCT C FLUX ...")
    cflux.prod.lu <- var.get.nc(nc,"acflux_products"
                                ,start=start.lu,count=count.lu
                                )
    print("reading LU_AREA ...")
    luarea.lu <- var.get.nc(nc,"lu_area"
                            ,start=start.lu,count=count.lu
                            )
    if (offset) {
      count.offset <- count.lu
      count.offset[length(count.offset)] <- 1
      print("reading C pool variables ...")
      littera.lu <- var.get.nc(nc,"littercarbon_ag",start=start.lu,count=count.offset)
      litterb.lu <- var.get.nc(nc,"littercarbon_bg",start=start.lu,count=count.offset)
      vegc.lu <- var.get.nc(nc,"vegcarbon",start=start.lu,count=count.offset)
      soilc.lu <- var.get.nc(nc,"soilcarbon",start=start.lu,count=count.offset)
      exuc.lu <- var.get.nc(nc,"exucarbon",start=start.lu,count=count.offset)
      prodc.lu <- var.get.nc(nc,"products",start=start.lu,count=count.offset)
      totc.lu <- littera.lu+litterb.lu+vegc.lu+soilc.lu+exuc.lu
      rm(littera.lu,litterb.lu,vegc.lu,soilc.lu,exuc.lu)
      gc()
    }
    close.nc(nc)
    gc() # garbage collection
    
    ## Read no-LU file
    print(paste("no-LU sim:",name.nolu))
    if (monthly) {
      nc <- open.nc( paste(dir,name.nolu,"_m.cdf",sep="") )
    } else {
      nc <- open.nc( paste(dir,name.nolu,".cdf",sep="") )
    }
    print("reading time ...")
    time <- var.get.nc(nc,"TIME")
    print("reading without-LU-variable (example) ...")
    nep.nolu <- var.get.nc(nc,"nep")

    ## determine whether no-landuse simulation has k dimension
    noludim <- FALSE
    if (length(dim(nep.nolu))==3){
      print("no LU dimension in file without LUC")
      noludim <- TRUE
      count.nolu[3] <- 1
      if (offset) { count.offset[3] <- 1 }
    } else {
      print("full LU dimension in file without LUC")
      llu.nolu <- dim(nep.nolu)[3]
    }

    ## Determine domain to be read (time)
    if (!is.na(tstart)&&!is.na(tend)){
      print(paste("cutting to years ",tstart,tend))
      start.nolu <- dim(nep.nolu)
      start.nolu[] <- 1
      start.nolu[length(start.nolu)] <- which.min(abs(time-tstart))
      count.nolu <- dim(nep.nolu)
      count.nolu[length(count.nolu)] <- which.min(abs(time-tend))-which.min(abs(time-tstart))+1
    } else {
      print(paste("reading all time steps "))
      start.nolu <- dim(nep.nolu)
      start.nolu[] <- 1
      count.nolu <- dim(nep.nolu)
    }
    rm(nep.nolu)
    gc() # garbage collection

    
    ## Cut time vector to desired length
    if (!is.na(tstart)&&!is.na(tend)){
      time <- time[which.min(abs(time-tstart)):which.min(abs(time-tend))]
    }
  
    ## Now reading variables from non-LU file
    print("reading NEP ...")
    nep.nolu <- var.get.nc(nc,"nep"
                           ,start=start.nolu,count=count.nolu
                           )
    print("reading LU_AREA ...")
    luarea.nolu <- var.get.nc(nc,"lu_area"
                              ,start=start.nolu,count=count.nolu
                              )
    if (offset) {
      ## print("start")
      ## print(start)
      ## print("count.offset")
      ## print(count.offset)
      print("reading C pool variables ...")
      littera.nolu <- var.get.nc(nc,"littercarbon_ag",start=start.nolu,count=count.offset)
      litterb.nolu <- var.get.nc(nc,"littercarbon_bg",start=start.nolu,count=count.offset)
      vegc.nolu <- var.get.nc(nc,"vegcarbon",start=start.nolu,count=count.offset)
      soilc.nolu <- var.get.nc(nc,"soilcarbon",start=start.nolu,count=count.offset)
      exuc.nolu <- var.get.nc(nc,"exucarbon",start=start.nolu,count=count.offset)
      totc.nolu <- littera.nolu+litterb.nolu+vegc.nolu+soilc.nolu+exuc.nolu
      rm(littera.nolu,litterb.nolu,vegc.nolu,soilc.nolu,exuc.nolu)
      gc()
    }
    close.nc(nc)
    gc() # garbage collection

    ## get gridcell area
    nc <- open.nc( paste(dir,name.nolu,".cdf",sep="") )
    area <- var.get.nc(nc,"area")
    close.nc(nc)
    gc() # garbage collection
    
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
    llon  <- dim(luarea.lu)[1]
    llat  <- dim(luarea.lu)[2]
    llu   <- dim(luarea.lu)[3]
    ltime <- dim(luarea.lu)[4]
            
    ## element-wise multiplication to weigh flux densities with
    ## land use area fractions [-> gC/m2/yr per gridcell and land use class]
    print("multiplying by luarea ...")
    lunep.nolu <- nep.nolu*luarea.nolu
    lunep.lu <- nep.lu*luarea.lu
    lucflux.prod.lu <- cflux.prod.lu
    if (offset){
      ## print(dim(totc.lu))
      ## print(dim(luarea.lu))
      lutotc.lu <- totc.lu*luarea.lu[,,,1]
      if (noludim){
        lutotc.nolu <- totc.nolu*luarea.nolu[,,1]
      } else {
        lutotc.nolu <- totc.nolu*luarea.nolu[,,,1]
      }
    }

    ## Reduce dimension: sum over land use class
    ## 'apply' function did not work for some reason [-> gC/m2/yr per gridcell]
    print("sum over luarea dimension ...")
    if (!noludim) {
      lunep.nolu.lusum <- array(NA,dim=c(dim(lunep.nolu)[1:2],dim(lunep.nolu)[4]))
      if (offset) {
        lutotc.nolu.lusum <- array(NA,dim=c(dim(lutotc.nolu)[1:2]))
      }
    } else {
      lunep.nolu.lusum <- lunep.nolu
      if (offset) {
        lutotc.nolu.lusum <- lutotc.nolu
      }
    }
    lunep.lu.lusum <- array(NA,dim=c(dim(lunep.lu)[1:2],dim(lunep.lu)[4]))
    lucflux.prod.lu.lusum <- array(NA,dim=c(dim(lucflux.prod.lu)[1:2],dim(lucflux.prod.lu)[4]))
    if (offset){
      lutotc.lu.lusum <- array(NA,dim=c(dim(lutotc.lu)[1:2]))
    }
    for (ilon in seq(llon)){
      for (ilat in seq(llat)){
        if (!is.na(luarea.lu[ilon,ilat,1,1])){
          for (itime in seq(ltime)){
            if (!noludim) {
              lunep.nolu.lusum[ilon,ilat,itime] <- sum(lunep.nolu[ilon,ilat,,itime], na.rm=TRUE)
              if (offset && itime==1){
                lutotc.nolu.lusum[ilon,ilat] <- sum(lutotc.nolu[ilon,ilat,], na.rm=TRUE )
              }
            }
            lunep.lu.lusum[ilon,ilat,itime] <- sum(lunep.lu[ilon,ilat,,itime], na.rm=TRUE)
            lucflux.prod.lu.lusum[ilon,ilat,itime] <- sum(lucflux.prod.lu[ilon,ilat,,itime], na.rm=TRUE)
            if (offset && itime==1){
              lutotc.lu.lusum[ilon,ilat] <- sum( lutotc.lu[ilon,ilat,]+prodc.lu[ilon,ilat,], na.rm=TRUE )
            }
          }
        }
      }
    }
    
    ## if (!noludim) {
    ##   lunep.nolu <- apply( lunep.nolu, c(1,2,4), FUN=sum )
    ## }
    ## lunep.lu <- apply( lunep.lu, c(1,2,4), FUN=sum )
    ## lucflux.prod.lu <- apply( lucflux.prod.lu, c(1,2,4), FUN=sum )
    
    ## Calculate land use C emissions: Difference of total terrestrial balance
    ## in the simulation with and without land use change activated.
    ## In the simulation with land use change, C flux from product decay has to
    ## be accounted for in addition of NEP.
    ## [-> gC/m2/yr per gridcell, LU flux]
    print("calculate fLUC as a field")
    f.luc.field <- lunep.nolu.lusum - (lunep.lu.lusum-lucflux.prod.lu.lusum )
    if (!is.na(mask.id)){
      print(paste("masking out for continent",mask.id))
      f.luc.field.out <- f.luc.field*NA
      for (itime in seq(ltime)){
        for (icont in seq(nconts)){
          if (mask.id=="EU"||mask.id=="eu"||mask.id=="europe") {
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,5]
          } else if (mask.id=="RU"||mask.id=="ru"||mask.id=="russia") {
            ## including former soviet union
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,8]
          } else if (mask.id=="LA"||mask.id=="la"||mask.id=="latin america") {
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,1]
          } else if (mask.id=="AF"||mask.id=="af"||mask.id=="africa") {
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,2]
          } else if (mask.id=="AS"||mask.id=="as"||mask.id=="asia") {
            ## Only south asia (without russia, former soviet union and china
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,3]
          } else if (mask.id=="CN"||mask.id=="cn"||mask.id=="china") {
            ## including japan
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,4]
          } else if (mask.id=="AU"||mask.id=="au"||mask.id=="australia") {
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,6]
          } else if (mask.id=="NA"||mask.id=="na"||mask.id=="North America") {
            f.luc.field.out[,,itime] <- f.luc.field[,,itime]*mask[,,7]
          } else {
            print("Mask identifier not valid")
          }
        }
      }
    } else {
      f.luc.field.out <- f.luc.field
    }
    if (offset) {
      offset.field <- lutotc.nolu.lusum - lutotc.lu.lusum
    }

    ## multiply with grid cell area to get absolute fluxes
    ## [-> gC/yr per gridcell]
    print("multiply field with grid cell area to get absolute values ...")
    f.luc.abs <- f.luc.field.out*NA
    for (itime in seq(ltime)){
      f.luc.abs[,,itime] <- f.luc.field.out[,,itime]*area
    }
    if (offset) { offset.abs <- offset.field*area }
    
    ## sum over all gridcells for each time step individually
    ## [-> PgC/yr global total]
    print("sum over all gridcells to get global total ...")
    ##f.luc <- apply( f.luc.abs, c(1), FUN=sum, na.rm=TRUE )
    f.luc <- rep(NA,ltime)
    for (itime in seq(ltime)){
      f.luc[itime] <- sum( f.luc.abs[,,itime], na.rm=TRUE )/1e15
    }
    if (offset) { offset.init <- sum( offset.abs, na.rm=TRUE )/1e15 }
    
    ## return land use flux time series 'f.luc.tseries', and spatially resolved field
    ## 'f.luc.field.out' attached to the list 'out.f.luc'.
    ## this tseries gives identical results as tseries when computed with ascii output
    ## (within numerical imprecision of ~1e-8 GtC/yr)
    print("attach output to list 'out.f.luc' ...")
    out.f.luc <- list()
    f.luc.tseries <- data.frame( time=time, f.luc=f.luc, cumf.luc=cumsum(f.luc) )
    out.f.luc$tseries <- f.luc.tseries
    out.f.luc$field <- f.luc.field.out
    out.f.luc$lon <- lon
    out.f.luc$lat <- lat
    if (offset){
      out.f.luc$offset.field <- offset.field
      out.f.luc$offset.init <- offset.init
    }
    
    return(out.f.luc)
        
  } else {
    ## ////////////////////////////////////////////
    ## use ASCII outputs
    ## --------------------------------------------
    print("reading LPX ASCII output ...")

    if (byflux) {
      ## --------------------------------------------
      ## compute by flux (default)
      ## --------------------------------------------
      
      tmp      <- read.table( paste( dir, "trans_", name.nolu, ".nep.out", sep="") )
      print(paste( dir, "trans_", name.nolu, ".nep.out", sep="" ))
      time.nolu<- tmp[,1]
      nep.nolu <- tmp[,2]
      if (grossfluxes) {
        tmp      <- read.table( paste( dir, "trans_", name.nolu, ".rh.out", sep="") )
        rh.nolu <- tmp[,2]
        tmp      <- read.table( paste( dir, "trans_", name.nolu, ".cflux_fire.out", sep="") )
        fire.nolu <- tmp[,2]
        tmp      <- read.table( paste( dir, "trans_", name.nolu, ".npp.out", sep="") )
        npp.nolu <- tmp[,2]
      }

      tmp      <- read.table( paste( dir, "trans_", name.lu,   ".nep.out", sep="") )
      nep.lu   <- tmp[,2]
      time.lu  <- tmp[,1]
      tmp      <- read.table( paste( dir, "trans_", name.lu, ".cflux_prod.out", sep="") )
      cflux.prod.lu <- tmp[,2]
      if (grossfluxes) {
        tmp      <- read.table( paste( dir, "trans_", name.lu, ".rh.out", sep="") )
        rh.lu    <- tmp[,2]
        tmp      <- read.table( paste( dir, "trans_", name.lu, ".cflux_fire.out", sep="") )
        fire.lu    <- tmp[,2]
        print(paste("opening ",dir, "trans_", name.lu, ".npp.out", sep=""))
        tmp      <- read.table( paste( dir, "trans_", name.lu, ".npp.out", sep="") )
        npp.lu    <- tmp[,2]
      }

    } else {
      ## --------------------------------------------
      ## compute by stocks (currently not possible:
      ## no ascii output for product pool)
      ## --------------------------------------------
      tmp      <- read.table( paste( dir, "trans_", name.nolu, ".totc.out", sep="") )
      time.nolu<- tmp[,1]
      nep.nolu <- tmp[,2]*NA
      nep.nolu[2:length(nep.nolu)] <- tmp[1:(length(nep.nolu)-1),2] - tmp[2:length(nep.nolu),2]
    
      tmp      <- read.table( paste( dir, "trans_", name.lu,   ".totc.out", sep="") )
      nep.lu   <- tmp[,2]*NA
      nep.lu[2:length(nep.nolu)] <-tmp[1:(length(nep.nolu)-1),2] - tmp[2:length(nep.nolu),2]
      time.lu  <- tmp[,1]
      tmp      <- read.table( paste( dir, "trans_", name.lu, ".cflux_prod.out", sep="") )
      cflux.prod.lu <- tmp[,2]      

    }

    if ( is.na(tstart) && is.na(tend) ){
      if (length(time.lu)!=length(time.nolu)) {
        print("selecting subset of time series available in both datasets")
        tstart.nolu <- head(time.nolu,n=1)
        tend.nolu <- tail(time.nolu,n=1)
        tstart.lu <- head(time.lu,n=1)
        tend.lu <- tail(time.lu,n=1)
        if (tstart.nolu < tstart.lu) {
          tstart <- tstart.lu
        } else {
          tstart <- tstart.nolu
        }
        if (tend.nolu < tend.lu){
          tend <- tend.nolu
        } else {
          tend <- tend.lu
        }
        nep.nolu  <- nep.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
        nep.lu <- nep.lu[which(time.lu==tstart):which(time.lu==tend)]
        cflux.prod.lu <- cflux.prod.lu[which(time.lu==tstart):which(time.lu==tend)]
        time.lu <- time.lu[which(time.lu==tstart):which(time.lu==tend)]
        if (grossfluxes) {
          rh.nolu  <- rh.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
          fire.nolu  <- fire.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
          npp.nolu  <- npp.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
          rh.lu <- rh.lu[which(time.lu==tstart):which(time.lu==tend)]
          fire.lu <- fire.lu[which(time.lu==tstart):which(time.lu==tend)]
          npp.lu <- npp.lu[which(time.lu==tstart):which(time.lu==tend)]
        }
        time.nolu <- time.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
        time.lu <- time.lu[which(time.lu==tstart):which(time.lu==tend)]
      }
    } else {
      nep.nolu  <- nep.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
      nep.lu    <- nep.lu[which(time.lu==tstart):which(time.lu==tend)]
      cflux.prod.lu <- cflux.prod.lu[which(time.lu==tstart):which(time.lu==tend)]
      if (grossfluxes) {
        rh.nolu <- rh.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
        fire.nolu <- fire.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
        npp.nolu <- npp.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
        rh.lu <- rh.lu[which(time.lu==tstart):which(time.lu==tend)]
        fire.lu <- fire.lu[which(time.lu==tstart):which(time.lu==tend)]
        npp.lu <- npp.lu[which(time.lu==tstart):which(time.lu==tend)]
      }
      time.nolu <- time.nolu[which(time.nolu==tstart):which(time.nolu==tend)]
      time.lu   <- time.lu[which(time.lu==tstart):which(time.lu==tend)]
    }
    
    out.f.luc <- list()
    f.luc <- data.frame( time=time.nolu, f.luc=(nep.nolu - (nep.lu-cflux.prod.lu)), nbp.nolu=nep.nolu, nbp.lu=(nep.lu-cflux.prod.lu) )
    out.f.luc$tseries <- f.luc
    if (grossfluxes) {
      deforflux <- data.frame( time=time.nolu, deforflux = (nep.nolu-npp.nolu - (nep.lu-npp.lu-cflux.prod.lu)) )
      out.f.luc$deforflux <- deforflux
      regrowth <- data.frame( time=time.nolu, regrowth = npp.lu - npp.nolu )
      out.f.luc$regrowth <- regrowth
    }
  }

  return(out.f.luc)
  
}
