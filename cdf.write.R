## /////////////////////////////////////////////////////////////////////////
## Function 'cdf.write' writes an array 'var' of 2, 3 or 4 dimensions,
## covering a geographical domain ('lon','lat') into a NetCDF file. 
##
## var:    Is an array of size length(lon) x length(lat) x length(z_dim) x length(time). The order of lon (1st dim), lat (2nd dim), time (last dim) is mandatory.
## varnam: Is a string defining the name of the variable that is written into the NetCDF output.
## lon:    Is a 1-dimensional array containing the longitudes of the grid cell's mid-points.
## lat:    Is a 1-dimensional array containing the latitudes of the grid cell's mid-points.
## filnam: Is a string defining the name of the NetCDF file to be written.
## z_dim:  Is an integer indicating the position at which the field 'var' is written in the 3rd dimension ("Z"). Defaults to NA (no z-dimension is written).
## time:   Is a vector containing the time coordinates. Defaults to NA (no time dim. is written).
## make.zdim : Write a z-dimension. If var is 2-dimensional, create an additional dimension of length 1.
## make.tdim : Write a time dimension. If var is 2-dimensional, create an additional dimension of length 1.
## nvars : Number of variables in output files. Maximum 4. Must all have same dimensions as first one. Defaults to 1.
## var2 (var3, var4): Array of output like var
## varnam2 (...): Name of output variable 2 (3,4)
## vartype: Data type of variables (not dimensions!) in output file. One of the set of predefined NetCDF external data types. The valid NetCDF external data types are NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, and NC_DOUBLE. Defaults to NC_FLOAT
##
## Beni Stocker, 17.06.2013
## -------------------------------------------------------------------------
cdf.write <- function(var,varnam,
                      lon,lat,
                      filnam,
                      z_dim=NA,time=NA,
                      make.zdim=FALSE,make.tdim=FALSE,
                      nvars=1,
                      var2=NA,varnam2=NA,
                      var3=NA,varnam3=NA,
                      var4=NA,varnam4=NA,
                      var5=NA,varnam5=NA,
                      vartype="NC_FLOAT",
                      verbose=FALSE
                      ){
  library(RNetCDF)
  
  ## Get the dimensionality of 'var' and whether the dimension vector provided fit.
  vardims <- dim(var)
  if (length(lon)!=vardims[1]) {print("longitude vector provided does not match data array");stop}
  if (length(lat)!=vardims[2]) {print("latitude vector provided does not match data array");stop}

  if (length(vardims)==4){
    make.zdim==TRUE
    make.tdim==TRUE
    if (is.na(make.zdim)){
      if (is.na(z_dim)) {print("no z_dim vector provided"); stop}
    }
    if (is.na(make.tdim)){
      if (is.na(time)) {print("no time vector provided"); stop}
    }
  } else {
    if (make.zdim) {
      if (length(z_dim)==1){
        if (is.na(z_dim)) {print("no z_dim vector provided"); stop}
      }
    }
    if (make.tdim) {
      if (length(time)==1){
        if (is.na(time)) {print("no time; vector provided"); stop}
      }
    }
  }

  if (length(vardims)==3 && is.na(z_dim) && is.na(time) ){
    print("array to be written is 3-dim.")
    print("provide z_dim or time")
    stop
  }

  if ((length(vardims)==2 && make.tdim)||(length(vardims)==2 && make.zdim)) {
    if (make.tdim) {
      if (length(time)!=1){
        print("time vector can only have length 1")
      } else {
        tmp <- array(NA,dim=c(dim(var),1))
        tmp[,,1] <- var
        var <- tmp
        if (nvars>1){
          tmp <- array(NA,dim=c(dim(var2),1))
          tmp[,,1] <- var2
          var2 <- tmp
        }
        if (nvars>2){
          tmp <- array(NA,dim=c(dim(var3),1))
          tmp[,,1] <- var3
          var3 <- tmp
        }
        if (nvars>3){
          tmp <- array(NA,dim=c(dim(var4),1))
          tmp[,,1] <- var4
          var4 <- tmp
        } 
        if (nvars>4){
          tmp <- array(NA,dim=c(dim(var5),1))
          tmp[,,1] <- var5
          var5 <- tmp
        } 
      }
    } else if (make.zdim) {
      if (length(z_dim)!=1){
        print("z_dim vector can only have length 1")
      } else {
        print("extending output array by one dimension of length 1")
        tmp <- array(NA,dim=c(dim(var),1))
        tmp[,,1] <- var
        var <- tmp
        if (nvars>1){
          tmp <- array(NA,dim=c(dim(var2),1))
          tmp[,,1] <- var2
          var2 <- tmp
        }
        if (nvars>2){
          tmp <- array(NA,dim=c(dim(var3),1))
          tmp[,,1] <- var3
          var3 <- tmp
        }
        if (nvars>3){
          tmp <- array(NA,dim=c(dim(var4),1))
          tmp[,,1] <- var4
          var4 <- tmp
        } 
        if (nvars>4){
          tmp <- array(NA,dim=c(dim(var5),1))
          tmp[,,1] <- var5
          var5 <- tmp
        } 
      }
    }  
  }
  
  ## Create a new NetCDF dataset and define dimensions
  if (verbose) {print(paste("Creating file ",filnam))}
  nc <- create.nc(filnam)

  if (verbose) {print("define dimensions...")}
  dim.def.nc(nc,"LONGITUDE",length(lon))
  dim.def.nc(nc,"LATITUDE",length(lat))
  if (make.zdim){
    dim.def.nc(nc,"Z",length(z_dim))
  }
  if (make.tdim){
    dim.def.nc(nc,"TIME",length(time))
  }
  
  ## Create variables and coordinate variables
  if (verbose) {print("define dimension variables...")}
  var.def.nc(nc,"LATITUDE","NC_FLOAT","LATITUDE")
  var.def.nc(nc,"LONGITUDE","NC_FLOAT","LONGITUDE")
  if (make.zdim){
    var.def.nc(nc,"Z","NC_INT","Z")
  }
  if (make.tdim){
    var.def.nc(nc,"TIME","NC_FLOAT","TIME")
  }
  
  ## The dimensions must be in the reverse order, as the fastest
  ## varying comes first and the slowest varying last.
  if (verbose) {print("define variables and their dimensions...")}
  if (make.tdim){
    if (make.zdim){
      var.def.nc(nc,varnam,vartype,c(0,1,2,3))
      if (nvars>1) {
        var.def.nc(nc,varnam2,vartype,c(0,1,2,3))
      }
      if (nvars>2) {
        var.def.nc(nc,varnam3,vartype,c(0,1,2,3))
      }
      if (nvars>3) {
        var.def.nc(nc,varnam4,vartype,c(0,1,2,3))
      }
      if (nvars>4) {
        var.def.nc(nc,varnam5,vartype,c(0,1,2,3))
      }
    }else{
      var.def.nc(nc,varnam,vartype,c(0,1,2))
      if (nvars>1) {
        var.def.nc(nc,varnam2,vartype,c(0,1,2))
      }
      if (nvars>2) {
        var.def.nc(nc,varnam3,vartype,c(0,1,2))
      }
      if (nvars>3) {
        var.def.nc(nc,varnam4,vartype,c(0,1,2))
      }
      if (nvars>4) {
        var.def.nc(nc,varnam5,vartype,c(0,1,2))
      }
    }
  }else{
    if (make.zdim){
      print("doin dis")
      var.def.nc(nc,varnam,vartype,c(0,1,2))
      if (nvars>1) {
        var.def.nc(nc,varnam2,vartype,c(0,1,2))
      }
      if (nvars>2) {
        var.def.nc(nc,varnam3,vartype,c(0,1,2))
      }
      if (nvars>3) {
        var.def.nc(nc,varnam4,vartype,c(0,1,2))
      }
      if (nvars>4) {
        var.def.nc(nc,varnam5,vartype,c(0,1,2))
      }
    }else{
      var.def.nc(nc,varnam,vartype,c(0,1))
      if (nvars>1) {
        var.def.nc(nc,varnam2,vartype,c(0,1))
      }
      if (nvars>2) {
        var.def.nc(nc,varnam3,vartype,c(0,1))
      }
      if (nvars>3) {
        var.def.nc(nc,varnam4,vartype,c(0,1))
      }
      if (nvars>4) {
        var.def.nc(nc,varnam5,vartype,c(0,1))
      }
    }
  }

  
  ## Add attributes
  if (verbose) {print("define attributes...")}
  att.put.nc(nc,varnam,"missing_value",vartype,-9999.)
  att.put.nc(nc,varnam,"_FillValue",vartype,-9999.)
  att.put.nc(nc,"NC_GLOBAL","title","NC_CHAR",varnam)
  if (nvars>1){
    att.put.nc(nc,varnam2,"missing_value",vartype,-9999.)
    att.put.nc(nc,varnam2,"_FillValue",vartype,-9999.)
    att.put.nc(nc,"NC_GLOBAL","title","NC_CHAR",varnam2)
  }
  if (nvars>2){
    att.put.nc(nc,varnam3,"missing_value",vartype,-9999.)
    att.put.nc(nc,varnam3,"_FillValue",vartype,-9999.)
    att.put.nc(nc,"NC_GLOBAL","title","NC_CHAR",varnam3)
  }
  if (nvars>3){
    att.put.nc(nc,varnam4,"missing_value",vartype,-9999.)
    att.put.nc(nc,varnam4,"_FillValue",vartype,-9999.)
    att.put.nc(nc,"NC_GLOBAL","title","NC_CHAR",varnam4)
  }
  if (nvars>4){
    att.put.nc(nc,varnam5,"missing_value",vartype,-9999.)
    att.put.nc(nc,varnam5,"_FillValue",vartype,-9999.)
    att.put.nc(nc,"NC_GLOBAL","title","NC_CHAR",varnam5)
  }
  att.put.nc(nc,"NC_GLOBAL","history","NC_CHAR", paste("Created on", date(), "by Beni Stocker"))
  att.put.nc(nc,"LONGITUDE","units","NC_CHAR","degrees_east")
  att.put.nc(nc,"LONGITUDE","long_name","NC_CHAR","LONGITUDE")
  att.put.nc(nc,"LONGITUDE","GridType","NC_CHAR","Cylindrical Equidistant projection Grid")
  att.put.nc(nc,"LATITUDE","units","NC_CHAR","degrees_north")
  att.put.nc(nc,"LATITUDE","long_name","NC_CHAR","LATITUDE")
  att.put.nc(nc,"LATITUDE","GridType","NC_CHAR","Cylindrical Equidistant projection Grid")
  if (make.tdim){
    att.put.nc(nc,"TIME","units","NC_CHAR","years")
    att.put.nc(nc,"TIME","long_name","NC_CHAR","TIME")
  }
  if (make.zdim){
    att.put.nc(nc,"Z","units","NC_CHAR","none")
    att.put.nc(nc,"Z","long_name","NC_CHAR","none")
  }
    
  ## Put the data
  if (verbose) {print("put the data...")}

  if (verbose) {print("longitude...")}
  var.put.nc(nc,"LONGITUDE",lon, NA, NA, na.mode=0)

  if (verbose) {print("latitude...")}
  var.put.nc(nc,"LATITUDE",lat, NA, NA, na.mode=0)

  if (make.tdim){
    if (verbose) {print("time...")}
    var.put.nc(nc,"TIME",time, NA, NA, na.mode=0)
  }

  if (make.zdim){
    if (verbose) {print("z-dimension...")}
    var.put.nc(nc,"Z",z_dim, NA, NA, na.mode=0)
  }

  if (verbose) {print("variable 1...")}
  var.put.nc(nc,varnam,var)

  if (nvars>1) {
    if (verbose) {print("variable 2...")}
    var.put.nc(nc,varnam2,var2)
  }

  if (nvars>2) {
    if (verbose) {print("variable 3...")}
    var.put.nc(nc,varnam3,var3)
  }

  if (nvars>3) {
    if (verbose) {print("variable 4...")}
    var.put.nc(nc,varnam4,var4)
  }

  if (nvars>4) {
    if (verbose) {print("variable 5...")}
    var.put.nc(nc,varnam5,var5)
  }
    
  sync.nc(nc)
  close.nc(nc)
}
