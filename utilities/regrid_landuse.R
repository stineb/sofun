####################################################################
## lu             3-dim array containing absolute landuse area on original grid (first two 
##		  dimensions). The third dimension holds the different landuse categories (crop, 
##		  pasture, ...) for which consistent fields are calculated (their sum cannot exceed 
##		  available land on output grid).
## land.avail     2-dim array containing the gridcell area fraction of available land on desired 
##                output grid
## fraction=FALSE defaults to FALSE, TRUE when 'lu' and 'land.avail' are fractions of grid cell area.
##                Note that sum of lu over categories can be greater than available land area used
##                in land model (e.g. LPX) when fraction=FALSE. This is due to inconsistency between
##                land mask on which land use data is provided and land mask of land model.
## mass=FALSE     defaults to FALSE, set to TRUE if used to regrid harvested mass (no "space" constraint).
##                In this case, also make sure that fraction is set to FALSE.
## harvest=FALSE  Can be set to true when regridding area-based harvest data: no aborting but reducing
##                areas when negative areas would occur.
## aligned=NA     Flag to declarea whether original and output grids are aligned (natural number of
##                original gridcells fit inside one output gridcell, AND grids have common lower
##                left corner (south-west corner)).
## lon=NA	  longitudes (grid cell center) of original grid. If NA, equally spaced grid points 
##                covering -180 to 180 degrees E are assumed, with length corresponding to number of 
##		  columns of 'lu'
## lat=NA	  latitudes (grid cell center) of original grid. If NA, equally spaced grid points 
##		  covering -90 to 90 degrees N are assumed, with length corresonding to number of 
##		  rows of 'lu'.
## dx=NA          Explicit declaration of original gridcell size in longitude. May be necessary
##                to avoid numerical imprecisions.
## dy=NA          Explicit declaration of original gridcell size in latitude. May be necessary
##                to avoid numerical imprecisions.
## lono=NA	  longitudes (grid cell center) of destination grid. If NA, equally spaced grid 
##		  points covering -180 to 180 degrees E are assumed, with length corresponding to 
##		  number of columns of land.avail.
## lato=NA	  latitudes (grid cell center) of destination grid. If NA, equally spaced grid 
##		  points covering -180 to 180 degrees E are assumed, with length corresponding to 
##		  number of columns of land.avail.
##
## Beni Stocker, 17.07.2013
## -----------------------------------------------------------------

regrid.landuse <- function(
                           lu,
                           land.avail,
                           fraction=FALSE,
                           mass=FALSE,
                           harvest=FALSE,
                           aligned=NA,
                           lon=NA,
                           lat=NA,
                           dx=NA,
                           dy=NA,
                           lono=NA,
                           lato=NA,
                           verbose=FALSE
                           ){
                               	                               	
  ## /////////////////////////////////////////////////////////////////
  ## INITIALISATIONS AND DEFINITIONS
  ## -----------------------------------------------------------------
  if (mass && fraction) {
    print("Set 'fraction' to FALSE when using 'mass'!")
    stop
  }
  
  ## Get total number of landuse categories for which consistent fields
  ## are required
  ncat <- dim(lu)[3]

  ## Input grid definition
  if (is.na(dx)){dx=lon[2]-lon[1]}
  if (is.na(dy)){dy=lat[2]-lat[1]}
  if (is.na(lon)){
    ## Interpret longitude information from 'lu'
    lon <- seq(from=-180+dx/2, to=180-dx/2, by=dx)
  } 
  if (is.na(lat)){
    ## Interpret latitude information from 'lu'
    lat <- seq(from=-90+dy/2, to=90-dy/2, by=dy)
  }

  ## Output grid definition	
  if (is.na(lono)){
    ## Interpret longitude information from 'land.avail'
    dxo <- 360/dim(land.avail)[1]
    lono <- seq(from=-180+dxo/2, to=180-dxo/2, by=dxo)
  } else {
    dxo <- lono[2]-lono[1]
  }
  if (is.na(lato)){
    ## Interpret latitude information from 'land.avail'
    dyo <- 180/dim(land.avail)[2]
    lato <- seq(from=-90+dyo/2, to=90-dyo/2, by=dyo)
  } else {
    dyo <- lato[2]-lato[1]
  }
  
  ## Define some auxiliary variables
  nlon <- length(lon)
  nlat <- length(lat)
  xl <- lon-dx/2
  xr <- lon+dx/2
  yb <- lat-dy/2
  yt <- lat+dy/2
  
  nlono <- length(lono)
  nlato <- length(lato)
  dxo <- lono[2]-lono[1]
  dyo <- lato[2]-lato[1]
  xol <- lono-dxo/2 #left edges of destination grid 
  xor <- lono+dxo/2 #right edges of dest. grid
  yob <- lato-dyo/2 #bottom edges of dest. grid
  yot <- lato+dyo/2 #top edges of dest. grid
    
  ## Initialize land use area on output grid
  aluo <- array( 0, dim=c(nlono,nlato,ncat) )

  ## Convert original landuse data to to absolute areas
  if (fraction) {
    if (verbose) { print("converting landuse fractions to absolute area (m2) ...") }
    for (i in seq(nlon)){
      for (j in seq(nlat)){
        lu[i,j,] <- lu[i,j,]*area(lat[j],dx,dy)/1e6                  # converted to km2
      }
    }
  }

  ## Define gridcell areas like they are calculated in LPX using function area().
  ## Use your own function if gridcell area is defined differently from LPX.
  if (!mass){
    if (verbose) {print("calculating destination gridcell area...")}
    area.lpx <- land.avail*NA
    for (k in seq(nlono)){
      for (l in seq(nlato)){
        if (fraction) {
          land.avail[k,l] <- land.avail[k,l]*area(lato[l],dxo,dyo)/1e6  # converted to km2
        }
        area.lpx[k,l] <- area(lato[l],dxo,dyo)/1e6                    # converted to km2
      }
    }
  }
  
  ## Check if original and destination grids are aligned (no overlapping cells)
  if (is.na(aligned)){
    if (verbose) {print("check if aligned...")}
    if (xl[1]==xol[1] && yb[1]==yob[1] &&
        dxo%%dx==0 && dyo%%dy==0 ) {aligned <- TRUE} else {aligned <- FALSE}
    if (verbose) {print(paste("...",aligned))}
  }
  
  ## Shift destination grid so that lower left corner of lower left grid cell
  ## lies in ursprung (0,0).
  xol_shifted <- xol-xl[1]  ## xl[1] is a negative number
  yob_shifted <- yob-yb[1]  ## yb[1] is a negative number

  if (!fraction){
    if (verbose) {
      print(paste(
                  "total landuse area in original file:",sum(lu[,,1],na.rm=TRUE)
                ))
    }
  }
 
  ## /////////////////////////////////////////////////////////////////
  ## AREA ALLOCATION TO DESTINATION GRID
  ## -----------------------------------------------------------------
  ## Loop over coarse grid.
  for (l in seq(nlato)){
    for (k in seq(nlono)){
      if (land.avail[k,l]>0){
      
        ## Find all (smaller) grid cells of original data that are either adjacent
        ## or fully inside the coarse grid cell.
        ## First, restrict search over (i,j) to limited domain given the geographical
        ## position of grid cell (k,l) and grid cell sizes (dx,dy, and dxo, dyo)
        starti <- (xol_shifted[k]-(xol_shifted[k] %% dx))/dx+1
        
        ##stopi <- starti+(dxo-(dxo %% dx))/dx-1
        stopi <- starti+(dxo/dx-(dx %% dxo)+1)  ## This should bi 8 in this case
        
        startj <- (yob_shifted[l]-(yob_shifted[l] %% dy))/dy+1
        
        ##stopj <- startj+(dyo-(dyo %% dy))/dy-1
        stopj <- startj+(dyo/dy-(dyo %% dy)+1)  ## This should bi 6 in this case
        
        ## print("===============================")
        ## print(paste("Cell",k,l,"extends over"))
        ## print(paste("x:",xol[k],":",xor[k]))
        ## print(paste("y:",yob[l],":",yot[l]))
        ## print("looking in:")
        ## print("-------------------------------")
        ## print(paste("Cells",starti,":",stopi,"/",startj,":",stopj,"extending over"))
        ## print(paste("in X:",xl[starti],":",xr[stopi]))
        ## print(paste("in Y:",yb[startj],":",yt[stopj]))
        
        for (i in starti:stopi){
          for (j in startj:stopj){
            
            ## XXX CAUTION: this is a hack to make life easier. Here, for HYDE data!
            ## Destination grid may extend beyond 180deg E, thus "wrap" search around globe.
            if (i > nlon) {
              i <- i %% nlon
              xl[i] <- (i-1)*dx
              xr[i] <- i*dx
            }

            ## ## Last elements in longitude and latitude are missing in HYDE data
            ## j <- min(nlat,j)
            ## i <- min(nlon,i)
            
            if (!is.na(lu[i,j,1])){
              
              if (aligned) {              
                
                a_to_dest <- 1	
                
              } else {
                
                ## Find all source grid cells that fall at least partially inside
                ## destination grid.
                isonleftmargin <- FALSE
                isonrightmargin <- FALSE
                isbetweenlrmargin <- FALSE
                isonlowermargin <- FALSE
                isontopmargin <- FALSE
                isbetweenbtmargin <- FALSE
                
                if(xol[k]>xl[i] && xol[k]<xr[i]){
                  isonleftmargin <- TRUE
                } else if (xor[k]>xl[i] && xor[k]<xr[i]){
                  isonrightmargin <- TRUE
                } else if (xol[k]<=xl[i] && xor[k]>=xr[i]){
                  isbetweenlrmargin <- TRUE
                }
                if(yob[l]>yb[j] && yob[l]<yt[j]){
                  isonlowermargin <- TRUE
                } else if(yob[l]<=yb[j] && yot[l]>=yt[j]){
                  isbetweenbtmargin <- TRUE
                } else if(yot[l]>yb[j] && yot[l]<yt[j]){
                  isontopmargin <- TRUE
                }
                
                ## Calculate area fraction from source grid cell (i,j) to
                ## be allocated to destination grid (k,l) (a_to_dest)
                ## according to fraction overlap.
                ## This assumes rectangular grid cells.
                scale_lowermargin <- (yt[j]-yob[l])/dy
                scale_topmargin <- (yot[l]-yb[j])/dy
                scale_leftmargin <- (xr[i]-xol[k])/dx
                scale_rightmargin <- (xor[k]-xl[i])/dx
                
                a_to_dest <- 0
                if(isonleftmargin){
                  if(isonlowermargin){
                    a_to_dest <- scale_leftmargin*scale_lowermargin
                  }else if(isbetweenbtmargin){
                    a_to_dest <- scale_leftmargin
                  }else if(isontopmargin){
                    a_to_dest <- scale_leftmargin*scale_topmargin
                  } 
                } else if(isbetweenlrmargin){
                  if(isonlowermargin){
                    a_to_dest <- scale_lowermargin
                  }else if(isbetweenbtmargin){
                    a_to_dest <- 1
                  }else if(isontopmargin){
                    a_to_dest <- scale_topmargin
                  } 
                } else if(isonrightmargin){
                  if(isonlowermargin){
                    a_to_dest <- scale_rightmargin*scale_lowermargin
                  }else if(isbetweenbtmargin){
                    a_to_dest <- scale_rightmargin
                  }else if(isontopmargin){
                    a_to_dest <- scale_rightmargin*scale_topmargin
                  } 
                }
              }
              
              ## Allocate land use area from source grid to destination grid
              ## and reduce land use area and land use area fraction on
              ## source grid, after it is allocated. Thus, if all land use
              ## areas can be allocated (no source land use area outside
              ## domain of destination grid), source land use areas should
              ## sum up to 0.
              for (m in seq(ncat)){
                if (mass){
                  ## No "space" constraint for mass-regridding (harvest)
                  fill <- lu[i,j,m] * a_to_dest
                  aluo[k,l,m] <- aluo[k,l,m] + fill
                  lu[i,j,m]   <- lu[i,j,m]   - fill
                } else {
                  ## Determine available land not yet used by other LU types
                  space <- land.avail[k,l]-sum(aluo[k,l,],na.rm=TRUE)
                  if (space>0){
                    fill <- min( space, lu[i,j,m]*a_to_dest )
                    aluo[k,l,m] <- aluo[k,l,m] + fill
                    lu[i,j,m]   <- lu[i,j,m]   - fill
                  }
                }
              }
            }
          } # end of startj:stopj loop
        } # end of starti:stopi loop
        
        
        ## Check if sum of land use areas fraction aluo[k,l,m] in destination grid cell
        ## is >1. This should not be the case, as lu[i,j,m] is always <1 and no other
        ## areas have been allocated to (k,l) at this point.
                                        #if (areao[k,l]>(land.avail[k,l])) {
        if (!mass && land.avail[k,l]-sum(aluo[k,l,])<(-1e-8)) {
          stop(print(paste(
                           "Too much land use area allocated to grid cell k=",
                           k,", l=",l, "landuse area =", sum(aluo[k,l,]),
                           "total available land area =",land.avail[k,l]
                           )))
        }                         
      }              
    }
  }
  
  ## Search for all cropland/pasture/urban areas lu(i,j,m) that have not
  ## yet been allocated to any destination grid cell (k,l). Such cells
  ## must lie outside the domain of the destination grid. Thus, they cannot
  ## be strictly associated to any destination grid cell, but are eqally
  ## distributed over all destination cells of the same or nearest latitude.
  ## print("checking if all LU areas have been allocated...")
  if (verbose) {print(paste(
              "total cropland in output file after initial round:",sum(aluo[,,1],na.rm=TRUE)
              ))}
  if (verbose) {print(paste(
              "total cropland area not yet allocated (km2):",
              sum(lu[,,1],na.rm=TRUE)
              ))}
  
  ## COMPUTATIONALLY EXPENSIVE! (but perfect)
  ileftover <- which(apply(lu,c(1),FUN=sum,na.rm=TRUE)>0)
  jleftover <- which(apply(lu,c(2),FUN=sum,na.rm=TRUE)>0)
  for (j in jleftover){
    for (i in ileftover){
      for (m in seq(ncat)){
      ## Check if there are land use areas not yet allocated here
        if (!is.na(lu[i,j,m])){
          if (lu[i,j,m]>0){
            ## Find nearest longitude and latitude in destination grid
            k <- which.min(abs(lono-lon[i]))
            l <- which.min(abs(lato-lat[j]))
            ## Search algorithm: jump from east to west of original position
            for (n in seq(2*nlono)){
              k_search <- (-1)^(n+1)*round((n+0.1)/2)+k
              if (k_search > nlono) {k_search <- k_search %% nlono} ## Wrap search around globe in latitudinal direction
              else if (k_search < 1) {k_search <- k_search + nlono}
              if (mass){
                ## No "space" constraint for mass-regridding (harvest)
                fill <- lu[i,j,m]
                if (!is.na(fill)) {
                  aluo[k_search,l,m] <- aluo[k_search,l,m] + fill
                  lu[i,j,m]          <- lu[i,j,m]          - fill
                } 
                if (lu[i,j,m]==0){break}
              } else {
                ## Determine available land not yet used
                space <- land.avail[k_search,l]-sum(aluo[k_search,l,],na.rm=TRUE)
                if (space>0){
                  fill <- min( space, lu[i,j,m] )
                  if (!is.na(fill)) {
                    aluo[k_search,l,m] <- aluo[k_search,l,m] + fill
                    lu[i,j,m]          <- lu[i,j,m]          - fill
                    if (!harvest && aluo[k_search,l,m]<0) {
                      stop("negative landuse area")
                    }
                  } 
                  if (lu[i,j,m]==0){break}
                }
              }
            }
          }
        }
      }
    }
  }
  if (verbose) {print(paste(
              "total cropland in output file in final round:     ",sum(aluo[,,1],na.rm=TRUE)
              ))}
  if (verbose) {print(paste(
              "total cropland area STILL not allocated (km2):    ",
              sum(lu[,,1],na.rm=T)
              ))}

  ## When fraction==TRUE, 'aluo' is already a gridcell fraction,
  ## otherwise, compute relative from absolute (aluo.rel)
  aluo.rel <- aluo
  if (!mass){
    for (k in seq(nlono)){
      for (l in seq(nlato)){
        for (m in seq(ncat)){
          aluo.rel[k,l,m] <- aluo[k,l,m]/area.lpx[k,l]
        }
      }
    } 	
  }
  
  ## Fill up all ocean cells with NA
  aluo[land.avail==0] <- NA
  aluo.rel[land.avail==0] <- NA

  ## Return a list of variables
  out.regrid.landuse <- list()
  out.regrid.landuse$lon <- lono
  out.regrid.landuse$lat <- lato
  out.regrid.landuse$lu  <- aluo
  if (!mass) {
    out.regrid.landuse$lu.rel  <- aluo.rel
    out.regrid.landuse$land.avail <- land.avail
  }
  return(out.regrid.landuse)
  
}
