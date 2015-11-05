## ----------------------------------------------------------------
## BIOME CLASSIFICATION ALGORITHM FOR LPX-Bern OUTPUT
## Using variables GDD, FPC_GRID, height
## after Prentice et al., 2011
## Adopted from Fortran subroutine 'findbiome.F' implemented within
## LPX-Bern.
## b.stocker@imperial.ac.uk
## ----------------------------------------------------------------
## NR.   BIOME
## 1     Tropical forest
## 2     Warm temperate forest
## 3     Temperate forest
## 4     Boreal forest
## 5     Tropical Savanna
## 6     Sclerophyll woodland
## 7     Temperate parkland
## 8     Boreal parkland
## 9     Desert
## 10    Dry grassland/shrubland
## 11    Shrub tundra
## 12    Tundra
## ----------------------------------------------------------------
biomisation <- function( gdd, fpc_grid, hgt, lunat=1.0 ) {

  ## ----------------------------------------------------------------
  ## Arguments
  ## ----------------------------------------------------------------
  ## gdd      : vector of annual growing degree days for an arbitrary number of years preceeding (and including) current year
  ## fpc_grid : vector of fractional plant cover of current year, vector of length 9 (for each LPX-PFT)
  ## hgt      : vector of vegetation height, vector of length 9 (for each LPX-PFT)
  ## lunat    : gridcell area fraction of natural land, scalar

  ## ----------------------------------------------------------------
  ## Parameters
  ## ----------------------------------------------------------------
  ## warm/cold cutoff (GDD, Prentice et al. 2011 used 350)
  gdd_cutoff  <- 1000.0
  fpc_cutoff1 <- c( 0.1, 0.1 )
  fpc_cutoff2 <- c( 9999.0, 0.6 )
  hgt_cutoff  <- c( 10.0, 6.0, 4.0, 4.0, 2.0 )

  ## ----------------------------------------------------------------
  ## Auxiliary variables
  ## ----------------------------------------------------------------
  ## GDD averaged over previous 31 years
  gdd_avg <- mean( gdd, na.rm=TRUE )
  if (gdd_avg > gdd_cutoff) { gdd_class <- 2 } else { gdd_class <- 1 }

  ## Total FPC, tree FPC
  fpc_tree  <- sum( fpc_grid[1:7], na.rm=TRUE )
  fpc_grass <- sum( fpc_grid[8:9], na.rm=TRUE )
  fpc_tot   <- sum( fpc_grid[1:9], na.rm=TRUE )

  ## Vegetation height, fpc_grid-weighted sum over trees
  if ( fpc_tree==0.0 ) {
    hgt_avg <- 0.0
  } else {
    hgt_avg <- sum( hgt[1:7] * fpc_grid[1:7], na.rm=TRUE ) / fpc_tree
  }

  ##--------------------------------------------------------
  ## Start biome classification algorithm
  ##--------------------------------------------------------
  ## Correct for case where natural land area is ~zero
  ## NOT ANYMORE: => take fpc_grid in cropland/pasture
  ## => classify to "closest" corresponding biome
  ##--------------------------------------------------------
  if (lunat<0.001){

    biome <- NA

    # useval <- max( sum( fpc_grid[12:13], na.rm=TRUE ), sum( fpc_grid[14:15], na.rm=TRUE ), na.rm=TRUE )
    # print( paste( "use value ", useval))
    # fpc_tot <- useval

  } else if (gdd_class==2) {
            
    ## Warm Biomes
    ##********************************************************
    
    if ( fpc_tot<fpc_cutoff1[ gdd_class ] ) {
      
      ## Desert
      ##--------------------------------------------------------
      biome <- 9
      
    } else if ( fpc_tot<fpc_cutoff2[ gdd_class ] ) {
      
      ## Dry grassland/shrubland
      ##--------------------------------------------------------
      biome <- 10
      
    } else {
      
      ## Forest/Savanna
      ##--------------------------------------------------------
      if ( max( fpc_grid[1:2] )>0.0) {
        ## tropical PFTs are present
        
        if (hgt_avg<hgt_cutoff[1] ) {
          ## Tropical savanna 
          biome <- 5
          
        } else {
          ## Tropical forest
          biome <- 1
          
        }
        
      } else if (fpc_grid[4]>fpc_grid[5] || fpc_grid[4]>fpc_grid[3]) {
        ## temp. broadl. evergr. PFT dominant
        
        if (hgt_avg<hgt_cutoff[2] ) {
          ## Sclerophyll woodland
          biome <- 6
        } else {
          ## Warm temperate forest
          biome <- 2
        }
        
      } else if (fpc_grid[3]>0.0 || fpc_grid[5]>0.0) {
        ## temp. needl. or temp. summergr. PFT present
        
        if (hgt_avg<hgt_cutoff[3] ) {
          ## Temperate parkland
          biome <- 7
        } else {
          ## Temperate forest
          biome <- 3
        }
        
      } else {
        
        ## This is in "warm biomes" because there may still be "warm biome areas" where 
        ## neither temperate nor (sub)tropical PFTs grow.

        if (hgt_avg<hgt_cutoff[4] ) {
          ## Boreal parkland
          biome <- 8
        } else {
          ## Boreal forest
          biome <- 4
        }
        
      }
      
    }
    
  } else {
    
    ## Cold Biomes
    ##********************************************************
    
    if (fpc_tot<fpc_cutoff1[gdd_class]) {
      
      ## Tundra
      ##--------------------------------------------------------
      biome <- 12
      
    } else if (hgt_avg<hgt_cutoff[5]) {
      
      ## Shrub tundra
      ##--------------------------------------------------------
      biome <- 11
    
    } else {
      
      ## This is in "cold biomes" because there may still be "cold biome areas"
      ## where no boreal PFTs grow.

      ## Forest/Savanna
      ##--------------------------------------------------------
      if (max(fpc_grid[1])>0.0) {
        ## tropical PFTs are present
      
        if (hgt_avg<hgt_cutoff[1] ) {
          ## Tropical savanna 
          biome <- 5
        } else {
          ## Tropical forest
          biome <- 1
        }

      } else if (fpc_grid[4]>fpc_grid[5] || fpc_grid[4]>fpc_grid[3]) {
        ## temp. broadl. evergr. PFT dominant

        if (hgt_avg<hgt_cutoff[2] ) {
          ## Sclerophyll woodland
          biome <- 6
        } else {
          ## Warm temperate forest
          biome <- 2
        }

      } else if (fpc_grid[3]>0.0 || fpc_grid[5]>0.0) {
        ## temp. needl. or temp. summergr. PFT present

        if (hgt_avg<hgt_cutoff[3] ) {
          ## Temperate parkland
          biome <- 7
        } else {
          ## Temperate forest
          biome <- 3
        }

      } else {

        if (hgt_avg<hgt_cutoff[4] ) {
          ## Boreal parkland
          biome <- 8
        } else {
          ## Boreal forest
          biome <- 4
        }

        }
      
    }

  }

  return( biome )

}