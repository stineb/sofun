plot_biome <- function( biome, time, lon, lat, outyear, plotfil ){
  ## ----------------------------------------------------------------
  ## PLOT THE NICE BIOME MAP
  ## Use in combination with biomisation algorithm 'biomisation.R'.
  ## Creates a plot ('plotfil').
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
  source("../utilities/mycolorbar.R")
  library(fields)
  library(sp)
  library(maps)

  if (dim(biome)[3]!=length(time)) {
    print( "length of time dimension of array 'biome' doesn't correspond to lenth of vector 'time'.")
  }

  iyr <- max( 1, which.min( abs( floor(time)-((outyear+0.5)-30) ) ) )

  print(paste("plot biomes for year", outyear, "..."))

  ## Replace NAs with 0 (ocean)
  biome[ is.na(biome) ] <- 0

    ## 1st Color key
  color <- c(
              "grey70",#            Ocean
              "darkgreen",#         Tropical forest
              "paleturquoise4",#    Warm temperate forest
              "limegreen",#         Temperate forest
              "royalblue4",#        Boreal forest
              "darkorange4",#       Tropical Savanna
              "darkolivegreen4",#   Sclerophyll woodland
              "palegreen",#         Temperate parkland
              "cadetblue1",#        Boreal parkland
              "yellow",#            Desert
              "darkorange",#        Dry grassland/shrubland
              "mediumpurple1",#     Shrub tundra
              "plum1"#              Tundra
              )

  label <- c(
             "Ocean/Ice", 
             "Tropical forest",
             "Warm temperate forest",
             "Temperate forest",
             "Boreal forest",
             "Tropical Savanna",
             "Sclerophyll woodland",
             "Temperate parkland",
             "Boreal parkland",
             "Desert",
             "Dry grassland/shrubland",
             "Shrub tundra",
             "Tundra"
             )

  ## lay out panels
  magn <- 3.5
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[2] <- 0.3*widths[2]
  heights <- rep(magn,nrows)
  order <- matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE)

  lev   <- (0:13)-0.5
  out.mycolorbar <- mycolorbar( color, lev, plot=FALSE )

  ## plot only latitudes ...
  ylim <- c( -60, 90 )

  pdf( plotfil, width=sum(widths), height=sum(heights) )

    panel <- layout(
                    order,
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
                    # layout.show(panel)

    ## map biomes
    par(mar=c(1,1,1,1),xaxs="i", yaxs="i",las=1)
    image(
          lon, lat,
          biome[,,iyr],
          ylim=ylim,
          yaxt="n", xaxt="n",
          col=out.mycolorbar$colors, breaks=out.mycolorbar$margins
          )

    map( add=TRUE, interior=FALSE )

    ## add legend on the right
    par( mar=c(1,0,1,0), xaxs="i", yaxs="i", las=1 )
    plot( (1:14)/14, 1:14, type="n", axes=FALSE )
    rect( rep(0,14), 0:13, rep(0.25,14), 1:14, col=color )
    text( rep(0.3,14), (0:13)+0.2, label, adj=c(0,0), cex=0.7 )

  dev.off()

}