      module outvars
      !////////////////////////////////////////////////////////////////
      !  Module contains all daily variables. These variables have no 
      !  spatial dimension. At the end of the daily loop, values may be 
      !  copied to spatial arrays.
      !----------------------------------------------------------------
        use params_core
        use classdefs

        implicit none

        !----------------------------------------------------------------
        ! Daily output variables
        !----------------------------------------------------------------
        real, dimension(npft,ndayyear,maxgrid) :: outdgpp
        real, dimension(npft,ndayyear,maxgrid) :: outdnpp
        real, dimension(npft,ndayyear,maxgrid) :: outdnup
        real, dimension(npft,ndayyear,maxgrid) :: outdcleaf
        real, dimension(npft,ndayyear,maxgrid) :: outdcex

        real, dimension(nlu,ndayyear,maxgrid)  :: outdnetmin
        real, dimension(nlu,ndayyear,maxgrid)  :: outdninorg

        !----------------------------------------------------------------
        ! Annual output variables
        !----------------------------------------------------------------
        real, dimension(npft,maxgrid) :: outagpp
        real, dimension(npft,maxgrid) :: outanpp
        real, dimension(npft,maxgrid) :: outanup
    
      end module outvars
