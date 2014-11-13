      module outvars
      !////////////////////////////////////////////////////////////////
      !  Module contains all daily variables. These variables have no 
      !  spatial dimension. At the end of the daily loop, values may be 
      !  copied to spatial arrays.
      !----------------------------------------------------------------
        use params_core
        use classdefs

        implicit none

        real, dimension(npft,maxgrid) :: outagpp
        type(carbon), dimension(npft,maxgrid)   :: outanpp
        type(nitrogen), dimension(npft,maxgrid) :: outanup


      end module outvars
