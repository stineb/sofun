module md_grid
  !////////////////////////////////////////////////////////////////
  ! Module contains global variables defining the model grid.
  ! A module 'gridvars_*' must contain the following subroutines:
  ! - getgrid
  !
  ! ... and define the following variables that are global within
  ! 'sofun' (but passed on to 'biosphere' as arguments).
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core
  use md_params_domain, only: type_params_domain

  implicit none

  private
  public gridtype, getgrid

  type gridtype
    real :: lon
    real :: lat
    real :: elv
    integer :: soilcode
  end type gridtype

contains

  function getgrid( spacetype, params_domain ) result( out_grid )
    !////////////////////////////////////////////////////////////////
    ! Defines grid variables
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal
    use md_params_domain, only: type_params_domain

    ! arguments
    character(len=*) :: spacetype
    type( type_params_domain ) :: params_domain

    ! function return variable
    type( gridtype ), dimension( maxgrid ) :: out_grid

    ! local variables
    integer :: jpngr, ilon, ilat

    if ( spacetype=='lonlat' ) then
      !----------------------------------------------------------------
      ! Spatial lon-lat simulations
      !----------------------------------------------------------------
      print*,'here, I should open the halfdeg landmask and topography files'
      print*,'landmask file: ', params_domain%filnam_landmask
      print*,'topography file: ', params_domain%filnam_topography
      print*,'soil file: ', params_domain%filnam_soil

      if ( params_domain%filnam_landmask=="landmaskfile_global_halfdeg.nc" ) then

        ! xxx test
        jpngr = 0
        do ilon=1,720
          do ilat=1,360
            jpngr=jpngr+1
            out_grid(jpngr)%lon = -179.75 + 0.25 * (ilon-1)
            out_grid(jpngr)%lat = -89.75 + 0.25 * (ilat-1)
          end do
        end do

        ! xxx cut to domain given by params_site%longitude_start and params_site%longitude_end

      else if ( params_domain%filnam_landmask=="landmaskfile_global_1x1deg.nc" ) then

        ! xxx test
        jpngr = 0
        do ilon=1,360
          do ilat=1,180
            jpngr=jpngr+1
            out_grid(jpngr)%lon = -179.5 + 0.5 * (ilon-1)
            out_grid(jpngr)%lat = -89.5 + 0.5 * (ilat-1)
          end do
        end do

        ! xxx cut to domain given by params_site%longitude_start and params_site%longitude_end

      end if

      out_grid(:)%elv = 100.0
      out_grid(:)%soilcode = 1

    else
      !----------------------------------------------------------------
      ! Site-scale simulations
      !----------------------------------------------------------------

      ! xxx try. otherwise loop over sites and allocate values for each 
      ! site into vectors containing all sites
      out_grid(:)%lon = params_domain%lon_site
      out_grid(:)%lat = params_domain%lat_site
      out_grid(:)%elv = params_domain%elv_site

    end if

  end function getgrid

end module md_grid

