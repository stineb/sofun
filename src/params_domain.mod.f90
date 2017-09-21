module md_params_domain
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_domain
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: maxgrid

  implicit none

  private 
  public paramstype_site, getpar_domain, type_params_domain

  type type_params_domain

    real :: lon_site
    real :: lat_site
    real :: elv_site

    character(len=256) :: filnam_landmask
    character(len=256) :: filnam_topography
    character(len=256) :: filnam_soil

  end type type_params_domain

  type( type_params_domain ) :: params_domain

contains

  function getpar_domain( sitename, spacetype ) result( out_params_domain )
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, getparlogical, getparint, getparstring

    ! arguments
    character(len=*) :: sitename
    character(len=*) :: spacetype

    ! function return variable
    type( type_params_domain ) :: out_params_domain
  
    write(0,*) 'reading site parameter file ', 'site_paramfils/'//trim(sitename)//'.parameter ...'

    if ( spacetype=='lonlat' ) then

      call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_landmask',   out_params_domain%filnam_landmask )
      call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_topography', out_params_domain%filnam_topography )
      call getparstring( 'site_paramfils/'//trim(sitename)//'.parameter', 'filnam_soil',       out_params_domain%filnam_soil )

    else

      out_params_domain%lon_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'longitude' )
      out_params_domain%lat_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'latitude' )
      out_params_domain%elv_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'altitude' )

      ! ! For C4 grasses, check if cultivar file is available for this site
      ! if (lGrC4) then
      !   inquire( file='site_paramfils/site_cultivars/site_cultivars_'//trim(sitename)//'.txt', exist=lcultfile )  
      ! end if

    end if 

  end function getpar_domain

end module md_params_domain
