module md_params_site
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_site
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: maxgrid

  implicit none

  private 
  public paramstype_site, getpar_site

  type paramstype_site
    real :: lon_site
    real :: lat_site
    real :: elv_site
    logical :: lTeBS
    logical :: lGrC3
    logical :: lGrC4
    integer :: soilcode
  end type paramstype_site

contains

  function getpar_site( sitename ) result( out_getpar_site )
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, getparlogical, getparint

    ! arguments
    character(len=*) :: sitename

    ! function return variable
    type( paramstype_site ) :: out_getpar_site
  
    write(0,*) 'reading site parameter file ', 'site_paramfils/'//trim(sitename)//'.parameter ...'

    ! out_getpar_site%lon_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'longitude' )
    ! out_getpar_site%lat_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'latitude' )
    ! out_getpar_site%elv_site = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'altitude' )

    ! order is important and must be equal to order in site parameter file
    ! as well as in getpar_modl()
    out_getpar_site%lTeBS = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lTeBS' )
    out_getpar_site%lGrC3 = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lGrC3' )
    out_getpar_site%lGrC4 = getparlogical( 'site_paramfils/'//trim(sitename)//'.parameter', 'lGrC4' )

    ! dimension length (maxgrid) is 1 for site-scale simulations
    out_getpar_site%soilcode = getparint( 'site_paramfils/'//trim(sitename)//'.parameter', 'soilcode' )

    ! ! For C4 grasses, check if cultivar file is available for this site
    ! if (lGrC4) then
    !   inquire( file='site_paramfils/site_cultivars/site_cultivars_'//trim(sitename)//'.txt', exist=lcultfile )  
    ! end if

  end function getpar_site

end module md_params_site
