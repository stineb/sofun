module md_params_domain
  !////////////////////////////////////////////////////////////////
  !  Module contains site-specific parameters read by getpar_domain
  ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
  ! contact: b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use md_params_core, only: maxgrid

  implicit none

  private 
  public getpar_domain, type_params_domain

  type type_params_domain

    ! For site scale simulations
    real :: lon_site
    real :: lat_site
    real :: elv_site
    real :: whc_site          ! water holding capacity of soil
    integer :: soilcode_site  ! soil code (WHSD)

    ! For global simulations
    character(len=256) :: filnam_landmask
    character(len=256) :: filnam_topography
    character(len=256) :: filnam_soil

    character(len=256) :: domain_name

  end type type_params_domain

  type( type_params_domain ) :: params_domain

contains

  function getpar_domain( sitename ) result( out_params_domain )
    !////////////////////////////////////////////////////////////////
    !  SR reads site-specific parameters from <sitename>.parameter
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal, getparlogical, getparint, getparstring

    ! arguments
    character(len=*) :: sitename

    ! function return variable
    type( type_params_domain ) :: out_params_domain
  
    print*, 'reading site parameter file ', 'site_paramfils/'//trim(sitename)//'.parameter ...'

    out_params_domain%lon_site      = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'longitude' )
    out_params_domain%lat_site      = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'latitude' )
    out_params_domain%elv_site      = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'altitude' )
    out_params_domain%whc_site      = getparreal( 'site_paramfils/'//trim(sitename)//'.parameter', 'whc' )
    out_params_domain%soilcode_site = getparint(  'site_paramfils/'//trim(sitename)//'.parameter', 'soilcode' )
    out_params_domain%domain_name   = sitename

  end function getpar_domain

end module md_params_domain
