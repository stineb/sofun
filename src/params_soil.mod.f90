module md_params_soil
  !////////////////////////////////////////////////////////////////
  ! Module containing soil parameters and functions to read them
  !----------------------------------------------------------------
  use md_grid, only: gridtype, domaininfo_type
  use netcdf
  use md_io_netcdf, only: check

  implicit none

  private
  public paramtype_soil, getsoil

  type paramtype_soil
    real :: whc
    real :: perc_k1
    real :: thdiff_wp
    real :: thdiff_whc15
    real :: thdiff_fc
    real :: forg
    real :: wbwp
    real :: por
    real :: fsand
    real :: fclay
    real :: fsilt
  end type

contains

  function getsoil( domaininfo, grid ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns the field of soil parameters
    !----------------------------------------------------------------
    ! arguments
    type( domaininfo_type ), intent(in) :: domaininfo
    type( gridtype ), dimension(domaininfo%maxgrid), intent(in) :: grid

    ! function return variable
    type(paramtype_soil) :: params_soil

    ! local variables
    type(paramtype_soil), dimension(9) :: soilparams_per_code

    integer :: isoilcode

    !----------------------------------------------------------------  
    ! Get soil parameters per code
    !----------------------------------------------------------------  
    do isoilcode=1,1
      soilparams_per_code(isoilcode) = get_soilparams_per_code( isoilcode )
    end do

    ! read from array to define grid type 
    params_soil = soilparams_per_code(domaininfo%soilcode)

    ! water holding capacity is read in separately from site parameter file
    params_soil%whc = domaininfo%whc

    return
    999  format (I2.2)

  end function getsoil


  function get_soilparams_per_code( soilcode ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns soil parameter values given soil code (except WHC)
    !----------------------------------------------------------------
    use md_sofunutils, only: getparreal

    ! arguments
    integer, intent(in) :: soilcode

    ! local variables
    character(len=2) :: soilcode_char

    ! function return variable
    type(paramtype_soil) :: params_soil

    write( soilcode_char, 999 ) soilcode

    params_soil%perc_k1      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'perc_k1' )
    params_soil%thdiff_wp    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_wp' )
    params_soil%thdiff_whc15 = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_whc15' )
    params_soil%thdiff_fc    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_fc' )
    params_soil%forg         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'forg' )
    params_soil%wbwp         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'wbwp' )
    params_soil%por          = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'por' )
    params_soil%fsand        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsand' )
    params_soil%fclay        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fclay' )
    params_soil%fsilt        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsilt' )

    return
    999  format (I2.2)

  end function get_soilparams_per_code

end module md_params_soil
