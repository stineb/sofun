module md_params_soil
  !////////////////////////////////////////////////////////////////
  ! Module containing soil parameters and functions to read them
  !----------------------------------------------------------------
  implicit none

  type paramtype_soil
    real :: perc_k1     
    real :: whc_eff     
    real :: thdiff_wp   
    real :: thdiff_whc15
    real :: thdiff_fc   
    real :: forg        
    real :: whc_blwp    
    real :: por         
    real :: fsand       
    real :: fclay       
    real :: fsilt       
  end type

contains

  function getsoil_field( soilcode_field ) result( params_soil_field )
    !////////////////////////////////////////////////////////////////
    ! Function returns array containing all soil parameter values
    !----------------------------------------------------------------
    use md_params_core, only: maxgrid

    ! arguments
    integer, dimension(maxgrid), intent(in) :: soilcode_field

    ! local variables
    integer :: jpngr

    ! function return variable
    type(paramtype_soil), dimension(maxgrid) :: params_soil_field

    do jpngr=1,maxgrid
      params_soil_field(jpngr) = getsoil( soilcode_field(jpngr) )
    end do

  end function getsoil_field


  function getsoil( soilcode ) result( params_soil )
    !////////////////////////////////////////////////////////////////
    ! Function returns soil parameter values given soil code
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
    params_soil%whc_eff      = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_eff' )
    params_soil%thdiff_wp    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_wp' )
    params_soil%thdiff_whc15 = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_whc15' )
    params_soil%thdiff_fc    = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'thdiff_fc' )
    params_soil%forg         = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'forg' )
    params_soil%whc_blwp     = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'whc_blwp' )
    params_soil%por          = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'por' )
    params_soil%fsand        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsand' )
    params_soil%fclay        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fclay' )
    params_soil%fsilt        = getparreal( trim('params/params_soil_sc'//soilcode_char//'.dat'), 'fsilt' )

    return
    999  format (I2.2)

  end function getsoil

end module md_params_soil

   
!       do jpngr = 1,maxgrid
! c     Take 'medium-coarse' soil type where map suggests 'organic' soil (aka peat).
! c     This is to avoid "circular reasoning" when determining distribution of peatlands
! c     online. Peatland LU class has still 'organic'-type soil parameters.
! c     Beni & Renato, May 2013.
!         if (soilcode(jpngr).eq.8) then
!           soilcode(jpngr) = 7
!         endif
!         if (maxgridcell(jpngr)) then
!           soilpar(15,jpngr) = store(2,soilcode(jpngr))      ! before termed 'whc'
!           soilpar(1,jpngr) = store(1,soilcode(jpngr))
!           soilpar(2,jpngr) = k2
!           soilpar(3,jpngr) = d1*soilpar(15,jpngr)
!           soilpar(4,jpngr) = d2*soilpar(15,jpngr)
!           soilpar(5,jpngr)=store(3,soilcode(jpngr))
!           soilpar(6,jpngr)=store(4,soilcode(jpngr))
!           soilpar(7,jpngr)=store(5,soilcode(jpngr))
! c     Mineral content = 1 - organic content - porosity
!           soilpar(8,jpngr)=1.0d0-store(6,soilcode(jpngr))
!      $         -store(8,soilcode(jpngr))
!           if (soilpar(8,jpngr).lt.1.0d-4) soilpar(8,jpngr) = 0.0d0
!           soilpar(9,jpngr)=store(6,soilcode(jpngr))
!           soilpar(10,jpngr)=store(7,soilcode(jpngr))
!           soilpar(11,jpngr)=store(8,soilcode(jpngr))        ! porosity
!           soilpar(12,jpngr)=store(2,soilcode(jpngr))+store(7,soilcode(jpngr))
!           soilpar(13,jpngr)=store(9,soilcode(jpngr))
!           soilpar(14,jpngr)=store(10,soilcode(jpngr))
!         endif
!       enddo
