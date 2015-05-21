module _outvars
!////////////////////////////////////////////////////////////////
!  Module contains all daily variables. These variables have no 
!  spatial dimension. At the end of the daily loop, values may be 
!  copied to spatial arrays.
!----------------------------------------------------------------
  use _params_core
  use _classdefs

  implicit none

  !----------------------------------------------------------------
  ! Daily output variables
  !----------------------------------------------------------------
  real, dimension(npft,ndayyear,maxgrid) :: outdgpp
  real, dimension(npft,ndayyear,maxgrid) :: outdnpp
  real, dimension(npft,ndayyear,maxgrid) :: outdnup
  real, dimension(npft,ndayyear,maxgrid) :: outdCleaf

  real, dimension(nlu,ndayyear,maxgrid)  :: outdninorg
  real, dimension(nlu,ndayyear,maxgrid)  :: outdtemp_soil
  
  !----------------------------------------------------------------
  ! Monthly output variables
  !----------------------------------------------------------------


  !----------------------------------------------------------------
  ! Annual output variables
  !----------------------------------------------------------------
  real, dimension(npft,maxgrid) :: outagpp
  real, dimension(npft,maxgrid) :: outanpp
  real, dimension(npft,maxgrid) :: outanup
  real, dimension(npft,maxgrid) :: outaCveg
  real, dimension(npft,maxgrid) :: outaCalloc
  real, dimension(npft,maxgrid) :: outaNalloc
  real, dimension(npft,maxgrid) :: outaCveg2lit
  real, dimension(npft,maxgrid) :: outaNveg2lit
  real, dimension(nlu, maxgrid) :: outaNinorg

contains

  subroutine initio()
    !////////////////////////////////////////////////////////////////
    ! Opens input/output files.
    !----------------------------------------------------------------
    use _params_siml, only: runname

    ! local variables
    character(len=256) :: prefix
    character(len=256) :: filnam

    prefix = "./output/"//trim(runname)

    !////////////////////////////////////////////////////////////////
    ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
    !----------------------------------------------------------------

    ! GPP
    filnam=trim(prefix)//'.d.gpp.out'
    open(101,file=filnam,err=999,status='unknown')

    ! NPP
    filnam=trim(prefix)//'.d.npp.out'
    open(102,file=filnam,err=999,status='unknown')

    ! LEAF C
    filnam=trim(prefix)//'.d.cleaf.out'
    open(103,file=filnam,err=999,status='unknown')

    ! N UPTAKE
    filnam=trim(prefix)//'.d.nup.out'
    open(104,file=filnam,err=999,status='unknown')

    ! INORGANIC N (NO3+NH4)
    filnam=trim(prefix)//'.d.ninorg.out'
    open(107,file=filnam,err=999,status='unknown')

    ! SOIL TEMPERATURE
    filnam=trim(prefix)//'.d.soiltemp.out'
    open(109,file=filnam,err=999,status='unknown')

    ! AIR TEMPERATURE
    filnam=trim(prefix)//'.d.temp.out'
    open(110,file=filnam,err=999,status='unknown')


    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------

    ! C ALLOCATED TO GROWTH 
    filnam=trim(prefix)//'.a.calloc.out'
    open(303,file=filnam,err=999,status='unknown')

    ! C VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.cveg2lit.out'
    open(307,file=filnam,err=999,status='unknown')

    ! N VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.nveg2lit.out'
    open(308,file=filnam,err=999,status='unknown')

    ! N ALLOCATED TO GROWTH 
    filnam=trim(prefix)//'.a.nalloc.out'
    open(309,file=filnam,err=999,status='unknown')

    ! GPP 
    filnam=trim(prefix)//'.a.gpp.out'
    open(310,file=filnam,err=999,status='unknown')

    ! NPP 
    filnam=trim(prefix)//'.a.npp.out'
    open(311,file=filnam,err=999,status='unknown')

    ! VEG C
    filnam=trim(prefix)//'.a.cveg.out'
    open(312,file=filnam,err=999,status='unknown')

    ! INORGANIC N (mean over days)
    filnam=trim(prefix)//'.a.ninorg.out'
    open(316,file=filnam,err=999,status='unknown')


    return


    999  stop 'INITIO: error opening output files'

  end subroutine initio


  subroutine getout_daily( jpngr, moy, doy )
    !////////////////////////////////////////////////////////////////
    ! SR called daily to sum up daily output variables.
    ! Note that output variables are collected only for those variables
    ! that are global anyway (e.g., outdcex). Others are not made 
    ! global just for this, but are collected inside the subroutine 
    ! where they are defined.
    !----------------------------------------------------------------
    use _params_core, only: ndayyear
    use _params_modl, only: npft, lu_category
    use _vars_core, only: dgpp, dnpp, dnup
    use _vars_core, only: pleaf, pninorg
    use _vars_core, only: dtemp_soil

    ! ARGUMENTS
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !print*,'netmin_litt ',netmin_litt(1)
    !print*,'netmin_soil ',netmin_soil(1)

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    outdgpp(:,doy,jpngr)   = dgpp(:)
    outdnpp(:,doy,jpngr)   = dnpp(:)%c12
    outdnup(:,doy,jpngr)   = dnup(:)%n14
    outdCleaf(:,doy,jpngr) = pleaf(:,jpngr)%c%c12

    outdninorg(:,doy,jpngr)   = pninorg(:,jpngr)%n14
    outdtemp_soil(:,doy,jpngr)= dtemp_soil(:,jpngr)

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    outagpp(:,jpngr)    = outagpp(:,jpngr) + dgpp(:)
    outanpp(:,jpngr)    = outanpp(:,jpngr) + dnpp(:)%c12
    outanup(:,jpngr)    = outanup(:,jpngr) + dnup(:)%n14
    outaninorg(:,jpngr) = outaNinorg(:,jpngr) + pninorg(:,jpngr)%n14 / ndayyear


  end subroutine getout_daily


  ! subroutine getout_monthly( jpngr, moy )
  !   !////////////////////////////////////////////////////////////////
  !   !  SR called once a year to gather monthly output variables.
  !   !----------------------------------------------------------------
  !   implicit none

  !   ! ARGUMENTS
  !   integer, intent(in) :: jpngr
  !   integer, intent(in) :: moy


  !   return

  ! end subroutine getout_monthly


  subroutine getout_annual( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use _params_modl, only: npft, lu_category
    use _vars_core, only: pleaf, proot

    implicit none

    ! ARGUMENTS
    integer, intent(in) :: jpngr

    outaCveg(:,jpngr)  = pleaf(:,jpngr)%c%c12 + proot(:,jpngr)%c%c12 !+ psapw(:,jpngr)%c%c12 + pwood(:,jpngr)%c%c12

    return

  end subroutine getout_annual


end module _outvars
