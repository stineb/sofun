module _outvars
!////////////////////////////////////////////////////////////////
!  Module contains all daily variables. These variables have no 
!  spatial dimension. At the end of the daily loop, values may be 
!  copied to spatial arrays.
! Copyright (C) 2015  Benjamin David Stocker
! contact: b.stocker@imperial.ac.uk
!----------------------------------------------------------------
  use _params_core
  use _classdefs

  implicit none

  !----------------------------------------------------------------
  ! Daily output variables
  !----------------------------------------------------------------
  real, dimension(npft,ndayyear,maxgrid) :: outdgpp
  real, dimension(npft,ndayyear,maxgrid) :: outdtransp
  real, dimension(npft,ndayyear,maxgrid) :: outdnpp
  real, dimension(npft,ndayyear,maxgrid) :: outdnup
  real, dimension(npft,ndayyear,maxgrid) :: outdCleaf
  real, dimension(npft,ndayyear,maxgrid) :: outdCroot
  real, dimension(npft,ndayyear,maxgrid) :: outdClabl
  real, dimension(npft,ndayyear,maxgrid) :: outdNlabl
  real, dimension(npft,ndayyear,maxgrid) :: outdClitt
  real, dimension(npft,ndayyear,maxgrid) :: outdNlitt
  real, dimension(nlu,ndayyear,maxgrid)  :: outdCsoil
  real, dimension(nlu,ndayyear,maxgrid)  :: outdNsoil

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
  real, dimension(npft,maxgrid) :: outaCveg2lit
  real, dimension(npft,maxgrid) :: outaNveg2lit
  real, dimension(nlu, maxgrid) :: outaNinorg

contains

  subroutine initoutput()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------

    ! daily output variables
    outdgpp(:,:,:)      = 0.0
    outdtransp(:,:,:)   = 0.0
    outdnpp(:,:,:)      = 0.0
    outdnup(:,:,:)      = 0.0
    outdCleaf(:,:,:)    = 0.0
    outdCroot(:,:,:)    = 0.0
    outdClabl(:,:,:)    = 0.0
    outdNlabl(:,:,:)    = 0.0
    outdClitt(:,:,:)    = 0.0
    outdCsoil(:,:,:)    = 0.0
    outdNlitt(:,:,:)    = 0.0
    outdNsoil(:,:,:)    = 0.0
    outdninorg(:,:,:)   = 0.0

    ! annual output variables
    outagpp(:,:)        = 0.0
    outanpp(:,:)        = 0.0
    outanup(:,:)        = 0.0
    outaCveg(:,:)       = 0.0
    outaCveg2lit(:,:)   = 0.0
    outaNveg2lit(:,:)   = 0.0
    outaninorg(:,:)     = 0.0


  end subroutine initoutput


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

    ! ROOT C
    filnam=trim(prefix)//'.d.croot.out'
    open(111,file=filnam,err=999,status='unknown')

    ! LABILE C
    filnam=trim(prefix)//'.d.clabl.out'
    open(112,file=filnam,err=999,status='unknown')

    ! LITTER C
    filnam=trim(prefix)//'.d.clitt.out'
    open(113,file=filnam,err=999,status='unknown')

    ! TRANSPIRATION
    filnam=trim(prefix)//'.d.transp.out'
    open(114,file=filnam,err=999,status='unknown')

    ! LABILE N
    filnam=trim(prefix)//'.d.nlabl.out'
    open(115,file=filnam,err=999,status='unknown')

    ! SOIL C
    filnam=trim(prefix)//'.d.csoil.out'
    open(118,file=filnam,err=999,status='unknown')

    ! LITTER N
    filnam=trim(prefix)//'.d.nlitt.out'
    open(119,file=filnam,err=999,status='unknown')

    ! SOIL N
    filnam=trim(prefix)//'.d.nsoil.out'
    open(120,file=filnam,err=999,status='unknown')

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------

    ! C VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.cveg2lit.out'
    open(307,file=filnam,err=999,status='unknown')

    ! N VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.nveg2lit.out'
    open(308,file=filnam,err=999,status='unknown')

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
    use _vars_core, only: dgpp, dnpp, dnup, pleaf, pninorg, proot, &
      plabl, plitt_bg, plitt_af, plitt_as, dtemp_soil, dtransp, &
      psoil_sl, psoil_fs

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
    outdgpp(:,doy,jpngr)        = dgpp(:)
    outdtransp(:,doy,jpngr)     = dtransp(:)
    outdnpp(:,doy,jpngr)        = dnpp(:)%c12
    outdnup(:,doy,jpngr)        = dnup(:)%n14
    outdCleaf(:,doy,jpngr)      = pleaf(:,jpngr)%c%c12
    outdCroot(:,doy,jpngr)      = proot(:,jpngr)%c%c12
    outdClabl(:,doy,jpngr)      = plabl(:,jpngr)%c%c12
    outdClabl(:,doy,jpngr)      = plabl(:,jpngr)%n%n14
    outdClitt(:,doy,jpngr)      = plitt_af(:,jpngr)%c%c12 + plitt_as(:,jpngr)%c%c12 + plitt_bg(:,jpngr)%c%c12
    outdNlitt(:,doy,jpngr)      = plitt_af(:,jpngr)%n%n14 + plitt_as(:,jpngr)%n%n14 + plitt_bg(:,jpngr)%n%n14
    outdCsoil(:,doy,jpngr)      = psoil_sl(:,jpngr)%c%c12 + psoil_fs(:,jpngr)%c%c12
    outdNsoil(:,doy,jpngr)      = psoil_sl(:,jpngr)%n%n14 + psoil_fs(:,jpngr)%n%n14
    outdninorg(:,doy,jpngr)     = pninorg(:,jpngr)%n14
    outdtemp_soil(:,doy,jpngr)  = dtemp_soil(:,jpngr)

    ! write(0,*) 'outdCsoil', outdCsoil(:,doy,jpngr)
    ! write(0,*) 'psoil_sl(:,jpngr)%c%c12', psoil_sl(:,jpngr)%c%c12
    ! write(0,*) 'psoil_fs(:,jpngr)%c%c12', psoil_fs(:,jpngr)%c%c12
    ! if (doy==ndayyear) stop

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

    ! arguments
    integer, intent(in) :: jpngr

    outaCveg(:,jpngr)  = pleaf(:,jpngr)%c%c12 + proot(:,jpngr)%c%c12 !+ psapw(:,jpngr)%c%c12 + pwood(:,jpngr)%c%c12

  end subroutine getout_annual


end module _outvars
