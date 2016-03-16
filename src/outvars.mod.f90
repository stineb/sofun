module _outvars
!////////////////////////////////////////////////////////////////
!  Module contains all daily variables. These variables have no 
!  spatial dimension. At the end of the daily loop, values may be 
!  copied to spatial arrays.
! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
! contact: b.stocker@imperial.ac.uk
!----------------------------------------------------------------
  use _params_core
  use _classdefs
  use _params_siml, only: & 
      loutdnpp       &
    , loutdnup       &
    , loutdCleaf     &
    , loutdCroot     &
    , loutdClabl     &
    , loutdNlabl     &
    , loutdClitt     &
    , loutdNlitt     &
    , loutdCsoil     &
    , loutdNsoil     &
    , loutdlai       &
    , loutdninorg    &
    , loutdtemp_soil

  implicit none

  !----------------------------------------------------------------
  ! Daily output variables
  !----------------------------------------------------------------
  real, allocatable, dimension(:,:,:) :: outdnpp
  real, allocatable, dimension(:,:,:) :: outdnup
  real, allocatable, dimension(:,:,:) :: outdCleaf
  real, allocatable, dimension(:,:,:) :: outdCroot
  real, allocatable, dimension(:,:,:) :: outdClabl
  real, allocatable, dimension(:,:,:) :: outdNlabl
  real, allocatable, dimension(:,:,:) :: outdClitt
  real, allocatable, dimension(:,:,:) :: outdNlitt
  real, allocatable, dimension(:,:,:) :: outdCsoil
  real, allocatable, dimension(:,:,:) :: outdNsoil
  real, allocatable, dimension(:,:,:) :: outdlai
  real, allocatable, dimension(:,:,:) :: outdninorg

  ! These are stored as dayly variables for annual output
  ! at day of year when LAI is at its maximum.
  real, dimension(npft,ndayyear) :: dnarea_mb
  real, dimension(npft,ndayyear) :: dnarea_cw
  real, dimension(npft,ndayyear) :: dlma
  real, dimension(npft,ndayyear) :: dcton_lm

  !----------------------------------------------------------------
  ! Annual output variables
  !----------------------------------------------------------------
  real, dimension(npft,maxgrid) :: outanpp
  real, dimension(npft,maxgrid) :: outanup
  real, dimension(npft,maxgrid) :: outaCveg
  real, dimension(npft,maxgrid) :: outaCveg2lit
  real, dimension(npft,maxgrid) :: outaNveg2lit
  real, dimension(nlu, maxgrid) :: outaNinorg
  real, dimension(npft,maxgrid) :: outanarea_mb
  real, dimension(npft,maxgrid) :: outanarea_cw
  real, dimension(npft,maxgrid) :: outalai
  real, dimension(npft,maxgrid) :: outalma
  real, dimension(npft,maxgrid) :: outacton_lm

contains

  subroutine initoutput()
    !////////////////////////////////////////////////////////////////
    ! Initialises all daily variables with zero.
    ! Called at the beginning of each year by 'biosphere'.
    !----------------------------------------------------------------
    if (loutdnpp      ) allocate( outdnpp      (npft,ndayyear,maxgrid) )
    if (loutdnup      ) allocate( outdnup      (npft,ndayyear,maxgrid) )
    if (loutdCleaf    ) allocate( outdCleaf    (npft,ndayyear,maxgrid) )
    if (loutdCroot    ) allocate( outdCroot    (npft,ndayyear,maxgrid) )
    if (loutdClabl    ) allocate( outdClabl    (npft,ndayyear,maxgrid) )
    if (loutdNlabl    ) allocate( outdNlabl    (npft,ndayyear,maxgrid) )
    if (loutdClitt    ) allocate( outdClitt    (npft,ndayyear,maxgrid) )
    if (loutdNlitt    ) allocate( outdNlitt    (npft,ndayyear,maxgrid) )
    if (loutdCsoil    ) allocate( outdCsoil    (nlu,ndayyear,maxgrid)  )
    if (loutdNsoil    ) allocate( outdNsoil    (nlu,ndayyear,maxgrid)  )
    if (loutdlai      ) allocate( outdlai      (npft,ndayyear,maxgrid) )
    if (loutdninorg   ) allocate( outdninorg   (nlu,ndayyear,maxgrid)  )

    ! annual output variables
    outanpp(:,:)        = 0.0
    outanup(:,:)        = 0.0
    outaCveg(:,:)       = 0.0
    outaCveg2lit(:,:)   = 0.0
    outaNveg2lit(:,:)   = 0.0
    outaninorg(:,:)     = 0.0
    outanarea_mb(:,:)   = 0.0
    outanarea_cw(:,:)   = 0.0
    outalai     (:,:)   = 0.0
    outalma     (:,:)   = 0.0
    outacton_lm (:,:)   = 0.0

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
    ! NPP
    if (loutdnpp) then 
      filnam=trim(prefix)//'.d.npp.out'
      open(102,file=filnam,err=999,status='unknown')
    end if

    ! LEAF C
    if (loutdCleaf    ) then
      filnam=trim(prefix)//'.d.cleaf.out'
      open(103,file=filnam,err=999,status='unknown')
    end if 

    ! N UPTAKE
    if (loutdnup      ) then
      filnam=trim(prefix)//'.d.nup.out'
      open(104,file=filnam,err=999,status='unknown')
    end if

    ! INORGANIC N (NO3+NH4)
    if (loutdninorg   ) then
      filnam=trim(prefix)//'.d.ninorg.out'
      open(107,file=filnam,err=999,status='unknown')
    end if

    ! AIR TEMPERATURE
    filnam=trim(prefix)//'.d.temp.out'
    open(110,file=filnam,err=999,status='unknown')

    ! ROOT C
    if (loutdCroot    ) then
      filnam=trim(prefix)//'.d.croot.out'
      open(111,file=filnam,err=999,status='unknown')
    end if

    ! LABILE C
    if (loutdClabl    ) then
      filnam=trim(prefix)//'.d.clabl.out'
      open(112,file=filnam,err=999,status='unknown')
    end if

    ! LITTER C
    if (loutdClitt    ) then
      filnam=trim(prefix)//'.d.clitt.out'
      open(113,file=filnam,err=999,status='unknown')
    end if

    ! LABILE N
    if (loutdNlabl    ) then
      filnam=trim(prefix)//'.d.nlabl.out'
      open(115,file=filnam,err=999,status='unknown')
    end if

    ! SOIL C
    if (loutdCsoil    ) then
      filnam=trim(prefix)//'.d.csoil.out'
      open(118,file=filnam,err=999,status='unknown')
    end if

    ! LITTER N
    if (loutdNlitt    ) then
      filnam=trim(prefix)//'.d.nlitt.out'
      open(119,file=filnam,err=999,status='unknown')
    end if

    ! SOIL N
    if (loutdNsoil    ) then
      filnam=trim(prefix)//'.d.nsoil.out'
      open(120,file=filnam,err=999,status='unknown')
    end if

    ! LAI
    if (loutdlai      ) then
      filnam=trim(prefix)//'.d.lai.out'
      open(121,file=filnam,err=999,status='unknown')
    end if

    !////////////////////////////////////////////////////////////////
    ! ANNUAL OUTPUT: OPEN ASCII OUTPUT FILES
    !----------------------------------------------------------------

    ! C VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.cveg2lit.out'
    open(307,file=filnam,err=999,status='unknown')

    ! N VEGETATION -> LITTER TRANSFER
    filnam=trim(prefix)//'.a.nveg2lit.out'
    open(308,file=filnam,err=999,status='unknown')

    ! NPP 
    filnam=trim(prefix)//'.a.npp.out'
    open(311,file=filnam,err=999,status='unknown')

    ! VEG C
    filnam=trim(prefix)//'.a.cveg.out'
    open(312,file=filnam,err=999,status='unknown')

    ! INORGANIC N (mean over days)
    filnam=trim(prefix)//'.a.ninorg.out'
    open(316,file=filnam,err=999,status='unknown')

    ! N UPTAKE
    filnam=trim(prefix)//'.a.nup.out'
    open(317,file=filnam,err=999,status='unknown')

    ! LAI (ANNUAL MAXIMUM)
    filnam=trim(prefix)//'.a.lai.out'
    open(318,file=filnam,err=999,status='unknown')

    ! METABOLIC NAREA (AT ANNUAL LAI MAXIMUM)
    filnam=trim(prefix)//'.a.narea_mb.out'
    open(319,file=filnam,err=999,status='unknown')

    ! CELL WALL NAREA (AT ANNUAL LAI MAXIMUM)
    filnam=trim(prefix)//'.a.narea_cw.out'
    open(320,file=filnam,err=999,status='unknown')

    ! LEAF C:N RATIO (AT ANNUAL LAI MAXIMUM)
    filnam=trim(prefix)//'.a.cton_lm.out'
    open(321,file=filnam,err=999,status='unknown')

    ! LMA (AT ANNUAL LAI MAXIMUM)
    filnam=trim(prefix)//'.a.lma.out'
    open(322,file=filnam,err=999,status='unknown')

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
    use _params_core, only: ndayyear, npft
    use _vars_core, only: dnpp, dnup, pleaf, pninorg, proot, plabl, &
      plitt_bg, plitt_af, plitt_as, psoil_sl, psoil_fs, lai_ind, &
      narea_metabolic, narea_structural, r_cton_leaf, lma

    ! arguments
    integer, intent(in) :: jpngr
    integer, intent(in) :: moy
    integer, intent(in) :: doy

    ! LOCAL VARIABLES
    integer :: pft

    !----------------------------------------------------------------
    ! DAILY
    ! Collect daily output variables
    ! so far not implemented for isotopes
    !----------------------------------------------------------------
    if (loutdnpp      ) outdnpp(:,doy,jpngr)       = dnpp(:)%c12
    if (loutdnup      ) outdnup(:,doy,jpngr)       = dnup(:)%n14
    if (loutdCleaf    ) outdCleaf(:,doy,jpngr)     = pleaf(:,jpngr)%c%c12
    if (loutdCroot    ) outdCroot(:,doy,jpngr)     = proot(:,jpngr)%c%c12
    if (loutdClabl    ) outdClabl(:,doy,jpngr)     = plabl(:,jpngr)%c%c12

    if (loutdNlabl    ) outdNlabl(:,doy,jpngr)     = plabl(:,jpngr)%n%n14
    if (loutdClitt    ) outdClitt(:,doy,jpngr)     = plitt_af(:,jpngr)%c%c12 + plitt_as(:,jpngr)%c%c12 + plitt_bg(:,jpngr)%c%c12
    if (loutdNlitt    ) outdNlitt(:,doy,jpngr)     = plitt_af(:,jpngr)%n%n14 + plitt_as(:,jpngr)%n%n14 + plitt_bg(:,jpngr)%n%n14
    if (loutdCsoil    ) outdCsoil(:,doy,jpngr)     = psoil_sl(:,jpngr)%c%c12 + psoil_fs(:,jpngr)%c%c12
    if (loutdNsoil    ) outdNsoil(:,doy,jpngr)     = psoil_sl(:,jpngr)%n%n14 + psoil_fs(:,jpngr)%n%n14
    if (loutdninorg   ) outdninorg(:,doy,jpngr)    = pninorg(:,jpngr)%n14
    if (loutdlai      ) outdlai(:,doy,jpngr)       = lai_ind(:,jpngr)
    
    dnarea_mb(:,doy)           = narea_metabolic(:)  
    dnarea_cw(:,doy)           = narea_structural(:)
    dcton_lm(:,doy)            = r_cton_leaf(:,jpngr)
    dlma(:,doy)                = lma(:,jpngr)

    ! write(0,*) 'psoil_sl(:,jpngr)%c%c12', psoil_sl(:,jpngr)%c%c12
    ! write(0,*) 'psoil_fs(:,jpngr)%c%c12', psoil_fs(:,jpngr)%c%c12
    ! if (doy==ndayyear) stop

    !----------------------------------------------------------------
    ! ANNUAL SUM OVER DAILY VALUES
    ! Collect annual output variables as sum of daily values
    !----------------------------------------------------------------
    outanpp(:,jpngr)    = outanpp(:,jpngr) + dnpp(:)%c12
    outanup(:,jpngr)    = outanup(:,jpngr) + dnup(:)%n14
    outaninorg(:,jpngr) = outaNinorg(:,jpngr) + pninorg(:,jpngr)%n14 / ndayyear

  end subroutine getout_daily


  subroutine getout_annual( jpngr )
    !////////////////////////////////////////////////////////////////
    !  SR called once a year to gather annual output variables.
    !----------------------------------------------------------------
    use _params_core, only: ndayyear, npft
    use _vars_core, only: pleaf, proot

    ! arguments
    integer, intent(in) :: jpngr

    ! local variables
    integer                  :: pft
    integer                  :: doy
    integer, dimension(npft) :: maxdoy
    ! intrinsic                :: maxloc

    ! Output annual value at day of peak LAI
    do pft=1,npft
      maxdoy(pft)=1
      do doy=2,ndayyear
        if ( outdlai(pft,doy,jpngr) > outdlai(pft,(doy-1),jpngr) ) then
          maxdoy(pft) = doy
        end if
      end do

      outaCveg    (pft,jpngr) = outdCleaf(pft,maxdoy(pft),jpngr) + outdCroot(pft,maxdoy(pft),jpngr) 
      outanarea_mb(pft,jpngr) = dnarea_mb(pft,maxdoy(pft))
      outanarea_cw(pft,jpngr) = dnarea_cw(pft,maxdoy(pft))
      outalai     (pft,jpngr) = outdlai(pft,maxdoy(pft),jpngr)
      outalma     (pft,jpngr) = dlma(pft,maxdoy(pft))
      outacton_lm (pft,jpngr) = dcton_lm(pft,maxdoy(pft))
    end do

    open(unit=666,file='cton_vs_lai.log',status='unknown')
    do doy=1,ndayyear
      write(666,*) outdlai(1,doy,1), ",", dcton_lm(1,doy)
    end do
    close(unit=666)


  end subroutine getout_annual


end module _outvars
