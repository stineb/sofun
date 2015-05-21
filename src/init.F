subroutine initglobal
  !////////////////////////////////////////////////////////////////
  !  Initialisation of all _pools on all gridcells at the beginning
  !  of the simulation.
  !  June 2014
  !  b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _classdefs
  use _vars_core, only: pleaf, proot, psapw, pwood, plabl, pexud, &
    plitt_af, plitt_as, plitt_bg, psoil_fs, psoil_sl, dwtot, pninorg, &
    lma, sla, lai_ind, fpc_ind
  use _vegdynamics, only: ispresent, height, fpc_grid, &
    nind, crownarea
  use _params_modl, only: tree
  use _params_core, only: npft
  use _ntransform, only: init_global_ntransform

  implicit none

  ! LOCAL VARIABLES
  integer :: pft

  ! initialise all _pools with zero
  pleaf(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))    ! leaves   , organic pool [gC/m2/ind. and gN/m2/ind.]
  proot(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))    ! roots    , organic pool [gC/m2/ind. and gN/m2/ind.]
  psapw(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))    ! sapwood  , organic pool [gC/m2/ind. and gN/m2/ind.]
  pwood(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))    ! heartwood, organic pool [gC/m2/ind. and gN/m2/ind.]
  plabl(:,:)    = orgpool(carbon(0.0),nitrogen(0.0))    ! labile   , organic pool [gC/m2/ind. and gN/m2/ind.]

  pexud(:,:)    = carbon(0.0)                           ! exudates in soil, carbon pool [gC/m2]

  plitt_af(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! above-ground fine   litter, organic pool [gC/m2 and gN/m2]
  plitt_as(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! above-ground coarse litter, organic pool [gC/m2 and gN/m2]
  plitt_bg(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! below-ground fine   litter, organic pool [gC/m2 and gN/m2]

  psoil_sl(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! fast decomposing soil organic matter, organic pool [gC/m2 and gN/m2]
  psoil_fs(:,:) = orgpool(carbon(0.0),nitrogen(0.0))    ! recalcitrant soil organic matter    , organic pool [gC/m2 and gN/m2]

  ! xxx try
  dwtot(:,:)    = 50.0  ! total soil water content (mm); corresponds to 'dwn' in STASH

  !! xxx try
  !! initialise _pools with non-zero values
  !! here example gridcell 9.5 E, 47.5 N
  ! psoil_fs(:,:) = orgpool( carbon(4000.0), nitrogen(150.0) )
  ! psoil_sl(:,:) = orgpool( carbon(10000.0), nitrogen(550.0) )
  
  ! xxx debug: change back to 0.0
  pninorg(:,:)  = nitrogen(10.0)

  ! Module specific initialisations
  call init_global_ntransform

  ! initialise other properties
  ispresent(:,:) = .false.
  fpc_grid(:,:)  = 0.0
  lai_ind(:,:)   = 0.0
  fpc_ind(:,:)   = 0.0
  height(:,:)    = 0.0

  do pft=1,npft
    if (tree(pft)) then
      nind(pft,:)      = 0.0
      crownarea(pft,:) = 0.0
    else
      nind(pft,:)      = 1.0
      crownarea(pft,:) = 1.0
    end if
  end do

  lma(:,:) = 0.0
  sla(:,:) = 0.0

end subroutine initglobal


subroutine initpft( pft, jpngr )
  !////////////////////////////////////////////////////////////////
  !  Initialisation of specified PFT on specified gridcell
  !  June 2014
  !  b.stocker@imperial.ac.uk
  !----------------------------------------------------------------
  use _classdefs
  use _vars_core, only: pleaf, proot, plabl, psapw, pwood, lai_ind, fpc_ind
  use _params_modl, only: tree, grass
  use _vegdynamics, only: ispresent, height, fpc_grid, nind, &
    crownarea

  implicit none
  integer, intent(in) :: pft
  integer, intent(in) :: jpngr

  ! initialise all _pools with zero
  pleaf(pft,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))
  proot(pft,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))
  plabl(pft,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))
  if (tree(pft)) then
    psapw(pft,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))
    pwood(pft,jpngr) = orgpool(carbon(0.0),nitrogen(0.0))
  endif

  ! initialise other properties
  fpc_grid(pft,jpngr)  = 0.0
  lai_ind(pft,jpngr)   = 0.0
  fpc_ind(pft,jpngr)   = 0.0
  height(pft,jpngr)    = 0.0
  if (tree(pft)) then
    nind(pft,jpngr)      = 0.0
    crownarea(pft,jpngr) = 0.0
  else
    nind(pft,jpngr)      = 1.0
    crownarea(pft,jpngr) = 1.0
  endif

  !if (.not.grass(pft)) then
  !sm_ind(pft,jpngr,1)=orgpool(carbon(0.0),nitrogen(0.0))
  !hm_ind(pft,jpngr,1)=orgpool(carbon(0.0),nitrogen(0.0))
  ! end if

  ! xxx put this somewhere else
  !if(.not.tree(pft)) crownarea(pft,jpngr)=1.0d0
  !leafon(pft,jpngr)=.true.
  !leafondays(pft,jpngr)=0.0d0
  !leafoffdays(pft,jpngr)=0.0d0


end subroutine initpft


subroutine initannual
  !////////////////////////////////////////////////////////////////
  !  Initialises all annually updated variables with zero.
  !----------------------------------------------------------------
  use _classdefs
  use _vars_core, only: aestab

  implicit none

  aestab(:) = orgpool(carbon(0.0),nitrogen(0.0))

end subroutine initannual


subroutine initdaily
  !////////////////////////////////////////////////////////////////
  !  Initialises all daily variables with zero.
  !----------------------------------------------------------------
  use _classdefs
  use _vars_core, only: dnpp, drsoil, drhet, dnep, dcex, dgpp, &
    drauto, drleaf, drroot, drgrow, dnup
  use _vars_core, only: ddoc
  use _nuptake, only: initdaily_nuptake
  use _waterbal, only: initdaily_waterbal
  use _gpp, only: initdaily_gpp

  implicit none

  dnpp(:)        = carbon(0.0)
  drsoil(:)      = carbon(0.0)
  drhet(:)       = carbon(0.0)
  dnep(:)        = carbon(0.0)
  dcex(:)        = 0.0
  dgpp(:)        = 0.0
  drauto(:)      = 0.0
  drleaf(:)      = 0.0
  drroot(:)      = 0.0
  drgrow(:)      = 0.0
  dnup(:)        = nitrogen(0.0)
  ddoc(:)        = 0.0

  ! initialise module-specific output variables
  call initdaily_nuptake()
  call initdaily_waterbal()
  call initdaily_gpp()

  !! xxx try
  !! hold soil carbon and litter _pools at a constant size
  !! here example gridcell 9.5 E, 47.5 N
  !psoil_fs = orgpool( carbon(4000.0), nitrogen(150.0) )
  !psoil_sl = orgpool( carbon(10000.0), nitrogen(550.0) )

  !plitt_af = orgpool( carbon(500.0), nitrogen(10.0))
  !plitt_as = orgpool( carbon(1500.0), nitrogen(30.0))
  !plitt_bg = orgpool( carbon(400.0), nitrogen(8.))

end subroutine initdaily


subroutine initoutput
  !////////////////////////////////////////////////////////////////
  ! Initialises all daily variables with zero.
  ! Called at the beginning of each year by 'biosphere'.
  !----------------------------------------------------------------
  use _outvars
  use _waterbal, only: initoutput_waterbal
  use _nuptake, only: initoutput_nuptake
  use _ntransform, only: initoutput_ntransform
  use _littersom, only: initoutput_littersom
  
  implicit none

  ! daily output variables
  outdgpp(:,:,:)      = 0.0
  outdnpp(:,:,:)      = 0.0
  outdnup(:,:,:)      = 0.0
  outdCleaf(:,:,:)    = 0.0
  outdninorg(:,:,:)   = 0.0

  ! annual output variables
  outagpp(:,:)        = 0.0
  outanpp(:,:)        = 0.0
  outanup(:,:)        = 0.0
  outaCveg(:,:)       = 0.0
  outaCalloc(:,:)     = 0.0
  outaNalloc(:,:)     = 0.0
  outaCveg2lit(:,:)   = 0.0
  outaNveg2lit(:,:)   = 0.0
  outaninorg(:,:)     = 0.0

  ! initialise module-specific output variables
  call initoutput_waterbal
  call initoutput_nuptake
  call initoutput_ntransform  ! nothing yet
  call initoutput_littersom

end subroutine initoutput
