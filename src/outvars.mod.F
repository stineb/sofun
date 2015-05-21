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


end module _outvars
