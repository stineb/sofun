subroutine initio
  !////////////////////////////////////////////////////////////////
  ! Opens input/output files.
  !----------------------------------------------------------------
  use params_siml, only: runname
  use waterbal_stash, only: initio_waterbal
  use nuptake_fun, only: initio_nuptake

  implicit none

  ! local variables
  character(len=256) :: prefix
  character(len=256) :: filnam


  !////////////////////////////////////////////////////////////////
  ! OPEN ASCII OUTPUT FILES FOR DAILY OUTPUT
  !----------------------------------------------------------------
  prefix = "./output/"//trim(runname)

  ! GPP
  filnam=trim(prefix)//'.d.gpp.dat'
  open(101,file=filnam,err=999,status='unknown')

  ! NPP
  filnam=trim(prefix)//'.d.npp.dat'
  open(102,file=filnam,err=999,status='unknown')

  ! LEAF C
  filnam=trim(prefix)//'.d.cleaf.dat'
  open(103,file=filnam,err=999,status='unknown')

  ! N UPTAKE
  filnam=trim(prefix)//'.d.nup.dat'
  open(104,file=filnam,err=999,status='unknown')

  ! C EXUDATION
  filnam=trim(prefix)//'.d.cex.dat'
  open(105,file=filnam,err=999,status='unknown')

  ! NET N MINERALISATION
  filnam=trim(prefix)//'.d.netmin.dat'
  open(106,file=filnam,err=999,status='unknown')

  ! INORGANIC N (NO3+NH4)
  filnam=trim(prefix)//'.d.ninorg.dat'
  open(107,file=filnam,err=999,status='unknown')

  ! BIOLOGICAL NITROGEN FIXATION OF FREE-LIVING  ORGANISMS
  filnam=trim(prefix)//'.d.nfixfree.dat'
  open(108,file=filnam,err=999,status='unknown')


  !////////////////////////////////////////////////////////////////
  ! OPEN ASCII OUTPUT FILES FOR ANNUAL OUTPUT
  !----------------------------------------------------------------

  ! SOIL C
  filnam=trim(prefix)//'.a.csoil.dat'
  open(301,file=filnam,err=999,status='unknown')

  ! LITTER C
  filnam=trim(prefix)//'.a.clitt.dat'
  open(302,file=filnam,err=999,status='unknown')

  ! C ALLOCATED TO GROWTH 
  filnam=trim(prefix)//'.a.calloc.dat'
  open(303,file=filnam,err=999,status='unknown')


  !----------------------------------------------------------------
  ! MODULE-SPECIFIC OUTPUT FILES
  !----------------------------------------------------------------
  call initio_waterbal( prefix )
  call initio_nuptake( prefix )


  return


  999  stop 'INITIO: error opening output files'

end subroutine initio
