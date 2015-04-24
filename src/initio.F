subroutine initio
  !////////////////////////////////////////////////////////////////
  ! Opens input/output files.
  !----------------------------------------------------------------
  use params_siml, only: runname
  use waterbal_stash, only: initio_waterbal
  use nuptake_fun, only: initio_nuptake
  use ntransform_xuri, only: initio_ntransform
  use littersom_lpj, only: initio_littersom

  implicit none

  ! local variables
  character(len=256) :: prefix
  character(len=256) :: filnam


  !////////////////////////////////////////////////////////////////
  ! DAILY OUTPUT: OPEN ASCII OUTPUT FILES 
  !----------------------------------------------------------------
  prefix = "./output/"//trim(runname)

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


  !----------------------------------------------------------------
  ! MODULE-SPECIFIC OUTPUT FILES
  !----------------------------------------------------------------
  call initio_waterbal( prefix )
  call initio_nuptake( prefix )
  call initio_ntransform( prefix ) ! nothing yet xxx
  call initio_littersom( prefix )


  return


  999  stop 'INITIO: error opening output files'

end subroutine initio
