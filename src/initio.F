      subroutine initio
      !////////////////////////////////////////////////////////////////
      ! Opens input/output files.
      !----------------------------------------------------------------
      use params_siml, only: runname

      implicit none

      ! local variables
      character(len=256) :: prefix
      character(len=256) :: filnam

      !//////////////////////////////////////////////////////////////////////////////
      ! OPEN ASCII OUTPUT FILES FOR DAILY OUTPUT
      !------------------------------------------------------------------------------
      prefix = "./output/"//trim(runname)

      filnam=trim(prefix)//'.d.gpp.dat'
      open(101,file=filnam,err=999,status='unknown')

      filnam=trim(prefix)//'.d.npp.dat'
      open(102,file=filnam,err=999,status='unknown')

      filnam=trim(prefix)//'.d.cleaf.dat'
      open(103,file=filnam,err=999,status='unknown')

      return

 999  stop 'INITIO: error opening output files'

      end subroutine initio
