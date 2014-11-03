      module params_siml
      !////////////////////////////////////////////////////////////////
      !  Module contains simulation parameters read in by getpar_siml
      !----------------------------------------------------------------
      use params_core, only: recycle

      implicit none
      integer :: runyears                    ! number of years of entire simulation (spinup+transient)
      integer :: spinupyears                 ! number of spinup years
      integer :: nyeartrend                  ! number of transient years
      integer :: firstyeartrend              ! year AD of first transient year
      logical :: do_spinup                   ! whether this simulation does spinup 
      logical :: init                        ! is true in first simulation year
      logical :: spinup                      ! is true during spinup
      integer :: transyear                   ! year in transient simulation
      character(len=256) :: runname
      character(len=6)   :: sitename
      character(len=256) :: input_dir
      character(len=256) :: co2_forcing_file

      contains

      subroutine getsteering( year, realyear)
      !////////////////////////////////////////////////////////////////
      !  SR defines variables used for steering simulation for each 
      !  simulation year (setting booleans for opening files, doing 
      !  spinup etc.)
      !----------------------------------------------------------------

      implicit none
      integer, intent(in) :: year
      integer, intent(out) :: realyear

      if (do_spinup) then
        realyear = max( firstyeartrend, year - spinupyears + firstyeartrend - 1 )
        if (year<=spinupyears) then
          spinup = .true.
        else  
          spinup = .false.
        endif
      else
        realyear = year + firstyeartrend
      endif
      
      if (year==1) then
        init = .true.
      else
        init = .false.
      endif 

      if (spinup) then
        transyear = mod(year,recycle)
      else
        transyear = realyear-firstyeartrend+1
      endif

      end subroutine getsteering
      
      subroutine getpar_siml( runname )
      !////////////////////////////////////////////////////////////////
      !  SR for reading and defining simulation parameters from file 
      !  <runname>.sofun.parameter. Only once at start of simulation.
      !----------------------------------------------------------------

      implicit none
      character(len=*), intent(in) :: runname
      
      ! Read in main model parameters
      write(0,*) 'reading parameter file ', runname//".sofun.parameter ..."

      call getstring(runname//".sofun.parameter",'sitename',sitename)
      call getlogical(runname//".sofun.parameter",'spinup',do_spinup)
      call getint(runname//".sofun.parameter",'spinupyears',spinupyears)
      call getint(runname//".sofun.parameter",'firstyeartrend',
     $     firstyeartrend)
      call getint(runname//".sofun.parameter",'nyeartrend',
     $     nyeartrend)
      
      call getstring(runname//".sofun.parameter",'input_dir',input_dir)
      call getstring(runname//".sofun.parameter",'co2_forcing_file',co2_forcing_file)

      if (do_spinup) then
        runyears = nyeartrend + spinupyears
      else
        runyears = nyeartrend
        spinupyears = 0
      endif
      
      end subroutine getpar_siml

      
      subroutine getreal( filename, paraname, paravalue )
      !////////////////////////////////////////////////////////////////
      !  "Low-level" SR for reading parameter values from text file
      !----------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: filename, paraname
      real,intent(out) :: paravalue

      integer :: filehandle
      character(len=40) :: readname, readvalue
      
      filehandle = 111
      open(filehandle,status='old',err=19,file=filename)
 9    read(filehandle,12,end=10)readname,readvalue
      if (trim(readname).eq.paraname) then
        read(readvalue,*) paravalue
        goto 11
      else
        goto 9
      endif
 10   continue
      write(0,*) 'GETPAR: '//paraname//' of type real not found'
      stop
      
 11   continue
 12   format(2a40)

      close(filehandle)

      return

 19   write(0,*) 'GETPAR: '//filename//' not found!'
      stop
      
      end

      subroutine getint( filename, paraname, paravalue )
      !////////////////////////////////////////////////////////////////
      ! "Low-level" SR for reading parameter values from text file
      !----------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: filename, paraname      
      integer, intent(out) :: paravalue

      integer :: filehandle
      character(len=40) :: readname, readvalue
      
      filehandle = 111
      open(filehandle,status='old',err=19,file=filename)
 9    read(filehandle,12,end=10)readname,readvalue
      if (trim(readname).eq.paraname) then
        read(readvalue,*) paravalue
        goto 11
      else
        goto 9
      endif
 10   continue
      write(0,*) 'GETPAR: in file '//filename//':'
      write(0,*) 'GETPAR: '//paraname//' of type integer not found'
      stop
      
 11   continue
 12   format(2a40)
      
      close(filehandle)
      
      return
      
 19   write(0,*) 'GETPAR: file '//filename//' not found:'
      write(0,*) filename
      stop
      
      end subroutine getint
      
      
      subroutine getlogical( filename, paraname, paravalue )
      !////////////////////////////////////////////////////////////////
      ! "Low-level" SR for reading parameter values from text file
      !----------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: filename, paraname
      logical, intent(out) :: paravalue

      integer :: filehandle
      character(len=40) :: readname, readvalue

      filehandle = 111
      open(filehandle,status='old',err=19,file=filename)
 9    read(filehandle,12,end=10)readname,readvalue
      if (trim(readname).eq.paraname) then
        read(readvalue,*) paravalue
        goto 11
      else
        goto 9
      endif
 10   continue
      write(0,*) 'GETPAR: '//paraname//' of type logical not found'
      stop
      
 11   continue
 12   format(2a40)

      close(filehandle)

      return

 19   write(0,*) 'GETPAR: '//filename//' not found!'
      stop
      
      end

      subroutine getstring( filename, paraname, paravalue )
      !////////////////////////////////////////////////////////////////
      ! "Low-level" SR for reading parameter values from text file
      !----------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: filename, paraname
      character(len=*), intent(out) :: paravalue

      integer :: filehandle,i
      character(len=40) :: readname
      character(len=1024) :: readvalue

      filehandle = 111
      open(filehandle,status='old',err=19,file=filename)
 9    read(filehandle,12,end=10)readname,readvalue
      if (trim(readname).eq.paraname) then
        ! Strip leading and trailing whitespace from readvalue
        i=1
        do while (readvalue(i:i).eq.' ' .and. i.lt.len(readvalue))
          i=i+1
        enddo
        paravalue = trim(readvalue(i:))
        goto 11
      else
        goto 9
      endif
 10   continue
      write(0,*) 'GETPAR: '//paraname//' of type string not found'
      stop
      
 11   continue
 12   format(a40,a)

      close(filehandle)

      return

 19   write(0,*) 'GETPAR: '//filename//' not found!'
      stop
      
      end
      

      subroutine getchar( filename, paraname, paravalue )
      !////////////////////////////////////////////////////////////////
      ! "Low-level" SR for reading parameter values from text file
      !----------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: filename, paraname
      character(len=1), intent(out) :: paravalue

      integer :: filehandle
      character(len=40) :: readname, readvalue

      filehandle = 111
      open(filehandle,status='old',err=19,file=filename)
 9    read(filehandle,12,end=10)readname,readvalue
      if (trim(readname).eq.paraname) then
        read(readvalue,*) paravalue
        goto 11
      else
        goto 9
      endif
 10   continue
      write(0,*) 'GETPAR: '//paraname//' of type char not found'
      stop
      
 11   continue
 12   format(2a40)

      close(filehandle)

      return

 19   write(0,*) 'GETPAR: '//filename//' not found!'
      stop

      end

c$$$  subroutine setreal(filename,paraname,default,paravalue)
c$$$          !////////////////////////////////////////////////////////////////
c$$$          ! "Low-level" SR for reading/defining parameter values from text file
c$$$          !----------------------------------------------------------------
c$$$
c$$$          implicit none
c$$$
c$$$          integer filehandle
c$$$          character(len=40) :: readname, readvalue
c$$$          character(len=*), intent(in) :: filename, paraname
c$$$
c$$$          real,intent(out) :: paravalue
c$$$          real,intent(in)  :: default
c$$$
c$$$          filehandle = 111
c$$$          open(filehandle,status='old',err=19,file=filename)
c$$$ 9        read(filehandle,12,end=10)readname,readvalue
c$$$          if (trim(readname).eq.paraname) then
c$$$             read(readvalue,*) paravalue
c$$$             goto 11
c$$$          else
c$$$             goto 9
c$$$          endif
c$$$ 10       continue
c$$$          paravalue = default
c$$$
c$$$ 11       continue
c$$$ 12       format(2a40)
c$$$
c$$$          close(filehandle)
c$$$
c$$$          return
c$$$
c$$$ 19       stop 'GETPAR: '//filename//' not found!'
c$$$
c$$$          end
c$$$
c$$$          subroutine setint(filename,paraname,default,paravalue)
c$$$          !////////////////////////////////////////////////////////////////
c$$$          ! "Low-level" SR for reading/defining parameter values from text file
c$$$          !----------------------------------------------------------------
c$$$
c$$$          implicit none
c$$$
c$$$          integer filehandle
c$$$          character(len=40) :: readname, readvalue
c$$$          character(len=*), intent(in) :: filename, paraname
c$$$
c$$$          integer, intent(out) :: paravalue
c$$$          integer :: default
c$$$
c$$$          filehandle = 111
c$$$          open(filehandle,status='old',err=19,file=filename)
c$$$ 9        read(filehandle,12,end=10)readname,readvalue
c$$$          if (trim(readname).eq.paraname) then
c$$$             read(readvalue,*) paravalue
c$$$             goto 11
c$$$          else
c$$$             goto 9
c$$$          endif
c$$$ 10       continue
c$$$          paravalue = default
c$$$
c$$$ 11       continue
c$$$ 12       format(2a40)
c$$$
c$$$          close(filehandle)
c$$$
c$$$          return
c$$$
c$$$ 19       stop 'GETPAR: '//filename//' not found!'
c$$$
c$$$          end
c$$$
c$$$          subroutine setlogical(filename,paraname,default,paravalue)
c$$$          !////////////////////////////////////////////////////////////////
c$$$          ! "Low-level" SR for reading/defining parameter values from text file
c$$$          !----------------------------------------------------------------
c$$$
c$$$          implicit none
c$$$
c$$$          character(len=*), intent(in) :: filename, paraname
c$$$          logical, intent(out) :: paravalue
c$$$
c$$$          integer filehandle
c$$$          character(len=40) :: readname, readvalue
c$$$          logical :: default
c$$$
c$$$          filehandle = 111
c$$$          open(filehandle,status='old',err=19,file=filename)
c$$$ 9        read(filehandle,12,end=10)readname,readvalue
c$$$          if (trim(readname).eq.paraname) then
c$$$             read(readvalue,*) paravalue
c$$$             goto 11
c$$$          else
c$$$             goto 9
c$$$          endif
c$$$ 10       continue
c$$$          paravalue = default
c$$$
c$$$ 11       continue
c$$$ 12       format(2a40)
c$$$
c$$$          close(filehandle)
c$$$
c$$$          return
c$$$
c$$$ 19       stop 'GETPAR: '//filename//' not found!'
c$$$
c$$$          end

c$$$
c$$$          subroutine setstring(filename,paraname,default,paravalue)
c$$$          !////////////////////////////////////////////////////////////////
c$$$          ! "Low-level" SR for reading/defining parameter values from text file
c$$$          !----------------------------------------------------------------
c$$$
c$$$          implicit none
c$$$
c$$$          ! Arguments
c$$$          character(len=*), intent(in) :: filename, paraname
c$$$          character(len=*), intent(out) :: paravalue, default
c$$$
c$$$          ! Local variables
c$$$          integer filehandle,i
c$$$          character readname*40,readvalue*1024
c$$$
c$$$          filehandle = 111
c$$$          open(filehandle,status='old',err=19,file=filename)
c$$$ 9        read(filehandle,12,end=10)readname,readvalue
c$$$          if (trim(readname).eq.paraname) then
c$$$             ! Strip leading and trailing whitespace from readvalue
c$$$             i=1
c$$$             do while (readvalue(i:i).eq.' ' .and. i.lt.len(readvalue))
c$$$                i=i+1
c$$$             enddo
c$$$             paravalue = trim(readvalue(i:))
c$$$             goto 11
c$$$          else
c$$$             goto 9
c$$$          endif
c$$$ 10       continue
c$$$          paravalue = default
c$$$
c$$$ 11       continue
c$$$
c$$$ 12       format(a40,a)
c$$$
c$$$          close(filehandle)
c$$$
c$$$          return
c$$$
c$$$ 19       stop 'GETPAR: '//filename//' not found!'
c$$$
c$$$          end

       end module params_siml

