      subroutine writeout_ascii( year )
      !/////////////////////////////////////////////////////////////////////////
      ! Write daily ASCII output
      !-------------------------------------------------------------------------
      use params_siml, only: firstyeartrend, spinupyears
      use params_core, only: ndayyear
      use params_modl, only: lu_category
      use statevars
      use pools
      use outvars

      ! Arguments
      integer, intent(in) :: year       ! simulation year

      ! Local variables
      real :: itime
      integer :: day


      ! xxx implement this: sum over gridcells? single output per gridcell?
      jpngr = 1

      !/////////////////////////////////////////////////////////////////////////
      ! Collect variables to output variables
      !-------------------------------------------------------------------------
      !do lu=1,nlu
      !  do pft=1,npft
      !    if (lu==lu_category(pft)) then
      !    end if
      !  end do
      !end do

      !/////////////////////////////////////////////////////////////////////////
      ! DAILY OUTPUT
      ! Write daily value, summed over all PFTs
      ! xxx implement taking sum over PFTs (and gridcells) in this land use category
      !-------------------------------------------------------------------------
      do day=1,ndayyear

        ! Define 'itime' as a decimal number corresponding to day in the year + year
        itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

        write(101,999) itime, sum(outdgpp(:,day,jpngr))
        write(102,999) itime, sum(outdnpp(:,day,jpngr))
        write(103,999) itime, sum(outdcleaf(:,day,jpngr))
        write(104,999) itime, sum(outdnup(:,day,jpngr))
        write(105,999) itime, sum(outdcex(:,day,jpngr))
        write(106,999) itime, sum(outdnetmin(:,day,jpngr))
        write(107,999) itime, sum(outdninorg(:,day,jpngr))

      end do

      return

999   format (F20.8,F20.8)

      end subroutine writeout_ascii
