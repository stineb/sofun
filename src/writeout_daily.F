      subroutine writeout_daily(year,day)
      !/////////////////////////////////////////////////////////////////////////
      ! Write daily ASCII output
      !-------------------------------------------------------------------------
      use params_siml, only: firstyeartrend, spinupyears
      use params_core, only: ndayyear
      use params_modl, only: lu_category
      use statevars
      use pools

      ! Arguments
      integer, intent(in) :: year
      integer, intent(in) :: day

      ! Local variables
      real :: itime
      integer :: lu, pft

      !/////////////////////////////////////////////////////////////////////////
      ! Collect variables to output variables
      !-------------------------------------------------------------------------
      do lu=1,nlu
        dreco(lu) = cplus( dreco(lu), drhet(lu) )
        do pft=1,npft
          if (lu==lu_category(pft)) then
            dreco(lu) = cplus( dreco(lu), drsoil(pft), carbon(drauto(pft)) )
          end if
        end do
      end do

      !/////////////////////////////////////////////////////////////////////////
      ! Write daily value, summed over all PFTs
      ! xxx implement taking sum over PFTs (and gridcells) in this land use category
      !-------------------------------------------------------------------------
      ! Define 'itime' as a decimal number corresponding to day in the year + year
      itime = real(year) + real(firstyeartrend) - real(spinupyears) + real(day-1)/real(ndayyear)

      write(101,999) itime, dgpp
      write(102,999) itime, dnpp%c12
      write(103,999) itime, pleaf%c%c12

      return

999   format (F20.8,F20.8)

      end subroutine writeout_daily
