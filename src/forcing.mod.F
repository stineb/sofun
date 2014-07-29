       module forcing
       !////////////////////////////////////////////////////////////////
       !  Module contains forcing variables (climate, co2, ...), and
       !  subroutines used to read forcing input files
       !----------------------------------------------------------------
         use params_core
         implicit none
         !real, allocatable :: in_co2(:)

         contains

         subroutine getco2( runname, realyear, pco2 )
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's atmospheric CO2 from input
         !----------------------------------------------------------------
         use params_core
         use params_siml

         implicit none
         character(len=*), intent(in) :: runname
         integer, intent(in) :: realyear
         real, intent(out) :: pco2

         pco2 = getvalreal( trim(input_dir)//trim(co2_forcing_file), realyear )

         end subroutine getco2

         subroutine getclimate( runname, realyear, dtemp, dprec )
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's daily temperature and precipitation.
         !----------------------------------------------------------------
         use params_core
         implicit none

         ! arguments
         character(len=*), intent(in) :: runname
         integer, intent(in) :: realyear
         real, dimension(ndayyear), intent(out) :: dtemp, dprec
           
         dtemp(:) = 12.
         dprec(:) = 2.7  ! ~1000 mm / 365 d

         end subroutine getclimate

 
         subroutine getndep( runname, realyear, dprec, dndep )
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's annual ndeposition and distributes it
         !  over days according to daily precipitation.
         !----------------------------------------------------------------
         use params_core
         implicit none
         
         ! arguments
         character(len=*), intent(in) :: runname
         integer, intent(in) :: realyear
         real, dimension(ndayyear), intent(out) :: dprec, dndep

         ! local variables
         real :: andep
         real, dimension(ndayyear) :: dprec_rel

         ! xxx testing
         dprec(:) = 3.

         ! Get annual N-deposition           
         andep = 1.  

         ! Distribute annual Ndep to days by daily precipitation
         dprec_rel(:) = dprec(:)/sum(dprec(:))
         dndep = andep * dprec_rel

         end subroutine getndep


         subroutine getlanduse( runname, realyear, landuse )
         !////////////////////////////////////////////////////////////////
         !  SR reads this year's annual landuse state
         !----------------------------------------------------------------
         use params_core
         use gridvars
         implicit none

         ! arguments
         character(len=*), intent(in)      :: runname
         integer, intent(in)               :: realyear
         real, dimension(3), intent(in)    :: landuse  ! array of cropland/pasture/built-up, to be "translated" into 'lu_area' inside this SR

         lu_area(lunat,:) = 1.
         lu_area(lucrop,:) = 0.

         end subroutine getlanduse

         !===========================LOW-LEVEL============================

         function getvalreal( filename, realyear, day, dm, mo )
         !////////////////////////////////////////////////////////////////
         !  SR reads one (annual) value corresponding to the given year 
         !  from a time series ascii file. 
         !----------------------------------------------------------------

         implicit none
         ! arguments
         character(len=*), intent(in) :: filename
         integer, intent(in) :: realyear
         integer, intent(in), optional :: day ! day in year (1:365)
         integer, intent(in), optional :: dm  ! day in month (1:31)
         integer, intent(in), optional :: mo  ! month in year (1:12)

         ! function return value
         real :: getvalreal

         ! local variables
         integer :: l
         real :: tmp(3) ! 3 so that an additional value for this year could be read
         real :: realyear_decimal 

         if (present(day)) then
           ! convert day number into decimal number
           realyear_decimal = real(realyear) + real(day)/real(ndayyear)
         endif

         open(20, file=filename, status='old',  form='formatted', err=888)

         if (present(day)) then
           ! find corresponding day in first column and read 3 values on this line
           read(20, 100, err=999) (tmp(l), l=1,3)  
           do while (abs(realyear_decimal-tmp(1)).gt.1.0d-8)
             read(20, 100, err=999) (tmp(l), l=1,3)
           enddo
  
         else
           ! find corresponding year in first column and read 3 values on this line
           read(20, 100, err=999) (tmp(l), l=1,3)  
           do while (abs(realyear-tmp(1)).gt.1.0d-8)
             read(20, 100, err=999) (tmp(l), l=1,3)
           enddo

         endif

         getvalreal = tmp(2) 

 100     format (30d16.8)
         close(20)

         return

 888     write(0,*) 'GETVALREAL: error opening file '//trim(filename)//'. Abort. '
         stop
 999     write(0,*) 'GETVALREAL: error reading file '//trim(filename)//'. Abort. '
         stop 

         end function getvalreal



         function getvalreal_STANDARD( filename, realyear, mo, dm )
         !////////////////////////////////////////////////////////////////
         !  SR reads one (annual) value corresponding to the given year 
         !  from a time series ascii file. 
         !----------------------------------------------------------------

         implicit none
         ! arguments
         character(len=*), intent(in) :: filename
         integer, intent(in) :: realyear
         integer, intent(in), optional :: mo  ! month in year (1:12)
         integer, intent(in), optional :: dm  ! day in month (1:31)

         ! function return value
         real :: getvalreal_STANDARD

         ! local variables
         integer :: inyear
         integer :: inmo
         integer :: indm
         real    :: inval1
         real    :: inval2

         !print*,'looking for realyear, mo, dm',realyear,mo,dm

         ! open file
         open(20, file='./input/'//filename, status='old', form='formatted', err=888)

         if (present(mo)) then
             if (present(dm)) then
                 ! read 2 daily values
                 read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
                 do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 .or. (dm-indm).ne.0 )
                   read(20, 100, err=999) inyear, inmo, indm, inval1, inval2
                 enddo

             else           
                 ! read 2 monthly values
                 read(20, 200, err=999) inyear, inmo, inval1, inval2
                 do while ( (realyear-inyear).ne.0 .or. (mo-inmo).ne.0 )
                   read(20, 200, err=999) inyear, inmo, inval1, inval2
                 enddo

             end if
         else
             ! read 2 annual values
             read(20, 300, err=999) inyear, inval1, inval2
             do while ( (realyear-inyear).ne.0 )
               read(20, 300, err=999) inyear, inval1, inval2
             enddo

         end if

         !print*,'found realyear, mo, dm      ',inyear,inmo,indm,inval1

         getvalreal_STANDARD = inval1

 100     format (I4,I3,I3,F9.7,F9.7)
 200     format (I4,I3,F9.7,F9.7)
 300     format (I4,F9.7,F9.7)
         close(20)

         return

 888     write(0,*) 'GETVALREAL_STANDARD: error opening file '//trim(filename)//'. Abort. '
         stop
 999     write(0,*) 'GETVALREAL_STANDARD: error reading file '//trim(filename)//'. Abort. '
         stop 

         end function getvalreal_STANDARD

       end module forcing
 