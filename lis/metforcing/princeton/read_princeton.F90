!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_princeton 
!  \label{read_princeton}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Beaudoing; Initial Specification adopted from 
!                                 LIS/retberg.F90
!  16 Feb 2016: Hiroko Beaudoing; Fixed indexes for precip and pressure fields 
!                                 so budget-bilinear interpolation applies to 
!                                 correct individual fields.
!  15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
!                           that is in 3D array (4D in version 1 & 2).
!  22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
!
! !INTERFACE:
subroutine read_princeton( order, n, findex, yr, mon, da, hr, ferror )

! !USES:
  use LIS_coreMod,          only : LIS_rc,LIS_domain, LIS_localPet, &
                                   LIS_masterproc
  use LIS_metforcingMod,    only : LIS_forc
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN
  use princeton_forcingMod, only : princeton_struc
  use LIS_forecastMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order     ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n         ! nest
  integer, intent(in)    :: findex    ! forcing index
  integer, intent(in)    :: yr,mon,da,hr     ! data and hour (multiple of 3)
  integer, intent(inout) :: ferror           ! set to zero if there's an error
!
! !DESCRIPTION:
!  For the given time, reads the parameters from 1 degree
!  PRINCETON data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  PRINCETON variables used to force LIS are:
!  mean values starting at timestep, available every 3 hours \newline
!
!  NOTE-1: be aware that PRINCETON has only total precipitation.
!  NOTE 2: only one wind component, it is magnitude \newline
!
!  PRINCETON FORCING VARIABLES: \newline
!  1. TAS       Air Temperature [K] \newline
!  2. SHUM      Specific humidity [kg kg-1] \newline
!  3. DSWRF     Downward shortwave radiation [W m-2] \newline
!  4. DLWRF     Downward longwave radiation [W m-2] \newline
!  5. WIND      Wind speed [m s-1] \newline
!  6. PRES      Surface pressure [Pa] \newline
!  7. PRCP      Precipitation [kg m-2 s-1] \newline
!
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    3 hourly instance, order=2, read the next 3 hourly instance)
!  \item[n]
!    index of the nest
!  \item[yr]
!    current year
!  \item[mon]
!    current month
!  \item[da]
!    current day of the year
!  \item[hr]
!    current hour of day
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[princetongrid\_2\_lisgrid](\ref{princetongrid_2_lisgrid}) \newline
!    transform the PRINCETON data to the LIS grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!
!EOP

! Specify Reanalysis PRINCETON forcing parameters & file paths
  integer,parameter ::  nvars =  5          ! number of variables: last
                                             ! variable is of interest
  integer,parameter ::  ndims =  4          ! number of dimensions
  integer,parameter ::  z =      1
  real,   parameter ::  no_data_value = 2.e+20

  integer, parameter :: N_PF = 7 ! # of PRINCETON forcing variables
  integer, parameter :: NF = 9   ! # of GLDAS forcing variables

  character(10), dimension(N_PF), parameter :: princeton_fv = (/  &
       'tas       ',    &
       'shum      ',    & 
       'dswrf     ',    &
       'dlwrf     ',    &
       'wind      ',    &
       'pres      ',    &
       'prcp      '     /)

  integer, dimension(N_PF), parameter :: fnum = (/ &
       31, 33, 34, 35, 36, 37, 38 /) 

  character(len=LIS_CONST_PATH_LEN) :: infile

! netcdf variables
  integer :: ncid, varid, status

  integer :: mo
  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  character :: cyr*4
  integer :: timestep

  integer :: i, j, v, ii, r, k, eindex, x, y
  integer :: ios                   ! set to non-zero if there's an error
  integer :: gldas, nprinceton     ! Size of I/O 1D fields
  integer :: iret, c
  real    :: gridDesco(50)         ! Input,output grid info arrays

  integer :: kk                    ! forecast index
  integer :: num_met               ! number of fields, V3 will have 6, 2 and 2.2 have 7 (precip)

  real,allocatable :: datain(:,:)  ! input data (longitude,latitude)
  real,allocatable :: temp2princeton(:,:,:)
  real,allocatable :: templdas(:,:,:)
  real,allocatable :: f(:)         ! 1D in fields
  real,allocatable :: go(:)        ! 1D out fields
  real,allocatable :: tg(:,:)      ! Interpolated 2D data field
  logical*1,allocatable :: lb(:)   ! input bitmap
  logical*1,allocatable :: lo(:)   ! output bitmaps
! following variables commented out since geogfill routine is not used
!  logical*1 :: geogmask(LIS_rc%lnc(n),LIS_rc%lnr(n))! 2D output bitmap
!  logical*1 :: tmask(LIS_rc%lnc(n),LIS_rc%lnr(n))   ! 2d valid temperature bitmap
! ______________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! If a problem, ferror is set to zero
   ferror = 1
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
   num_met = 7

   ! Allocate memory
   allocate(datain(princeton_struc(n)%ncold,princeton_struc(n)%nrold))
   allocate(temp2princeton(princeton_struc(n)%ncold,princeton_struc(n)%nrold,NF))
   allocate(f(princeton_struc(n)%ncold*princeton_struc(n)%nrold))
   allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   allocate(lb(princeton_struc(n)%ncold*princeton_struc(n)%nrold))
   allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   allocate(tg(LIS_rc%lnc(n),LIS_rc%lnr(n)))  
   allocate(templdas(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nf), stat=ios)
 
   if(ios .ne.0) then 
     write(LIS_logunit,*) '[ERR] Error allocating templdas,',LIS_localPet
     stop 344
   endif

   temp2princeton = 0.0     ! initialize

  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------
   mn=LIS_rc%mn ! Time of input file
   ss=0
   ts=0
   call LIS_tick(timenow,doy,gmt,yr,mon,da,hr,mn,ss,real(ts))

   ! One file per year--use day of year to get to the time record
   timestep = 8*(doy - 1) + (1 + hr/3)

   if(LIS_masterproc) then
      write(LIS_logunit,*)'[INFO] Order and month-timestep', order, timestep, doy, hr
   endif

!=== Open PRINCETON forcing files ===

! If version 3 is being used, set num_met to six because
!    the 3-hr precip fields are not available. Turn off
!    if future versions include 3-hr precip fields.
  if(princeton_struc(n)%version == "3") then
    num_met = 6
  endif

  ! Loop over forecast index:
  do kk= princeton_struc(n)%st_iterid, princeton_struc(n)%en_iterid

     ! Modified by KRA to implement forecast mode:
     if(LIS_rc%forecastMode.eq.0) then ! hindcast run
        write(cyr, '(i4.4)') yr

     else !forecast mode
        !sample yr, mo, da
        call LIS_sample_forecastDate(n, kk, findex, yr, mo, da) 
        write(cyr, '(i4.4)') yr
     endif

     ! Forcing variable loop:
     do v = 1, num_met  ! Number of met fields in Princeton data

       ! File name for data year/variable(v)_3hourly_year-year.nc
       infile=trim(princeton_struc(n)%princetondir)//'/'//cyr//'/'//&
              trim(princeton_fv(v))//'_3hourly_'//cyr//'-'//cyr//'.nc'

       ! Open netCDF file.
       status = nf90_open(infile, nf90_NoWrite, ncid)
       status = nf90_inq_varid(ncid, trim(princeton_fv(v)), varid)

       if(status/=0) then
         if(LIS_masterproc) then 
            write(LIS_logunit,*)'[ERR] Problem opening file: ',infile,status
            write(LIS_logunit,*)'[ERR]  Stopping...'
            call LIS_endrun
         endif
         call LIS_endrun
       else
         if(LIS_masterproc) write(LIS_logunit,*)'[INFO] Opened file: ',infile
       end if

      if (princeton_struc(n)%version=="2.2" .OR. princeton_struc(n)%version=="3")then
       status = nf90_get_var(ncid, varid, datain, &
                                     start=(/1,1,timestep/), &
       count=(/princeton_struc(n)%ncold,princeton_struc(n)%nrold,1/))
      else
       status = nf90_get_var(ncid, varid, datain, &
                                     start=(/1,1,1,timestep/), &
       count=(/princeton_struc(n)%ncold,princeton_struc(n)%nrold,1,1/))
      endif

       ! Close netCDF file.
       status=nf90_close(ncid)

      !----------------------------------------------------------------
      ! Change data from PRINCETON grid convention to LIS-domain
      ! Shift longitudes by 180deg. 
      !----------------------------------------------------------------
       call princetongrid_2_lisgrid(princeton_struc(n)%ncold,&
                                    princeton_struc(n)%nrold,datain)

      !-----------------------------------------------------------------
      ! Filter out any unrealistic forcing values.
      ! Transfer PRINCETON forcing fields to LIS format
      !-----------------------------------------------------------------
       do j=1,princeton_struc(n)%nrold
         do i=1,princeton_struc(n)%ncold

           select case (v)
            case (1)! tair
              temp2princeton(i,j,1) = datain(i,j)
            case (2)! qair
              temp2princeton(i,j,2) = datain(i,j)
            case (3)! Shortwave
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2princeton(i,j,3) = datain(i,j)
            case (4)! Longwave
              temp2princeton(i,j,4) = datain(i,j)
            case (5)! Wind
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2princeton(i,j,5) = datain(i,j)  !Since absolute wind speed 
                                                   ! let U=WIND and V=0.0
            case (6)! pressure
              temp2princeton(i,j,7) = datain(i,j)    
            case (7)! Total precipitation: convective precp=0
              if (datain(i,j) < 0.0) then
                datain(i,j) = 0.0
              endif
              temp2princeton(i,j,8) = datain(i,j)
           end select
         enddo
       enddo
     end do !v or variable loop

     !--------------------------------------------------------------
     ! Interpolate each forcing variable to LIS domain
     !--------------------------------------------------------------

     !=== Initialize input & output grid arrays
     gridDesco = 0
           
     !=== Set input & output grid array values (reanlPRINCETON to LIS)
     gridDesco = LIS_rc%gridDesc(n,:)

     !=== Define input & output data bitmaps
     nprinceton = princeton_struc(n)%ncold*princeton_struc(n)%nrold
     gldas  = LIS_rc%lnc(n)*LIS_rc%lnr(n)

     !== valid value over land and ocean for Princeton data
     lb = .true.
     lo = .false.
   !  tmask = .false.

     templdas(:,:,9) = 0.0

     do v=1,NF-1      ! do not process convective precip field, which is set to 0
        if (v .ne. 6) then ! not the v-wind component, which is set to zero.
      
          !=== Transferring current data to 1-D array for interpolation
          c=0
          do i=1,princeton_struc(n)%nrold
             do j=1,princeton_struc(n)%ncold
                c = c + 1
                f(c) = temp2princeton(j,i,v)
!                if(v.eq.1 .and. f(c).ne.-9999.0) write(LIS_logunit,*) c,f(c)
             enddo
          enddo

          !=== Interpolate if forcing and model grids are not both one deg.
          if( LIS_rc%gridDesc(n,9).ne.1.0 ) then

           select case( LIS_rc%met_interp(findex) )

             case( "bilinear" )
               call bilinear_interp(gridDesco,lb,f,lo,go,&
                  princeton_struc(n)%mi,mo, & 
                  LIS_domain(n)%lat, LIS_domain(n)%lon,&
                  princeton_struc(n)%w111,princeton_struc(n)%w121,&
                  princeton_struc(n)%w211,princeton_struc(n)%w221,& 
                  princeton_struc(n)%n111,princeton_struc(n)%n121,&
                  princeton_struc(n)%n211,princeton_struc(n)%n221,&
                  LIS_rc%udef, iret)

             case( "budget-bilinear" )
               if(v.eq.8) then 
                 call conserv_interp(gridDesco,lb,f,lo,go,&
                     princeton_struc(n)%mi,mo, & 
                     LIS_domain(n)%lat, LIS_domain(n)%lon,&
                     princeton_struc(n)%w112,princeton_struc(n)%w122,&
                     princeton_struc(n)%w212,princeton_struc(n)%w222,& 
                     princeton_struc(n)%n112,princeton_struc(n)%n122,&
                     princeton_struc(n)%n212,princeton_struc(n)%n222,&
                     LIS_rc%udef,iret)
               else
                 call bilinear_interp(gridDesco,lb,f,lo,go,princeton_struc(n)%mi,mo, & 
                     LIS_domain(n)%lat, LIS_domain(n)%lon,&
                     princeton_struc(n)%w111,princeton_struc(n)%w121,&
                     princeton_struc(n)%w211,princeton_struc(n)%w221, & 
                     princeton_struc(n)%n111,princeton_struc(n)%n121,&
                     princeton_struc(n)%n211,princeton_struc(n)%n221,LIS_rc%udef, iret)
               endif

           end select

         else ! forcing and model grids both one degree
            k = 0
            do r=1,LIS_rc%lnr(n)
               y = r + (LIS_rc%gridDesc(n,4) + 89.50) / 1.000
               do c=1,LIS_rc%lnc(n)
                  x = c + (LIS_rc%gridDesc(n,5) + 179.50) / 1.000
                  k = k + 1
                  eindex = ((y - 1) * 360) + x
                  go(k) = f(eindex)
                  lo(k) = lb(eindex)
               end do
            end do

         end if ! LIS_rc%domain 
    
         !=== Convert data to original 3D array & a 2D array to 
         !=== fill in of missing points due to geography difference  
!         tg = -9999.0
         tg = 0.0 
         c = 0
         do j = 1, LIS_rc%lnr(n)
            do i = 1, LIS_rc%lnc(n)
!              write(LIS_logunit,*) i,j,gindex(i,j),lo(i+c)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then 
!                geogmask(i,j) = lo(i+c)
                tg(i,j) = go(i+c)
              endif
            enddo
            c = c + LIS_rc%lnc(n)
         enddo

         !== no need to fill in for Princeton
         ! call geogfill2(n, LIS_rc%lnc(n),LIS_rc%lnr(n),geogmask,tg,v,tmask)
         !       write(LIS_logunit,*) gindex(21,3),tg(21,3),v

         do j = 1, LIS_rc%lnr(n)
           do i = 1, LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then 
                if ((tg(i,j) .lt. 0) .and. (tg(i,j) .ne. LIS_rc%udef)) then
                   !For version 3, replace the value with an undefined value instead of ending the run.
                   !For all other versions, keep old call to end the program.
                   if(princeton_struc(n)%version == "3") then
                      write(LIS_logunit,*)'[WARN] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
                      write(LIS_logunit,*)'[WARN] Check output to make sure the missing neighbor'
                      write(LIS_logunit,*)'is isolated and does not distort the output.'      
                      tg(i,j) = LIS_rc%udef ! New code to change data to undefined
                   else
                      write(LIS_logunit,*)'[WARN] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
                      call LIS_endrun()  ! Old code to call a fatal error
                   endif
                endif
                templdas(i,j,v) = tg(i,j)
              endif
           end do ! c
        enddo    ! r

      else ! v==6, v-wind component, always zero
         templdas(:,:,6) = 0.0
      endif
    enddo  ! end v=variable loop

    do v= 1,NF
      ! Fill in undefined and ocean points
      do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
          if (LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1)then
              princeton_struc(n)%metdata1(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
           else
              princeton_struc(n)%metdata2(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
           endif
          endif
        enddo !c
      enddo   !r
    enddo     !v

  end do  ! End forecast member index loop

  ! Deallocate local interp-based variables:
  deallocate(datain)
  deallocate(temp2princeton)
  deallocate(f)
  deallocate(go)
  deallocate(lb)
  deallocate(lo)
  deallocate(tg)
  deallocate(templdas)
#endif

end subroutine read_princeton


!BOP
! 
! !ROUTINE: princetongrid_2_lisgrid
! \label{princetongrid_2_lisgrid}
!
! !REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from 
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
!  29 Jan 2007: Hiroko Kato; modified berggrid_2_lisgrid for Princeton data
!
! !INTERFACE:
subroutine princetongrid_2_lisgrid( nx, ny, grid_data )
  use LIS_logMod,           only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS:   
  integer, intent(in)                   :: nx, ny
  real, intent(inout), dimension(nx,ny) :: grid_data
!
! !DESCRIPTION:
! Changes grid data from PRINCETON data convention to LIS convention
!
! PRINCETON: Sorth-to-Nouth around Greenwich Meridian
! Global grid. Data are written in NetCDF from ``lower left to 
! upper right'' starting at 0.5-degree grid point center coordinates: 
! 0.5E,89.5S and going to 0.5W,89.5N.
! 
! LIS: South-to-North around Date Line
! Full global grid.  Starts at the southernmost latitude and date line, 
! going east and then north.
!
!EOP
  
  integer :: i, j, m
  real    :: tmp, tmp_data1(nx)
  
  ! ------------------------------------------------------------------
  ! Some checks
  ! ------------------------------------------------------------------
  if( ((nx /= 360) .or. (ny /= 180)) .and. ((nx /= 1440) .or. (ny /= 600)) ) then
     write(LIS_logunit,*) '[ERR] Rprincetongrid_2_gldasgrid(): This routine has only been'
     write(LIS_logunit,*) '  checked for nx=360 and ny=180 (version 2 and 2.2) and for '
     write(LIS_logunit,*) '  nx=1440 and ny=600 (version 3). ' 
     write(LIS_logunit,*) '  Make sure you know what you are doing. STOPPING.'
     call LIS_endrun()
  end if  
  if ((mod(nx,2) /= 0) .or. (mod(ny,2) /= 0)) then
     write(LIS_logunit,*) '[ERR] Rprincetongrid_2_gldasgrid(): This routine can only work'
     write(LIS_logunit,*) '  for even nx and ny. Make sure you know'
     write(LIS_logunit,*) '  what you are doing. STOPPING.'
     call LIS_endrun()
  end if
  
  !-------------------------------------------------------------------
  
  do j=1,ny
     
     do i=1,nx
     tmp_data1(i)    = grid_data(i,j)
     end do
     
     do i=1,nx/2

        ! shift longitudes (wrapping around Greenwhich Meridian becomes
        !  wrapping around Date Line)
        m = i + nx/2
        tmp          = tmp_data1(i)
        tmp_data1(i) = tmp_data1(m)
        tmp_data1(m) = tmp
        
     end do
     
     do i=1,nx
     grid_data(i,j) = tmp_data1(i)
     end do
     
  end do

end subroutine Princetongrid_2_lisgrid



