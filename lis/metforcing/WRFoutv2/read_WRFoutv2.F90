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
! !ROUTINE: read_WRFoutv2 
!  \label{read_WRFoutv2}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Beaudoing; Initial specification
!  14 Mar 2013: Sujay Kumar; Initial Specification in LIS
!  20 Nov 2020: K.R. Arsenault; Updated for different WRF output files
!
! !INTERFACE:
subroutine read_WRFoutv2( order, n, findex, yr, mon, da, hr, ferror )

! !USES:
  use LIS_coreMod,          only : LIS_rc,LIS_domain, LIS_localPet, &
                                   LIS_masterproc
  use LIS_metforcingMod,    only : LIS_forc
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_logMod,           only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN
  use WRFoutv2_forcingMod,  only : WRFoutv2_struc
  use LIS_forecastMod
  use LIS_mpiMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order         ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n             ! nest
  integer, intent(in)    :: findex        ! forcing index
  integer, intent(in)    :: yr,mon,da,hr  ! date and hour (multiple of 3)
  integer, intent(inout) :: ferror        ! set to zero if there's an error
!
! !DESCRIPTION:
!  For the given time, reads the parameters from 4-km
!  WRF CONUS domain data, transforms into 8 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  WRF variables used to force LIS are:
!  instanteous values starting at each timestep, available every 1-hour, \newline
!  except precipitation is accumulation over the previous 1-hour \newline
!
!  WRFOUT FORCING VARIABLES: \newline
!  1. T2          Air Temperature [K] \newline
!  2. Q2          Specific humidity [kg kg-1] \newline
!  3. SWDNB       Downward shortwave radiation [W m-2] \newline
!  4. LWDNB       Downward longwave radiation [W m-2] \newline
!  5. U10         10m U-dir wind vector [m s-1] \newline
!  6. V10         10m V-dir wind vector [m s-1] \newline
!  7. PSFC        Surface pressure [Pa] \newline
!  8. PREC_ACC_NC Accumulated Precipitation [mm] \newline
!
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
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
!  \item[WRFoutv2grid\_2\_lisgrid](\ref{WRFoutv2grid_2_lisgrid}) \newline
!    transform the WRF out data to the LIS grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor](\ref{conserv_interp}) \newline
!    select nearest neighbor value in the forcing data 
!  \end{description}
!
!EOP

  integer, parameter :: NF = 8   ! # of S forcing variables
  integer, parameter :: N_PF = 8 ! # of WRF out forcing variables

  character(11), dimension(N_PF), parameter :: WRFoutv2_fv = (/  &
       'T2         ',  &
       'Q2         ',  & 
       'SWDNB      ',  &
       'LWDNB      ',  &
       'U10        ',  &
       'V10        ',  &
       'PSFC       ',  &
       'PREC_ACC_NC'   /)

  character(len=LIS_CONST_PATH_LEN) :: infile

! netcdf variables
  integer :: ncid, varid, status
  integer :: latid, lonid

  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  integer :: timestep
  character :: cyr*4

  integer :: i, j, v, ii, r, k, eindex, x, y
  integer :: ios            ! set to non-zero if there's an error
  integer :: ierr
  integer :: nWRFoutv2      ! Size of I/O 1D fields
  integer :: iret, c
  real    :: gridDesco(50)         ! Input,output grid info arrays

  integer :: kk                    ! forecast index
  integer :: mo 

  real,allocatable :: lat(:,:)     ! input data (longitude,latitude)
  real,allocatable :: lon(:,:)     ! input data (longitude,latitude)

  real,allocatable :: datain(:,:)  ! input data (longitude,latitude)
  real,allocatable :: temp2WRFoutv2(:,:,:)
  real,allocatable :: templdas(:,:,:)
  real,allocatable :: f(:)         ! 1D in fields
  real,allocatable :: go(:)        ! 1D out fields
  real,allocatable :: tg(:,:)      ! Interpolated 2D data field
  logical*1,allocatable :: lb(:)   ! input bitmap
  logical*1,allocatable :: lo(:)   ! output bitmaps
! __________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! If a problem, ferror is set to zero
   ferror = 1
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

   ! Allocate memory
   allocate(lat(WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr))
   allocate(lon(WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr))
   allocate(datain(WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr))
   allocate(temp2WRFoutv2(WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr,NF))

   allocate(f(WRFoutv2_struc(n)%nc*WRFoutv2_struc(n)%nr))
   allocate(go(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   allocate(lb(WRFoutv2_struc(n)%nc*WRFoutv2_struc(n)%nr))
   allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
   allocate(tg(LIS_rc%lnc(n),LIS_rc%lnr(n)))  

   allocate(templdas(LIS_rc%lnc(n),LIS_rc%lnr(n),LIS_rc%nf), stat=ios)
   if(ios.ne.0) then 
     write(LIS_logunit,*) '[ERR] Error allocating templdas,',LIS_localPet
     call LIS_endrun
   endif

   temp2WRFoutv2 = 0.0     ! initialize

  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------
   mn=LIS_rc%mn ! Time of input file
   ss=0
   ts=0
   call LIS_tick(timenow,doy,gmt,yr,mon,da,hr,mn,ss,real(ts))

 ! One file per year--use day of year to get to the time record
   ! First year of data available begins in Oct 1, 2000:
   if( yr == 2000 ) then 
     timestep = (24*(doy - 1) + (1 + hr)) - 6576 ! Value represents: 24 * DOY of 09-30-2000
   else
     ! All subsequent years:
     timestep = (24*(doy - 1) + (1 + hr))  
   endif

   if(LIS_masterproc) then
      write(LIS_logunit,*)'[INFO] Order, timestep, doy, hr ::', &
            order, timestep, doy, hr
   endif

!=== Open WRF out forcing files ===

  ! Loop over forecast index:
  do kk= WRFoutv2_struc(n)%st_iterid, WRFoutv2_struc(n)%en_iterid

     ! Hindcast | OL run
     if(LIS_rc%forecastMode.eq.0) then 
        write(cyr, '(i4.4)') yr

     else ! Forecast mode (e.g., ESP)
        !sample yr, mo, da
        call LIS_sample_forecastDate(n, kk, findex, yr, mo, da) 
        write(cyr, '(i4.4)') yr
     endif

     ! Forcing variable loop:
     do v = 1, LIS_rc%met_nf(findex)  ! Number of met fields in WRFoutv2 data

       ! File name for data wrf_annual_CTRL_variable(v)_year.nc
       infile=trim(WRFoutv2_struc(n)%WRFoutv2dir)//&
              '/wrf_annual_CTRL_'//trim(WRFoutv2_fv(v))//'_'//cyr//'.nc'

       ! Open netCDF file.
       status = nf90_open(trim(infile), nf90_NoWrite, ncid)
       if(status/=0) then
         if(LIS_masterproc) then 
            write(LIS_logunit,*)'[ERR] Problem opening file: ',trim(infile),status
            call LIS_endrun
         endif
       else
         if(LIS_masterproc) then
           write(LIS_logunit,*)'[INFO] Opened file: ',trim(infile)
         endif
       end if

       datain = LIS_rc%udef
       status = nf90_inq_varid(ncid, trim(WRFoutv2_fv(v)), varid)
       status = nf90_get_var(ncid, varid, datain, &
                     start=(/1,1,timestep/), &
                     count=(/WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr,1/))
! KRA TO BE UPDATED TO SUPPORT PARALLEL READING IN OF FILES ...
!       status = nf90_get_var(ncid, varid, datain, &
!                     start=(/1,1,timestep/), &
!                     count=(/WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr,1/))
! KRA TO BE UPDATED TO SUPPORT PARALLEL READING IN OF FILES ...

       if( trim(WRFoutv2_fv(v)) == "T2" .and. timestep == 1 ) then
         ! Read in lat and lon info from the file:
         lat = LIS_rc%udef
         status = nf90_inq_varid(ncid, "XLAT", latid)
         status = nf90_get_var(ncid, latid, lat, &
                       start=(/1,1/), &
                       count=(/WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr/))

         lon = LIS_rc%udef
         status = nf90_inq_varid(ncid, "XLONG", lonid)
         status = nf90_get_var(ncid, lonid, lon, &
                       start=(/1,1/), &
                       count=(/WRFoutv2_struc(n)%nc,WRFoutv2_struc(n)%nr/))
       endif

       ! Close netCDF file.
       status=nf90_close(ncid)

      !-----------------------------------------------------------------
      ! Filter out any unrealistic forcing values.
      ! Transfer WRF out forcing fields to LIS format
      !-----------------------------------------------------------------
       do j=1,WRFoutv2_struc(n)%nr
         do i=1,WRFoutv2_struc(n)%nc

           select case (v)
            case (1) ! Tair
              temp2WRFoutv2(i,j,1) = datain(i,j)

            case (2) ! Qair
              temp2WRFoutv2(i,j,2) = datain(i,j)

            case (3) ! Shortwave
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2WRFoutv2(i,j,3) = datain(i,j)

            case (4) ! Longwave
              temp2WRFoutv2(i,j,4) = datain(i,j)

            case (5) ! U-vector wind
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2WRFoutv2(i,j,5) = datain(i,j)   ! Rotated wind

            case (6) ! V-vector wind
              if (datain(i,j) < 0.0001) then
                datain(i,j) = 0.0001
              endif
              temp2WRFoutv2(i,j,6) = datain(i,j)   ! Rotated wind

            case (7) ! Pressure
              temp2WRFoutv2(i,j,7) = datain(i,j)    

            case (8) ! Total precipitation: convective precp=0
              if (datain(i,j) < 0.0) then
                datain(i,j) = 0.0
              endif
              temp2WRFoutv2(i,j,8) = datain(i,j)
           end select
         enddo
       enddo
     end do !v or variable loop

     !--------------------------------------------------------------
     ! Interpolate each forcing variable to LIS domain
     !--------------------------------------------------------------

     !=== Initialize input & output grid arrays
     gridDesco = 0
           
     !=== Set input & output grid array values (WRF out to LIS)
     gridDesco = LIS_rc%gridDesc(n,:)

     !=== Define input & output data bitmaps
     nWRFoutv2 = WRFoutv2_struc(n)%nc*WRFoutv2_struc(n)%nr

     !== valid value over land and ocean for WRFoutv2 data
     lb = .true.
     lo = .false.

     templdas(:,:,LIS_rc%nf) = 0.0

     do v=1,NF      ! do not process convective precip field, which is set to 0
      
        !== Transferring current data to 1-D array for interpolation
        c=0
        do j=1,WRFoutv2_struc(n)%nr
           do i=1,WRFoutv2_struc(n)%nc
              c = c + 1
              f(c) = temp2WRFoutv2(i,j,v)
           enddo
        enddo

        !== Interpolate forcing to model grids
        select case( LIS_rc%met_interp(findex) )

          case( "bilinear" )
            call bilinear_interp(gridDesco,lb,f,lo,go,&
                 WRFoutv2_struc(n)%mi,mo, & 
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 WRFoutv2_struc(n)%w111,WRFoutv2_struc(n)%w121,&
                 WRFoutv2_struc(n)%w211,WRFoutv2_struc(n)%w221,& 
                 WRFoutv2_struc(n)%n111,WRFoutv2_struc(n)%n121,&
                 WRFoutv2_struc(n)%n211,WRFoutv2_struc(n)%n221,&
                 LIS_rc%udef, iret)

         case( "budget-bilinear" )
          if(v.eq.8) then 
            call conserv_interp(gridDesco,lb,f,lo,go,&
                 WRFoutv2_struc(n)%mi,mo, & 
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 WRFoutv2_struc(n)%w112,WRFoutv2_struc(n)%w122,&
                 WRFoutv2_struc(n)%w212,WRFoutv2_struc(n)%w222,& 
                 WRFoutv2_struc(n)%n112,WRFoutv2_struc(n)%n122,&
                 WRFoutv2_struc(n)%n212,WRFoutv2_struc(n)%n222,&
                 LIS_rc%udef,iret)
          else
            call bilinear_interp(gridDesco,lb,f,lo,go,WRFoutv2_struc(n)%mi,mo, & 
                 LIS_domain(n)%lat, LIS_domain(n)%lon,&
                 WRFoutv2_struc(n)%w111,WRFoutv2_struc(n)%w121,&
                 WRFoutv2_struc(n)%w211,WRFoutv2_struc(n)%w221,& 
                 WRFoutv2_struc(n)%n111,WRFoutv2_struc(n)%n121,&
                 WRFoutv2_struc(n)%n211,WRFoutv2_struc(n)%n221,LIS_rc%udef, iret)
          endif

         case( "neighbor" )
           call neighbor_interp(gridDesco,lb,f,lo,go,&
                WRFoutv2_struc(n)%mi,mo,&
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                WRFoutv2_struc(n)%n113,LIS_rc%udef,iret)

        end select

    
        !== Convert data to original 3D array & a 2D array to 
        !== Fill in of missing points due to geography difference  
        tg = 0.0 
        c = 0
        do j = 1, LIS_rc%lnr(n)
            do i = 1, LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then 
                tg(i,j) = go(i+c)
              endif
            enddo
            c = c + LIS_rc%lnc(n)
         enddo

         !== No need to fill in for WRFoutv2

         do j = 1, LIS_rc%lnr(n)
           do i = 1, LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(i,j).ne.-1) then 
                if ((tg(i,j) .lt. 0) .and. (tg(i,j) .ne. LIS_rc%udef)) then
                  write(LIS_logunit,*)'[ERR] No nearest neighbors, v, i, j',v,i,j,tg(i,j)
                  call LIS_endrun()  
                endif
                templdas(i,j,v) = tg(i,j)
              endif
           end do ! c
        enddo     ! r

    enddo  ! end v=variable loop

    do v= 1, NF
      ! Fill in undefined and ocean points
      do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
          if (LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1)then
              WRFoutv2_struc(n)%metdata1(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
           else
              WRFoutv2_struc(n)%metdata2(kk,v,LIS_domain(n)%gindex(c,r))=templdas(c,r,v)
           endif
          endif
        enddo !c
      enddo   !r
    enddo     !v

  end do  ! End forecast member index loop

  ! Deallocate local interp-based variables:
  deallocate(datain)
  deallocate(temp2WRFoutv2)
  deallocate(f)
  deallocate(go)
  deallocate(lb)
  deallocate(lo)
  deallocate(tg)
  deallocate(templdas)
#endif

end subroutine read_WRFoutv2


