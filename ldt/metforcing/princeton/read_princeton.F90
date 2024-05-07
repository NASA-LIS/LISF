!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_princeton 
!  \label{read_princeton}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LDT/retberg.F90
!
! !INTERFACE:
subroutine read_princeton( order, n, findex, yr, mon, da, hr, ferror )

! !USES:
  use LDT_coreMod,          only : LDT_rc,LDT_domain, LDT_localPet, &
                                   LDT_masterproc
  use LDT_metforcingMod,    only : LDT_forc
  use LDT_timeMgrMod,       only : LDT_tick
  use LDT_logMod,           only : LDT_logunit, LDT_endrun
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use princeton_forcingMod, only : princeton_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n      ! nest
  integer, intent(in)    :: order  ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: findex
  integer, intent(in)    :: yr,mon,da,hr     ! data and hour (multiple of 3)
  integer, intent(inout) :: ferror           ! set to zero if there's an error

! !DESCRIPTION:
!  For the given time, reads the parameters from 1 degree
!  PRINCETON data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain.

!  PRINCETON variables used to force LDT are:
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
!  \item[princetongrid\_2\_ldtgrid](\ref{princetongrid_2_ldtgrid}) \newline
!    transform the PRINCETON data to the LDT grid 
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP

! Specify Reanalysis PRINCETON forcing parameters & file paths
  integer,parameter::      nvars =    5          ! number of variables: last
                                                 ! variable is of interest
  integer,parameter::      ndims =    4          ! number of dimensions
  integer,parameter::      z =      1
  real, parameter :: no_data_value = 2.e+20

  integer, parameter :: N_PF=7 ! # of PRINCETON forcing variables
  integer, parameter :: NF=9   ! # of GLDAS forcing variables

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

  character(len=LDT_CONST_PATH_LEN) :: infile

! netcdf variables
  integer :: ncid, varid, status

  integer :: mo
  real*8  :: timenow
  real    :: gmt
  integer :: doy,mn,ss,ts
  character :: cyr*4, cmo*2 
  integer :: timestep

  integer :: i, j, v, ii, r, kk, eindex, x, y
  integer :: ios                         ! set to non-zero if there's an error


  integer :: gldas,nprinceton                 ! Size of I/O 1D fields
  integer :: iret,c
  real :: gridDesco(20)                       ! Input,output grid info arrays

  real,allocatable :: datain(:,:)  ! input data (longitude,latitude)
  real,allocatable :: temp2princeton(:,:,:)
  real,allocatable :: templdas(:,:,:)
  real,allocatable :: f(:)! 1D in fields
  real,allocatable :: go(:)  ! 1D out fields
  real,allocatable :: tg(:,:)  ! Interpolated 2D data field
  logical*1,allocatable :: lb(:) ! input bitmap
  logical*1,allocatable :: lo(:)      ! output bitmaps
! following variables commented out since geogfill routine is not used
!  logical*1 :: geogmask(LDT_rc%lnc(n),LDT_rc%lnr(n))! 2D output bitmap
!  logical*1 :: tmask(LDT_rc%lnc(n),LDT_rc%lnr(n))   ! 2d valid temperature bitmap

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  ! if a problem, ferror is set to zero
  ferror = 1
  mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
  !----------------------------------------------------------------
  ! Establish which file timestep the date & hour correspond to
  ! and which file to open.
  !----------------------------------------------------------------

  mn=LDT_rc%mn ! Time of input file
  ss=0
  ts=0

  call LDT_tick(timenow,doy,gmt,yr,mon,da,hr,mn,ss,real(ts))
!  print*,'timenow, doy:',timenow, doy

!== One file per year--use day of year to get to the time record
  timestep = 8*(doy - 1) + (1 + hr/3)

  if(LDT_masterproc) &
       write(LDT_logunit,*) 'order and month-timestep', order, timestep, doy, hr
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mon

!== allocate memory
  allocate(datain(princeton_struc(n)%nc,princeton_struc(n)%nr))
  allocate(temp2princeton(princeton_struc(n)%nc,princeton_struc(n)%nr,NF))
  allocate(f(princeton_struc(n)%nc*princeton_struc(n)%nr))
  allocate(go(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
  allocate(lb(princeton_struc(n)%nc*princeton_struc(n)%nr))
  allocate(lo(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
  allocate(tg(LDT_rc%lnc(n),LDT_rc%lnr(n)))  
  allocate(templdas(LDT_rc%lnc(n),LDT_rc%lnr(n),LDT_rc%nf), stat=ios)

  if(ios .ne.0) then 
     write(LDT_logunit,*) 'Error allocating templdas.',LDT_localPet
     stop 344
  endif

  temp2princeton = 0.0                 ! initialize

!=== Open PRINCETON forcing files

     do v = 1, 7 !N_PF
!     File name for data year/variable(v)_3hourly_year-year.nc
      infile=trim(princeton_struc(n)%princetondir)//'/'//cyr//'/'//&
             trim(princeton_fv(v))//'_3hourly_'//cyr//'-'//cyr//'.nc'

!     Open netCDF file.
      status = nf90_open(infile, nf90_NoWrite, ncid)
      status = nf90_inq_varid(ncid, trim(princeton_fv(v)), varid)

      if (status/=0) then
         if(LDT_masterproc) then 
            write(LDT_logunit,*) 'Problem opening file: ', trim(infile),status
            write(LDT_logunit,*) 'Stopping...'
         endif
         call LDT_endrun
      else
         if(LDT_masterproc) write(LDT_logunit,*) 'Opened file: ', trim(infile)
      end if

      status = nf90_get_var(ncid, varid, datain, &
                                     start=(/1,1,1,timestep/), &
      count=(/princeton_struc(n)%nc,princeton_struc(n)%nr,1,1/))

!     Close netCDF file.
      status=nf90_close(ncid)

     !----------------------------------------------------------------
     ! Change data from PRINCETON grid convention to GLDAS one 
     ! Shift longitudes by 180deg. 
     !----------------------------------------------------------------
     call princetongrid_2_ldtgrid(princeton_struc(n)%nc,&
                                  princeton_struc(n)%nr,datain)

  !-----------------------------------------------------------------
  ! Filter out any unrealdttic forcing values.
  ! Transfer PRINCETON forcing fields to GLDAS format
  !-----------------------------------------------------------------
      do j=1,princeton_struc(n)%nr
       do i=1,princeton_struc(n)%nc

    select case (v)
    case (1)! tair
        temp2princeton(i,j,1) = datain(i,j)
    case (2)! qair
        temp2princeton(i,j,2) = datain(i,j)
    case (3)! Shortwave
        if (datain(i,j) < 0.0001) &
             datain(i,j) = 0.0001
        temp2princeton(i,j,3) = datain(i,j)
    case (4)! Longwave
        temp2princeton(i,j,4) = datain(i,j)
    case (5)! Wind
        if (datain(i,j) < 0.0001) &
             datain(i,j) = 0.0001
        temp2princeton(i,j,5) = datain(i,j)    !Since absolute wind speed 
                                               ! let U=WIND and V=0.0
    case (6)! Total precipitation: convective precp=0
        if (datain(i,j) < 0.0) &
             datain(i,j) = 0.0
        temp2princeton(i,j,7) = datain(i,j)    
    case (7)! pressure
        temp2princeton(i,j,8) = datain(i,j)
    end select

       end do
      enddo

     end do !v
  !-----------------------------------------------------------------
  ! Interpolate each forcing variable to GLDAS domain
  !-----------------------------------------------------------------

  !=== Initialize input & output grid arrays
  gridDesco = 0
           
  !=== Set input & output grid array values (reanlPRINCETON to GLDAS)
        
  gridDesco = LDT_rc%gridDesc(n,:)

  !=== Define input & output data bitmaps
  nprinceton = princeton_struc(n)%nc*princeton_struc(n)%nr
  gldas  = LDT_rc%lnc(n)*LDT_rc%lnr(n)
!== valid value over land and ocean for Princeton data
  lb = .true.
  lo = .false.
!  tmask = .false.

  templdas(:,:,9) = 0.0

  do v=1,NF-1      ! do not process convective precip field, which is set to 0
     if (v .ne. 6) then ! not the v-wind component, which is set to zero.
     
      !=== Transferring current data to 1-D array for interpolation
      c=0
      do i=1,princeton_struc(n)%nr
          do j=1,princeton_struc(n)%nc
             c = c + 1
             f(c) = temp2princeton(j,i,v)
!             if(v.eq.1 .and. f(c).ne.-9999.0) write(LDT_logunit,*) c,f(c)
          enddo
      enddo

      !=== Interpolate if forcing and model grids are not both one deg.
      if (LDT_rc%gridDesc(n,9).ne.1.0) then

        select case( LDT_rc%met_gridtransform(findex) )

          case( "bilinear" )
             call bilinear_interp(gridDesco,lb,f,lo,go,&
                  princeton_struc(n)%mi,mo, & 
                  LDT_domain(n)%lat, LDT_domain(n)%lon,&
                  princeton_struc(n)%w111,princeton_struc(n)%w121,&
                  princeton_struc(n)%w211,princeton_struc(n)%w221,& 
                  princeton_struc(n)%n111,princeton_struc(n)%n121,&
                  princeton_struc(n)%n211,princeton_struc(n)%n221,&
                  LDT_rc%udef, iret)

          case( "budget-bilinear" )
             if(v.eq.7) then 
                call conserv_interp(gridDesco,lb,f,lo,go,&
                     princeton_struc(n)%mi,mo, & 
                     LDT_domain(n)%lat, LDT_domain(n)%lon,&
                     princeton_struc(n)%w112,princeton_struc(n)%w122,&
                     princeton_struc(n)%w212,princeton_struc(n)%w222,& 
                     princeton_struc(n)%n112,princeton_struc(n)%n122,&
                     princeton_struc(n)%n212,princeton_struc(n)%n222,&
                     LDT_rc%udef,iret)
             else
                call bilinear_interp(gridDesco,lb,f,lo,go,princeton_struc(n)%mi,mo, & 
                     LDT_domain(n)%lat, LDT_domain(n)%lon,&
                     princeton_struc(n)%w111,princeton_struc(n)%w121,&
                     princeton_struc(n)%w211,princeton_struc(n)%w221, & 
                     princeton_struc(n)%n111,princeton_struc(n)%n121,&
                     princeton_struc(n)%n211,princeton_struc(n)%n221,LDT_rc%udef, iret)
             endif
          end select

      else ! forcing and model grids both one degree
        kk = 0
        do r=1,LDT_rc%lnr(n)
          y = r + (LDT_rc%gridDesc(n,4) + 89.50) / 1.000
          do c=1,LDT_rc%lnc(n)
            x = c + (LDT_rc%gridDesc(n,5) + 179.50) / 1.000
            kk = kk + 1
            eindex = ((y - 1) * 360) + x
            go(kk) = f(eindex)
            lo(kk) = lb(eindex)
          end do
        end do

      end if ! LDT_rc%domain 
    
       !=== Convert data to original 3D array & a 2D array to 
       !=== fill in of missing points due to geography difference  
!      tg = -9999.0
      tg = 0.0 
      c = 0
      do j = 1, LDT_rc%lnr(n)
         do i = 1, LDT_rc%lnc(n)
!            write(LDT_logunit,*) i,j,gindex(i,j),lo(i+c)
            if(LDT_domain(n)%gindex(i,j).ne.-1) then 
!               geogmask(i,j) = lo(i+c)
               tg(i,j) = go(i+c)
            endif
         enddo
         c = c + LDT_rc%lnc(n)
      enddo

!== no need to fill in for Princeton
!      call geogfill2(n, LDT_rc%lnc(n),LDT_rc%lnr(n),geogmask,tg,v,tmask)
      !       write(LDT_logunit,*) gindex(21,3),tg(21,3),v

      do j = 1, LDT_rc%lnr(n)
         do i = 1, LDT_rc%lnc(n)
            if(LDT_domain(n)%gindex(i,j).ne.-1) then 
               if ((tg(i,j) .lt. 0) .and. (tg(i,j) .ne. LDT_rc%udef)) then
                  write(LDT_logunit,*) 'No nearest neighbours, v, i, j',v,i,j,tg(i,j)
                  stop
               endif
               templdas(i,j,v) = tg(i,j)
            endif
         end do !c
       enddo ! r

    else ! v==6, v-wind component, always zero
       templdas(:,:,6) = 0.0
    endif
    
 enddo !v

  do v= 1,NF
    !=== Fill in undefined and ocean points
      do r = 1,LDT_rc%lnr(n)
       do c = 1,LDT_rc%lnc(n)
        if (LDT_domain(n)%gindex(c,r).ne.-1) then
        if(order.eq.1)then
           LDT_forc(n,findex)%metdata1(v,LDT_domain(n)%gindex(c,r))=templdas(c,r,v)
        else
           LDT_forc(n,findex)%metdata2(v,LDT_domain(n)%gindex(c,r))=templdas(c,r,v)
        endif
       endif
      enddo !c
     enddo !r
  enddo !v
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
! !ROUTINE: princetongrid_2_ldtgrid
! \label{princetongrid_2_ldtgrid}
!
! !REVISION HISTORY:
!  10 Apr 2002: Urszula Jambor;  Code adapted from 
!               ecmwfgrid_2_grid2catgrid, by R. Reichle
!  29 Jan 2007: Hiroko Kato; modified berggrid_2_ldtgrid for Princeton data
! !INTERFACE:
subroutine princetongrid_2_ldtgrid( nx, ny, grid_data )

  implicit none
! !ARGUMENTS:   
  integer, intent(in)                   :: nx, ny
  real, intent(inout), dimension(nx,ny) :: grid_data
!
! !DESCRIPTION:
! Changes grid data from PRINCETON data convention to LDT convention
!
! PRINCETON: Sorth-to-Nouth around Greenwich Meridian
! Global grid. Data are written in NetCDF from ``lower left to 
! upper right'' starting at 0.5-degree grid point center coordinates: 
! 0.5E,89.5S and going to 0.5W,89.5N.
! 
! LDT: South-to-North around Date Line
! Full global grid.  Starts at the southernmost latitude and date line, 
! going east and then north.
!EOP
  
  integer :: i, j, m
  real :: tmp, tmp_data1(nx)
  
  ! ------------------------------------------------------------------
  ! some checks
  
  if ((nx /= 360) .or. (ny /= 180)) then
     write (*,*) 'Rprincetongrid_2_gldasgrid(): This routine has only been'
     write (*,*) 'checked for nx=360 and ny=180. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
  end if  
  if ((mod(nx,2) /= 0) .or. (mod(ny,2) /= 0)) then
     write (*,*) 'Rprincetongrid_2_gldasgrid(): This routine can only work'
     write (*,*) 'for even nx and ny. Make sure you know'
     write (*,*) 'what you are doing. STOPPING.'
     stop
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

end subroutine Princetongrid_2_ldtgrid



