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
! !ROUTINE: read_AWRAL
! \label{read_AWRAL}
!
! !REVISION HISTORY:
!  30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:
subroutine read_AWRAL( order, n, findex, year, doy, ferror_AWRAL )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_metforcingMod,  only : LIS_forc
  use AWRAL_forcingMod,    only : AWRAL_struc
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order     ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n         ! nest
  integer, intent(in)    :: findex    ! forcing index
  integer, intent(in) :: year,doy
  character(4) :: cyear
  integer             :: ferror_AWRAL

! !DESCRIPTION:
!  For the given time, reads parameters from AWRAL datasets
!  and interpolates to a designated user-domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ferror\_AWRAL]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_AWRAL](\ref{interp_AWRAL}) \newline
!    spatially interpolates the STAGE4 data
!  \end{description}
!EOP

  integer            :: i, j, t, iret ! Loop indicies and error flags
  integer            :: ftn, ndata, index1,v
  integer, parameter :: N_AF = 6 ! # of AWRAL forcing variables
  character(10), dimension(N_AF), parameter :: awral_fv = (/  &
       'tat       ',    &
       'avpt      ',    &
       'rgt       ',    &
       'radcskyt  ',    &
       'u2t       ',    &
       'pt        '     /)


  character(len=LIS_CONST_PATH_LEN) :: var_fname
  real,allocatable  :: datain(:,:) ! input data (lat,lon)
  logical            :: file_exists            
! netcdf variables
  integer :: ncid, varid, status 
  integer            :: timestep
  real               :: missingValue, latt, lont
  integer            :: pds5_val, pds7_val

  integer            :: c1,r1,c,r,rowt,colt
  real,allocatable   :: temp2awral(:,:,:)

  allocate(datain(AWRAL_struc(n)%ncol,AWRAL_struc(n)%nrow))
  allocate(temp2awral(AWRAL_struc(n)%ncol,AWRAL_struc(n)%nrow,N_AF))

!=== End Variable Definition =======================
  if(order.eq.1) then 
     AWRAL_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     AWRAL_struc(n)%metdata2 = LIS_rc%udef
  endif


!-- Set necessary parameters for call to interp_AWRAL   
  write(cyear, '(i4.4)') year

!-- Check initially if file exists:
  do v = 1, N_AF  ! N_AF
    var_fname = trim(AWRAL_struc(n)%AWRALdir)//'/'//trim(awral_fv(v))//'_'//trim(cyear)//'.nc'
    inquire (file=trim(var_fname), exist=file_exists ) ! Check if file exists
     if (.not. file_exists)  then 
       write(LIS_logunit,*)"[ERR] Missing AWRAL file: ", trim(var_fname)
       ferror_AWRAL = 1
       return
    endif
  enddo

!-- Get timestep from cdoy --!
  timestep = doy

!-- Check forcing is on the same grid as the model --!
  if((AWRAL_struc(n)%nrow .ne. LIS_rc%gnr(n)) .or. (AWRAL_struc(n)%ncol .ne. LIS_rc%gnc(n))) then
     if(LIS_masterproc) then
          write(LIS_logunit,*)'[ERR] Problem using AWRAL forcing: Forcing must be on the same grid as the model'
          write(LIS_logunit,*)'[ERR] Remapping not implemented. Stopping...'
          call LIS_endrun
     endif
  endif
  

  do v = 1, N_AF ! N_AF
    ndata = AWRAL_struc(n)%ncol * AWRAL_struc(n)%nrow
    datain = 0.0
    var_fname = trim(AWRAL_struc(n)%AWRALdir)//'/'//trim(awral_fv(v))//'_'//trim(cyear)//'.nc'
    write(LIS_logunit,*)"[INFO] Attempting to read file: ", trim(var_fname)
    !-- netcdf reader --!
    ! Open netCDF file.
    status = nf90_open(var_fname, nf90_NoWrite, ncid)
    status = nf90_inq_varid(ncid, trim(awral_fv(v)), varid)
  
    if(status/=0) then
       if(LIS_masterproc) then
            write(LIS_logunit,*)'[ERR] Problem opening file: ',trim(var_fname),status
            write(LIS_logunit,*)'[ERR]  Stopping...'
            call LIS_endrun
       endif
         call LIS_endrun
    else
       if(LIS_masterproc) write(LIS_logunit,*)'[INFO] Opened file: ',trim(var_fname)
    endif

    status = nf90_get_var(ncid, varid, datain, &
             start=(/1,1,timestep/), &
             count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/))
 
    status=nf90_close(ncid)
    do j=1,AWRAL_struc(n)%nrow
         do i=1,AWRAL_struc(n)%ncol
            temp2awral(i,j,v) = datain(i,j)           
         enddo
    enddo
  enddo


  do v = 1, N_AF ! N_AF   
    do j = 1, LIS_rc%lnr(n)
        do i = 1, LIS_rc%lnc(n)
           index1 = LIS_domain(n)%gindex(i,j)
           if(index1 .ne. -1) then
             latt = LIS_domain(n)%grid(index1)%lat
             lont = LIS_domain(n)%grid(index1)%lon
             call awrallatlon_2_globalgrid(latt, lont, AWRAL_struc(n)%gridDesci(4), AWRAL_struc(n)%gridDesci(5), &
                                          AWRAL_struc(n)%gridDesci(9), AWRAL_struc(n)%gridDesci(10), rowt, colt)
             if ( temp2awral(colt,rowt,v) .ne. LIS_rc%udef ) then
               if(order.eq.1) then 
                 AWRAL_struc(n)%metdata1(1,v,index1) = temp2awral(colt,rowt,v)
               elseif(order.eq.2) then 
                 AWRAL_struc(n)%metdata2(1,v,index1) = temp2awral(colt,rowt,v)
               endif
             endif
           endif
        enddo
    enddo
  enddo
  deallocate(datain)
  deallocate(temp2awral)
end subroutine read_AWRAL

!BOP
!
! !ROUTINE: awrallatlon_2_globalgrid
! \label{awrallatlon_2_globalgrid}
!
! !REVISION HISTORY:
!  05.06.2019: Wendy Sharples
!
! !INTERFACE:
subroutine awrallatlon_2_globalgrid(lat, lon, latllcnr, lonllcnr, deltalat, deltalon, row, col)
  use LIS_logMod,           only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS:
  real, intent(in)                     :: lat, lon, latllcnr, lonllcnr, deltalat, deltalon
  integer, intent(inout)               :: row, col
!
! !DESCRIPTION:
! Gets global i,j from lat,lon
!
!
!EOP

  real                                :: row_real, col_real
  
  row_real = (deltalat + ABS(latllcnr - lat))/deltalat
  col_real = (deltalon + ABS(lonllcnr - lon))/deltalon
  row = NINT(row_real)
  col = NINT(col_real) 
end subroutine awrallatlon_2_globalgrid

