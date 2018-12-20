!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
subroutine read_AWRAL( n, fname, index, order, year, day, ferror_AWRAL )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use AWRAL_forcingMod,    only : AWRAL_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: order     ! lower(1) or upper(2) time interval bdry
  integer, intent(in)    :: n         ! nest
  integer, intent(in)    :: index    ! forcing index
  character(len=80)   :: fname          
  integer, intent(in) :: year
  integer, intent(in) :: day
  character(4) :: cyear
  character(3) :: cday

  integer             :: ferror_AWRAL

! !DESCRIPTION:
!  For the given time, reads parameters from AWRAL datasets
!  and interpolates to a designated user-domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the AWRAL file : file naming convention is currently: AWRALdir/AWRAL_calendaryear_calendarday/var.grb'
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

  character(10) :: var_fname
  real               :: regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real               :: AWRALin(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)
  logical*1          :: lb(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)
  logical            :: file_exists            
  integer            :: igrib
! netcdf variables
  integer :: ncid, varid, status 
  integer            :: timestep
  real               :: missingValue
  integer            :: pds5_val, pds7_val

  integer            :: c1,r1,c,r
  real               :: AWRALin1(AWRAL_struc(n)%ncol*AWRAL_struc(n)%nrow)


  write ( cyear, '(i4)' ) year
  write ( cday, '(i3)' ) day

!=== End Variable Definition =======================
  if(order.eq.1) then 
     AWRAL_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     AWRAL_struc(n)%metdata2 = LIS_rc%udef
  endif

!-- Set necessary parameters for call to interp_AWRAL   

!-- Check initially if file exists:
  do v = 1,  AWRAL_struc(n)%N_AF  ! N_AF
    var_fname = fname//'/'//AWRAL_struc(n)%awral_fv(v)//'_'//cyear//'.nc'
    inquire (file=trim(var_fname), exist=file_exists ) ! Check if file exists
     if (.not. file_exists)  then 
       write(LIS_logunit,*)"[ERR] Missing AWRAL file: ", fname
       ferror_AWRAL = 1
       return
    endif
  enddo

!-- Get timestep from cdoy --!
  timestep = day - 1

  do v = 1, AWRAL_struc(n)%N_AF ! N_AF
    regrid = LIS_rc%udef
    ndata = AWRAL_struc(n)%ncol * AWRAL_struc(n)%nrow
    AWRALin = 0.0
    var_fname = fname//'/'//AWRAL_struc(n)%awral_fv(v)//'_'//cyear//'.nc'
     
    !-- netcdf reader --!
    ! Open netCDF file.
    status = nf90_open(var_fname, nf90_NoWrite, ncid)
    status = nf90_inq_varid(ncid, trim(AWRAL_struc(n)%awral_fv(v)), varid)
  
    if(status/=0) then
       if(LIS_masterproc) then
            write(LIS_logunit,*)'[ERR] Problem opening file: ',var_fname,status
            write(LIS_logunit,*)'[ERR]  Stopping...'
            call LIS_endrun
       endif
         call LIS_endrun
    else
       if(LIS_masterproc) write(LIS_logunit,*)'[INFO] Opened file: ',var_fname
    end if

    status = nf90_get_var(ncid, varid, AWRALin, &
                                     start=(/1,1,timestep/), &
    count=(/AWRAL_struc(n)%ncol,AWRAL_struc(n)%nrow,1/))

    do r=1, AWRAL_struc(n)%nrow
        do c=1,AWRAL_struc(n)%ncol
           c1 = c
           r1 = AWRAL_struc(n)%nrow-r+1
           AWRALin1(c1+(r1-1)*AWRAL_struc(n)%ncol) = &
                AWRALin(c+(r-1)*AWRAL_struc(n)%ncol)
        enddo
    enddo

    lb = .false. 
    do t=1,ndata
        if ( AWRALin1(t) .ne. missingvalue ) then
           lb(t) = .true. 
        endif
    enddo

    call interp_AWRAL( n, index, ndata, AWRALin1, lb, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n), LIS_rc%lnr(n), regrid )

    do j = 1, LIS_rc%lnr(n)
        do i = 1, LIS_rc%lnc(n)
           
           !if ( precip_regrid(i,j) .ne. -1.0 ) then
           if ( regrid(i,j) .ne. LIS_rc%udef ) then
              index1 = LIS_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    AWRAL_struc(n)%metdata1(1,v,index1) = regrid(i,j) 
                 elseif(order.eq.2) then 
                    AWRAL_struc(n)%metdata2(1,v,index1) = regrid(i,j) 
                 endif
              endif
           endif
           
        enddo
    enddo

  ! Close netCDF file.
    status=nf90_close(ncid)
  enddo
     
end subroutine read_AWRAL

