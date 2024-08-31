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
! !ROUTINE: read_stg4
! \label{read_stg4}
!
! !REVISION HISTORY:
!  25 May 2006: Kristi Arsenault;  Data and code implementation
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!
! !INTERFACE:
subroutine read_stg4( n, fname,findex,order, ferror_stg4 )

! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_logMod,         only : LDT_verify, LDT_logunit, LDT_endrun
  use LDT_metforcingMod,  only : LDT_forc
  use stg4_forcingMod,    only : stg4_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=80)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_stg4

! !DESCRIPTION:
!  For the given time, reads parameters from STAGE4 datasets
!  and interpolates to a designated user-domain.
!  NOTE:: These subroutines use the READ\_GRIB routines for
!         for opening and reading the STAGE IV grib files.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the hourly STAGE4 file
!  \item[ferror\_stg4]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_stg4](\ref{interp_stg4}) \newline
!    spatially interpolates the STAGE4 data
!  \end{description}
!EOP

  integer         :: i, j, t, iret ! Loop indicies and error flags
  integer         :: ftn, ndata, index1
  real            :: precip_regrid(LDT_rc%lnc(n),LDT_rc%lnr(n))     
  real            :: stg4in(stg4_struc(n)%ncol*stg4_struc(n)%nrow)
  logical*1       :: lb(stg4_struc(n)%ncol*stg4_struc(n)%nrow)
  logical         :: file_exists            ! Check Stage IV file status 
  integer         :: igrib
  real            :: missingValue
  integer         :: pds5_val, pds7_val

!=== End Variable Definition =======================
  if(order.eq.1) then 
     LDT_forc(n,findex)%metdata1 = LDT_rc%udef
  elseif(order.eq.2) then 
     LDT_forc(n,findex)%metdata2 = LDT_rc%udef
  endif

!-- Set necessary parameters for call to interp_stg4   
  !precip_regrid = -1.0
  precip_regrid = LDT_rc%udef

!-- Check initially if file exists:
  inquire (file=fname, exist=file_exists ) ! Check if file exists
  if (.not. file_exists)  then 
     write(LDT_logunit,*)"** Missing STAGE IV precipitation file: ", fname
     ferror_stg4 = 1
     return
  endif

#if(defined USE_GRIBAPI) 

  ndata = stg4_struc(n)%ncol * stg4_struc(n)%nrow
  stg4in = 0.0

  call grib_open_file(ftn,trim(fname),'r',iret)
  call LDT_verify(iret,'error grib_open_file in read_stg4')
  
  call grib_new_from_file(ftn,igrib,iret)
  call LDT_verify(iret,'error in grib_new_from_file in read_stg4')
  
  call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
  call LDT_verify(iret, 'error in grib_get: indicatorOfParameter in read_stg4')
  
  call grib_get(igrib,'level',pds7_val,iret)
  call LDT_verify(iret, 'error in grib_get: level in read_stg4')

  if(pds5_val.eq.61.and.pds7_val.eq.0) then 
     call grib_get(igrib,'values',stg4in,iret)
     call LDT_verify(iret, 'error in grib_get:values in read_stg4')
     
     call grib_get(igrib,'missingValue',missingValue,iret)
     call LDT_verify(iret, 'error in grib_get:missingValue in read_stg4')
     
     lb = .false. 
     do t=1,ndata
        if ( stg4in(t) .ne. missingvalue ) then
           lb(t) = .true. 
        endif
     enddo

     call interp_stg4( n, findex, ndata, stg4in, lb, LDT_rc%gridDesc(n,:), &
          LDT_rc%lnc(n), LDT_rc%lnr(n), precip_regrid )
        
     do j = 1, LDT_rc%lnr(n)
        do i = 1, LDT_rc%lnc(n)
           
           !if ( precip_regrid(i,j) .ne. -1.0 ) then
           if ( precip_regrid(i,j) .ne. LDT_rc%udef ) then
              index1 = LDT_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    LDT_forc(n,findex)%metdata1(1,index1) = precip_regrid(i,j) 
                 elseif(order.eq.2) then 
                    LDT_forc(n,findex)%metdata2(1,index1) = precip_regrid(i,j) 
                 endif
              endif
           endif
           
        enddo
     enddo

     call grib_release(igrib,iret)
     call LDT_verify(iret,'error in grib_release in read_stg4')

  else
     write(LDT_logunit,*) 'Could not retrieve entries in file: ',trim(fname)
     ferror_stg4 = 1
     return
  endif

  call grib_close_file(ftn)

#endif
     
end subroutine read_stg4

