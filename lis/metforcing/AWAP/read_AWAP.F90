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
! !ROUTINE: read_AWAP
! \label{read_AWAP}
!
! !REVISION HISTORY:
!  30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:
subroutine read_AWAP( n, fname,findex,order, ferror_AWAP )

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use AWAP_forcingMod,    only : AWAP_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_AWAP

! !DESCRIPTION:
!  For the given time, reads parameters from AWAP datasets
!  and interpolates to a designated user-domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the AWAP file
!  \item[ferror\_AWAP]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_AWAP](\ref{interp_AWAP}) \newline
!    spatially interpolates the STAGE4 data
!  \end{description}
!EOP

  integer            :: i, j, t, iret ! Loop indicies and error flags
  integer            :: ftn, ndata, index1
  real               :: precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))     
  real               :: AWAPin(AWAP_struc(n)%ncol*AWAP_struc(n)%nrow)
  logical*1          :: lb(AWAP_struc(n)%ncol*AWAP_struc(n)%nrow)
  logical            :: file_exists            
  integer            :: igrib
  real               :: missingValue
  integer            :: pds5_val, pds7_val

  integer            :: c1,r1,c,r
  real               :: AWAPin1(AWAP_struc(n)%ncol*AWAP_struc(n)%nrow)

!=== End Variable Definition =======================
  if(order.eq.1) then 
     AWAP_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     AWAP_struc(n)%metdata2 = LIS_rc%udef
  endif

!-- Set necessary parameters for call to interp_AWAP   
  !precip_regrid = -1.0
  precip_regrid = LIS_rc%udef

!-- Check initially if file exists:
  inquire (file=fname, exist=file_exists ) ! Check if file exists
  if (.not. file_exists)  then 
     write(LIS_logunit,*)"[ERR] Missing AWAP precipitation file: ", fname
     ferror_AWAP = 1
     return
  endif

#if(defined USE_GRIBAPI) 

  ndata = AWAP_struc(n)%ncol * AWAP_struc(n)%nrow
  AWAPin = 0.0

  call grib_open_file(ftn,trim(fname),'r',iret)
  call LIS_verify(iret,'error grib_open_file in read_AWAP')
  
  call grib_new_from_file(ftn,igrib,iret)
  call LIS_verify(iret,'error in grib_new_from_file in read_AWAP')
  
  call grib_get(igrib,'indicatorOfParameter',pds5_val,iret)
  call LIS_verify(iret, 'error in grib_get: indicatorOfParameter in read_AWAP')
  
  call grib_get(igrib,'level',pds7_val,iret)
  call LIS_verify(iret, 'error in grib_get: level in read_AWAP')

  if(pds5_val.eq.20.and.pds7_val.eq.0) then 
     call grib_get(igrib,'values',AWAPin,iret)
     call LIS_verify(iret, 'error in grib_get:values in read_AWAP')
     
     call grib_get(igrib,'missingValue',missingValue,iret)
     call LIS_verify(iret, 'error in grib_get:missingValue in read_AWAP')
     
     do r=1, AWAP_struc(n)%nrow
        do c=1,AWAP_struc(n)%ncol
           c1 = c
           r1 = AWAP_struc(n)%nrow-r+1
           AWAPin1(c1+(r1-1)*AWAP_struc(n)%ncol) = &
                AWAPin(c+(r-1)*AWAP_struc(n)%ncol)
        enddo
     enddo

     lb = .false. 
     do t=1,ndata
        if ( AWAPin1(t) .ne. missingvalue ) then
           lb(t) = .true. 
        endif
     enddo

     call interp_AWAP( n, findex, ndata, AWAPin1, lb, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n), LIS_rc%lnr(n), precip_regrid )

     do j = 1, LIS_rc%lnr(n)
        do i = 1, LIS_rc%lnc(n)
           
           !if ( precip_regrid(i,j) .ne. -1.0 ) then
           if ( precip_regrid(i,j) .ne. LIS_rc%udef ) then
              index1 = LIS_domain(n)%gindex(i,j)
              if(index1 .ne. -1) then
                 if(order.eq.1) then 
                    AWAP_struc(n)%metdata1(1,index1) = precip_regrid(i,j) 
                 elseif(order.eq.2) then 
                    AWAP_struc(n)%metdata2(1,index1) = precip_regrid(i,j) 
                 endif
              endif
           endif
           
        enddo
     enddo

     call grib_release(igrib,iret)
     call LIS_verify(iret,'error in grib_release in read_AWAP')

  else
     write(LIS_logunit,*) 'Could not retrieve entries in file: ',trim(fname)
     ferror_AWAP = 1
     return
  endif

  call grib_close_file(ftn)

#endif
     
end subroutine read_AWAP

