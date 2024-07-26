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
! !ROUTINE: read_cmap
! \label{read_cmap}
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  02 Dec 2014: KR Arsenault; Updated CMAP reader
!  
! !INTERFACE:
subroutine read_cmap( n, fname, findex, order, ferror_cmap, filehr)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_warning, LIS_verify
  use LIS_metforcingMod, only : LIS_forc
  use cmap_forcingMod,   only : cmap_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_cmap
  integer             :: filehr
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CMAP data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the 6 hour CMAP file
!  \item[ferror\_cmap]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_cmap](\ref{interp_cmap}) \newline
!    spatially interpolates the CMAP data
!  \end{description}
!EOP

  integer                :: ftn
  real                   :: precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                :: ncmap
  real, allocatable      :: cmapin(:)
  logical*1,allocatable  :: lb(:)
  integer                :: index1
  integer                :: i,j,t
  logical                :: file_exists
  integer                :: igrib
  integer                :: iv, iv_total
  integer                :: pds5_val, pds7_val
  integer                :: pds5(1), pds7(1)
  logical                :: var_status(1)
  logical                :: var_found
  integer                :: kk,var_index
  real                   :: missingValue
  integer                :: nvars
  integer                :: iret,rc,status
  
!=== End Variable Definition =======================

#if(defined USE_GRIBAPI) 
  pds5 = (/ 059 /) !parameter
  pds7 = (/ 000 /) !htlev2

  iv_total = 1
  inquire (file=trim(fname), exist=file_exists)
! File exists:
  if (file_exists) then   
     
     ncmap = cmap_struc(n)%ncold*cmap_struc(n)%nrold
     allocate(cmapin(ncmap))
     allocate(lb(ncmap)) 

   ! Open the CMAP grib file:
     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             'Could not open file: ',trim(fname)
        ferror_cmap = 0
        return
     endif
     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_cmap')
     
   ! Search for appropiate CMAP PPT variable:
     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call LIS_warning(iret,'error in grib_new_from_file in read_cmap')
        
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_cmap = 0
           deallocate(lb)
           deallocate(cmapin)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_cmap')
        
        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_cmap')
        
        var_found = .false. 
        
        do iv=1,iv_total
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv))) then 
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo

        call grib_get(igrib,'values',cmapin,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_cmap')
        
        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_cmap = 0
           deallocate(lb)
           deallocate(cmapin)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_cmap')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_cmap')
        
        if(var_found) then 
           lb = .false.
           do t=1,ncmap
              if(cmapin(t).ne.missingValue) lb(t) = .true. 
           enddo

         ! Spatially interpolate the CMAP field to target LIS grid:
           call interp_cmap(n,findex,ncmap,cmapin,lb,LIS_rc%gridDesc, &
                            LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)

           do j = 1,LIS_rc%lnr(n)
              do i = 1,LIS_rc%lnc(n)
                 if (precip_regrid(i,j) .ge. 0.0) then
                    index1 = LIS_domain(n)%gindex(i,j)
                    if(index1 .ne. -1) then
                       if(order.eq.2) then 
                          cmap_struc(n)%metdata2(1,index1) = &
                               precip_regrid(i,j) !*3600.0
                       endif
                    endif
                 endif
              enddo
           enddo
        endif
     end do
  
     call grib_close_file(ftn)
     
     deallocate(lb)
     deallocate(cmapin)
     
     do kk=1,iv_total
        if(.not.var_status(kk)) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_cmap = 0
           
           return
        endif
     enddo
  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(fname)
     ferror_cmap = 0
  endif
#endif     
     
end subroutine read_cmap





