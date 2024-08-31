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
! !ROUTINE: read_cmap
! \label{read_cmap}
!
! !REVISION HISTORY:
!  17 Jul 2001: Jon Gottschalck; Initial code
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  10 Oct 2014: KR Arsenault: Added to LDT
!  
! !INTERFACE:
subroutine read_cmap( n, filename, findex, order, ferror_cmap, filehr)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_warning, LDT_verify
  use LDT_metforcingMod, only : LDT_forc
  use cmap_forcingMod,   only : cmap_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)   :: filename          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_cmap
  integer             :: filehr
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  CMAP data and interpolates to the LDT domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[filename]
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
  real                   :: precip_regrid(LDT_rc%lnc(n),LDT_rc%lnr(n))
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
  inquire (file=trim(filename), exist=file_exists)

! File exists:
  if (file_exists) then   
     
     ncmap = cmap_struc(n)%nc * cmap_struc(n)%nr
     allocate(cmapin(ncmap))
     allocate(lb(ncmap)) 

   ! Open the CMAP grib file:
     call grib_open_file(ftn,trim(filename),'r',iret)
     if(iret.ne.0) then 
        write(LDT_logunit,*) &
            "Could not open file: ",trim(filename)
        ferror_cmap = 0
        return
     endif
     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, "error in grib_count_in_file in read_cmap")
     
   ! Search for appropiate CMAP PPT variable:
     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call LDT_warning(iret,'error in grib_new_from_file in read_cmap')
        
        if(iret.ne.0) then 
           write(LDT_logunit,*) &
                "Could not retrieve entries in file: ",trim(filename)
           ferror_cmap = 0
           deallocate(lb)
           deallocate(cmapin)
           return           
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LDT_verify(rc, "error in grib_get: indicatorOfParameter in read_cmap")
        
        call grib_get(igrib,'level',pds7_val,rc)
        call LDT_verify(rc, "error in grib_get: level in read_cmap")
        
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

      ! Get CMAP PPT values from the native grib file:
        call grib_get(igrib,'values',cmapin,rc)
        call LDT_warning(rc, 'error in grib_get:values in read_cmap')
        
        if(rc.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(filename)
           ferror_cmap = 0
           deallocate(lb)
           deallocate(cmapin)
           return           
        endif

      ! OBtain what the the designated missingValue stored is:
        call grib_get(igrib,'missingValue',missingValue,rc)
        call LDT_verify(rc, "error in grib_get:missingValue in read_cmap")

        call grib_release(igrib,rc)
        call LDT_verify(rc, "error in grib_release in read_cmap")
        
      ! If PPT data present, map logical values to True:
        if( var_found ) then 
           lb = .false.
           do t=1,ncmap
              if(cmapin(t).ne.missingValue) lb(t) = .true. 
           enddo
           
         ! Spatially interpolate native CMAP file to LIS run domain:
           call interp_cmap(n,ncmap,cmapin,lb,LDT_rc%gridDesc, &
                       LDT_rc%lnc(n),LDT_rc%lnr(n),precip_regrid)

         ! Transfer spatially interpolated 1D PPT field to 2D:
           do j = 1,LDT_rc%lnr(n)
              do i = 1,LDT_rc%lnc(n)
                 if (precip_regrid(i,j) .ge. 0.0) then
                    index1 = LDT_domain(n)%gindex(i,j)
                    if(index1 .ne. -1) then
                       if(order.eq.2) then 
                          LDT_forc(n,findex)%metdata2(1,index1) = &
                               precip_regrid(i,j) ! *3600.0 - note: CMAP/GDAS already rate
                       endif
                    endif
                 endif
              enddo
           enddo

        endif  ! If variable found ...
     end do    ! Search for appropiate variable
  
     call grib_close_file(ftn)
     
     deallocate(lb)
     deallocate(cmapin)
     
     do kk=1,iv_total
        if(.not.var_status(kk)) then 
           write(LDT_logunit,*) &
                "Could not retrieve entries in file: ",trim(filename)
           ferror_cmap = 0
           
           return
        endif
     enddo

! File does not exist:
  else
     write(LDT_logunit,*) "Could not find file: ",trim(filename)
     ferror_cmap = 0
  endif

#endif     
     
end subroutine read_cmap

