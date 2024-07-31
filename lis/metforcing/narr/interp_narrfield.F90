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
! !ROUTINE: interp_narrfield
! \label{interp_narrfield}
!
! !REVISION HISTORY:
!  30 APR 2009: Sujay Kumar; Initial Specification
!  25 Jan 2012: Sujay Kumar; switched to the use of grib_api
!
! !INTERFACE:
subroutine interp_narrfield(n, narrfile, iv,pds5,pds6,pds7,metdata)
! !USES:  
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_logMod,         only : LIS_logunit,LIS_warning, LIS_verify
  use narr_forcingMod,    only : narr_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif
  implicit none
  
! !ARGUMENTS:
!
  integer,          intent(in) :: n
  character(len=*), intent(in) :: narrfile
  integer,          intent(in) :: iv
  integer,          intent(in) :: pds5
  integer,          intent(in) :: pds6
  integer,          intent(in) :: pds7
  real                         :: metdata(LIS_rc%ngrid(n))

! !DESCRIPTION:
!
!EOP
  integer               :: ftn
  integer               :: nnarr
  integer               :: iret, rc,status
  integer               :: igrib,kk,t,nvars
  logical               :: file_exists, var_found
  real                  :: missingValue
  integer               :: pds5_val,pds6_val,pds7_val
  logical*1, allocatable    :: lb(:)
  real,      allocatable    :: gi(:)
  real                  :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer               :: c,r

#if (defined USE_GRIBAPI) 

  nnarr = narr_struc(n)%nc*narr_struc(n)%nr

  inquire(file=trim(narrfile),exist=file_exists)
  if(file_exists) then 
     
     call grib_open_file(ftn,trim(narrfile),'r',iret)
     if(iret.ne.0) then
        write(LIS_logunit,*) 'Could not open file ',trim(narrfile)
        return
     endif
     
     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret,'error in grib_count_in_file in interp_narr')
     
     allocate(lb(nnarr))
     allocate(gi(nnarr))
     
     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call LIS_warning(iret, 'error in grib_new_from_file in interp_narr')
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(narrfile)
           deallocate(lb)
           deallocate(gi)
           return
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in interp_narr')

        call grib_get(igrib,'indicatorOfTypeOfLevel',pds6_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfTypeOfLevel in interp_narr')

        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in interp_narr')
        
        var_found = .false. 

        if((pds5_val.eq.pds5).and.&
             (pds6_val.eq.pds6).and.&
             (pds7_val.eq.pds7)) then 

           var_found = .true. 
           exit
        endif

        gi = -9999.0
        
        call grib_get(igrib, 'values',gi,rc)
        call LIS_warning(rc, 'error in grib_get:values in interp_narr')
        
        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(narrfile)
           deallocate(lb)
           deallocate(gi)
           return           
        endif
        
        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in interp_narr')
        
        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in interp_narr')

        if(var_found) then 
           
           lb = .false. 
           do t=1,nnarr
              if(gi(t).ne.missingValue) lb(t) = .true. 
           enddo
           
           call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
                lo, go, nnarr,LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                LIS_domain(n)%lat, LIS_domain(n)%lon,&
                narr_struc(n)%w111, narr_struc(n)%w121, &
                narr_struc(n)%w211, narr_struc(n)%w221, &
                narr_struc(n)%n111, narr_struc(n)%n121, &
                narr_struc(n)%n211, narr_struc(n)%n221, &
                LIS_rc%udef, iret)
           
           do r=1,LIS_rc%lnr(n)
              do c=1,LIS_rc%lnc(n)
                 if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                    metdata(LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              enddo
           enddo
        endif
     enddo
     
     call grib_close_file(ftn)
     deallocate(gi)
     deallocate(lb)

     if(.not.var_found) then 
        write(LIS_logunit,*) &
             'Could not retrieve entries in file: ',trim(narrfile)
        return
     endif
  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(narrfile)
  endif
#endif   
end subroutine interp_narrfield

