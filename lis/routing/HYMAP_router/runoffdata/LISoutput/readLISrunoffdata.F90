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
! !DESCRIPTION: 
! 
! !REVISION HISTORY: 
! 7 Jan 2016: Sujay Kumar, Initial implementation
! 
! !USES: 
subroutine readLISrunoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_logMod
  use LISrunoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))

  !caveats 
  ! 1) assumes the LIS outputs are the same output interval as that of
  ! the HYMAP model timestep. No temporal aggregation is done. 
  !
  ! 3) Assumes that the units of Qs and Qsb are in kg/m2s
  !
  ! 4) Assumes that the LIS outputs are in the same grid/map projection/
  ! resolution.
  ! 
  ! 5) LIS outputs are in NETCDF format. 
  !
  integer                       :: c,r,t
  integer, allocatable      :: nqs(:,:)
  real,   allocatable       :: qs(:,:),qs_t(:)
  real,   allocatable       :: qsb(:,:),qsb_t(:)
  integer                   :: ios, nid,qsid,qsbid
  character(len=LIS_CONST_PATH_LEN) :: filename
  logical               :: file_exists
  logical               :: check_Flag
  !create LIS filename

  call LIS_create_output_filename(n, &
       filename, check_flag, &
       model_name='SURFACEMODEL', &
       odir=LISrunoffdata_struc(n)%odir,&
       writeint=LISrunoffdata_struc(n)%outInterval)
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  
  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)
     allocate(qs(LIS_rc%gnc(n),LIS_rc%gnr(n)))
     allocate(qsb(LIS_rc%gnc(n),LIS_rc%gnr(n)))
     
     ios = nf90_open(path=filename,&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in readLISrunoffdata')
     
     ios = nf90_inq_varid(nid,'Qs_tavg',qsid)
     call LIS_verify(ios,'failed to read Qs_tavg field in readLISrunoffdata')
     
     ios = nf90_inq_varid(nid,'Qsb_tavg',qsbid)
     call LIS_verify(ios,'failed to read Qsb_tavg field in readLISrunoffdata')
     
     if(LIS_rc%wopt.eq."2d gridspace") then 
        ios = nf90_get_var(nid,qsid,qs)
        call LIS_verify(ios, 'failed to read Qs_tavg field in readLISrunoffdata')
        
        ios = nf90_get_var(nid,qsbid,qsb)
        call LIS_verify(ios, 'failed to read Qsb_tavg field in readLISrunoffdata')
     elseif(LIS_rc%wopt.eq."1d tilespace") then 
        
        allocate(qs_t(LIS_rc%glbntiles(n)))
        allocate(qsb_t(LIS_rc%glbntiles(n)))
        allocate(nqs(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        
        ios = nf90_get_var(nid,qsid,qs_t)
        call LIS_verify(ios, 'failed to read Qs_tavg field in readLISrunoffdata')
        
        ios = nf90_get_var(nid,qsbid,qsb_t)
        call LIS_verify(ios, 'failed to read Qsb_tavg field in readLISrunoffdata')
        
        do t=1,LIS_rc%glbntiles(n)
           c = LIS_domain(n)%tile(t)%col
           r = LIS_domain(n)%tile(t)%row
           
           qs(c,r) = qs(c,r) + qs_t(t)
           qsb(c,r) = qsb(c,r) + qsb_t(t)
           nqs(c,r) = nqs(c,r) + 1
        enddo
        
        do r=1,LIS_rc%gnr(n)
           do c=1,LIS_rc%gnc(n)
              if(nqs(c,r).gt.0) then 
                 qs(c,r) = qs(c,r)/nqs(c,r)
                 qsb(c,r) = qsb(c,r)/nqs(c,r)
              else
                 qs(c,r) = LIS_rc%udef
                 qsb(c,r) = LIS_rc%udef
              endif
           enddo
        enddo
        deallocate(qs_t)
        deallocate(qsb_t)
        deallocate(nqs)
     endif
     
     do r=1,LIS_rc%gnr(n)
        do c=1,LIS_rc%gnc(n)
           if(qs(c,r).ne.-9999.0) then 
!                          surface_runoff(c,r) = qs(c,r)/(LIS_CONST_RHOFW)
!                          baseflow(c,r) = qsb(c,r)/(LIS_CONST_RHOFW)
              surface_runoff(c,r) = qs(c,r)
              baseflow(c,r) = qsb(c,r)
           else
              surface_runoff(c,r) = 0.0
              baseflow(c,r) = 0.0
           endif
        enddo
     enddo
     
     deallocate(qs)
     deallocate(qsb)
     call LIS_verify(nf90_close(nid))
  else
     write(LIS_logunit,*) 'Failed to find '//trim(filename)
     call LIS_endrun()
  endif
#endif


end subroutine readLISrunoffdata
