!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! 
! !DESCRIPTION: 
! 
! !REVISION HISTORY: 
!  7 Jan 2016: Sujay Kumar, Initial implementation
! 17 Mar 2016: Augusto Getirana, changes in input file name generation and surface runoff and baseflow variables - this will reduce the number of times input files are read
! 
! !USES: 
subroutine readLISrunoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_logMod
  use LISrunoffdataMod
  use LIS_fileIOMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: total_evapotranspiration(LIS_rc%gnc(n),LIS_rc%gnr(n))

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
  !Added total evapotranspiration (Evap)
  integer                       :: c,r,t
  integer, allocatable      :: nqs(:,:)

  !ag - 17Mar2016
  !real,   allocatable       :: qs(:,:),qs_t(:)
  !real,   allocatable       :: qsb(:,:),qsb_t(:)
  real,   allocatable       :: qs_t(:),qsb_t(:),evap_t(:)
  
  integer                   :: ios, nid,qsid,qsbid,evapid
  character*100         :: filename
  logical               :: file_exists
  logical               :: check_Flag
  !create LIS filename

  call LIS_create_output_filename(n, &
       filename, check_flag, &
       model_name='SURFACEMODEL', &
       odir=LISrunoffdata_struc(n)%odir,&
       writeint=LISrunoffdata_struc(n)%outInterval)
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  
  if(trim(LISrunoffdata_struc(n)%previous_filename)/=trim(filename))then
  
    LISrunoffdata_struc(n)%previous_filename=filename
    
    inquire(file=filename, exist=file_exists)
    if(file_exists) then 
      write(LIS_logunit,*) 'Reading '//trim(filename)
      !ag - 17Mar2016
      !allocate(qs(LIS_rc%gnc(n),LIS_rc%gnr(n)))
      !allocate(qsb(LIS_rc%gnc(n),LIS_rc%gnr(n)))
     
      ios = nf90_open(path=filename,&
           mode=NF90_NOWRITE,ncid=nid)
      call LIS_verify(ios,'Error in readLISrunoffdata')
     
      ios = nf90_inq_varid(nid,'Qs_tavg',qsid)
      call LIS_verify(ios,'failed to read Qs_tavg field in readLISrunoffdata')
     
      ios = nf90_inq_varid(nid,'Qsb_tavg',qsbid)
      call LIS_verify(ios,'failed to read Qsb_tavg field in readLISrunoffdata')
     
      ios = nf90_inq_varid(nid,'Evap_tavg',evapid)
      call LIS_verify(ios,'failed to read Evap_tavg field in readLISrunoffdata')
     
      if(LIS_rc%wopt.eq."2d gridspace") then 
        ios = nf90_get_var(nid,qsid,LISrunoffdata_struc(n)%qs)
        call LIS_verify(ios, 'failed to read Qs_tavg field in readLISrunoffdata')
        
        ios = nf90_get_var(nid,qsbid,LISrunoffdata_struc(n)%qsb)
        call LIS_verify(ios, 'failed to read Qsb_tavg field in readLISrunoffdata')

        ios = nf90_get_var(nid,evapid,LISrunoffdata_struc(n)%evap)
        call LIS_verify(ios, 'failed to read Evap_tavg field in readLISrunoffdata')
      elseif(LIS_rc%wopt.eq."1d tilespace") then 
        
        allocate(qs_t(LIS_rc%glbntiles(n)))
        allocate(qsb_t(LIS_rc%glbntiles(n)))
        allocate(evap_t(LIS_rc%glbntiles(n)))
        allocate(nqs(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        
        ios = nf90_get_var(nid,qsid,qs_t)
        call LIS_verify(ios, 'failed to read Qs_tavg field in readLISrunoffdata')
        
        ios = nf90_get_var(nid,qsbid,qsb_t)
        call LIS_verify(ios, 'failed to read Qsb_tavg field in readLISrunoffdata')

        ios = nf90_get_var(nid,evapid,evap_t)
        call LIS_verify(ios, 'failed to read Evap_tavg field in readLISrunoffdata')
        
        do t=1,LIS_rc%glbntiles(n)
          c = LIS_domain(n)%tile(t)%col
          r = LIS_domain(n)%tile(t)%row
           
          LISrunoffdata_struc(n)%qs(c,r) = LISrunoffdata_struc(n)%qs(c,r) + qs_t(t)
          LISrunoffdata_struc(n)%qsb(c,r) = LISrunoffdata_struc(n)%qsb(c,r) + qsb_t(t)
          LISrunoffdata_struc(n)%evap(c,r) = LISrunoffdata_struc(n)%evap(c,r) + evap_t(t)
          nqs(c,r) = nqs(c,r) + 1
        enddo
        
        do r=1,LIS_rc%gnr(n)
          do c=1,LIS_rc%gnc(n)
            if(nqs(c,r).gt.0) then 
              LISrunoffdata_struc(n)%qs(c,r) = LISrunoffdata_struc(n)%qs(c,r)/nqs(c,r)
              LISrunoffdata_struc(n)%qsb(c,r) = LISrunoffdata_struc(n)%qsb(c,r)/nqs(c,r)
              LISrunoffdata_struc(n)%evap(c,r) = LISrunoffdata_struc(n)%evap(c,r)/nqs(c,r)
            else
              LISrunoffdata_struc(n)%qs(c,r) = LIS_rc%udef
              LISrunoffdata_struc(n)%qsb(c,r) = LIS_rc%udef
              LISrunoffdata_struc(n)%evap(c,r) = LIS_rc%udef
            endif
          enddo
        enddo
        deallocate(qs_t)
        deallocate(qsb_t)
        deallocate(evap_t)
        deallocate(nqs)
      endif

      call LIS_verify(nf90_close(nid))
    else
      write(LIS_logunit,*) 'Failed to find '//trim(filename)
      call LIS_endrun()
    endif
  endif  

#endif

  !ag - 17Mar2016
  where(LISrunoffdata_struc(n)%qs/=LIS_rc%udef)
    surface_runoff = LISrunoffdata_struc(n)%qs
    baseflow = LISrunoffdata_struc(n)%qsb
    total_evapotranspiration = LISrunoffdata_struc(n)%evap
  else where
    surface_runoff = 0.0
    baseflow = 0.0
    total_evapotranspiration = 0.0
  end where

end subroutine readLISrunoffdata
