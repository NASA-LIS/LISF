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
! !REVISION HISTORY: 
! 7 Jan 2016: Sujay Kumar, Initial implementation
! 
! !USES: 
subroutine readERAILandrunoffdata(n,surface_runoff, baseflow)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use ERAILandrunoffdataMod
  use LIS_fileIOMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none

! !ARGUMENTS: 
  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%gnc(n),LIS_rc%gnr(n))
  real                         :: baseflow(LIS_rc%gnc(n),LIS_rc%gnr(n))

!
! !DESCRIPTION: 
!
!  This subroutine reads the runoff fields from the ERA interim land data, 
!  subsets to the current model time and conducts the geospatial transform
!  to the LIS grid. 
!  
!  The arguments are: 
!  \begin{description}
!   \item[n]               index of the nest
!   \item[surface_runoff]  surface runoff field generated from the ERA interim
!                          land data
!   \item[baseflow]        baseflow field generated from the ERA interim land
!                          data 
!  \end{description}
!EOP 

  integer                       :: nc,nr
  integer                       :: c,r,t
  real,   allocatable           :: qs(:,:)
  real,   allocatable           :: qsb(:,:)
  integer                       :: ios
  integer                       :: ftn
  integer                       :: qsid, qsbid
  integer                       :: tindex
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time1
  real                          :: gmt
  logical                       :: file_exists
  type(ESMF_TimeInterval)       :: deltaT
  integer                       :: time(1), timeid
  type(ESMF_Time)               :: currTime, startYrTime
  integer                       :: hr_elapsed
  real                          :: scale_qs, scale_qsb
  real                          :: offset_qs, offset_qsb

  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr =ERAILandrunoffdata_struc(n)%outInterval*((LIS_rc%hr)/&
       ERAILandrunoffdata_struc(n)%outInterval)
  mn =0
  ss =0
  ts =ERAILandrunoffdata_struc(n)%outInterval

  call LIS_tick(time1,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_ERAILand_filename(ERAILandrunoffdata_struc(n)%odir,&
       yr, filename)

  surface_runoff = 0.0
  baseflow       = 0.0

  nc = ERAILandrunoffdata_struc(n)%nc
  nr = ERAILandrunoffdata_struc(n)%nr

  tindex = hr+1

  allocate(qs(nc,nr))
  allocate(qsb(nc,nr))
  
  qs = LIS_rc%udef
  qsb = LIS_rc%udef
  

  inquire(file=filename, exist=file_exists)
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading '//trim(filename)

       ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(ios.ne.0) then
        write(LIS_logunit,*) '[ERR] Opening ERAILand file ',trim(filename),&
             'failed'
        call LIS_endrun()
     endif

     call ESMF_TimeSet(currTime, yy=yr, &
          mm = mo, dd = da, h = 0, m = 0, calendar = LIS_calendar, &
          rc=ios)
     call ESMF_TimeSet(startYrTime, yy=yr, &
          mm = 1, dd = 2, h = 0 , m = 0, calendar=LIS_calendar, &
          rc=ios)

     deltaT = startYrTime - ERAILandrunoffdata_struc(n)%startTime
     call ESMF_TimeIntervalGet(deltaT, h=hr_elapsed)

     call LIS_verify(nf90_inq_varid(ftn,"time",timeid),&
          'nf90_inq_varid failed for time')
     call LIS_verify(nf90_get_var(ftn,timeid, time, &
          start=(/1/), count=(/1/)),&
          'Error in nf90_get_var time')
     if(hr_elapsed.ne.time(1)) then 
        write(LIS_logunit,*) '[ERR] The starting time in the ERA interim fcst file'
        write(LIS_logunit,*) '[ERR] do not match the computed starting time'
        call LIS_endrun()
     endif

     call ESMF_TimeSet(startYrTime, yy=yr, &
          mm = 1, dd = 1, h = 0 , m = 0, calendar=LIS_calendar, &
          rc=ios)

     deltaT = currTime - startYrTime
     call ESMF_TimeIntervalGet(deltaT, h=tindex)
     
     if(mod(tindex,24).eq.0) then 
        write(LIS_logunit,*) '[INFO] Reading ERA Interim Land file ',&
             trim(filename)
        tindex = tindex/24 + 1

        call LIS_verify(nf90_inq_varid(ftn,"sro",qsid),&
             'nf90_inq_varid failed for sro')
        call LIS_verify(nf90_inq_varid(ftn,"ssro",qsbid),&
             'nf90_inq_varid failed for ssro')        
     
        call LIS_verify(nf90_get_var(ftn,qsid,qs,&
             start=(/1,1,tindex/),count=(/nc,nr,1/)),&
             'Error in nf90_get_var sro')
        call LIS_verify(nf90_get_var(ftn,qsid,qsb,&
             start=(/1,1,tindex/),count=(/nc,nr,1/)),&
             'Error in nf90_get_var ssro')
     
        call LIS_verify(nf90_get_att(ftn,qsid, 'scale_factor',&
             scale_qs), &
             'Error in nf90_get_att scale_factor')
        call LIS_verify(nf90_get_att(ftn,qsid, 'add_offset',&
             offset_qs), &
             'Error in nf90_get_att add_offset')

        call LIS_verify(nf90_get_att(ftn,qsbid, 'scale_factor',&
             scale_qsb), &
             'Error in nf90_get_att scale_factor')
        call LIS_verify(nf90_get_att(ftn,qsbid, 'add_offset',&
             offset_qsb), &
             'Error in nf90_get_att add_offset')

     endif
     call LIS_verify(nf90_close(ftn))
  endif

  call interp_ERAILandrunoffdata(n,nc,nr,qs,scale_qs,offset_qs, surface_runoff)
  call interp_ERAILandrunoffdata(n,nc,nr,qsb,scale_qsb, offset_qsb, baseflow)

  
  deallocate(qs)
  deallocate(qsb)

    
end subroutine readERAILandrunoffdata

!BOP
! 
! !ROUTINE: create_ERAILand_filename
! \label{create_ERAILand_filename}
!
! !INTERFACE: 
subroutine create_ERAILand_filename(odir, yr, filename)

  use LIS_logMod

! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the ERA Interim Land 
! runoff data based on the given date 
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            ERA interim land base directory
!   \item[yr]              year of data
!   \item[filename]        Name of the ERA inteirm land file
!  \end{description}
! 
!EOP

  character*4             :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/era-interim_stream_fcst2_'//trim(fyr)//'.nc'
  

end subroutine create_ERAILand_filename




!BOP
!
! !ROUTINE: interp_ERAILandrunoffdata
!  \label{interp_ERAILandrunoffdata}
!
! !INTERFACE:
  subroutine interp_ERAILandrunoffdata(n, nc,nr,var_input,scale_f, offset_f,&
       var_output)
!
! !USES:
    use LIS_coreMod
    use LIS_constantsMod
    use ERAILandrunoffdataMod
      
    implicit none
! !ARGUMENTS:
    integer            :: n
    integer            :: nc
    integer            :: nr
    real               :: var_input(nc,nr)
    real               :: scale_f, offset_f
    real               :: var_output(LIS_rc%lnc(n), LIS_rc%lnr(n))
!
! !DESCRIPTION:
!   This subroutine spatially transforms the ERA interim land variable to the
!   model (LIS) output grid and resolution, using a neighbor interpolation
!   or upscaling approach. 
!
!   The arguments are:
!   \begin{description}
!    \item[n]       index of the nest
!    \item[nc]      number of columns in the input grid
!    \item[nr]      number of rows in the input grid
!    \item[var_input] input variable to be transformed
!    \item[scale_f]   scale factor associated with the variable
!    \item[offset_f]   offset factor associated with the variable
!    \item[var_output] resulting transformed field
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
    
    real               :: var_input_1d(nc*nr)
    logical*1          :: lb(nc*nr)
    integer            :: ios
    integer            :: c,r,c1
    logical*1          :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    real               :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

    var_output = LIS_rc%udef
    lb = .false.
    do r = 1,nr
       do c = 1,nc
          if(c.ge.1.and.c.le.240) then 
             c1 = c + 240
          else
             c1 = c - 240
          endif
          
          var_input(c,r) = var_input(c,r)*scale_f + offset_f
          if (var_input(c,r).ge.0) then
             lb(c1+((nr-r+1)-1)*nc) = .true.
             var_input_1d(c1+((nr-r+1)-1)*nc) = var_input(c,r) 
          else
             var_input_1d(c1+((nr-r+1)-1)*nc) = LIS_rc%udef
          endif
       enddo
    enddo
    
    if(LIS_isAtAfinerResolution(n,0.75)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input_1d,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            ERAILandrunoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            ERAILandrunoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))*LIS_CONST_RHOFW/&
               (24*3600.0)
       enddo
    enddo

  end subroutine interp_ERAILandrunoffdata


