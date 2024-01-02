!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readERA5Obs
! \label{readERA5Obs}
!
! !INTERFACE: 
subroutine readERA5Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use ERA5obsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This plugin processes the ER5 surface meteorology data
!   (currently only precip and near surface air temperature 
!    variables are supported)
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  5 Dec 2020: Sujay Kumar, Initial Specification
! 
!EOP

  real                    :: timenow
  logical                 :: alarmCheck
  integer                 :: c,r, k,nc,nr
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real                    :: gmt
  integer                 :: t
  type(ESMF_Time)         :: era5time1, era5time2, initTime
  type(ESMF_TimeInterval) :: dayInterval
  character(len=100)      :: var_name
  real*8                  :: lis_prevtime
  real                    :: prcp(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: tair(LVT_rc%lnc,LVT_rc%lnr)
  integer                 :: status

  prcp = 0.0
  tair = 0.0

  nc = ERA5obs(source)%nc
  nr = ERA5obs(source)%nr
  
  if(ERA5obs(source)%mo.ne.LVT_rc%dmo(source)) then 

     ERA5obs(source)%mo = LVT_rc%dmo(source)

     call ESMF_TimeSet(era5obs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readEra5obs')

     yr = LVT_rc%dyr(source)
     mo = LVT_rc%dmo(source)
     da = LVT_rc%dda(source)
     hr = LVT_rc%dhr(source)
     mn = LVT_rc%dmn(source)
     ss = LVT_rc%dss(source)
     
     call process_ERA5data(source, yr, mo, da)
     
  endif
  
  call ESMF_TimeSet(era5time1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readERA5obs')
  

  t = nint((era5time1 - era5obs(source)%starttime)/&
       era5obs(source)%ts) + 1
  
  call select_timeslice_era5var(source, t, prcp,  &
       era5obs(source)%prcp)
  call select_timeslice_era5var(source, t, tair,  &
       era5obs(source)%tair)

  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source, prcp,&
       vlevel=1,units="kg/m2s")
  call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source, tair,&
       vlevel=1,units="K")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(prcp(c,r).ne.LVT_rc%udef) then 
           prcp(c,r) = prcp(c,r) * 3600.0
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_TOTALPRECIP,source, prcp,&
       vlevel=1,units="kg/m2")

end subroutine readERA5Obs

!BOP
!
! !ROUTINE: process_ERA5data
! \label{process_ERA5data}
!
! !INTERFACE: 
subroutine process_ERA5data(source, yr, mo, da)
! !USES: 
  use LVT_coreMod
  use LVT_logMod
  use ERA5obsMod
! 
! !DESCRIPTION: 
!  This subroutine reads and spatially interpolates
!  data corresponding to one ERA5 datafile
!EOP
        
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  integer                :: yr
  integer                :: mo
  integer                :: da

  integer                :: ftn, ftn_flx
  character*100          :: fname
  logical                :: file_exists
  integer                :: prcpid, tairid
  real, allocatable      :: prcp(:,:)
  real, allocatable      :: tair(:,:)

  integer                :: rec_size 
  integer                :: k,iret
  integer                :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/


#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  call create_ERA5_filename(ERA5obs(source)%odir,&
       yr, mo, da, fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading ERA5 file ',trim(fname)
     
    if((mod(yr,4) .eq. 0 .and. mod(yr, 100).ne.0) &!leap year
          .or.(mod(yr,400) .eq.0)) then 
        days(2) = 29
     else 
        days(2) = 28
     endif

     rec_size = days(mo)*24+1

     era5obs(source)%prcp = LVT_rc%udef
     era5obs(source)%tair = LVT_rc%udef
     allocate(prcp(era5obs(source)%npts,rec_size))
     allocate(tair(era5obs(source)%npts,rec_size))
     
     prcp = LVT_rc%udef
     tair = LVT_rc%udef

     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 

        call LVT_verify(nf90_inq_varid(ftn,"Rainf",prcpid),&
             'nf90_inq_varid failed for Rainf')
        call LVT_verify(nf90_inq_varid(ftn,"Tair",tairid),&
             'nf90_inq_varid failed for Tair')
        
        call LVT_verify(nf90_get_var(ftn,prcpid,prcp),&
             'Error in nf90_get_var Rainf')
        call LVT_verify(nf90_get_var(ftn,tairid,tair),&
             'Error in nf90_get_var TSURF')
        
        call LVT_verify(nf90_close(ftn))

        do k=1,rec_size
           call interp_era5var2d(source,prcp(:,k),&
                era5obs(source)%prcp(:,:,k))
           call interp_era5var2d(source,tair(:,k),&
                era5obs(source)%tair(:,:,k))
        enddo
     endif
     deallocate(prcp)
     deallocate(tair)
  end if


#endif
end subroutine process_ERA5data

!BOP
!
! !ROUTINE: select_timeslice_era5var
! \label{select_timeslice_era5var}
!
! !INTERFACE: 
subroutine select_timeslice_era5var(source, t, var, era5var)
! !USES: 
  use LVT_coreMod
  use ERA5obsMod
! 
! !DESCRIPTION: 
!  This routine selects the ERA5 data corresponding to 
!  a particular instance. 
!EOP

  integer       :: source
  integer       :: t
  real          :: var(LVT_rc%lnc, LVT_rc%lnr)
  real          :: era5var(LVT_rc%lnc, LVT_rc%lnr, era5obs(source)%ntimes)
  integer       :: c,r

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var(c,r).ne.LVT_rc%udef) then 
           var(c,r) = era5var(c,r,t)
        endif
     enddo
  enddo
  
end subroutine select_timeslice_era5var


!BOP
!
! !ROUTINE: interp_era5var2d
! \label{interp_era5var2d}
!
! !INTERFACE: 
subroutine interp_era5var2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use ERA5obsMod
  use LVT_logMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)       
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(era5obs(source)%npts)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the ERA5 fields to the 
!  target LVT domain
!
!EOP
  real              :: var_inp_1d(era5obs(source)%nc*era5obs(source)%nr)
  logical*1         :: input_bitmap(era5obs(source)%nc*era5obs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r,k
  real              :: pcp_tmp
  integer           :: iret

  nc = era5obs(source)%nc
  nr = era5obs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        k = c+(r-1)*era5obs(source)%nc
        if(era5obs(source)%g2p(c,r).gt.0) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(era5obs(source)%g2p(c,r))
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp_1d(c+(r-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo
  
    if(LVT_isAtAfinerResolution(era5obs(source)%datares)) then
       call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d, &
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, &
            era5obs(source)%rlat, & 
            era5obs(source)%rlon, &
            era5obs(source)%n11, &
            LVT_rc%udef, iret)
       
    else
       call upscaleByAveraging(&
            nc*nr, &
            LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
            era5obs(source)%n11, input_bitmap, &
            var_inp_1d, output_bitmap, var_out_1d)
       
    endif

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
             var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
          endif
       enddo
    enddo
    
     
end subroutine interp_era5var2d


!BOP
! 
! !ROUTINE: create_ERA5_filename
! \label{create_ERA5_filename}
!
! !INTERFACE: 
subroutine create_ERA5_filename(odir,yr,mo,da,filename)
! 
! !USES:   
  use LVT_logMod

  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the ERA5 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            ERA5 base directory
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[da]              day of data
!   \item[filename]        Name of the ERA5 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo

  filename = trim(odir)//'/FORCING_'//trim(fyr)//trim(fmo)//'.nc'
 
end subroutine create_ERA5_filename


