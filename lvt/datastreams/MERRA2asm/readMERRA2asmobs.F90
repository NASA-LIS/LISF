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
! !ROUTINE: readMERRA2asmObs
! \label{readMERRA2asmObs}
!
! !INTERFACE: 
subroutine readMERRA2asmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use MERRA2asmobsMod
          
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
!   This plugin processes MERRA2 single-leved diagnostic data available from 
!   NASA GES-DISC.
!   
!  NOTES: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  3 June 2017: Eric Kemp, Initial Specification
! 
!EOP

  real                    :: timenow
  logical                 :: alarmCheck
  integer                 :: c,r, k,nc,nr
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real                    :: gmt
  integer                 :: t
  type(ESMF_Time)         :: merra2time1, merra2time2, initTime
  type(ESMF_TimeInterval) :: dayInterval
  character(len=100)      :: var_name
  real*8                  :: lis_prevtime
  real                    :: t2m(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: qv2m(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: u10m(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: v10m(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: spd10m(LVT_rc%lnc,LVT_rc%lnr)
 
  integer                 :: status

  t2m = 0.0
  qv2m = 0.0
  u10m = 0.0
  v10m = 0.0
  spd10m = 0.0

  nc = MERRA2asmobs(source)%nc
  nr = MERRA2asmobs(source)%nr
  
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) &
       + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(alarmCheck.or.(LVT_rc%dda(source).ne.&
       MERRA2asmobs(source)%da).or.LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     call ESMF_TimeSet(merra2asmobs(source)%starttime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readMerra2asmobs')
     MERRA2asmobs(source)%da = LVT_rc%dda(source)
     
     yr = LVT_rc%dyr(source)
     mo = LVT_rc%dmo(source)
     da = LVT_rc%dda(source)
     hr = LVT_rc%dhr(source)
     mn = LVT_rc%dmn(source)
     ss = LVT_rc%dss(source)
     
     call process_MERRA2asmdata(source, yr, mo, da)
     
  endif
  
  call ESMF_TimeSet(merra2time1,yy=LVT_rc%dyr(source), &
       mm = LVT_rc%dmo(source), &
       dd = LVT_rc%dda(source), &
       h = LVT_rc%dhr(source), &
       m = LVT_rc%dmn(source), &
       calendar = LVT_calendar, &
       rc=status)
  call LVT_verify(status, 'error in timeset: readMERRA2asmobs')
  

  t = nint((merra2time1 - merra2asmobs(source)%starttime)/&
       merra2asmobs(source)%ts) + 1
  
  call aggregate_merra2asmvar(source, t, t2m, &
       merra2asmobs(source)%t2m)
  call aggregate_merra2asmvar(source, t, qv2m, &
       merra2asmobs(source)%qv2m)
  call aggregate_merra2asmvar(source, t, u10m, &
       merra2asmobs(source)%u10m)
  call aggregate_merra2asmvar(source, t, v10m, &
       merra2asmobs(source)%v10m)

  call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,source,t2m,&
       vlevel=1,units="K")
  call LVT_logSingleDataStreamVar(LVT_MOC_QAIRFORC,source,qv2m,&
       vlevel=1,units="kg/kg")
  do r = 1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if ( u10m(c,r) .ne. LVT_rc%udef .and. &
             v10m(c,r) .ne. LVT_rc%udef) then
           spd10m(c,r) = sqrt(u10m(c,r)*u10m(c,r) + v10m(c,r)*v10m(c,r))
        end if
     end do ! c
  end do ! r
  call LVT_logSingleDataStreamVar(LVT_MOC_WINDFORC,source,spd10m,&
       vlevel=1,units="m/s")
   do r=1,LVT_rc%lnr
      do c=1,LVT_rc%lnc
         if(spd10m(c,r).ne.LVT_rc%udef) then 
            spd10m(c,r) = spd10m(c,r)*0.001*86400. ! m/s to km/day
         endif
      enddo ! c
   enddo ! r
  call LVT_logSingleDataStreamVar(LVT_MOC_WINDFORC,source,spd10m,&
       vlevel=1,units="km/day")

end subroutine readMERRA2asmObs

subroutine process_MERRA2asmdata(source, yr, mo, da)

  use LVT_coreMod
  use LVT_logMod
  use MERRA2asmobsMod
          
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  
  integer                :: source
  integer                :: yr
  integer                :: mo
  integer                :: da

  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: t2mid,qv2mid,u10mid,v10mid
  real         :: t2m(merra2asmobs(source)%nc, merra2asmobs(source)%nr,24)
  real         :: qv2m(merra2asmobs(source)%nc, merra2asmobs(source)%nr,24)
  real         :: u10m(merra2asmobs(source)%nc, merra2asmobs(source)%nr,24)
  real         :: v10m(merra2asmobs(source)%nc, merra2asmobs(source)%nr,24)

  integer                :: k,iret
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  call create_MERRA2asm_filename(MERRA2asmobs(source)%odir,&
       yr, mo, da, fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  if(file_exists) then 
     write(LVT_logunit,*) '[INFO] Reading MERRA2 asm file ',trim(fname)
     
     iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn)
     if(iret.eq.0) then 
        
        call LVT_verify(nf90_inq_varid(ftn,"T2M",t2mid),&
             'nf90_inq_varid failed for T2M')
        call LVT_verify(nf90_inq_varid(ftn,"QV2M",qv2mid),&
             'nf90_inq_varid failed for QV2M')
        call LVT_verify(nf90_inq_varid(ftn,"U10M",u10mid),&
             'nf90_inq_varid failed for U10M')
        call LVT_verify(nf90_inq_varid(ftn,"V10M",v10mid),&
             'nf90_inq_varid failed for V10M')
        
        call LVT_verify(nf90_get_var(ftn,t2mid,t2m),&
             'Error in nf90_get_var T2M')
        call LVT_verify(nf90_get_var(ftn,qv2mid,qv2m),&
             'Error in nf90_get_var QV2M')
        call LVT_verify(nf90_get_var(ftn,u10mid,u10m),&
             'Error in nf90_get_var U10M')
        call LVT_verify(nf90_get_var(ftn,v10mid,v10m),&
             'Error in nf90_get_var V10M')

        call LVT_verify(nf90_close(ftn))

        do k=1,24
           call interp_merra2asmvar2d(source,t2m(:,:,k),&
                merra2asmobs(source)%t2m(:,:,k))
           call interp_merra2asmvar2d(source,qv2m(:,:,k),&
                merra2asmobs(source)%qv2m(:,:,k))
           call interp_merra2asmvar2d(source,u10m(:,:,k),&
                merra2asmobs(source)%u10m(:,:,k))
           call interp_merra2asmvar2d(source,v10m(:,:,k),&
                merra2asmobs(source)%v10m(:,:,k))
        enddo
     endif
  else
     write(LVT_logunit,*) '[WARN] Cannot find MERRA2 asm file ', &
          trim(fname)     
  end if
#endif
end subroutine process_MERRA2asmdata

subroutine aggregate_merra2asmvar(source, t, var, merra2var)

  use LVT_coreMod

  integer       :: source
  integer       :: t
  real          :: var(LVT_rc%lnc, LVT_rc%lnr)
  real          :: merra2var(LVT_rc%lnc, LVT_rc%lnr, 24)
  integer       :: c,r

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(var(c,r).ne.LVT_rc%udef) then 
           var(c,r) = merra2var(c,r,t)
        endif
     enddo
  enddo
  
end subroutine aggregate_merra2asmvar



!BOP
!
! !ROUTINE: interp_merra2asmvar2d
! \label{interp_merra2asmvar2d}
!
! !INTERFACE: 
subroutine interp_merra2asmvar2d(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use MERRA2asmobsMod
! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(merra2asmobs(source)%nc,merra2asmobs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the MERRA2 fields to the 
!  target LVT domain
!
!EOP

  real              :: var_inp_1d(merra2asmobs(source)%nc*merra2asmobs(source)%nr)
  logical*1         :: input_bitmap(merra2asmobs(source)%nc*merra2asmobs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r
  integer           :: iret


  nc = merra2asmobs(source)%nc
  nr = merra2asmobs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        if(var_inp(c,r).ne.1e15) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp(c,r) = LVT_rc%udef
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
        endif
     enddo
  enddo
  
  if(LVT_isAtAfinerResolution(merra2asmobs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          merra2asmobs(source)%rlat, & 
          merra2asmobs(source)%rlon, &
          merra2asmobs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          merra2asmobs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo

end subroutine interp_merra2asmvar2d


!BOP
! 
! !ROUTINE: create_MERRA2asm_filename
! \label{create_MERRA2asm_filename}
!
! !INTERFACE: 
subroutine create_MERRA2asm_filename(odir,yr,mo,da,filename)
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
! This routine creates a timestamped filename for the MERRA2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            MERRA2 base directory
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[da]              day of data
!   \item[filename]        Name of the MERRA2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  character*50            :: prefix

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  if (yr==1979 .and. mo>=2) then
     prefix = 'MERRA2_100'
  elseif (yr>1979 .and. yr<=1991) then
     prefix = 'MERRA2_100'
  elseif ( yr >= 1992 .and. yr <= 2000 ) then ! Since 2000 is last full year
     prefix = 'MERRA2_200'
  elseif ( yr >= 2001 .and. yr <= 2009 ) then ! Since 2009 is last full year
     prefix = 'MERRA2_300'
  elseif ( yr >= 2010 ) then
     prefix = 'MERRA2_400'
  else
!     write(LVT_logunit,*) '[ERR] merra2files: date out of range'
!     write(LVT_logunit,*) '[ERR] Supported years are from 1979-2-1 through ...'
!     call LVT_endrun()	
     filename = "none"
  endif
      
  filename = trim(odir)//'/'//trim(prefix)//'/Y'//trim(fyr)//&
       '/M'//trim(fmo)//'/'//trim(prefix)//'.inst1_2d_asm_Nx.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'.nc4'
  
end subroutine create_MERRA2asm_filename


