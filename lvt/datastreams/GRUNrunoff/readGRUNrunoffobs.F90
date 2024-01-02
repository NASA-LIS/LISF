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
! !ROUTINE: readGRUNrunoffobs
! \label{readGRUNrunoffobs}
!
! !INTERFACE: 
subroutine readGRUNrunoffobs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use LVT_timeMgrMod
  use GRUNrunoff_obsMod
          
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
!   This plugin processes the GRUN runoff data at a monthly
!   interval
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  30 Jul 2022   Sujay Kumar  Initial Specification
! 
!EOP

  real                    :: timenow
  logical                 :: alarmCheck
  integer                 :: c,r, k,nc,nr
  integer                 :: yr, mo, da, hr, mn, ss, doy
  real                    :: gmt
  integer                 :: t
  integer                 :: iret
  integer                 :: runoffid
  integer                 :: ftn
  type(ESMF_Time)         :: gruntime
  type(ESMF_TimeInterval) :: offset
  integer                 :: off_month
  character(len=100)      :: var_name,fname
  logical                 :: file_exists
  real*8                  :: lis_prevtime
  real                    :: runoff(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: runoff_in(GRUNobs(source)%nc,GRUNobs(source)%nr)
  integer                 :: status

  runoff = LVT_rc%udef
  
  nc = GRUNobs(source)%nc
  nr = GRUNobs(source)%nr

  if(GRUNobs(source)%mo.ne.LVT_rc%dmo(source).or.&
       GRUNobs(source)%startFlag) then 
     GRUNobs(source)%startFlag = .false. 
     GRUNobs(source)%mo = LVT_rc%dmo(source)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  
     inquire(file=trim(GRUNobs(source)%odir),exist=file_exists) 

     fname = GRUNobs(source)%odir
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading GRUN file ',trim(fname)
        
        call ESMF_TimeSet(gruntime,yy=LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h = LVT_rc%dhr(source), &
             m = LVT_rc%dmn(source), &
             calendar = LVT_calendar, &
             rc=status)
        call LVT_verify(status, 'error in timeset: readGRUNrunoffobs')
       
        off_month = (LVT_rc%dyr(source)-1902)*12+LVT_rc%dmo(source)
        
        iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
             ncid = ftn)
        if(iret.eq.0) then 

           call LVT_verify(nf90_inq_varid(ftn,"Runoff",runoffid),&
                'nf90_inq_varid failed for Runoff')
           call LVT_verify(nf90_get_var(ftn,runoffid,runoff_in,&
                start = (/1,1,off_month/),&
                count = (/GRUNobs(source)%nc,GRUNobs(source)%nr,1/)),&
                'Error in nf90_get_var Runoff')
           
           call LVT_verify(nf90_close(ftn))

           call interp_grunndata(source, runoff_in, runoff)           

        endif
     endif
  endif
#endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(runoff(c,r).ne.-9999.0) then
           runoff(c,r) = runoff(c,r)/86400.0 !to kg/m2s
        endif
     enddo
  enddo
    
  call LVT_logSingleDataStreamVar(LVT_MOC_RUNOFF,source, runoff,&
       vlevel=1,units="kg/m2s")

 end subroutine readGRUNrunoffobs



!BOP
!
! !ROUTINE: interp_grunndata
! \label{interp_grunndata}
!
! !INTERFACE: 
subroutine interp_grunndata(source, var_inp,var_out)
! !USES: 
  use LVT_coreMod
  use GRUNrunoff_obsMod
  use LVT_logMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)       
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer           :: source
  real              :: var_inp(GRUNobs(source)%nc,GRUNobs(source)%nr)
  real              :: var_out(LVT_rc%lnc,LVT_rc%lnr)
! 
! !DESCRIPTION: 
!  This routine interpolates/upscales the GRUN runoff data to the 
!  target LVT domain
!
!EOP
  real              :: var_inp_1d(GRUNobs(source)%nc*GRUNobs(source)%nr)
  logical*1         :: input_bitmap(GRUNobs(source)%nc*GRUNobs(source)%nr)
  real              :: var_out_1d(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: output_bitmap(LVT_rc%lnc*LVT_rc%lnr)
  integer           :: nc, nr, c,r,k
  real              :: pcp_tmp
  integer           :: iret

  nc = GRUNobs(source)%nc
  nr = GRUNobs(source)%nr
  
  input_bitmap = .false. 
  do r=1,nr
     do c=1,nc
        k = c+(r-1)*GRUNobs(source)%nc
        if(var_inp(c,r).ge.0) then 
           var_inp_1d(c+(r-1)*nc) = var_inp(c,r)
           input_bitmap(c+(r-1)*nc) = .true. 
        else
           var_inp_1d(c+(r-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo
  
  if(LVT_isAtAfinerResolution(GRUNobs(source)%datares)) then
     call neighbor_interp(LVT_rc%gridDesc, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d, &
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          GRUNobs(source)%rlat, & 
          GRUNobs(source)%rlon, &
          GRUNobs(source)%n11, &
          LVT_rc%udef, iret)
     
  else
     call upscaleByAveraging(&
          nc*nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          GRUNobs(source)%n11, input_bitmap, &
          var_inp_1d, output_bitmap, var_out_1d)
     
  endif
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(output_bitmap(c+(r-1)*LVT_rc%lnc)) then 
           var_out(c,r) = var_out_1d(c+(r-1)*LVT_rc%lnc)
        endif
     enddo
  enddo
  
end subroutine interp_grunndata






