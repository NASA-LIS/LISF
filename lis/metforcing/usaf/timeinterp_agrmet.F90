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
! !ROUTINE: timeinterp_agrmet
! \label{timeinterp_agrmet}
!
! !REVISION HISTORY:
!   25 Jul 2005: Sujay Kumar; Initial Specification
!   10 Oct 2006: Sujay Kumar: Switched to using ESMF_State for storing
!               forcing data. 
!   30 Apr 2012: Chris Franks: If end-of-run set wt1=0 & wt2=1 for
!                temporal interpolation
!   18 Feb 2020: Eric Kemp: Removed convective precip.
!   16 Feb 2022: K. Arsenault: Commented out "btime" calls, since not used
!
! !INTERFACE:
subroutine timeinterp_agrmet(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_endofrun
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only  : LIS_FORC_Base_State
  use LIS_timeMgrMod, only     : LIS_time2date
  use LIS_logMod,     only     : LIS_verify
  use LIS_constantsMod, only   : LIS_MS2KMDAY
  use LIS_pluginIndices
  use AGRMET_forcingMod, only : agrmet_struc
!EOP
  implicit none
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Precipitation and radiation fields
!  are not temporally interpolated, and the previous hourly value
!  is used. All other variables are linearly interpolated between 
!  the hourly blocks. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!  
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!  \end{description}
!EOP
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8  :: btime
  real    :: wt1,wt2,gmt1,gmt2
  integer :: index1, t
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)

  if(LIS_rc%runmode.eq.LIS_agrmetrunId) then 
     if(LIS_rc%run_model) then 
        btime=agrmet_struc(n)%agrmettime1
        call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
        btime=agrmet_struc(n)%agrmettime2
        call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
        
        call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
             rc=status)
        call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
                
        call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
        call LIS_verify(status,'Error in FieldGet for tmp')
        
        call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
        call LIS_verify(status,'Error in FieldGet for q2')
        
        call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
        call LIS_verify(status,'Error in FieldGet for swd')
        
        call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
        call LIS_verify(status,'Error in FieldGet for lwd')
        
        call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
        call LIS_verify(status,'Error in FieldGet for uwind')
        
        call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
        call LIS_verify(status,'Error in FieldGet for vwind')
        
        call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
        call LIS_verify(status,'Error in FieldGet for psurf')
        
        call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
        call LIS_verify(status,'Error in FieldGet for pcp')
        
        !-----------------------------------------------------------------------
        !  Interpolate Data in Time
        !-----------------------------------------------------------------------
        wt1=(agrmet_struc(n)%agrmettime2-LIS_rc%time)/ & 
             (agrmet_struc(n)%agrmettime2-agrmet_struc(n)%agrmettime1)
        wt2=1.0-wt1
        
        do t=1,LIS_rc%ntiles(n)
           index1    = LIS_domain(n)%tile(t)%index
           if(agrmet_struc(n)%metdata2(1,index1).ne.-9999.0) then 
              tmp(t)   = agrmet_struc(n)%metdata2(1,index1)
           else
              tmp(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(2,index1).ne.-9999.0) then 
              q2(t)    = agrmet_struc(n)%metdata2(2,index1)
           else
              q2(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(3,index1).ne.-9999.0) then 
              swd(t)   = agrmet_struc(n)%metdata2(3,index1)
           else
              swd(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(4,index1).ne.-9999.0) then 
              lwd(t)   = agrmet_struc(n)%metdata2(4,index1)
           else
              lwd(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(5,index1).ne.-9999.0) then 
              uwind(t) = agrmet_struc(n)%metdata2(5,index1)
           else
              uwind(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(6,index1).ne.-9999.0) then            
              vwind(t) = agrmet_struc(n)%metdata2(6,index1)
           else
              vwind(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(7,index1).ne.-9999.0) then 
              psurf(t) = agrmet_struc(n)%metdata2(7,index1)
           else
              psurf(t) = -9999.0
           endif
           if(agrmet_struc(n)%metdata2(8,index1).ne.-9999.0) then 
              pcp(t)  = agrmet_struc(n)%metdata2(8,index1)           
           else
              pcp(t) = -9999.0
           endif
        enddo
     endif

  ! Retrospective USAF/NAFPA forcing (offline):
  else
      ! Note below code is not used in this routine:
!     btime=agrmet_struc(n)%agrmettime1
!     call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
!     btime=agrmet_struc(n)%agrmettime2
!     call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')
     
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')
        
     call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
     call LIS_verify(status)
     
     call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
     call LIS_verify(status)
     
     !-----------------------------------------------------------------------
     !  Interpolate Data in Time (the order of variables below is important)
     !-----------------------------------------------------------------------
     if (LIS_endofrun()) then
        wt1 = 0
        wt2 = 1
     else
        wt1=(agrmet_struc(n)%agrmettime2-LIS_rc%time)/ & 
            (agrmet_struc(n)%agrmettime2-agrmet_struc(n)%agrmettime1)
        wt2=1.0-wt1
     endif
     
     do t=1,LIS_rc%ntiles(n)
        index1    = LIS_domain(n)%tile(t)%index
        if(agrmet_struc(n)%metdata1(1,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(1,index1).ne.-9999.0) then 
           tmp(t)   = wt1*agrmet_struc(n)%metdata1(1,index1)+&
                wt2*agrmet_struc(n)%metdata2(1,index1)
         ! Kludge fix for some bad data points over the high lat. regions
           if(tmp(t).gt.360.) then  ! In Kelvins (~87 deg C)
              tmp(t) = 273.1        ! Assign to Tfrz
           endif
        else
           tmp(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(2,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(2,index1).ne.-9999.0) then 
           q2(t)    = wt1*agrmet_struc(n)%metdata1(2,index1)+&
                wt2*agrmet_struc(n)%metdata2(2,index1)
        else
           q2(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(3,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(3,index1).ne.-9999.0) then            
           swd(t)   = wt1*agrmet_struc(n)%metdata1(3,index1)+&
                wt2*agrmet_struc(n)%metdata2(3,index1)
        else
           swd(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(4,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(4,index1).ne.-9999.0) then 
           lwd(t)   = wt1*agrmet_struc(n)%metdata1(4,index1)+&
                wt2*agrmet_struc(n)%metdata2(4,index1)
        else
           lwd(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(5,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(5,index1).ne.-9999.0) then 
           uwind(t) = (wt1*agrmet_struc(n)%metdata1(5,index1)+&
                wt2*agrmet_struc(n)%metdata2(5,index1))/LIS_MS2KMDAY
           vwind(t) = 0.0
        else
           uwind(t) = -9999.0
           vwind(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(6,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(6,index1).ne.-9999.0) then 
           psurf(t) = wt1*agrmet_struc(n)%metdata1(6,index1)+&
                wt2*agrmet_struc(n)%metdata2(6,index1)
        else
           psurf(t) = -9999.0
        endif
        if(agrmet_struc(n)%metdata1(7,index1).ne.-9999.0.and.&
             agrmet_struc(n)%metdata2(7,index1).ne.-9999.0) then 
           pcp(t)  = agrmet_struc(n)%metdata2(7,index1)/(3*3600) !Convert accum input to rate
        else
           pcp(t) = -9999.0
        endif
     enddo
  endif

end subroutine timeinterp_agrmet

