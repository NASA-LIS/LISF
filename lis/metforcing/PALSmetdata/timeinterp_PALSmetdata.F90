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
! !ROUTINE: timeinterp_PALSmetdata
! \label{timeinterp_PALSmetdata}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
!
! !INTERFACE:
subroutine timeinterp_PALSmetdata(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_constantsMod, only  : LIS_CONST_SOLAR
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod, only : LIS_tick, LIS_time2date
  use LIS_logMod, only :LIS_logunit, LIS_verify, LIS_endrun
  use PALSmetdata_forcingMod, only : PALSmetdata_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the PALS forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous 3 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 3 hourly blocks. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    advances or retracts time by the specified amount
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: zdoy
  real :: zw1, zw2
  real :: czm, cze, czb
  real :: wt1, wt2,swt1,swt2
  real :: gmt1, gmt2, tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)

  btime=PALSmetdata_struc(n)%fcsttime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  
  tempbdoy=bdoy
  tempgmt1=gmt1
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=PALSmetdata_struc(n)%fcsttime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt2=gmt2
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(PALSmetdata_struc(n)%fcsttime2-LIS_rc%time)/ & 
       (PALSmetdata_struc(n)%fcsttime2-PALSmetdata_struc(n)%fcsttime1)
  wt2=1.0-wt1
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1


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

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     zdoy=LIS_rc%doy
     !  compute and apply zenith angle weights
     call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
          LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
     
     if(PALSmetdata_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) =  PALSmetdata_struc(n)%metdata2(3,index1) ! * zw1 !Grey Nearing
        if (swd(t) < 0) then
           write(LIS_logunit,*) '2 warning!!!  SW radiation is negative!!'
           write(LIS_logunit,*) 'sw=', swd(t), '... negative'
           write(LIS_logunit,*) 'GEOS5=', PALSmetdata_struc(n)%metdata2(3,index1)
           call LIS_endrun
        end if
        
        if (swd(t).gt.LIS_CONST_SOLAR) then
           swd(t)=PALSmetdata_struc(n)%metdata2(3,index1)
        endif
     endif
  enddo

  !do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if(PALSmetdata_struc(n)%metdata2(8,index1).ne.LIS_rc%udef) then 
        pcp(t)=PALSmetdata_struc(n)%metdata2(8,index1)	
        pcp(t)  = pcp(t)
     endif
  enddo

  !linearly interpolate everything else

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(1,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then 
        tmp(t) =PALSmetdata_struc(n)%metdata1(1,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(1,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(2,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then 
        q2(t) =PALSmetdata_struc(n)%metdata1(2,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(2,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(4,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(4,index1).ne.LIS_rc%udef)) then 
        lwd(t) =PALSmetdata_struc(n)%metdata1(4,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(4,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(5,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(5,index1).ne.LIS_rc%udef)) then 
        uwind(t) =PALSmetdata_struc(n)%metdata1(5,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(5,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(6,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(6,index1).ne.LIS_rc%udef)) then 
        vwind(t) =PALSmetdata_struc(n)%metdata1(6,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(6,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((PALSmetdata_struc(n)%metdata1(7,index1).ne.LIS_rc%udef).and.&
          (PALSmetdata_struc(n)%metdata2(7,index1).ne.LIS_rc%udef)) then 
        psurf(t) =PALSmetdata_struc(n)%metdata1(7,index1)*wt1+ & 
             PALSmetdata_struc(n)%metdata2(7,index1)*wt2
     endif
  enddo

end subroutine timeinterp_PALSmetdata
 
