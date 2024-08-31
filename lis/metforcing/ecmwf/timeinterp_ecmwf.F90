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
! !ROUTINE: timeinterp_ecmwf
! \label{timeinterp_ecmwf}
!
! !REVISION HISTORY:
!   5 Nov 2001: Urszula Jambor; Initial Specification
!  10 Oct 2006: Sujay Kumar: Switched to using ESMF_State for storing
!               forcing data. 
! !INTERFACE:
subroutine timeinterp_ecmwf(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,            only : LIS_rc, LIS_domain
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,      only : LIS_forc, LIS_FORC_Base_State
  use LIS_constantsMod,       only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,         only : LIS_time2date, LIS_tick
  use LIS_logMod,             only : LIS_logunit, LIS_verify
  use ecmwf_forcingMod,       only : ecmwf_struc
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
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
  integer :: t,zdoy,idoy,iyr,imo,ida,ihr,imn,iss
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8  :: btime,inittime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2,igmt
  real    :: zw1,zw2,its
  integer :: tid,tid1,tid2
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  !=== Reset GMT times
  btime=ecmwf_struc(n)%ecmwftime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=ecmwf_struc(n)%ecmwftime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  !=== Need lower SW time boundary for call to zterp based on time2
  inittime=btime
  call LIS_time2date( inittime, idoy, igmt, iyr, imo, ida, ihr, imn )
  select case (ihr)
  case(00,12,24)
     its = -12*60*60
  case(03,15)
     its =  -3*60*60
  case(06,18)
     its =  -6*60*60
  case(09,21)
     its =  -9*60*60
  end select
  iss = 0 
  call LIS_tick( inittime, idoy, igmt, iyr, imo, ida, ihr, imn, iss, its )

  !===  Interpolate Data in Time
  wt1=(ecmwf_struc(n)%ecmwftime2-LIS_rc%time)/ & 
       (ecmwf_struc(n)%ecmwftime2-ecmwf_struc(n)%ecmwftime1)
  wt2=1.0-wt1

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

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')


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

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

!only the particular ensemble member is initialized. 
  do t= 1, LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy = LIS_rc%doy
     call zterp(0,LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon, &
          igmt,gmt2,LIS_rc%gmt,zdoy, &
          zw1,zw2,czb,cze,czm,LIS_rc)
     swd(t) = zw1 * ecmwf_struc(n)%metdata2(3,index1)
     if ((swd(t).ne.LIS_rc%udef).and.  &
          (swd(t).lt.0)                ) then
        if (swd(t) > -0.00001) then !arbitrary!!
           swd(t) = 0.0             !threshold!!
        else
           write(LIS_logunit,*) 'Stopping because SW forcing not udef but lt 0, '
           write(LIS_logunit,*) t,swd(t)
!           stop
        endif
     endif
     if (swd(t).gt.LIS_CONST_SOLAR) then
        swd(t) = ecmwf_struc(n)%metdata2(3,index1)
     endif
  enddo
  do t = 1,LIS_rc%ntiles(n)
     index1   = LIS_domain(n)%tile(t)%index
     lwd(t)   = ecmwf_struc(n)%metdata2(4,index1)
     pcp(t)   = ecmwf_struc(n)%metdata2(8,index1)
     cpcp(t)  = ecmwf_struc(n)%metdata2(9,index1)
     tmp(t)   = ecmwf_struc(n)%metdata1(1,index1)*wt1 + &
                ecmwf_struc(n)%metdata2(1,index1)*wt2
     q2(t)    = ecmwf_struc(n)%metdata1(2,index1)*wt1 + &
                ecmwf_struc(n)%metdata2(2,index1)*wt2
     uwind(t) = ecmwf_struc(n)%metdata1(5,index1)*wt1 + &
                ecmwf_struc(n)%metdata2(5,index1)*wt2
     vwind(t) = ecmwf_struc(n)%metdata1(6,index1)*wt1 + &
                ecmwf_struc(n)%metdata2(6,index1)*wt2
     psurf(t) = ecmwf_struc(n)%metdata1(7,index1)*wt1 + &
                ecmwf_struc(n)%metdata2(7,index1)*wt2
  enddo
end subroutine timeinterp_ecmwf
