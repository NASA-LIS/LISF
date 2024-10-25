!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_gfs
! \label{timeinterp_gfs}
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!  10 Oct 2006: Sujay Kumar: Switched to using ESMF_State for storing
!               forcing data. 
! !INTERFACE:
subroutine timeinterp_gfs(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only  : LIS_rc, LIS_domain
  use LIS_FORC_AttributesMod 
  use LIS_timeMgrMod,     only  : LIS_time2date, LIS_tick
  use LIS_constantsMod,   only  : LIS_CONST_SOLAR
  use LIS_metforcingMod, only  : LIS_forc, LIS_FORC_Base_State
  use LIS_logMod,         only  : LIS_logunit, LIS_verify, LIS_endrun
  use gfs_forcingMod,    only  : gfs_struc

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
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
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
  integer :: t,zdoy,idoy,iyr,imo,ida,ihr,imn,its,iss
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime,inittime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2,igmt
  real :: zw1,zw2
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  btime=gfs_struc(n)%gfstime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gfs_struc(n)%gfstime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  
  wt1 = (gfs_struc(n)%gfstime2-LIS_rc%time) / & 
       (gfs_struc(n)%gfstime2-gfs_struc(n)%gfstime1)
  wt2 = 1.0 - wt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')

  
  call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
  call LIS_verify(status)
        
  call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(tmpField,localDE=0, farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  pcp  = 0.0
  cpcp = 0.0
  
  do t = 1, LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy = LIS_rc%doy
     call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
          LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
     swd(t) = zw1 * gfs_struc(n)%metdata2(3,index1)
     if (swd(t) < 0) then
        write(LIS_logunit,*) '2 warning!!!  SW radiation is negative!!'
        write(LIS_logunit,*) 'sw=', swd(t), '... negative'
        write(LIS_logunit,*) 'gfs2=', gfs_struc(n)%metdata2(3,index1)
        call LIS_endrun
     end if
     
     if (swd(t).gt.LIS_CONST_SOLAR) then
        swd(t)=gfs_struc(n)%metdata2(3,index1)
     endif
  end do
  do t = 1, LIS_rc%ntiles(n)     
     index1 = LIS_domain(n)%tile(t)%index
     lwd(t) = gfs_struc(n)%metdata2(4,index1)     
     pcp(t) = gfs_struc(n)%metdata2(8,index1)    
     cpcp(t) = gfs_struc(n)%metdata2(9,index1)     
     tmp(t) = wt1 * gfs_struc(n)%metdata1(1,index1) & 
          + wt2 *gfs_struc(n)%metdata2(1,index1)
     q2(t) = wt1 * gfs_struc(n)%metdata1(2,index1) & 
          + wt2 *gfs_struc(n)%metdata2(2,index1)
     uwind(t) = wt1 * gfs_struc(n)%metdata1(5,index1) & 
          + wt2 *gfs_struc(n)%metdata2(5,index1)
     vwind(t) = wt1 * gfs_struc(n)%metdata1(6,index1) & 
          + wt2 *gfs_struc(n)%metdata2(6,index1)
     psurf(t) = wt1 * gfs_struc(n)%metdata1(7,index1) & 
          + wt2 *gfs_struc(n)%metdata2(7,index1)
  end do

end subroutine timeinterp_gfs
