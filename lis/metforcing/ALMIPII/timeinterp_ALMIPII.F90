!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_ALMIPII
! \label{timeinterp_ALMIPII}
!
! !REVISION HISTORY:
!
! 02Feb2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine timeinterp_ALMIPII(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_constantsMod, only  : LIS_CONST_SOLAR
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use LIS_timeMgrMod, only : LIS_calendar
  use LIS_logMod, only :LIS_logunit, LIS_verify
  use ALMIPII_forcingMod, only : ALMIPII_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!
!EOP
  integer :: zdoy
  real :: zw1, zw2
  real :: czm, cze, czb
  real :: wt1, wt2
  real :: gmt1, gmt2
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss, tempbts
  integer            :: status
  type(ESMF_Time)    :: currTime
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),snowf(:)
  

  call ESMF_TimeSet(currTime, yy = LIS_rc%yr, &
       mm = LIS_rc%mo, &
       dd = LIS_rc%da, &
       h  = LIS_rc%hr, &
       m  = LIS_rc%mn, & 
       s  = LIS_rc%ss, &
       calendar = LIS_calendar, & 
       rc = status)  
!=== Interpolate Data in time      
  wt1=(ALMIPII_struc(n)%time2-currTime)/ & 
       (ALMIPII_struc(n)%time2-ALMIPII_struc(n)%time1)
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

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Snowf%varname(1),snowfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Snowf in the forcing variables list')
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     
     if(ALMIPII_struc(n)%metdata1(3,index1).ne.LIS_rc%udef.and.&
          ALMIPII_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) =ALMIPII_struc(n)%metdata1(3,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(3,index1)*wt2
     endif
     
  enddo

  !do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if(ALMIPII_struc(n)%metdata2(8,index1).ne.LIS_rc%udef) then 
        pcp(t)=ALMIPII_struc(n)%metdata2(8,index1)
        pcp(t)  = pcp(t)
     endif
  enddo

  call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if(ALMIPII_struc(n)%metdata2(9,index1).ne.LIS_rc%udef) then 
        snowf(t)=ALMIPII_struc(n)%metdata2(9,index1)	
        snowf(t) = snowf(t)
     endif
  enddo
  !linearly interpolate everything else
  
  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(1,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then 
        tmp(t) =ALMIPII_struc(n)%metdata1(1,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(1,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(2,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then 
        q2(t) =ALMIPII_struc(n)%metdata1(2,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(2,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(4,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(4,index1).ne.LIS_rc%udef)) then 
        lwd(t) =ALMIPII_struc(n)%metdata1(4,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(4,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(5,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(5,index1).ne.LIS_rc%udef)) then 
        uwind(t) =ALMIPII_struc(n)%metdata1(5,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(5,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(6,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(6,index1).ne.LIS_rc%udef)) then 
        vwind(t) =ALMIPII_struc(n)%metdata1(6,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(6,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((ALMIPII_struc(n)%metdata1(7,index1).ne.LIS_rc%udef).and.&
          (ALMIPII_struc(n)%metdata2(7,index1).ne.LIS_rc%udef)) then 
        psurf(t) =ALMIPII_struc(n)%metdata1(7,index1)*wt1+ & 
             ALMIPII_struc(n)%metdata2(7,index1)*wt2
     endif
  enddo
end subroutine timeinterp_ALMIPII
 
