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
! !ROUTINE: timeinterp_gswp1
! \label{timeinterp_gswp1}
!
! !REVISION HISTORY:
! 11 Dec 2003: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine timeinterp_gswp1(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod, only        : LIS_get_nstep, LIS_time2date
  use LIS_constantsMod, only       : LIS_CONST_SOLAR
  use gswp1_forcingMod, only : gswp1_struc
  use LIS_logMod, only         : LIS_logunit, LIS_verify, LIS_endrun
! !ARGUMENTS: 
  implicit none

  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous hourly value
!  is used. All other variables are linearly interpolated between 
!  the hourly blocks. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP

!==== Local Variables=======================

  integer :: c,f,zdoy
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  integer :: nforce, index1
  real*8 :: btime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real :: zw1,zw2
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  nforce = LIS_rc%met_nf(findex)

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

  btime=gswp1_struc(n)%gswp1time1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gswp1_struc(n)%gswp1time2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------
  wt1 = (gswp1_struc(n)%gswp1time2-LIS_rc%time)/ &
       (gswp1_struc(n)%gswp1time2-gswp1_struc(n)%gswp1time1)
  wt2 = 1.0-wt1
  
  do f = 1,nforce
     if (f.eq.3) then
        if (LIS_rc%shortflag.eq.2) then
!-----------------------------------------------------------------------
! Got Time Averaged SW
!-----------------------------------------------------------------------
           do c = 1,LIS_rc%ntiles(n)
              zdoy=LIS_rc%doy
              index1 = LIS_domain(n)%tile(c)%index
              call zterp(1,LIS_domain(n)%grid(index1)%lat,&
                   LIS_domain(n)%grid(index1)%lon, &
                   gmt1,gmt2,LIS_rc%gmt,zdoy, &
                   zw1,zw2,czb,cze,czm,LIS_rc)
              swd(c) = gswp1_struc(n)%metdata1(f,index1)*wt1+&
                   gswp1_struc(n)%metdata2(f,index1)*wt2
              
              if ((swd(c).ne.LIS_rc%udef).and. &
                   (swd(c).lt.0)) then
                 if (swd(c).gt.-0.00001) then
                    swd(c) = 0.0
                 else
                    write(LIS_logunit,*)'ERR: timeinterp_gswp1 -- Stopping ', &
                         'because forcing not udef but lt0,'
                    call LIS_endrun
                 endif
              endif
              if (swd(c).gt.LIS_CONST_SOLAR) then
                 swd(c)=gswp1_struc(n)%metdata2(3,index1)
              endif
              
           enddo
        endif
        
     else
        
        do c = 1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(c)%index
           tmp(c) =gswp1_struc(n)%metdata1(1,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(1,index1)*wt2
           q2(c) =gswp1_struc(n)%metdata1(2,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(2,index1)*wt2
           lwd(c) =gswp1_struc(n)%metdata1(4,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(4,index1)*wt2
           uwind(c) =gswp1_struc(n)%metdata1(5,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(5,index1)*wt2
           vwind(c) =gswp1_struc(n)%metdata1(6,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(6,index1)*wt2
           psurf(c) =gswp1_struc(n)%metdata1(7,index1)*wt1+ & 
                gswp1_struc(n)%metdata2(7,index1)*wt2
           pcp(c) =gswp1_struc(n)%metdata2(8,index1)
           cpcp(c) =gswp1_struc(n)%metdata2(9,index1)
        enddo
     endif
  enddo
84 format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)

  return

end subroutine timeinterp_gswp1

