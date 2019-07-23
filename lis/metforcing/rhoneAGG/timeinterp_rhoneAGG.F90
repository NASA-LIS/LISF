!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_rhoneAGG
! \label{timeinterp_rhoneAGG}
!
! !REVISION HISTORY:
!  5 Nov 2003: Dave Mocko, Initial Specification
! !INTERFACE:
subroutine timeinterp_rhoneAGG(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only :LIS_rc, LIS_domain, LIS_localPet
  use LIS_constantsMod, only  : LIS_CONST_SOLAR
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only: LIS_forc, LIS_FORC_Base_State
  use LIS_timeMgrMod,     only: LIS_time2date
  use LIS_logMod, only : LIS_logunit, LIS_verify, LIS_endrun
  use rhoneAGG_forcingMod, only : rhoneAGG_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)     :: n 
  integer, intent(in)     :: findex
! !DESCRIPTION:
!  Opens, reads, and interpolates RHONEAGG forcing.
!
!    TIME1 = most recent past data
!    TIME2 = nearest future data
!
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day.
!
! \subsection{Core Functions of time$_-$interp$_-$rhoneAGG}
!  \begin{description}
!  \item[zterp]
!   Performs zenith angle-based temporal interpolation
!  \end{description}
!
!EOP

  integer :: t,zdoy
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime
  real :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real :: zw1,zw2
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  
  btime=rhoneAGG_struc(n)%rhoneAGGtime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=rhoneAGG_struc(n)%rhoneAGGtime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------
  wt1=(rhoneAGG_struc(n)%rhoneAGGtime2-LIS_rc%time)/ & 
       (rhoneAGG_struc(n)%rhoneAGGtime2-rhoneAGG_struc(n)%rhoneAGGtime1)
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
  
  if (LIS_rc%shortflag.eq.2) then
!-----------------------------------------------------------------------
! Got Time Averaged SW
!-----------------------------------------------------------------------
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        zdoy=LIS_rc%doy
        
        call zterp(0,LIS_domain(n)%grid(index1)%lat,&
             LIS_domain(n)%grid(index1)%lon, & 
             gmt1,gmt2,LIS_rc%gmt,zdoy, & 
             zw1,zw2,czb,cze,czm,LIS_rc)
        swd(t) =rhoneAGG_struc(n)%metdata2(3,index1)*zw1
        if ((swd(t).ne.LIS_rc%udef).and. & 
             (swd(t).lt.0) ) then
           if (swd(t) > -0.00001) then 
              swd(t) = 0.0 
           else
              write(LIS_logunit,*)'ERR: timeinterp_rhoneAGG -- Stopping because ', & 
                      'forcing not udef but lt0,'
              write(LIS_logunit,*)'ERR: timeinterp_rhoneAGG -- ', & 
                   t,swd(t),rhoneAGG_struc(n)%metdata2(3,index1), & 
                   ' (',LIS_localPet,')'
              call LIS_endrun
           end if
        endif
        
        if (swd(t).gt.LIS_CONST_SOLAR) then
           swd(t)=rhoneAGG_struc(n)%metdata2(3,index1)
        endif
     enddo
  endif
!-----------------------------------------------------------------------
! precip variable Block Interpolation
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     pcp(t) = rhoneAGG_struc(n)%metdata2(8,index1)
     cpcp(t) = rhoneAGG_struc(n)%metdata2(9,index1)
  enddo

  if (LIS_rc%longflag.eq.1) then 
!-----------------------------------------------------------------------
!    Got Instantaneous LW
!-----------------------------------------------------------------------
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        lwd(t) = rhoneAGG_struc(n)%metdata1(4,index1)*wt1+ & 
             rhoneAGG_struc(n)%metdata2(4,index1)*wt2 
     enddo
  endif

  if (LIS_rc%longflag.eq.2) then 
!-----------------------------------------------------------------------
!    Got Time Averaged LW
!-----------------------------------------------------------------------
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        lwd(t) = rhoneAGG_struc(n)%metdata2(4,index1)
     enddo
  endif

!-----------------------------------------------------------------------
!     Linearly interpolate everything else
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     tmp(t) =rhoneAGG_struc(n)%metdata1(1,index1)*wt1+ & 
          rhoneAGG_struc(n)%metdata2(1,index1)*wt2
     
     q2(t) =rhoneAGG_struc(n)%metdata1(2,index1)*wt1+ & 
          rhoneAGG_struc(n)%metdata2(2,index1)*wt2
     
     uwind(t) =rhoneAGG_struc(n)%metdata1(5,index1)*wt1+ & 
          rhoneAGG_struc(n)%metdata2(5,index1)*wt2
     
     vwind(t) =rhoneAGG_struc(n)%metdata1(6,index1)*wt1+ & 
          rhoneAGG_struc(n)%metdata2(6,index1)*wt2
     
     psurf(t) =rhoneAGG_struc(n)%metdata1(7,index1)*wt1+ & 
          rhoneAGG_struc(n)%metdata2(7,index1)*wt2
  enddo
end subroutine timeinterp_rhoneAGG

