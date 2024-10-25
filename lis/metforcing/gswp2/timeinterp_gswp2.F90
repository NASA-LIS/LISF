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
! !ROUTINE: timeinterp_gswp2
! \label{timeinterp_gswp2}
!
! !REVISION HISTORY:
!  20 Feb 2004; Sujay Kumar : Initial Specification
!  10 Oct 2006: Sujay Kumar: Switched to using ESMF_State for storing
!               forcing data. 
! !INTERFACE:
subroutine timeinterp_gswp2(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_constantsMod,   only : LIS_CONST_SOLAR
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod,     only : LIS_time2date
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use gswp2_forcingMod,    only : gswp2_struc

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
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!   \item[finterp\_gswp2](\ref{finterp_gswp2}) \newline
!    routine to compute temporal disaggregation of forcing
!  \end{description}
!EOP
  integer :: t,zdoy
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  real*8 :: btime
  real :: wt1,wt2,gmt1,gmt2
  integer :: madtt
  integer :: sub_dt = 1
  character(len=1), dimension(10) :: trp_flag = &
       (/'I','I','L','L','I','X','I','L','L','X'/)
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  zdoy=LIS_rc%doy
  btime=gswp2_struc(n)%gswp2time1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  btime=gswp2_struc(n)%gswp2time2
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

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Snowf%varname(1),snowfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Snowf in the forcing variables list')

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

  call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

  if ( trim(LIS_rc%met_tinterp(findex)) == "linear") then

   !-----------------------------------------------------------------------
   !  Interpolate Data in Time
   !-----------------------------------------------------------------------
     wt1=(gswp2_struc(n)%gswp2time2-LIS_rc%time)/ & 
          (gswp2_struc(n)%gswp2time2-gswp2_struc(n)%gswp2time1)
     wt2=1.0-wt1
     
     if (LIS_rc%shortflag.eq.2) then
   !-----------------------------------------------------------------------
   ! Got Time Averaged SW
   !-----------------------------------------------------------------------
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           swd(t)=gswp2_struc(n)%metdata1(3,index1)*wt1+ & 
                gswp2_struc(n)%metdata2(3,index1)*wt2 
           if ((swd(t).ne.LIS_rc%udef).and. & 
                (swd(t).lt.0) ) then
              if (swd(t) > -0.00001) then 
                 swd(t) = 0.0 
              else
                 write(LIS_logunit,*)'ERR: timeinterp_gswp2 -- Stopping because ', & 
                      'forcing not udef but lt0,'
                 write(LIS_logunit,*)'ERR: timeinterp_gswp2 -- ', & 
                      t,swd(t),gswp2_struc(n)%metdata2(3,index1)
                 call LIS_endrun
              end if
           endif
           
           if (swd(t).gt.LIS_CONST_SOLAR) then
              swd(t)=gswp2_struc(n)%metdata2(3,index1)
           endif
        enddo
     endif
   !-----------------------------------------------------------------------
   ! precip variable Block Interpolation
   !-----------------------------------------------------------------------
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        pcp(t) = gswp2_struc(n)%metdata2(8,index1)
        cpcp(t) = gswp2_struc(n)%metdata2(9,index1)
     enddo
     if (LIS_rc%longflag.eq.1) then 
   !-----------------------------------------------------------------------
   !    Got Instantaneous LW
   !-----------------------------------------------------------------------
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           lwd(t)=gswp2_struc(n)%metdata1(4,index1)*wt1+ & 
                gswp2_struc(n)%metdata2(4,index1)*wt2 
        enddo
     endif
     if (LIS_rc%longflag.eq.2) then 
   !-----------------------------------------------------------------------
   !    Got Time Averaged LW
   !-----------------------------------------------------------------------
        do t=1,LIS_rc%ntiles(n)
           index1 = LIS_domain(n)%tile(t)%index
           lwd(t) =gswp2_struc(n)%metdata2(4,index1)
        enddo
     endif

   !-----------------------------------------------------------------------
   !     Linearly interpolate everything else
   !-----------------------------------------------------------------------
     do t=1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(t)%index
        tmp(t) =gswp2_struc(n)%metdata1(1,index1)*wt1+ & 
             gswp2_struc(n)%metdata2(1,index1)*wt2
        q2(t) =gswp2_struc(n)%metdata1(2,index1)*wt1+ & 
             gswp2_struc(n)%metdata2(2,index1)*wt2
        uwind(t) =gswp2_struc(n)%metdata1(5,index1)*wt1+ & 
             gswp2_struc(n)%metdata2(5,index1)*wt2
        vwind(t) =gswp2_struc(n)%metdata1(6,index1)*wt1+ & 
             gswp2_struc(n)%metdata2(6,index1)*wt2
        psurf(t) =gswp2_struc(n)%metdata1(7,index1)*wt1+ & 
             gswp2_struc(n)%metdata2(7,index1)*wt2
     enddo
  else
     ! Note:  the GSWP2 temporal interpolation routine, finterp_gswp2, returns
     ! an array of length madtt of interpolated values.
     ! Helin Wei's version only returns the value corresponding to the
     ! current time-step, sub_dt, w.r.t. the given 3-hour forcing interval.

     madtt = 3600 / LIS_rc%ts * 3 ! number of time-steps in a 
                                 ! 3-hourly forcing interval
     do t = 1, LIS_rc%ntiles(n)
        ! Note that some grid-cells may have undefined forcing.
        ! When these are temporally interpolated, they get a 
        ! value of 0.0.  Trap these points and set the temporal
        ! interpolation result to undefined.
        ! This check belongs higher up in the code.  This case could
        ! happen for other domains / base forcing sets.
        index1 = LIS_domain(n)%tile(t)%index
        if ( gswp2_struc(n)%metdata1(1,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(1,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(1,index1) == LIS_rc%udef ) then
           tmp(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(1), &
                0.0, gswp2_struc(n)%metdata1(1,index1), &
                gswp2_struc(n)%metdata2(1,index1), &
                gswp2_struc(n)%metdata3(1,index1), &
                madtt, sub_dt,         &
                tmp(t))
        endif
        
        if ( gswp2_struc(n)%metdata1(2,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(2,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(2,index1) == LIS_rc%udef ) then
           q2(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(2), &
                0.0, gswp2_struc(n)%metdata1(2,index1), &
                gswp2_struc(n)%metdata2(2,index1), &
                gswp2_struc(n)%metdata3(2,index1), &
                madtt, sub_dt,         &
                q2(t))
        endif
        if ( gswp2_struc(n)%metdata1(3,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(3,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(3,index1) == LIS_rc%udef ) then
           swd(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(3), &
                0.0, gswp2_struc(n)%metdata1(3,index1), &
                gswp2_struc(n)%metdata2(3,index1), &
                gswp2_struc(n)%metdata3(3,index1), &
                madtt, sub_dt,         &
                swd(t))
        endif
        
        if ( gswp2_struc(n)%metdata1(4,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(4,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(4,index1) == LIS_rc%udef ) then
           lwd(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(4), &
                0.0, gswp2_struc(n)%metdata1(4,index1), &
                gswp2_struc(n)%metdata2(4,index1), &
                gswp2_struc(n)%metdata3(4,index1), &
                madtt, sub_dt,         &
                lwd(t))
        endif
        if ( gswp2_struc(n)%metdata1(5,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(5,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(5,index1) == LIS_rc%udef ) then
           uwind(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(5), &
                0.0, gswp2_struc(n)%metdata1(5,index1), &
                gswp2_struc(n)%metdata2(5,index1), &
                gswp2_struc(n)%metdata3(5,index1), &
                madtt, sub_dt,         &
                uwind(t))
        endif
           
        if ( gswp2_struc(n)%metdata1(6,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(6,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(6,index1) == LIS_rc%udef ) then
           vwind(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(6), &
                0.0, gswp2_struc(n)%metdata1(6,index1), &
                gswp2_struc(n)%metdata2(6,index1), &
                gswp2_struc(n)%metdata3(6,index1), &
                madtt, sub_dt,         &
                vwind(t))
        endif
        if ( gswp2_struc(n)%metdata1(7,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(7,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(7,index1) == LIS_rc%udef ) then
           psurf(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(7), &
                0.0, gswp2_struc(n)%metdata1(7,index1), &
                gswp2_struc(n)%metdata2(7,index1), &
                gswp2_struc(n)%metdata3(7,index1), &
                madtt, sub_dt,         &
                psurf(t))
        endif
        if ( gswp2_struc(n)%metdata1(8,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(8,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(8,index1) == LIS_rc%udef ) then
           pcp(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(8), &
                0.0, gswp2_struc(n)%metdata1(8,index1), &
                gswp2_struc(n)%metdata2(8,index1), &
                gswp2_struc(n)%metdata3(8,index1), &
                madtt, sub_dt,         &
                pcp(t))
        endif
        if ( gswp2_struc(n)%metdata1(9,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(9,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(9,index1) == LIS_rc%udef ) then
           snowf(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(9), &
                0.0, gswp2_struc(n)%metdata1(9,index1), &
                gswp2_struc(n)%metdata2(9,index1), &
                gswp2_struc(n)%metdata3(9,index1), &
                madtt, sub_dt,         &
                snowf(t))
        endif
        if ( gswp2_struc(n)%metdata1(10,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata2(10,index1) == LIS_rc%udef .or. &
             gswp2_struc(n)%metdata3(10,index1) == LIS_rc%udef ) then
           cpcp(t) = LIS_rc%udef
        else
           call finterp_gswp2(0, trp_flag(10), &
                0.0, gswp2_struc(n)%metdata1(10,index1), &
                gswp2_struc(n)%metdata2(10,index1), &
                gswp2_struc(n)%metdata3(10,index1), &
                madtt, sub_dt,         &
                cpcp(t))
        endif
     enddo
     
     sub_dt = sub_dt + 1
     if ( sub_dt > madtt ) then
        sub_dt = 1
     endif
84   format('now',i4,4i3,2x,'pvt ',a22,' nxt ',a22)
  end if
end subroutine timeinterp_gswp2

