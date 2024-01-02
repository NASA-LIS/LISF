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
! !ROUTINE: timeinterp_galwemge
! \label{timeinterp_galwemge}
!
! !REVISION HISTORY:
!
! 09 May 2022: Yeosang Yoon, Initial specification
!
! !INTERFACE:
subroutine timeinterp_galwemge(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod
  use LIS_metforcingMod
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod
  use LIS_logMod
  use galwemge_forcingMod
  use LIS_forecastMod

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
  integer :: zdoy
  real    :: zw1, zw2
  real    :: czm, cze, czb
  real    :: wt1, wt2, swt1, swt2
  real    :: gmt1, gmt2
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2,tempbts
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn,tempbss
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)
  integer            :: mfactor, m, k, tid, fcsthr_intv
! ________________________________________

  btime=galwemge_struc(n)%fcsttime1
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

  btime=galwemge_struc(n)%fcsttime2
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

  !Interpolate Data in Time
  wt1=(galwemge_struc(n)%fcsttime2-LIS_rc%time)/ &
       (galwemge_struc(n)%fcsttime2-galwemge_struc(n)%fcsttime1)
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

  ! Metforcing ensemble member count factor
  mfactor = LIS_rc%nensem(n)/galwemge_struc(n)%max_ens_members

  ! Downward shortwave radiation (average):
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,galwemge_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index
           zdoy=LIS_rc%doy

           ! Compute and apply zenith angle weights
           call zterp(0,LIS_domain(n)%grid(index1)%lat,&
                LIS_domain(n)%grid(index1)%lon,&
                gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

           if(galwemge_struc(n)%metdata1(3,m,index1).ne.LIS_rc%udef.and.&
              galwemge_struc(n)%metdata2(3,m,index1).ne.LIS_rc%udef) then
              !swd(tid) = galwemge_struc(n)%metdata1(3,m,index1)*zw1+&
              !         galwemge_struc(n)%metdata2(3,m,index1)*zw2

              swd(tid) = zw1 * galwemge_struc(n)%metdata2(3,m,index1)
              ! In cases of small cos(zenith) angles, use linear weighting
              !  to avoid overly large weights
              !if((swd(tid).gt.galwemge_struc(n)%metdata1(3,m,index1).and. &
              !   swd(tid).gt.galwemge_struc(n)%metdata2(3,m,index1)).and. &
              !   (czb.lt.0.1.or.cze.lt.0.1))then
              !   swd(tid) = galwemge_struc(n)%metdata1(3,m,index1)*swt1+ &
              !            galwemge_struc(n)%metdata2(3,m,index1)*swt2
              !endif

              if(swd(tid).lt.0.0) then
                 write(unit=LIS_logunit,fmt=*)'[ERR] SW radiation is unphysical'
                 write(unit=LIS_logunit,fmt=*)'[ERR] it is', LIS_localPet, tid, swd(tid)
                 write(unit=LIS_logunit,fmt=*)'[ERR] data1=',galwemge_struc(n)%metdata1(3,m,index1)
                 write(unit=LIS_logunit,fmt=*)'[ERR] data2=',galwemge_struc(n)%metdata2(3,m,index1)
                 write(unit=LIS_logunit,fmt=*)'[ERR] wts=',wt1,wt2
                 call LIS_endrun
              endif

              if(swd(tid).gt.LIS_CONST_SOLAR) then
                !swd(tid)=galwemge_struc(n)%metdata2(3,m,index1)
                !In cases of small cos(zenith) angles, use linear weighting
                ! to avoid overly large weights
                swd(tid) = galwemge_struc(n)%metdata1(3,m,index1)*wt1+ &
                          galwemge_struc(n)%metdata2(3,m,index1)*wt2
              endif

           endif
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! precip variable Block Interpolation
!-----------------------------------------------------------------------
  ! Total precipitation field (accumulated):
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,galwemge_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index

           if(galwemge_struc(n)%metdata2(8,m,index1).ne.LIS_rc%udef) then
              ! account for the accum fields
              pcp(tid)=galwemge_struc(n)%metdata2(8,m,index1)/(3600*galwemge_struc(n)%fcst_hour)
              if(pcp(tid).lt.0) then
                 pcp(tid) = 0.0
              endif
           endif
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! Linearly interpolate everything else
!-----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,galwemge_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index
           
           ! 2-meter air temp
           if((galwemge_struc(n)%metdata1(1,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(1,m,index1).ne.LIS_rc%udef)) then
              tmp(tid) = galwemge_struc(n)%metdata1(1,m,index1)*wt1 + &
                         galwemge_struc(n)%metdata2(1,m,index1)*wt2
           endif
           ! Specific humidity
           if((galwemge_struc(n)%metdata1(2,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(2,m,index1).ne.LIS_rc%udef)) then
              q2(tid)  = galwemge_struc(n)%metdata1(2,m,index1)*wt1 + &
                         galwemge_struc(n)%metdata2(2,m,index1)*wt2
           endif
           ! Downward longwave field
           if((galwemge_struc(n)%metdata1(4,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(4,m,index1).ne.LIS_rc%udef)) then
              lwd(tid)  = galwemge_struc(n)%metdata1(4,m,index1)*wt1 + &
                          galwemge_struc(n)%metdata2(4,m,index1)*wt2
           endif
           ! U-wind component
           if((galwemge_struc(n)%metdata1(5,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(5,m,index1).ne.LIS_rc%udef)) then
              uwind(tid) = galwemge_struc(n)%metdata1(5,m,index1)*wt1+&
                           galwemge_struc(n)%metdata2(5,m,index1)*wt2
           endif
           ! V-wind component
           if((galwemge_struc(n)%metdata1(6,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(6,m,index1).ne.LIS_rc%udef)) then
              vwind(tid) = galwemge_struc(n)%metdata1(6,m,index1)*wt1 + &
                           galwemge_struc(n)%metdata2(6,m,index1)*wt2
           endif
           ! Surface pressure field
           if((galwemge_struc(n)%metdata1(7,m,index1).ne.LIS_rc%udef).and.&
                (galwemge_struc(n)%metdata2(7,m,index1).ne.LIS_rc%udef)) then
              psurf(tid) = galwemge_struc(n)%metdata1(7,m,index1)*wt1 + &
                           galwemge_struc(n)%metdata2(7,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

end subroutine timeinterp_galwemge
