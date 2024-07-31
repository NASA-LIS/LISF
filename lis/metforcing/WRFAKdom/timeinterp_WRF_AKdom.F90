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
!
! !ROUTINE: timeinterp_WRF_AKdom
!  \label{timeinterp_WRF_AKdom}
! 
! !REVISION HISTORY: 
! 21 Jun 2021; K.R. Arsenault, Updated for different WRF AK files
! 
! !INTERFACE:
subroutine timeinterp_WRF_AKdom(n, findex)

! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use LIS_constantsMod,    only : LIS_CONST_SOLAR
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,   only : LIS_forc, LIS_FORC_Base_State
  use LIS_timeMgrMod,      only : LIS_time2date, LIS_tick
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
  use WRF_AKdom_forcingMod, only : WRFAK_struc
  use LIS_forecastMod,     only : LIS_get_iteration_index

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. All native input fields are 
!  instantaneous, so for the "zterp" option of 1 is selected.
!
!  Precipitation represents accumulation over the previous hour,
!   and converted here to rate.
!  All other variables are linearly interpolated between 
!  the hourly blocks. 
!  Only total precipitation is available (no convective).
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  real    :: wt1,wt2,swt1,swt2
  real    :: zw1,zw2,czb,cze,czm,gmt1,gmt2
  integer :: zdoy, t
  integer :: index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime, newtime1, newtime2
  real    :: tempgmt1,tempgmt2, tempbts
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss

  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  integer            :: mfactor, m, k, kk
! ____________________________________________________________
  
! Time-1
  btime=WRFAK_struc(n)%WRFouttime1
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

! Time-2
  btime=WRFAK_struc(n)%WRFouttime2
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

  !===  Interpolate Data in Time
  wt1=(WRFAK_struc(n)%WRFouttime2-LIS_rc%time)/&
       (WRFAK_struc(n)%WRFouttime2-WRFAK_struc(n)%WRFouttime1)
  wt2=1.0-wt1

  ! Solar zenith angle weights for SWdown:
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1

  !===
 
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

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
          rc=status)
     call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')
  endif

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)
        
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  endif

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

   pcp  = 0.0

   if (LIS_FORC_CRainf%selectOpt.eq.1) then
     cpcp = 0.0
   endif

   ! Multi-factor for metforcing ensemble number:
   mfactor = LIS_rc%nensem(n)/WRFAK_struc(n)%nIter

   ! SWdown radiation block:
   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)

       zdoy = LIS_rc%doy
       ! Compute and apply zenith angle weights
       call zterp( 1, LIS_domain(n)%grid(index1)%lat,&
            LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
            LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

       if(WRFAK_struc(n)%metdata1(kk,3,index1).ne.LIS_rc%udef.and.&
         WRFAK_struc(n)%metdata2(kk,3,index1).ne.LIS_rc%udef) then

         swd(t) = WRFAK_struc(n)%metdata1(kk,3,index1)*zw1+&
                  WRFAK_struc(n)%metdata2(kk,3,index1)*zw2

         ! In cases of small cos(zenith) angles, use linear weighting
         !  to avoid overly large weights
         if( swd(t) > WRFAK_struc(n)%metdata1(kk,3,index1).and. &
             swd(t) > WRFAK_struc(n)%metdata2(kk,3,index1).and. &
            (czb.lt.0.1.or.cze.lt.0.1)) then

           swd(t) = WRFAK_struc(n)%metdata1(kk,3,index1)*swt1+&
                    WRFAK_struc(n)%metdata2(kk,3,index1)*swt2
         endif
       end if

       ! Check range of SWDown values:
       if (swd(t) < 0) then
          write(LIS_logunit,*) '[ERR] WRFout V2 warning! SW radiation is negative!!'
          write(LIS_logunit,*) 'sw=', swd(t), '... negative'
          write(LIS_logunit,*) 'WRFout1 =', WRFAK_struc(n)%metdata1(kk,3,index1)
          write(LIS_logunit,*) 'WRFout2 =', WRFAK_struc(n)%metdata2(kk,3,index1)
          call LIS_endrun
       end if
       if (swd(t).gt.LIS_CONST_SOLAR) then
          write(LIS_logunit,fmt=*)'[WARN] SW radiation too high!!'
          write(LIS_logunit,fmt=*)'[WARN] it is', swd(t)
          write(LIS_logunit,fmt=*)'[WARN] data1=',WRFAK_struc(n)%metdata1(kk,3,index1)
          write(LIS_logunit,fmt=*)'[WARN] data2=',WRFAK_struc(n)%metdata2(kk,3,index1)
          write(LIS_logunit,fmt=*)'[WARN] zw1=',zw1,'zw2=',zw2
          write(LIS_logunit,fmt=*)'[WARN] swt1=',swt1,'swt2=',swt2
          write(LIS_logunit,fmt=*)'[WARN] Setting SWD to time 1 value'
          swd(t)=WRFAK_struc(n)%metdata1(kk,3,index1)
       endif

     end do
  enddo   ! End for SWdown block

  ! Time Averaged Longwave, Block Interpolation
   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)
 
       ! Linearly interpolate instantaneous fields:
       tmp(t) = wt1 * WRFAK_struc(n)%metdata1(kk,1,index1) & 
              + wt2 * WRFAK_struc(n)%metdata2(kk,1,index1)
       q2(t) = wt1 * WRFAK_struc(n)%metdata1(kk,2,index1) & 
             + wt2 * WRFAK_struc(n)%metdata2(kk,2,index1)
       lwd(t) = wt1 * WRFAK_struc(n)%metdata1(kk,4,index1) &
              + wt2 * WRFAK_struc(n)%metdata2(kk,4,index1)
       uwind(t) = wt1 * WRFAK_struc(n)%metdata1(kk,5,index1) & 
                + wt2 * WRFAK_struc(n)%metdata2(kk,5,index1)
       vwind(t) = wt1 * WRFAK_struc(n)%metdata1(kk,6,index1) & 
                + wt2 * WRFAK_struc(n)%metdata2(kk,6,index1)
       psurf(t) = wt1 * WRFAK_struc(n)%metdata1(kk,7,index1) & 
                + wt2 * WRFAK_struc(n)%metdata2(kk,7,index1)
       ! Handle accum fields:
       pcp(t) = WRFAK_struc(n)%metdata1(kk,8,index1)/3600.  ! Convert to mm/sec

     enddo
  enddo

end subroutine timeinterp_WRF_AKdom

