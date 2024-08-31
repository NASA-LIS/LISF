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
! !ROUTINE: timeinterp_gefs
! \label{timeinterp_gefs}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
!
! !INTERFACE:
subroutine timeinterp_gefs(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc,LIS_domain
  use LIS_constantsMod,     only : LIS_CONST_SOLAR
  use LIS_metforcingMod,    only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod,       only : LIS_tick, LIS_time2date
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_endrun
  use gefs_forcingMod,      only : gefs_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the GEFS forecast forcing data to the current model 
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
  integer          :: zdoy,k,m,tid,mfactor
  real             :: zw1,zw2
  real             :: czm,cze,czb
  real             :: wt1,wt2
  real             :: gmt1,gmt2
  integer          :: t,index1
  integer          :: bdoy,byr,bmo,bda,bhr,bmn
  real*8           :: btime
  integer          :: status

  type(ESMF_Field) :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field) :: psurfField,pcpField
  real,pointer     :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer     :: swd(:),lwd(:),psurf(:),pcp(:)

! _____________________________________________________________

  ! Times for downward SW radiation time interp
  btime=gefs_struc(n)%fcsttime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)

  btime=gefs_struc(n)%fcsttime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)

  
 !== Interpolate data in time      

  ! Check if bookend times differ, else stop ...
  if( (gefs_struc(n)%fcsttime2-gefs_struc(n)%fcsttime1)==0 ) then
     write(LIS_logunit,*) " GEFS forecast times are same (in time interp) ... "
     call LIS_endrun
  endif

  wt1=(gefs_struc(n)%fcsttime2-LIS_rc%time)/ & 
       (gefs_struc(n)%fcsttime2-&
        gefs_struc(n)%fcsttime1)
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


  ! Metforcing ensemble member count factor
  mfactor = LIS_rc%nensem(n)/gefs_struc(n)%max_ens_members

  ! Downward shortwave radiation (average):
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  zdoy=LIS_rc%doy
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 

           ! Compute and apply zenith angle weights
           call zterp( 0, LIS_domain(n)%grid(index1)%lat,   &
                LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
                LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc )

           if( gefs_struc(n)%metdata2(3,m,index1).ne.LIS_rc%udef .and. &
               gefs_struc(n)%metdata1(3,m,index1).ne.LIS_rc%udef ) then

              swd(tid) = zw1 * gefs_struc(n)%metdata2(3,m,index1)

             if (swd(tid) < 0) then
                write(LIS_logunit,*) '[ERR] SW radiation is negative!'
                write(LIS_logunit,*) '[ERR] sw =', swd(tid)
                write(LIS_logunit,*) '[ERR] zw1=', zw1
                write(LIS_logunit,*) '[ERR] GEFS data =', gefs_struc(n)%metdata2(3,m,index1)
                call LIS_endrun
             end if
             if( swd(tid).gt.LIS_CONST_SOLAR ) then 
!                swd(tid)=gefs_struc(n)%metdata2(3,m,index1)
               ! In cases of small cos(zenith) angles, use linear weighting
               !  to avoid overly large weights
               swd(tid) = gefs_struc(n)%metdata1(3,m,index1)*wt1+ &
                          gefs_struc(n)%metdata2(3,m,index1)*wt2
             endif

           endif

        enddo
     enddo
  enddo
  
  ! 2-meter air temp (instantaneous):
  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(1,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(1,m,index1).ne.LIS_rc%udef)) then 
              tmp(tid) = gefs_struc(n)%metdata1(1,m,index1)*wt1+ & 
                         gefs_struc(n)%metdata2(1,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  ! Specific humidity (instantaneous):
  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(2,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(2,m,index1).ne.LIS_rc%udef)) then 
              q2(tid) = gefs_struc(n)%metdata1(2,m,index1)*wt1+ & 
                        gefs_struc(n)%metdata2(2,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  ! Downward longwave field (average):
  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(4,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(4,m,index1).ne.LIS_rc%udef)) then 
              lwd(tid) = gefs_struc(n)%metdata1(4,m,index1)*wt1+ & 
                         gefs_struc(n)%metdata2(4,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  ! U-wind component (instantaneous):
  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(5,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(5,m,index1).ne.LIS_rc%udef)) then 
              uwind(tid) = gefs_struc(n)%metdata1(5,m,index1)*wt1+ & 
                           gefs_struc(n)%metdata2(5,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  ! V-wind component (instantaneous):
  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(6,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(6,m,index1).ne.LIS_rc%udef)) then 
              vwind(tid) = gefs_struc(n)%metdata1(6,m,index1)*wt1+ & 
                           gefs_struc(n)%metdata2(6,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  ! Surface pressure field (instantaneous):
  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((gefs_struc(n)%metdata1(7,m,index1).ne.LIS_rc%udef).and.&
                (gefs_struc(n)%metdata2(7,m,index1).ne.LIS_rc%udef)) then 
              psurf(tid) = gefs_struc(n)%metdata1(7,m,index1)*wt1+ & 
                           gefs_struc(n)%metdata2(7,m,index1)*wt2
           endif
        enddo
     enddo
  enddo
  
  ! Total precipitation field (accumulated):
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,gefs_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if(gefs_struc(n)%metdata2(8,m,index1).ne.LIS_rc%udef) then
              if( gefs_struc(n)%gefs_fcsttype .eq. "Reforecast2" ) then
                 pcp(tid)=gefs_struc(n)%metdata2(8,m,index1)/(3600*6)
              else
                 pcp(tid)=gefs_struc(n)%metdata2(8,m,index1)/(3600*3)
              endif
              ! Eventually account for the 3-hour accum fields
              !  which alternate with the 6-hour up to hour 69 
              !  (need to subtract off the alternating 3-hour from
              !   6-hour accumulated precip to get at having 3-hourly)
           endif
        enddo
     enddo
  enddo

end subroutine timeinterp_gefs
 
