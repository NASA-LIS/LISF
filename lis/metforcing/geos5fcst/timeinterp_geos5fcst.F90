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
! !ROUTINE: timeinterp_geos5fcst
! \label{timeinterp_geos5fcst}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
!
! !INTERFACE:
subroutine timeinterp_geos5fcst(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc,LIS_domain
  use LIS_constantsMod,     only : LIS_CONST_SOLAR
  use LIS_metforcingMod,    only : LIS_forc, LIS_FORC_Base_State
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod,       only : LIS_tick, LIS_time2date
  use LIS_logMod,           only : LIS_logunit, LIS_verify, LIS_endrun
  use geos5fcst_forcingMod, only : geos5fcst_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the GEOS5 forecast forcing data to the current model 
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
  real             :: zw1, zw2
  real             :: czm, cze, czb
  real             :: wt1, wt2,swt1,swt2
  real             :: gmt1, gmt2, tempbts
  integer          :: t,index1
  integer          :: bdoy,byr,bmo,bda,bhr,bmn
  real*8           :: btime,newtime1,newtime2
  real             :: tempgmt1,tempgmt2
  integer          :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer          :: tempbss
  integer          :: status
  type(ESMF_Field) :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field) :: psurfField,pcpField,cpcpField,fhgtField
  type(ESMF_Field) :: pardrField,pardfField,swlandField,snowfField
  real,pointer     :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer     :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer     :: fheight(:),pardr(:),pardf(:)
  real,pointer     :: swland(:),snowf(:)

  logical          :: forcing_z, forcing_ch
  integer :: kk 

  btime=geos5fcst_struc(n)%fcsttime1
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
  
  btime=geos5fcst_struc(n)%fcsttime2
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
  wt1=(geos5fcst_struc(n)%fcsttime2-LIS_rc%time)/ & 
       (geos5fcst_struc(n)%fcsttime2-&
       geos5fcst_struc(n)%fcsttime1)
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

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')

  if (LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Snowf%varname(1),snowfField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable Snowf in the forcing variables list')
  endif


  if (LIS_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Pardr%varname(1),pardrField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable PARDR in the forcing variables list')
  endif

  if (LIS_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Pardf%varname(1),pardfField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable PARDF in the forcing variables list')
  endif

  if (LIS_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
          LIS_FORC_Swnet%varname(1),swlandField,&
          rc=status)
     call LIS_verify(status, &
          'Error: Enable SWnet in the forcing variables list')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Forc_Hgt%varname(1),&
          fhgtField, rc=status)
     forcing_z = .true.
  else
     forcing_z = .false.
  endif

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  mfactor = LIS_rc%nensem(n)/geos5fcst_struc(n)%max_ens_members

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           zdoy=LIS_rc%doy
     !  compute and apply zenith angle weights
           call zterp( 0, LIS_domain(n)%grid(index1)%lat,&
                LIS_domain(n)%grid(index1)%lon, gmt1, gmt2, & 
                LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

           if(geos5fcst_struc(n)%metdata2(3,m,index1).ne.LIS_rc%udef) then 
              swd(tid) = zw1 * geos5fcst_struc(n)%metdata2(3,m,index1)

              if (swd(tid) < 0) then
                 write(LIS_logunit,*) '[ERR] 2 warning!!!  SW radiation is negative!!'
                 write(LIS_logunit,*) '[ERR] sw=', swd(tid), '... negative'
                 write(LIS_logunit,*) '[ERR] GEOS5=', geos5fcst_struc(n)%metdata2(3,m,index1)
                 call LIS_endrun
              end if
              
              if (swd(tid).gt.LIS_CONST_SOLAR) then
                 swd(tid)=geos5fcst_struc(n)%metdata2(3,m,index1)
              endif
           endif
        enddo
     enddo
  enddo
  
  !linearly interpolate everything else

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(1,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(1,m,index1).ne.LIS_rc%udef)) then 
              tmp(tid) =geos5fcst_struc(n)%metdata1(1,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(1,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(2,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(2,m,index1).ne.LIS_rc%udef)) then 
              q2(tid) =geos5fcst_struc(n)%metdata1(2,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(2,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(4,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(4,m,index1).ne.LIS_rc%udef)) then 
              lwd(tid) =geos5fcst_struc(n)%metdata1(4,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(4,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(5,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(5,m,index1).ne.LIS_rc%udef)) then 
              uwind(tid) =geos5fcst_struc(n)%metdata1(5,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(5,m,index1)*wt2
           endif
        enddo
     enddo
  enddo
  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(6,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(6,m,index1).ne.LIS_rc%udef)) then 
              vwind(tid) =geos5fcst_struc(n)%metdata1(6,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(6,m,index1)*wt2
           endif
        enddo
     enddo
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if((geos5fcst_struc(n)%metdata1(7,m,index1).ne.LIS_rc%udef).and.&
                (geos5fcst_struc(n)%metdata2(7,m,index1).ne.LIS_rc%udef)) then 
              psurf(tid) =geos5fcst_struc(n)%metdata1(7,m,index1)*wt1+ & 
                   geos5fcst_struc(n)%metdata2(7,m,index1)*wt2
           endif
        enddo
     enddo
  enddo
  
  !do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if(geos5fcst_struc(n)%metdata2(8,m,index1).ne.LIS_rc%udef) then 
              pcp(tid)=geos5fcst_struc(n)%metdata2(8,m,index1)	
           endif
        enddo
     enddo
  enddo

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status)
     
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,geos5fcst_struc(n)%max_ens_members
           do k=1,mfactor
              tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
              index1 = LIS_domain(n)%tile(tid)%index 
              if(geos5fcst_struc(n)%metdata2(9,m,index1).ne.LIS_rc%udef) then 
                 snowf(tid)=geos5fcst_struc(n)%metdata2(9,m,index1)	
              endif
           enddo
        enddo
     enddo
  endif

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,geos5fcst_struc(n)%max_ens_members
        do k=1,mfactor
           tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
           index1 = LIS_domain(n)%tile(tid)%index 
           if(geos5fcst_struc(n)%metdata2(10,m,index1).ne.LIS_rc%udef) then 
              cpcp(tid)=geos5fcst_struc(n)%metdata2(10,m,index1)	
              ! INPUT is actually convective precip fraction
              ! Calc actual CPC below  
           endif
        enddo
     end do
  enddo

  if ( forcing_z ) then
     call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,geos5fcst_struc(n)%max_ens_members
           do k=1,mfactor
              tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
              index1 = LIS_domain(n)%tile(tid)%index 
              if((geos5fcst_struc(n)%metdata1(11,m,index1).ne.LIS_rc%udef).and.&
                   (geos5fcst_struc(n)%metdata2(11,m,index1).ne.LIS_rc%udef)) then 
                 fheight(tid) =geos5fcst_struc(n)%metdata1(11,m,index1)*wt1+ & 
                      geos5fcst_struc(n)%metdata2(11,m,index1)*wt2
              endif
           enddo
        enddo
     enddo
  endif


  if (LIS_FORC_Pardr%selectOpt.eq.1) then
     call ESMF_FieldGet(pardrField,localDE=0,farrayPtr=pardr,rc=status)
     call LIS_verify(status)
     
     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,geos5fcst_struc(n)%max_ens_members
           do k=1,mfactor
              tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
              index1 = LIS_domain(n)%tile(tid)%index 
              if((geos5fcst_struc(n)%metdata1(12,m,index1).ne.LIS_rc%udef).and.&
                   (geos5fcst_struc(n)%metdata2(12,m,index1).ne.LIS_rc%udef)) then 
                 pardr(tid) =geos5fcst_struc(n)%metdata1(12,m,index1)*wt1+ & 
                      geos5fcst_struc(n)%metdata2(12,m,index1)*wt2
              endif
           enddo
        enddo
     enddo
  endif
  if (LIS_FORC_Pardf%selectOpt.eq.1) then
     call ESMF_FieldGet(pardfField,localDE=0,farrayPtr=pardf,rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,geos5fcst_struc(n)%max_ens_members
           do k=1,mfactor
              tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
              index1 = LIS_domain(n)%tile(tid)%index 
              if((geos5fcst_struc(n)%metdata1(13,m,index1).ne.LIS_rc%udef).and.&
                   (geos5fcst_struc(n)%metdata2(13,m,index1).ne.LIS_rc%udef)) then 
                 pardf(tid) =geos5fcst_struc(n)%metdata1(13,m,index1)*wt1+ & 
                      geos5fcst_struc(n)%metdata2(13,m,index1)*wt2
              endif
           enddo
        enddo
     enddo
  endif


  if (LIS_FORC_Swnet%selectOpt.eq.1) then
     call ESMF_FieldGet(swlandField,localDE=0,farrayPtr=swland,rc=status)
     call LIS_verify(status)

     do t=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
        do m=1,geos5fcst_struc(n)%max_ens_members
           do k=1,mfactor
              tid = (t-1)*LIS_rc%nensem(n)+(m-1)*mfactor+k
              index1 = LIS_domain(n)%tile(tid)%index 
              if((geos5fcst_struc(n)%metdata1(14,m,index1).ne.LIS_rc%udef).and.&
                   (geos5fcst_struc(n)%metdata2(14,m,index1).ne.LIS_rc%udef)) then 
                 swland(tid) =geos5fcst_struc(n)%metdata1(14,m,index1)*wt1+ & 
                      geos5fcst_struc(n)%metdata2(14,m,index1)*wt2
              endif
           enddo
        enddo
     enddo
  endif

end subroutine timeinterp_geos5fcst
 
