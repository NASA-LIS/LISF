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
! !ROUTINE: timeinterp_era5
! \label{timeinterp_era5}
!
! !REVISION HISTORY:
! 23 dec 2019: Sujay Kumar, initial code 
!
! !INTERFACE:

subroutine timeinterp_era5(n,findex)

! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_metforcingMod
  use LIS_FORC_AttributesMod
  use LIS_constantsMod
  use LIS_timeMgrMod
  use LIS_logMod
  use era5_forcingMod
  use LIS_forecastMod
  use LIS_ran2_gasdev


  implicit none

! !ARGUMENTS:
  integer, intent(in):: n
  integer, intent(in):: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation 
!  is not temporally interpolated, and the previous 1 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 1 hourly blocks. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: k,kk,zdoy
  real    :: zw1, zw2
  real    :: czm, cze, czb
  real    :: wt1, wt2,swt1,swt2
  real    :: gmt1, gmt2, tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss
  integer            :: status
  integer            :: mfactor,m
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)


  btime=era5_struc(n)%era5time1
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
  
  btime=era5_struc(n)%era5time2
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
  wt1=(era5_struc(n)%era5time2-LIS_rc%time)/ & 
       (era5_struc(n)%era5time2-era5_struc(n)%era5time1)
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


  mfactor = LIS_rc%nensem(n)/era5_struc(n)%nIter

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        zdoy=LIS_rc%doy
        call zterp(1,LIS_domain(n)%grid(index1)%lat,&
             LIS_domain(n)%grid(index1)%lon,&
             gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

        kk = LIS_get_iteration_index(n, k, index1, mfactor)

        if (era5_struc(n)%metdata1(kk,3,index1).ne.LIS_rc%udef.and.&
             era5_struc(n)%metdata2(kk,3,index1).ne.LIS_rc%udef) then 
           swd(t) = era5_struc(n)%metdata1(kk,3,index1)*zw1 + &
                era5_struc(n)%metdata2(kk,3,index1)*zw2

           ! In cases of small cos(zenith) angles, use linear weighting
           !  to avoid overly large weights

           if((swd(t).gt.era5_struc(n)%metdata1(kk,3,index1).and. & 
                swd(t).gt.era5_struc(n)%metdata2(kk,3,index1)).and. & 
                (czb.lt.0.1.or.cze.lt.0.1))then
              swd(t) = era5_struc(n)%metdata1(kk,3,index1)*swt1+ & 
                   era5_struc(n)%metdata2(kk,3,index1)*swt2
           endif
           
           if (swd(t).gt.LIS_CONST_SOLAR) then
              write(unit=LIS_logunit,fmt=*) &
                   '[WARN] sw radiation too high in ERA5!!'
              write(unit=LIS_logunit,fmt=*)'[WARN] it is',swd(t)
              write(unit=LIS_logunit,fmt=*)'[WARN] era5data1=',&
                   era5_struc(n)%metdata1(kk,3,index1)
              write(unit=LIS_logunit,fmt=*)'[WARN] era5data2=',&
                   era5_struc(n)%metdata2(kk,3,index1)
              write(unit=LIS_logunit,fmt=*)'[WARN] zw1=',zw1,'zw2=',zw2
              swd(t) = LIS_CONST_SOLAR
              write(unit=LIS_logunit,fmt=*)'[WARN] forcing set to ',swd(t) 
           endif
        endif

        if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
           if (swd(t).gt.-0.00001) then 
              swd(t) = 0.0 
           else
              write(LIS_logunit,*) &
                   '[ERR] timeinterp_era5 -- Stopping because ', & 
                   'forcing not udef but lt0,'
              write(LIS_logunit,*)'[ERR] timeinterp_era5 -- ', & 
                   t,swd(t),era5_struc(n)%metdata2(kk,3,index1), & 
                   ' (',LIS_localPet,')'
              call LIS_endrun
           endif
        endif
        
        if (swd(t).gt.LIS_CONST_SOLAR) then
           swd(t)=era5_struc(n)%metdata2(kk,3,index1)
        endif
     enddo
  enddo
 
!-----------------------------------------------------------------------
! precip variable - constant rate over the ERA5 hour
!-----------------------------------------------------------------------  
  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)
        pcp(t) = era5_struc(n)%metdata1(kk,7,index1)

        if ( pcp(t) < 0 ) then
           pcp(t) = 0
        endif
     enddo
  enddo

!-----------------------------------------------------------------------
! Linearly interpolate everything else 
!-----------------------------------------------------------------------

  do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)

        tmp(t) = era5_struc(n)%metdata1(kk,1,index1)*wt1 + &
             era5_struc(n)%metdata2(kk,1,index1)*wt2
        q2(t)  = era5_struc(n)%metdata1(kk,2,index1)*wt1 + &
             era5_struc(n)%metdata2(kk,2,index1)*wt2
        lwd(t)  = era5_struc(n)%metdata1(kk,4,index1)*wt1 + &
             era5_struc(n)%metdata2(kk,4,index1)*wt2

        uwind(t) = era5_struc(n)%metdata1(kk,5,index1)*wt1+&
             era5_struc(n)%metdata2(kk,5,index1)*wt2
        vwind(t) = 0.0
        psurf(t) = era5_struc(n)%metdata1(kk,6,index1)*wt1 + &
             era5_struc(n)%metdata2(kk,6,index1)*wt2

     enddo
  enddo

end subroutine timeinterp_era5
  
