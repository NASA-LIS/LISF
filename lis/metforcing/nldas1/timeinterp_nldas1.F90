!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_nldas1
! \label{timeinterp_nldas1}
!
! !REVISION HISTORY:
!
! 02Feb2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine timeinterp_nldas1(n, findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_constantsMod, only  : LIS_CONST_SOLAR
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use LIS_timeMgrMod, only : LIS_tick, LIS_time2date
  use LIS_logMod, only :LIS_logunit, LIS_verify
  use nldas1_forcingMod, only : nldas1_struc
 
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
  real :: zw1, zw2
  real :: czm, cze, czb
  real :: wt1, wt2,swt1,swt2
  real :: gmt1, gmt2, tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  
  btime=nldas1_struc(n)%nldas1time1
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
  tempbts=900
  call LIS_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=nldas1_struc(n)%nldas1time2
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
  tempbts=900
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(nldas1_struc(n)%nldas1time2-LIS_rc%time)/ & 
       (nldas1_struc(n)%nldas1time2-nldas1_struc(n)%nldas1time1)
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
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     zdoy=LIS_rc%doy
     !  compute and apply zenith angle weights
     call zterp(1,LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)
     
     if(nldas1_struc(n)%metdata1(3,index1).ne.LIS_rc%udef.and.&
          nldas1_struc(n)%metdata2(3,index1).ne.LIS_rc%udef) then 
        swd(t) = nldas1_struc(n)%metdata1(3,index1)*zw1+&
             nldas1_struc(n)%metdata2(3,index1)*zw2
        
           !       In cases of small cos(zenith) angles, use linear weighting
           !       to avoid overly large weights
           
        if((swd(t).gt.nldas1_struc(n)%metdata1(3,index1).and. & 
             swd(t).gt.nldas1_struc(n)%metdata2(3,index1)).and. & 
             (czb.lt.0.1.or.cze.lt.0.1))then
           swd(t) =nldas1_struc(n)%metdata1(3,index1)*swt1+ & 
                nldas1_struc(n)%metdata2(3,index1)*swt2
        endif
     endif
     if (swd(t).gt.LIS_CONST_SOLAR) then
        write(unit=LIS_logunit,fmt=*)'warning, sw radiation too high!!'
        write(unit=LIS_logunit,fmt=*)'it is', swd(t)
        write(unit=LIS_logunit,fmt=*)'ncepdata1=',nldas1_struc(n)%metdata1(3,index1)
        write(unit=LIS_logunit,fmt=*)'ncepdata2=',nldas1_struc(n)%metdata2(3,index1)
        write(unit=LIS_logunit,fmt=*)'zw1=',zw1,'zw2=',zw2
        write(unit=LIS_logunit,fmt=*)'swt1=',swt1,'swt2=',swt2
        
        swd(t) =nldas1_struc(n)%metdata1(3,index1)*swt1+ & 
             nldas1_struc(n)%metdata2(3,index1)*swt2
        write(unit=LIS_logunit,fmt=*)'forcing set to ',swd(t) 
     endif
  enddo

  !do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if(nldas1_struc(n)%metdata2(8,index1).ne.LIS_rc%udef) then 
        pcp(t)=nldas1_struc(n)%metdata2(8,index1)
        pcp(t)  = pcp(t)/(60.0*60.0)
     endif
  enddo

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if(nldas1_struc(n)%metdata2(9,index1).ne.LIS_rc%udef) then 
        cpcp(t)=nldas1_struc(n)%metdata2(9,index1)	
        cpcp(t) = cpcp(t)/(60.0*60.0)	
     endif
  enddo

  !linearly interpolate everything else
  
  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(1,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then 
        tmp(t) =nldas1_struc(n)%metdata1(1,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(1,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(2,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then 
        q2(t) =nldas1_struc(n)%metdata1(2,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(2,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(4,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(4,index1).ne.LIS_rc%udef)) then 
        lwd(t) =nldas1_struc(n)%metdata1(4,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(4,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(5,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(5,index1).ne.LIS_rc%udef)) then 
        uwind(t) =nldas1_struc(n)%metdata1(5,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(5,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(6,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(6,index1).ne.LIS_rc%udef)) then 
        vwind(t) =nldas1_struc(n)%metdata1(6,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(6,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((nldas1_struc(n)%metdata1(7,index1).ne.LIS_rc%udef).and.&
          (nldas1_struc(n)%metdata2(7,index1).ne.LIS_rc%udef)) then 
        psurf(t) =nldas1_struc(n)%metdata1(7,index1)*wt1+ & 
             nldas1_struc(n)%metdata2(7,index1)*wt2
     endif
  enddo
end subroutine timeinterp_nldas1
 
