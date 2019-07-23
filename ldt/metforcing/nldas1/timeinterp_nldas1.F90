!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
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
  use LDT_coreMod,       only : LDT_rc,LDT_domain
  use LDT_constantsMod,  only : LDT_CONST_SOLAR
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_FORC_Base_State, LDT_forc
  use LDT_timeMgrMod,    only : LDT_tick, LDT_time2date
  use LDT_logMod,        only : LDT_logunit, LDT_verify
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
!   \item[LDT\_time2date](\ref{LDT_time2date}) \newline
!    converts the time to a date format
!   \item[LDT\_tick](\ref{LDT_tick}) \newline
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
  call LDT_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  
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
  call LDT_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=nldas1_struc(n)%nldas1time2
  call LDT_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
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
  call LDT_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(nldas1_struc(n)%nldas1time2-LDT_rc%time)/ & 
       (nldas1_struc(n)%nldas1time2-nldas1_struc(n)%nldas1time1)
  wt2=1.0-wt1
  swt1=(newtime2-LDT_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable CRainf in the forcing variables list')
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     zdoy=LDT_rc%doy
     !  compute and apply zenith angle weights
     call zterp(1,LDT_domain(n)%grid(index1)%lat,&
          LDT_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LDT_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LDT_rc)
     
     if(LDT_forc(n,findex)%metdata1(3,index1).ne.LDT_rc%udef.and.&
          LDT_forc(n,findex)%metdata2(3,index1).ne.LDT_rc%udef) then 
        swd(t) = LDT_forc(n,findex)%metdata1(3,index1)*zw1+&
             LDT_forc(n,findex)%metdata2(3,index1)*zw2
        
           !       In cases of small cos(zenith) angles, use linear weighting
           !       to avoid overly large weights
           
        if((swd(t).gt.LDT_forc(n,findex)%metdata1(3,index1).and. & 
             swd(t).gt.LDT_forc(n,findex)%metdata2(3,index1)).and. & 
             (czb.lt.0.1.or.cze.lt.0.1))then
           swd(t) =LDT_forc(n,findex)%metdata1(3,index1)*swt1+ & 
                LDT_forc(n,findex)%metdata2(3,index1)*swt2
        endif
     endif
     if (swd(t).gt.LDT_CONST_SOLAR) then
        write(unit=LDT_logunit,fmt=*)'warning, sw radiation too high!!'
        write(unit=LDT_logunit,fmt=*)'it is', swd(t)
        write(unit=LDT_logunit,fmt=*)'ncepdata1=',LDT_forc(n,findex)%metdata1(3,index1)
        write(unit=LDT_logunit,fmt=*)'ncepdata2=',LDT_forc(n,findex)%metdata2(3,index1)
        write(unit=LDT_logunit,fmt=*)'zw1=',zw1,'zw2=',zw2
        write(unit=LDT_logunit,fmt=*)'swt1=',swt1,'swt2=',swt2
        
        swd(t) =LDT_forc(n,findex)%metdata1(3,index1)*swt1+ & 
             LDT_forc(n,findex)%metdata2(3,index1)*swt2
        write(unit=LDT_logunit,fmt=*)'forcing set to ',swd(t) 
     endif
  enddo

  !do block precipitation interpolation
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LDT_verify(status)


  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if(LDT_forc(n,findex)%metdata2(8,index1).ne.LDT_rc%udef) then 
        pcp(t)=LDT_forc(n,findex)%metdata2(8,index1)
        pcp(t)  = pcp(t)/(60.0*60.0)
     endif
  enddo

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LDT_verify(status)
  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if(LDT_forc(n,findex)%metdata2(9,index1).ne.LDT_rc%udef) then 
        cpcp(t)=LDT_forc(n,findex)%metdata2(9,index1)	
        cpcp(t) = cpcp(t)/(60.0*60.0)	
     endif
  enddo

  !linearly interpolate everything else
  
  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LDT_verify(status)
  
  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(1,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(1,index1).ne.LDT_rc%udef)) then 
        tmp(t) =LDT_forc(n,findex)%metdata1(1,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(1,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LDT_verify(status)
  
  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(2,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(2,index1).ne.LDT_rc%udef)) then 
        q2(t) =LDT_forc(n,findex)%metdata1(2,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(2,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(4,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(4,index1).ne.LDT_rc%udef)) then 
        lwd(t) =LDT_forc(n,findex)%metdata1(4,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(4,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(5,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(5,index1).ne.LDT_rc%udef)) then 
        uwind(t) =LDT_forc(n,findex)%metdata1(5,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(5,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(6,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(6,index1).ne.LDT_rc%udef)) then 
        vwind(t) =LDT_forc(n,findex)%metdata1(6,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(6,index1)*wt2
     endif
  enddo

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index 
     if((LDT_forc(n,findex)%metdata1(7,index1).ne.LDT_rc%udef).and.&
          (LDT_forc(n,findex)%metdata2(7,index1).ne.LDT_rc%udef)) then 
        psurf(t) =LDT_forc(n,findex)%metdata1(7,index1)*wt1+ & 
             LDT_forc(n,findex)%metdata2(7,index1)*wt2
     endif
  enddo
end subroutine timeinterp_nldas1
 
